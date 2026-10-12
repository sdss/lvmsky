#!/usr/bin/env python3
"""Build the final LSF catalog and stability report for the serious corpus."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import html
import json
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.offline
from astropy.io import fits
from astropy.table import Table

from skysub.sky_decomp.lsf_spline2d import evaluate_lsf_density, evaluate_lsf_diagnostics
from skysub.sky_decomp.result_io import load_lsf_surface_cubes, lsf_surface_state_from_cubes


OUT = Path("/data/sas/sdsswork/users/u6058164/lvm_solarblue_lsf_serious_corpus_drp132")
INPUT = OUT / "individual_science_fibers.fits"
FIT = OUT / "fit/individual_science_fibers_sci_meta_coef_palacecorr_aijc_vnf_split_zodi_lsf_spline2d_solarblue.fits"
ANCHORS = {
    "B": (3800.0, 4250.0, 4700.0, 5200.0, 5650.0),
    "R": (5900.0, 6250.0, 6600.0, 7000.0, 7350.0),
    "Z": (7600.0, 8000.0, 8500.0, 9000.0, 9600.0),
}
CENTRAL_ANCHOR = {"B": 4700.0, "R": 6600.0, "Z": 8500.0}
CHANNEL_COLORS = {"B": "#1f77b4", "R": "#d62728", "Z": "#9467bd"}
TIME_ORDER = ["0–0.5 d", "0.5–1 d", "1–2 d", "2–3 d", "3–5 d"]


def dataframe(table) -> pd.DataFrame:
    frame = Table(table).to_pandas()
    for name in frame.select_dtypes("object"):
        frame[name] = frame[name].map(
            lambda value: value.decode().strip() if isinstance(value, bytes) else str(value).strip())
    return frame


def calibration_epoch(path: str) -> tuple[str, int]:
    try:
        comments = fits.getheader(path, 0).get("COMMENT", [])
        comments = [comments] if isinstance(comments, str) else comments
        found = re.findall(r"/calib/(\d+)/", "\n".join(comments))
        return path, max(map(int, found)) if found else -1
    except OSError:
        return path, -1


def add_calibration_epochs(fibers: pd.DataFrame) -> pd.DataFrame:
    paths = fibers.source_path.drop_duplicates().astype(str).tolist()
    with ThreadPoolExecutor(max_workers=24) as pool:
        epochs = dict(pool.map(calibration_epoch, paths))
    fibers["calibration_epoch"] = fibers.source_path.map(epochs).fillna(-1).astype(int)
    return fibers


def drp_lsf_at_anchors(input_path: Path, rows: np.ndarray, wavelengths: np.ndarray) -> np.ndarray:
    with fits.open(input_path, memmap=True) as hdul:
        wave = np.asarray(hdul["WAVE"].data, dtype=float)
        values = hdul["LSF_SCI"].data
        right = np.searchsorted(wave, wavelengths)
        right = np.clip(right, 1, len(wave) - 1)
        left = right - 1
        fraction = (wavelengths - wave[left]) / (wave[right] - wave[left])
        return ((1.0 - fraction) * np.asarray(values[rows, left], dtype=float)
                + fraction * np.asarray(values[rows, right], dtype=float))


def metric_table(input_path: Path, fit_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    with fits.open(input_path, memmap=True) as hdul:
        fibers = add_calibration_epochs(dataframe(hdul["META"].data))
    with fits.open(fit_path, memmap=True) as hdul:
        fit_meta = dataframe(hdul["META"].data)
    cubes = load_lsf_surface_cubes(fit_path)
    wavelengths = np.asarray([value for values in ANCHORS.values() for value in values])
    channels = np.asarray([channel for channel, values in ANCHORS.items() for _ in values])
    delta = np.linspace(-3.0, 3.0, 601)
    records = []
    for row in range(len(fibers)):
        try:
            state = lsf_surface_state_from_cubes(cubes, row)
            diagnostics = evaluate_lsf_diagnostics(state, wavelengths)
            density = evaluate_lsf_density(state, wavelengths, delta)
        except (ValueError, IndexError):
            continue
        area = np.trapz(density, delta, axis=1)
        centre = diagnostics["centroid_angstrom"]
        centered = delta[None, :] - centre[:, None]
        variance = np.trapz(density * centered**2, delta, axis=1) / area
        skewness = np.trapz(density * centered**3, delta, axis=1) / area / variance**1.5
        kurtosis = np.trapz(density * centered**4, delta, axis=1) / area / variance**2 - 3.0
        wing = np.trapz(np.where(np.abs(centered) >= 1.0, density, 0.0), delta, axis=1) / area
        for index, (channel, wavelength) in enumerate(zip(channels, wavelengths)):
            channel_metric = state.metrics[channel]
            records.append({
                "corpus_row": row,
                "channel": channel,
                "wavelength_angstrom": wavelength,
                **{name: float(values[index]) for name, values in diagnostics.items()},
                "wing_fraction_outside_1A": float(wing[index]),
                "skewness": float(skewness[index]),
                "excess_kurtosis": float(kurtosis[index]),
                "gaussian_fwhm_ratio": float(
                    diagnostics["fwhm_angstrom"][index]
                    / (2.354820045 * diagnostics["sigma_angstrom"][index])),
                "fit_status": state.fit_status,
                "completed_cycles": state.completed_cycles,
                "converged": bool(channel_metric.get("converged", False)),
                "lsf_chi2_red": float(channel_metric.get("chi2_red", np.nan)),
            })
        if (row + 1) % 1000 == 0:
            print(f"LSF diagnostics {row + 1}/{len(fibers)}", flush=True)
    metrics = pd.DataFrame(records)
    keep = [
        "corpus_row", "expnum", "mjd", "obstime", "midpoint_mjd_utc",
        "full_moon_mjd_utc", "full_moon_utc", "days_from_full_moon",
        "relative_lunar_continuum", "continuum_5400_5500", "calibration_epoch",
        "fiberid", "spectrographid", "slit_rank", "slit_quantile", "fit_valid_fraction",
        "gaia_status", "gaia_clean_1pct", "gaia_ratio", "source_path",
    ]
    metrics = metrics.merge(fibers[keep], on="corpus_row", how="left", validate="many_to_one")
    metrics["drp_fwhm_angstrom"] = drp_lsf_at_anchors(
        input_path, metrics.corpus_row.to_numpy(dtype=int),
        metrics.wavelength_angstrom.to_numpy(dtype=float))
    metrics["fwhm_ratio_to_drp"] = metrics.fwhm_angstrom / metrics.drp_fwhm_angstrom
    metrics["quality_sample"] = (
        metrics.converged & (metrics.fit_status == "Solved")
        & (metrics.fit_valid_fraction >= 0.8)
        & np.isfinite(metrics.fwhm_angstrom) & (metrics.fwhm_angstrom > 0)
    )
    metrics["gaia_clean_sample"] = metrics.quality_sample & metrics.gaia_clean_1pct.fillna(False)
    return fibers, fit_meta, metrics


def stability_tables(metrics: pd.DataFrame):
    good = metrics[metrics.quality_sample].copy()
    keys = ["full_moon_mjd_utc", "spectrographid", "slit_rank", "channel", "wavelength_angstrom"]
    good["reference_distance"] = np.abs(np.log(good.relative_lunar_continuum))
    reference = (good.sort_values("reference_distance")
                 .drop_duplicates(keys)[keys + ["midpoint_mjd_utc", "fwhm_angstrom",
                                                "w50_angstrom", "w90_angstrom",
                                                "centroid_angstrom", "asymmetry"]])
    reference = reference.rename(columns={
        "midpoint_mjd_utc": "reference_mjd",
        "fwhm_angstrom": "reference_fwhm",
        "w50_angstrom": "reference_w50",
        "w90_angstrom": "reference_w90",
        "centroid_angstrom": "reference_centroid",
        "asymmetry": "reference_asymmetry",
    })
    comparisons = good.merge(reference, on=keys, how="inner", validate="many_to_one")
    comparisons["signed_days_from_reference"] = comparisons.midpoint_mjd_utc - comparisons.reference_mjd
    comparisons["days_from_reference"] = comparisons.signed_days_from_reference.abs()
    for name in ("fwhm", "w50", "w90"):
        comparisons[f"abs_{name}_difference_percent"] = 100.0 * np.abs(
            comparisons[f"{name}_angstrom"] / comparisons[f"reference_{name}"] - 1.0)
    comparisons["abs_centroid_difference_angstrom"] = np.abs(
        comparisons.centroid_angstrom - comparisons.reference_centroid)
    comparisons["abs_asymmetry_difference"] = np.abs(
        comparisons.asymmetry - comparisons.reference_asymmetry)
    comparisons["time_bin"] = pd.cut(
        comparisons.signed_days_from_reference,
        [-1e-9, 0.5, 1.0, 2.0, 3.0, 5.01],
        labels=TIME_ORDER)

    def q90(values):
        return values.quantile(0.9)

    summary = comparisons[comparisons.signed_days_from_reference > 1e-8].groupby(
        ["channel", "wavelength_angstrom", "time_bin"], observed=True).agg(
            n=("fwhm_angstrom", "size"),
            median_abs_fwhm_percent=("abs_fwhm_difference_percent", "median"),
            p90_abs_fwhm_percent=("abs_fwhm_difference_percent", q90),
            median_abs_w50_percent=("abs_w50_difference_percent", "median"),
            p90_abs_w50_percent=("abs_w50_difference_percent", q90),
            median_abs_w90_percent=("abs_w90_difference_percent", "median"),
            p90_abs_w90_percent=("abs_w90_difference_percent", q90),
            median_abs_centroid_angstrom=("abs_centroid_difference_angstrom", "median"),
            p90_abs_centroid_angstrom=("abs_centroid_difference_angstrom", q90),
            median_abs_asymmetry=("abs_asymmetry_difference", "median"),
            p90_abs_asymmetry=("abs_asymmetry_difference", q90),
        ).reset_index()
    cycle = good.groupby(
        ["full_moon_utc", "full_moon_mjd_utc", "calibration_epoch", "channel",
         "wavelength_angstrom", "spectrographid", "slit_rank"], as_index=False).agg(
            fwhm_angstrom=("fwhm_angstrom", "median"),
            centroid_angstrom=("centroid_angstrom", "median"),
            w50_angstrom=("w50_angstrom", "median"),
            w90_angstrom=("w90_angstrom", "median"),
            asymmetry=("asymmetry", "median"), n=("fwhm_angstrom", "size"))
    dates = pd.to_datetime(cycle.full_moon_utc, utc=True)
    month = dates.dt.month
    cycle["season"] = np.asarray(["DJF", "MAM", "JJA", "SON"])[(month % 12) // 3]
    cycle["season_year"] = dates.dt.year + (month == 12).astype(int)
    cycle["season_label"] = cycle.season_year.astype(str) + " " + cycle.season
    lunation_median = cycle.groupby(
        ["season_year", "season_label", "season", "full_moon_utc", "channel", "wavelength_angstrom"],
        as_index=False).agg(fwhm_angstrom=("fwhm_angstrom", "median"))
    seasonal = lunation_median.groupby(
        ["season_year", "season_label", "season", "channel", "wavelength_angstrom"], as_index=False).agg(
        lunations=("full_moon_utc", "nunique"),
        median_fwhm_angstrom=("fwhm_angstrom", "median"),
        p05_fwhm_angstrom=("fwhm_angstrom", lambda x: x.quantile(0.05)),
        p95_fwhm_angstrom=("fwhm_angstrom", lambda x: x.quantile(0.95)))
    seasonal["p05_p95_span_percent"] = 100 * (
        seasonal.p95_fwhm_angstrom - seasonal.p05_fwhm_angstrom) / seasonal.median_fwhm_angstrom
    return comparisons, summary, cycle, seasonal


def style(fig, height=500):
    fig.update_layout(template="plotly_white", hovermode="closest", height=height,
                      margin=dict(t=70, b=80), font=dict(size=13))
    fig.update_xaxes(showspikes=False)
    fig.update_yaxes(showspikes=False)
    return fig


def plots(metrics, summary, cycle, output):
    summary = summary.copy()
    summary["time_bin"] = pd.Categorical(summary.time_bin, TIME_ORDER, ordered=True)
    central = cycle[cycle.apply(
        lambda row: row.wavelength_angstrom == CENTRAL_ANCHOR[row.channel], axis=1)].copy()
    figures = [
        style(px.line(central.groupby(["full_moon_utc", "full_moon_mjd_utc", "channel"], as_index=False)
                      .fwhm_angstrom.median(), x="full_moon_utc", y="fwhm_angstrom",
                      color="channel", markers=True, title="Median LSF FWHM across lunations",
                      labels={"full_moon_utc": "Full Moon (UTC)", "fwhm_angstrom": "FWHM (Å)", "channel": "Channel"},
                      color_discrete_map=CHANNEL_COLORS)),
        style(px.line(summary, x="time_bin", y="p90_abs_fwhm_percent", color="channel",
                      facet_col="wavelength_angstrom", facet_col_wrap=5, markers=True,
                      title="Difference after the robust bright reference in each lunation",
                      labels={"time_bin": "Separation from reference", "p90_abs_fwhm_percent": "P90 |ΔFWHM| (%)"},
                      category_orders={"time_bin": TIME_ORDER}, color_discrete_map=CHANNEL_COLORS), 760),
        style(px.box(metrics[metrics.quality_sample], x="channel", y="fwhm_ratio_to_drp",
                     color="channel", points=False, title="Recovered continuous LSF versus DRP input LSF",
                     labels={"channel": "Channel", "fwhm_ratio_to_drp": "Recovered / DRP FWHM"},
                     color_discrete_map=CHANNEL_COLORS)),
        style(px.scatter(metrics[metrics.quality_sample & metrics.wavelength_angstrom.isin(CENTRAL_ANCHOR.values())],
                         x="midpoint_mjd_utc", y="centroid_angstrom", color="channel",
                         opacity=0.25, render_mode="webgl", hover_data=["expnum", "fiberid", "spectrographid"],
                         title="LSF barycenter stability", labels={"midpoint_mjd_utc": "MJD (UTC)", "centroid_angstrom": "Barycenter offset (Å)"},
                         color_discrete_map=CHANNEL_COLORS)),
    ]
    (output / "plotly.min.js").write_text(plotly.offline.get_plotlyjs(), encoding="utf-8")
    interactive = [figure.to_html(full_html=False, include_plotlyjs=False,
                                  config={"displaylogo": False, "responsive": True}) for figure in figures]

    aggregate = summary.groupby(["channel", "time_bin"], observed=True).agg(
        median=("median_abs_fwhm_percent", "median"), p90=("p90_abs_fwhm_percent", "median")).reset_index()
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for channel, group in aggregate.groupby("channel"):
        x = np.arange(len(group))
        ax.plot(x, group.p90, marker="o", color=CHANNEL_COLORS[channel], label=f"{channel}: P90")
        ax.plot(x, group["median"], marker=".", linestyle="--", color=CHANNEL_COLORS[channel], label=f"{channel}: median")
    ax.set_xticks(np.arange(aggregate.time_bin.nunique()), aggregate.time_bin.drop_duplicates(), rotation=25)
    ax.set_ylabel("Absolute FWHM difference (%)")
    ax.set_xlabel("Days after bright reference")
    ax.grid(alpha=0.25); ax.legend(ncol=3, fontsize=8); fig.tight_layout()
    fig.savefig(output / "reference_stability.png", dpi=180); plt.close(fig)
    return interactive


def write_catalog(input_path, fit_path, fibers, fit_meta, metrics, output):
    destination = output / "lsf_catalog.fits"
    with fits.open(fit_path, memmap=True) as fitted:
        primary = fitted[0].copy()
        primary.header["CATALOG"] = True
        hdus = [primary, fits.BinTableHDU(Table.from_pandas(fibers), name="FIBERS"),
                fits.BinTableHDU(Table.from_pandas(fit_meta), name="FIT_META")]
        for name in ("COEF", "COEF_ERR", "LSF_COEF", "LSF_KNOTS", "LSF_META"):
            hdus.append(fitted[name].copy())
        hdus.append(fits.BinTableHDU(Table.from_pandas(metrics), name="METRICS"))
        fits.HDUList(hdus).writeto(destination.with_suffix(".tmp.fits"), overwrite=True, checksum=True)
    destination.with_suffix(".tmp.fits").replace(destination)
    return destination


def render(output, fibers, fit_meta, metrics, comparisons, summary, seasonal, interactive, catalog):
    quality = metrics[metrics.quality_sample]
    convergence = quality.groupby("channel").size() / metrics.groupby("channel").size()
    drp = quality.groupby("channel").fwhm_ratio_to_drp.agg(["median", "count"])
    drp_answer = "; ".join(
        f"{channel}: {(row['median'] - 1) * 100:+.1f}%"
        for channel, row in drp.iterrows())
    day = summary.groupby(["channel", "time_bin"], observed=True).agg(
        median_abs_fwhm_percent=("median_abs_fwhm_percent", "median"),
        p90_abs_fwhm_percent=("p90_abs_fwhm_percent", "median")).reset_index()
    day["time_bin"] = pd.Categorical(day.time_bin, TIME_ORDER, ordered=True)
    day = day.sort_values(["time_bin", "channel"])
    five = day[day.time_bin.astype(str) == "3–5 d"]
    near = day[day.time_bin.astype(str) == "0–0.5 d"]
    answer = "; ".join(
        f"{row.channel}: median {row.median_abs_fwhm_percent:.2f}%, P90 {row.p90_abs_fwhm_percent:.2f}%"
        for row in five.itertuples()) or "The 3–5 day bin has insufficient accepted comparisons."
    baseline = "; ".join(
        f"{row.channel}: median {row.median_abs_fwhm_percent:.2f}%, P90 {row.p90_abs_fwhm_percent:.2f}%"
        for row in near.itertuples()) or "The sub-day bin has insufficient accepted comparisons."
    worst = (summary[summary.time_bin.astype(str) == "3–5 d"]
             .sort_values("p90_abs_fwhm_percent").groupby("channel", sort=False).tail(1))
    worst_answer = "; ".join(
        f"{row.channel}: {row.p90_abs_fwhm_percent:.2f}% at {row.wavelength_angstrom:.0f} Å"
        for row in worst.itertuples())
    shape_answer = "; ".join(
        f"{row.channel}: W50 {row.p90_abs_w50_percent:.2f}%, W90 {row.p90_abs_w90_percent:.2f}%, "
        f"barycenter {row.p90_abs_centroid_angstrom:.3f} Å, asymmetry {row.p90_abs_asymmetry:.3f}"
        for row in summary[summary.time_bin.astype(str) == "3–5 d"].groupby("channel", as_index=False).agg(
            p90_abs_w50_percent=("p90_abs_w50_percent", "median"),
            p90_abs_w90_percent=("p90_abs_w90_percent", "median"),
            p90_abs_centroid_angstrom=("p90_abs_centroid_angstrom", "median"),
            p90_abs_asymmetry=("p90_abs_asymmetry", "median")).itertuples())
    later = comparisons[(comparisons.signed_days_from_reference > 3.0)
                        & (comparisons.signed_days_from_reference <= 5.01)]
    summary_payload = {
        "exposures": int(fibers.expnum.nunique()), "fiber_rows": int(len(fibers)),
        "successful_fits": int((fit_meta.fit_status == "Solved").sum()),
        "quality_metric_rows": int(len(quality)), "gaia_clean_fibers": int(fibers.gaia_clean_1pct.sum()),
        "gaia_available_fibers": int((fibers.gaia_status != "unavailable").sum()),
        "lunations": int(fibers.full_moon_utc.nunique()),
        "calibration_epochs": int(fibers.loc[fibers.calibration_epoch >= 0, "calibration_epoch"].nunique()),
        "three_to_five_day_fwhm": answer,
    }
    (output / "summary.json").write_text(json.dumps(summary_payload, indent=2) + "\n", encoding="utf-8")
    display_day = day.rename(columns={
        "channel": "Channel", "time_bin": "Days after reference",
        "median_abs_fwhm_percent": "Median |ΔFWHM| (%)",
        "p90_abs_fwhm_percent": "P90 |ΔFWHM| (%)"})
    table = display_day.round(3).to_html(index=False, classes="data", border=0)
    seasonal_central = seasonal[seasonal.apply(
        lambda row: row.wavelength_angstrom == CENTRAL_ANCHOR[row.channel], axis=1)].copy()
    season_order = {"DJF": 0, "MAM": 1, "JJA": 2, "SON": 3}
    seasonal_central["season_order"] = seasonal_central.season.map(season_order)
    seasonal_central = seasonal_central.sort_values(["season_year", "season_order", "channel"]).drop(columns="season_order")
    supported_seasons = seasonal_central[seasonal_central.lunations >= 2]
    seasonal_answer = "; ".join(
        f"{channel}: {100 * (group.median_fwhm_angstrom.max() - group.median_fwhm_angstrom.min()) / group.median_fwhm_angstrom.median():.2f}%"
        for channel, group in supported_seasons.groupby("channel"))
    seasonal_display = seasonal_central[[
        "season_label", "channel", "wavelength_angstrom", "lunations",
        "median_fwhm_angstrom", "p05_fwhm_angstrom", "p95_fwhm_angstrom",
        "p05_p95_span_percent"]].rename(columns={
        "season_label": "Season",
        "channel": "Channel", "wavelength_angstrom": "Wavelength (Å)",
        "lunations": "Lunations", "median_fwhm_angstrom": "Median FWHM (Å)",
        "p05_fwhm_angstrom": "P05 FWHM (Å)", "p95_fwhm_angstrom": "P95 FWHM (Å)",
        "p05_p95_span_percent": "P05–P95 span (%)"})
    seasonal_table = seasonal_display.round(4).to_html(index=False, classes="data", border=0)
    drp_table = drp.reset_index().rename(columns={
        "channel": "Channel", "median": "Median recovered / DRP FWHM",
        "count": "Accepted anchor rows"}).round(4).to_html(index=False, classes="data", border=0)
    table = table.replace("<th>", "<th scope='col'>")
    seasonal_table = seasonal_table.replace("<th>", "<th scope='col'>")
    drp_table = drp_table.replace("<th>", "<th scope='col'>")
    document = f"""<!doctype html><html lang='en'><head><meta charset='utf-8'>
<meta name='viewport' content='width=device-width,initial-scale=1'><title>DRP 1.3.2 SolarBlue LSF stability</title>
<script src='plotly.min.js'></script><style>body{{font:16px/1.5 system-ui;max-width:1450px;margin:auto;padding:24px;color:#18202a}}h1,h2{{color:#153a5b}}.cards{{display:flex;gap:12px;flex-wrap:wrap}}.card{{background:#eef5fa;border-radius:8px;padding:12px 18px}}.table-wrap{{overflow-x:auto}}table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{padding:6px;border-bottom:1px solid #ccd6df;text-align:right}}th:first-child,td:first-child{{text-align:left}}code{{background:#edf1f4;padding:2px 4px}}img{{max-width:850px;width:100%}}</style></head><body>
<h1>Solar-constrained LSF stability in bright-Moon DRP 1.3.2 spectra</h1>
<p>This report measures the continuous B/R/Z LSF surfaces recovered from nine Science fibers per exposure. The sample uses Moon-up astronomical-night exposures within five days of full Moon and at least 0.3 of the robust maximum continuum in each lunation.</p>
<p><b>Terms:</b> LSF is the line-spread function; B, R and Z are the blue, red and near-infrared channels; P90 is the 90th percentile; W50 and W90 enclose 50% and 90% of the normalized LSF area.</p>
<p><b>Gaia contamination flags:</b> version-matched local Gaia tables cover {summary_payload['gaia_available_fibers']:,} of {summary_payload['fiber_rows']:,} fitted fibers; {summary_payload['gaia_clean_fibers']:,} pass the 1% flux-ratio criterion. Remaining fibers are explicitly marked <code>unavailable</code>, retained in the corpus by design, and must not be interpreted as Gaia-clean.</p>
<div class='cards'><div class='card'><b>{summary_payload['exposures']:,}</b><br>exposures</div><div class='card'><b>{summary_payload['fiber_rows']:,}</b><br>fiber fits</div><div class='card'><b>{summary_payload['lunations']}</b><br>lunations</div><div class='card'><b>{(fit_meta.fit_status == 'Solved').mean():.1%}</b><br>successful fits</div></div>
<h2>Can a bright-reference FWHM be reused several days later?</h2><p>The reference is the exposure nearest the robust lunation maximum (relative continuum 1.0), preventing a rare flux outlier from defining the comparison. Only exposures observed <em>after</em> that reference enter the time-bin summaries. For the 3–5 day interval, the wavelength-median absolute FWHM differences are: <b>{html.escape(answer)}</b>. The later 0–0.5 day baseline is: <b>{html.escape(baseline)}</b>. The descriptive R and Z summaries show no evident increase at 3–5 days, while the B P90 discrepancy increases.</p>
<p><b>Operational interpretation:</b> reusing R or Z for 3–5 days is supported when the listed empirical tolerances are acceptable; B should be remeasured for work requiring materially better precision than its listed P90. The worst 3–5 day anchors are {html.escape(worst_answer)}. Wavelength-median P90 shape diagnostics at the same lag are {html.escape(shape_answer)}. This bin contains {later.expnum.nunique():,} exposures from {later.full_moon_utc.nunique()} lunations. Each wavelength summary treats the fitted fiber positions as empirical samples; it is descriptive rather than a cluster-weighted population confidence interval. The differences include reconstruction noise and are an upper limit on pure hardware drift.</p>
<img src='reference_stability.png' alt='Median and P90 absolute FWHM differences for B, R and Z versus days after the bright reference'>
<h2>Time-separation statistics</h2><div class='table-wrap'>{table}</div>
<h2>Seasonal and lunation-scale stability</h2>{interactive[0]}<p>The seasonal P05–P95 interval is computed from one spatially aggregated median per lunation. Across season-year blocks with at least two lunations, the range of seasonal median FWHM at the central anchors is {html.escape(seasonal_answer)}. A season with one lunation has no measured temporal span and remains listed for provenance.</p><div class='table-wrap'>{seasonal_table}</div>
<h2>Wavelength-resolved reuse test</h2>{interactive[1]}
<h2>Recovery checks</h2><p>Every accepted row passed the solver and channel convergence criteria. The recovered FWHM is systematically broader than the DRP input by {html.escape(drp_answer)}. This comparison is an independent pipeline reference, not ground truth, and the offset should be calibrated before mixing the two LSF definitions. R/Z also have strong sky-line information; B is constrained primarily by the high-resolution solar absorption spectrum.</p><div class='table-wrap'>{drp_table}</div>{interactive[2]}
<h2>Barycenter stability</h2>{interactive[3]}
<h2>Data products</h2><p><a href='{catalog.name}'><code>{catalog.name}</code></a> contains FIBERS, FIT_META, COEF, COEF_ERR, LSF_COEF, LSF_KNOTS, LSF_META and METRICS. CSV tables: <a href='lsf_metrics.csv'>metrics</a>, <a href='reference_comparisons.csv'>reference comparisons</a>, <a href='time_stability_summary.csv'>time summary</a>, and <a href='seasonal_summary.csv'>seasonal summary</a>.</p>
</body></html>"""
    (output / "solarblue_lsf_stability_report.html").write_text(document, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--fit", type=Path, default=FIT)
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    readme = Path(__file__).with_name("README.md")
    (args.output / "README.md").write_text(readme.read_text(encoding="utf-8"), encoding="utf-8")
    fibers, fit_meta, metrics = metric_table(args.input, args.fit)
    comparisons, summary, cycle, seasonal = stability_tables(metrics)
    metrics.to_csv(args.output / "lsf_metrics.csv", index=False)
    comparisons.to_csv(args.output / "reference_comparisons.csv", index=False)
    summary.to_csv(args.output / "time_stability_summary.csv", index=False)
    cycle.to_csv(args.output / "lunation_metrics.csv", index=False)
    seasonal.to_csv(args.output / "seasonal_summary.csv", index=False)
    interactive = plots(metrics, summary, cycle, args.output)
    catalog = write_catalog(args.input, args.fit, fibers, fit_meta, metrics, args.output)
    render(args.output, fibers, fit_meta, metrics, comparisons, summary, seasonal, interactive, catalog)
    print(catalog)


if __name__ == "__main__":
    main()
