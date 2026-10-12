"""Actual LSF surfaces and wavelength diagnostics for geometrically selected fibers."""
from pathlib import Path
import html
import json

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
import mistune

from skysub.sky_decomp.lsf_reference import load_reference_lsf
from skysub.sky_decomp.lsf_spline2d import evaluate_lsf_density, evaluate_lsf_diagnostics, mspline_bin_integrals, evaluate_bspline_basis

CHANNELS = ("B", "R", "Z")
METRICS = ("fwhm_angstrom", "w50_angstrom", "w90_angstrom", "centroid_angstrom", "asymmetry")
LABELS = ("FWHM (Å)", "W50 (Å)", "W90 (Å)", "Centroid (Å)", "Asymmetry")


def select_fibers(catalog, destination, records):
    """Nearest working fibers at the lower, middle and upper observed slit positions."""
    survey = catalog[catalog.scope == "survey"]
    with fits.open(destination / survey.file.iloc[0]) as hdus:
        source = str(hdus["SOURCES"].data[0]["path"]).strip()
    with fits.open(source) as hdus:
        geometry = Table(hdus["INPUT_META"].data)[["fiberid", "telescope", "spectrographid", "ypix_b"]].to_pandas()
    geometry = geometry.drop_duplicates("fiberid")
    complete = survey.groupby("fiberid").channel.nunique()
    geometry = geometry[geometry.fiberid.isin(complete[complete == 3].index)]
    rows = []
    for (telescope, sp), group in geometry.groupby(["telescope", "spectrographid"]):
        available = group.sort_values("ypix_b")
        if len(available) < 3 or not np.isfinite(available.ypix_b).all():
            raise ValueError("Three finite working slit positions are required per group")
        lower, upper = available.ypix_b.min(), available.ypix_b.max()
        used = set()
        for label, fraction in (("lower", 0.), ("middle", .5), ("upper", 1.)):
            remaining = available[~available.fiberid.isin(used)]
            row = remaining.loc[(remaining.ypix_b-(lower+fraction*(upper-lower))).abs().idxmin()]
            used.add(int(row.fiberid))
            rows.append(dict(telescope=str(telescope), spectrographid=int(sp), fiberid=int(row.fiberid),
                             slit_position=label, ypix_b=float(row.ypix_b), target_fraction=fraction))
    selected = pd.DataFrame(rows)
    assert len(selected) == 27 and not selected.fiberid.duplicated().any()
    selected.to_csv(records / "selected_fibers.csv", index=False)
    return selected


def wavelength_grid(state, channel):
    knots, degree = state.knot_vectors[channel], state.degrees[channel]
    return np.linspace(knots[degree], knots[-degree-1], 161)


def render_fiber(row, catalog, destination, views):
    fiber = int(row.fiberid)
    refs = catalog[catalog.fiberid == fiber]
    survey_path = destination / refs[refs.scope == "survey"].file.iloc[0]
    state = load_reference_lsf(survey_path, fiber, allow_epoch_pooling=True)
    waves = {c: wavelength_grid(state, c) for c in CHANNELS}
    delta = np.linspace(-2.5, 2.5, 161)
    density = {c: evaluate_lsf_density(state, waves[c], delta) for c in CHANNELS}
    assert all(np.isfinite(z).all() and z.min() >= 0 for z in density.values())
    bins = mspline_bin_integrals(np.array([-2.5]), np.array([2.5]), n_basis=11, half_width=2.5)[0]
    for coefficient in state.coefficients.values():
        np.testing.assert_allclose(bins @ coefficient, 1, atol=2e-12)
    vmax = max(z.max() for z in density.values())
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ci, channel in enumerate(CHANNELS):
        ax = axes[ci]
        image = ax.pcolormesh(waves[channel], delta, density[channel].T,
                              shading="auto", cmap="viridis", vmin=0, vmax=vmax)
        ax.set(xlabel="Wavelength (Å)", ylabel="Offset Δλ (Å)", title=channel)
        ax.xaxis.set_major_locator(MaxNLocator(4))
    label = f"{row.telescope} / SP{row.spectrographid} / fiber {fiber} / {row.slit_position} slit"
    fig.suptitle(label+"\nSurvey-pooled LSF: descriptive average across calibration epochs", fontsize=13)
    fig.subplots_adjust(left=.055, right=.88, bottom=.15, top=.77, wspace=.28)
    colorbar_ax = fig.add_axes([.91, .15, .012, .62])
    fig.colorbar(image, cax=colorbar_ax, label="LSF density (Å⁻¹)")
    colorbar_ax.tick_params(labelsize=9)
    surface_name = f"fiber_{fiber}_surface.png"
    fig.savefig(views / surface_name, dpi=135)
    plt.close(fig)
    fig, axes = plt.subplots(5, 3, figsize=(14, 13), sharex="col")
    periods = refs[refs.scope == "lunation"].groupby("file", sort=True)
    dates = refs[refs.scope == "lunation"].weighted_mjd.to_numpy()
    norm = Normalize(float(dates.min()), float(dates.max()))
    records = []
    cases = [(survey_path, "survey")] + [(destination / name, "lunation") for name, _ in periods]
    for path, scope in cases:
        selected = refs[refs.file == str(path.relative_to(destination))]
        if set(selected.channel) != set(CHANNELS):
            raise ValueError("Selected fiber requires B/R/Z coverage for wavelength views")
        current = state if scope == "survey" else load_reference_lsf(path, fiber)
        for ci, channel in enumerate(CHANNELS):
            wave = waves[channel]
            diagnostic = evaluate_lsf_diagnostics(current, wave)
            metadata = selected[selected.channel == channel].iloc[0]
            color = "black" if scope == "survey" else plt.get_cmap("viridis")(norm(float(metadata.weighted_mjd)))
            assert all(np.isfinite(diagnostic[m]).all() for m in METRICS)
            anchor = dict(B=4700., R=6600., Z=8500.)[channel]
            check = evaluate_lsf_diagnostics(current, np.array([anchor]))
            for metric in METRICS[:4]:
                np.testing.assert_allclose(check[metric][0], metadata[metric+"_reference"], atol=1e-9, rtol=1e-9)
            for mi, metric in enumerate(METRICS):
                ax = axes[mi, ci]
                ax.plot(wave, diagnostic[metric], color=color, lw=2 if scope == "survey" else .75,
                        alpha=1 if scope == "survey" else .7, zorder=5 if scope == "survey" else 1)
                ax.set_ylabel(LABELS[mi])
                ax.grid(alpha=.18)
            table = pd.DataFrame({m: diagnostic[m] for m in METRICS})
            table["wavelength_angstrom"], table["fiberid"], table["channel"] = wave, fiber, channel
            table["scope"], table["file"] = scope, str(path.relative_to(destination))
            table["calibration_epoch"], table["full_moon_utc"] = int(metadata.calibration_epoch), metadata.full_moon_utc
            table["n_exposures"], table["n_nights"], table["weighted_mjd"] = int(metadata.n_exposures), int(metadata.n_nights), float(metadata.weighted_mjd)
            records.append(table)
    for ci, channel in enumerate(CHANNELS):
        axes[0, ci].set_title(channel)
        axes[-1, ci].set_xlabel("Wavelength (Å)")
        axes[-1, ci].xaxis.set_major_locator(MaxNLocator(4))
    fig.suptitle(label, fontsize=13)
    fig.legend(handles=[Line2D([0], [0], color="black", lw=2, label="Whole-sample mean LSF"),
                        Line2D([0], [0], color=plt.get_cmap("viridis")(.6), lw=1, label="Epoch–lunation mean LSF (color = date)")],
               loc="upper center", bbox_to_anchor=(.47, .952), ncol=2, frameon=False)
    fig.subplots_adjust(left=.08, right=.86, bottom=.07, top=.88, wspace=.52, hspace=.2)
    colorbar_ax = fig.add_axes([.90, .16, .012, .62])
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="viridis"), cax=colorbar_ax,
                 label="Reference MJD (UTC)")
    colorbar_ax.tick_params(labelsize=9)
    curve_name = f"fiber_{fiber}_wavelength.png"
    fig.savefig(views / curve_name, dpi=125)
    plt.close(fig)
    print(f"Rendered {label}", flush=True)
    return surface_name, curve_name, pd.concat(records, ignore_index=True)


def render_metric_maps(catalog, destination, figures):
    """All available physical fibers; rows are categorical IDs, with no interpolation."""
    from analysis.lsf_epoch_coverage.analyze_stability import diagnostics
    survey = catalog[catalog.scope == "survey"]
    images = []
    with fits.open(destination / survey.file.iloc[0], memmap=True) as hdus:
        for (telescope, sp), rows in survey.groupby(["telescope", "spectrographid"]):
            fig, axes = plt.subplots(5, 3, figsize=(15, 14), sharex="col")
            for ci, channel in enumerate(CHANNELS):
                selected = rows[rows.channel == channel].sort_values("fiberid")
                ids = selected.reference_id.to_numpy(dtype=int)
                cf = hdus["LSF_CONFIG"].data[ci]
                knots = np.asarray(hdus["LSF_KNOTS"].data[ids[0]], dtype=float)
                wave = np.linspace(knots[int(cf["degree"])], knots[-int(cf["degree"])-1], 161)
                basis = evaluate_bspline_basis(wave, knots, int(cf["degree"]))
                coef = np.asarray(hdus["LSF_COEF"].data[ids], dtype=float)
                masses = np.einsum("wb,ftb->fwt", basis, coef)
                np.testing.assert_allclose(masses.sum(axis=2), 1, atol=2e-12)
                flat = masses.reshape(-1, 11)
                values = {m: [] for m in METRICS}
                for start in range(0, len(flat), 512):
                    current = diagnostics(flat[start:start+512], float(cf["offset_half_width_angstrom"]))
                    for metric in METRICS:
                        values[metric].append(current[metric])
                values = {m: np.concatenate(v).reshape(len(ids), len(wave)) for m, v in values.items()}
                # Verify the batched raster against the existing single-state evaluator.
                for index in sorted(set([0, len(ids)//2, len(ids)-1])):
                    state = load_reference_lsf(destination / survey.file.iloc[0], int(selected.fiberid.iloc[index]), allow_epoch_pooling=True)
                    reference = evaluate_lsf_diagnostics(state, wave[[0, 80, 160]])
                    for metric in METRICS:
                        np.testing.assert_allclose(values[metric][index, [0, 80, 160]], reference[metric], atol=1e-9, rtol=1e-9)
                positions = np.arange(len(ids))
                ticks = np.unique(np.linspace(0, len(ids)-1, min(6, len(ids))).round().astype(int))
                for mi, metric in enumerate(METRICS):
                    ax = axes[mi, ci]
                    image = ax.pcolormesh(wave, positions, values[metric], shading="nearest", cmap="viridis")
                    ax.set_ylabel("Physical fiber ID")
                    ax.set_yticks(ticks, selected.fiberid.to_numpy()[ticks])
                    fig.colorbar(image, ax=ax, pad=.02, label=LABELS[mi])
                axes[0, ci].set_title(channel)
                axes[-1, ci].set_xlabel("Wavelength (Å)")
                axes[-1, ci].xaxis.set_major_locator(MaxNLocator(4))
            fig.suptitle(f"{telescope} / SP{sp}: whole-sample mean LSF metrics\nEach row is one measured physical fiber; rows sorted by ID; no interpolation between fibers", fontsize=13)
            fig.subplots_adjust(left=.065, right=.97, bottom=.065, top=.92, wspace=.48, hspace=.24)
            name = f"metrics_{telescope}_sp{sp}.png"
            fig.savefig(figures / name, dpi=140)
            plt.close(fig)
            images.append((f"{telescope} / SP{sp}", name))
            print(f"Mapped all fibers: {telescope} / SP{sp}", flush=True)
    return images


def publish_reference_views(catalog, manifest, selection, destination):
    destination = Path(destination)
    views = destination / "figures"
    views.mkdir(exist_ok=True)
    records = destination.parent / "lsf_reference_audit" / "publication_records"
    records.mkdir(parents=True, exist_ok=True)
    loader = Path(__file__).resolve().parents[2] / "skysub/sky_decomp/lsf_reference.py"
    portable = loader.read_text().split("\ndef reference_covariance", 1)[0]
    portable = portable.replace("from .lsf_", "from skysub.sky_decomp.lsf_").replace("from .result_io", "from skysub.sky_decomp.result_io")
    (destination / "load_lsf.py").write_text(portable)
    selected = select_fibers(catalog, destination, records)
    sections, tables = [], []
    for (telescope, sp), fibers in selected.groupby(["telescope", "spectrographid"]):
        sections.append(f'<h2>{html.escape(telescope)} / SP{sp}</h2>')
        for row in fibers.itertuples(index=False):
            surface, curve, values = render_fiber(row, catalog, destination, views)
            tables.append(values)
            fiber_rows = catalog[catalog.fiberid == row.fiberid]
            coverage = []
            for channel in CHANNELS:
                metadata = fiber_rows[(fiber_rows.scope == "survey") & (fiber_rows.channel == channel)].iloc[0]
                periods = len(fiber_rows[(fiber_rows.scope == "lunation") & (fiber_rows.channel == channel)])
                coverage.append(f"{channel}: {int(metadata.n_exposures)} exposures, {int(metadata.n_nights)} nights, {periods} epoch–lunation curves")
            sections.append(f'<h3>Fiber {row.fiberid} — {row.slit_position} slit position</h3>'
                            f'<p>{"; ".join(coverage)}.</p>'
                            f'<img src="figures/{surface}?v=maps" alt="2D LSF map for fiber {row.fiberid}">'
                            f'<img src="figures/{curve}?v=maps" alt="Wavelength-dependent LSF diagnostics for fiber {row.fiberid}">')
    maps = render_metric_maps(catalog, destination, views)
    pd.concat(tables, ignore_index=True).to_parquet(records / "selected_wavelength_diagnostics.parquet", index=False)
    summary = catalog.groupby("scope").agg(surfaces=("reference_id", "size"), files=("file", "nunique")).reset_index()
    catalog.to_parquet(records / "reference_catalog.parquet", index=False)
    index_columns = ["scope", "fiberid", "telescope", "spectrographid", "channel", "calibration_epoch",
                     "full_moon_utc", "reference_id", "n_exposures", "n_nights", "file"]
    catalog[index_columns].to_csv(records / "reference_catalog.csv", index=False)
    summary.to_csv(records / "summary.csv", index=False)
    stats = dict(reference_surfaces=len(catalog), fits_files=int(catalog.file.nunique()), selected_fibers=len(selected),
                 rendered_surfaces=len(selected)*3, wavelength_figures=len(selected),
                 input_channel_measurements=len(selection), exposures=int(selection.expnum.nunique()),
                 fibers=int(selection.fiberid.nunique()), nights=int(selection.mjd.nunique()))
    old_summary = records / "summary.json"
    previous = json.loads(old_summary.read_text()) if old_summary.exists() else {}
    previous.update(stats)
    old_summary.write_text(json.dumps(previous, indent=2)+"\n")
    method = Path(__file__).with_name("LSF_REFERENCES_METHODS.md").read_text()
    counts = "\n".join(f"| {r.scope} | {r.surfaces:,} | {r.files} |" for r in summary.itertuples(index=False))
    readme = method.replace("{{SUMMARY_TABLE}}", "| Scope | Surfaces | FITS files |\n|---|---:|---:|\n"+counts)
    metric_rows = []
    survey = catalog[catalog.scope == "survey"]
    for channel, anchor in zip(CHANNELS, (4700, 6600, 8500)):
        rows = survey[survey.channel == channel]
        cells = []
        for metric in METRICS[:4]:
            lo, median, hi = np.percentile(rows[metric+"_reference"], [16, 50, 84])
            cells.append(f"{median:.3f} [{lo:.3f}, {hi:.3f}]")
        metric_rows.append(f"| {channel} | {anchor} | " + " | ".join(cells) + " |")
    readme = readme.replace("{{CHANNEL_METRICS}}", "| Channel | Wavelength (Å) | FWHM (Å) | W50 (Å) | W90 (Å) | Centroid (Å) |\n|---|---:|---:|---:|---:|---:|\n"+"\n".join(metric_rows))
    (destination / "README.md").write_text(readme)
    for scope in ("survey", "epoch", "lunation"):
        (destination / scope / "README.md").write_text(
            f"# {scope.capitalize()} LSF references\n\nCompact, self-describing FITS for physical fibers and B/R/Z channels. "
            "[Method and loading example](../README.md).\n")
    for name, rows in catalog.groupby("file"):
        file = destination / name
        epochs = sorted(rows.calibration_epoch.unique().tolist())
        moons = sorted(rows.full_moon_utc.unique().tolist())
        file.with_suffix(".README.md").write_text(
            f"# {file.name}\n\nScope: {rows.scope.iloc[0]}. Surfaces: {len(rows):,}. Epoch IDs: {epochs}. "
            f"Lunation labels: {moons}.\n\n"
            "Each row is one physical fiber and channel. Average = sum of scalar-weighted coefficient matrices; "
            "weights use coefficient MAD scales, robust RMS distances, and equal night mass before downweighting. "
            "Gaia aperture G-band flux / current FLUX+SKY synthetic G-band flux <0.1; the <0.01 flag is retained. "
            "Only the parent technically accepted fits contribute.\n\n"
            "[Full equations, criteria, FITS schema and loading](../README.md). "
            "The AUDIT header identifies the full calculation archive; loading LSF coefficients needs only this compact file.\n")
    prefix, files = readme.split("## Files and use", 1)
    markdown = mistune.create_markdown(plugins=["table", "math"])
    map_sections = "<h2>Mean-solution maps: wavelength × physical fiber</h2>" + "".join(
        f'<h3>{html.escape(label)}</h3><img src="figures/{name}?v=legend" alt="Mean LSF metrics versus wavelength and fiber ID">'
        for label, name in maps)
    report = r'''<!doctype html><html lang="en"><meta charset="utf-8"><title>LSF references: method, maps and coefficients</title>
<style>body{max-width:1450px;margin:24px auto;padding:0 20px;font:16px/1.5 sans-serif;color:#17202a}img{width:100%;height:auto}h2{border-top:1px solid #aaa;padding-top:20px}table{border-collapse:collapse}td,th{padding:5px 12px;border:1px solid #ccc}pre{background:#f4f5f6;padding:14px;overflow:auto}.math{overflow-x:auto}</style>
<script>window.MathJax={tex:{inlineMath:[['\\(','\\)']]},svg:{fontCache:'local'}};</script>
<script defer src="assets/mathjax-tex-svg-full.js"></script>
'''+markdown(prefix)+map_sections+'<h2>Selected physical fibers: density maps and wavelength curves</h2>'+"".join(sections)+markdown("## Files and use"+files)+'''<p><a href="README.md">README</a> · <a href="survey/">Survey coefficient FITS</a> · <a href="epoch/">Epoch coefficient FITS</a> · <a href="lunation/">Epoch–lunation coefficient FITS</a></p></html>'''
    (destination / "lsf_reference_surfaces.html").write_text(report)
    return stats
