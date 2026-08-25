"""Patch notebook_sky_interpolation_triplet_dual_encoder_group_mlp_split_zodi.ipynb

Turn the hard-coded moon_zodi_spline2/lvmsframe_median_stack_1.2.1_p40_p70 path
prefix into two module-level variables at the top of the notebook so switching
input corpora is a one-line change:

    DECOMP_DATA_ROOT = 'moon_zodi_spline2'
    DECOMP_STEM      = 'lvmsframe_median_stack_1.2.1_p40_p70'

All 15 downstream references are rewritten to f-string derivations.
"""

from __future__ import annotations

import json
from pathlib import Path

REPO = Path("/Users/droryn/prog/lvm/lvmsky")
NB = REPO / "skysub" / (
    "notebook_sky_interpolation_triplet_dual_encoder_group_mlp_split_zodi.ipynb"
)

OLD_ROOT = "moon_zodi_spline2"
OLD_STEM = "lvmsframe_median_stack_1.2.1_p40_p70"
OLD_PREFIX = f"{OLD_ROOT}/{OLD_STEM}"

# (needle-in-cell-source, replacement) pairs. Each needle is anchored to a
# single line so it works even after the .ipynb JSON has been split into a
# list of strings. Line-terminating "\n" is intentionally kept on both sides.
LINE_REPLACEMENTS = [
    (
        f"_DECOMP_STEM = '{OLD_PREFIX}'\n",
        "_DECOMP_STEM = f'{DECOMP_DATA_ROOT}/{DECOMP_STEM}'\n",
    ),
    (
        f"                          '{OLD_PREFIX}')\n",
        "                          f'{DECOMP_DATA_ROOT}/{DECOMP_STEM}')\n",
    ),
    (
        f"INPUT_FITS_FOR_BASIS = '{OLD_PREFIX}_every10.fits'\n",
        "INPUT_FITS_FOR_BASIS = f'{DECOMP_DATA_ROOT}/{DECOMP_STEM}_every10.fits'\n",
    ),
    (
        f"WAVELENGTH_CACHE = Path('{OLD_ROOT}/coef_wavelengths_basis_v4.npz')"
        "  # bumped to v4 for split_zodi decomposition (new basis: Zodi_bs"
        " coefficients + moon_albedo carrier) cache based on p70 improved"
        " continuumnuum basis\n",
        "WAVELENGTH_CACHE = Path(f'{DECOMP_DATA_ROOT}/coef_wavelengths_basis_v4.npz')"
        "  # bumped to v4 for split_zodi decomposition (new basis: Zodi_bs"
        " coefficients + moon_albedo carrier) cache based on p70 improved"
        " continuumnuum basis\n",
    ),
    (
        f"with fits.open('{OLD_PREFIX}_meta_only.fits') as hdul:\n",
        "with fits.open(f'{DECOMP_DATA_ROOT}/{DECOMP_STEM}_meta_only.fits') as hdul:\n",
    ),
    (
        f'_e10_stem = "{OLD_PREFIX}_every10"\n',
        '_e10_stem = f"{DECOMP_DATA_ROOT}/{DECOMP_STEM}_every10"\n',
    ),
    (
        f'    _e10_stem = "{OLD_PREFIX}_every10"\n',
        '    _e10_stem = f"{DECOMP_DATA_ROOT}/{DECOMP_STEM}_every10"\n',
    ),
    (
        f'_e10_stem_ab = "{OLD_PREFIX}_every10"\n',
        '_e10_stem_ab = f"{DECOMP_DATA_ROOT}/{DECOMP_STEM}_every10"\n',
    ),
    (
        f"    INPUT_FITS = '{OLD_PREFIX}.fits'\n",
        "    INPUT_FITS = f'{DECOMP_DATA_ROOT}/{DECOMP_STEM}.fits'\n",
    ),
    (
        f"    NEAR_FITS = '{OLD_PREFIX}_decomp_sky1_lsf_surface_iterative_split_zodi.fits'\n",
        "    NEAR_FITS = f'{DECOMP_DATA_ROOT}/{DECOMP_STEM}_decomp_sky1_lsf_surface_iterative_split_zodi.fits'\n",
    ),
    (
        f"    FAR_FITS = '{OLD_PREFIX}_decomp_sky2_lsf_surface_iterative_split_zodi.fits'\n",
        "    FAR_FITS = f'{DECOMP_DATA_ROOT}/{DECOMP_STEM}_decomp_sky2_lsf_surface_iterative_split_zodi.fits'\n",
    ),
    (
        f"    SCI_FITS = '{OLD_PREFIX}_decomp_sci_lsf_surface_iterative_split_zodi.fits'\n",
        "    SCI_FITS = f'{DECOMP_DATA_ROOT}/{DECOMP_STEM}_decomp_sci_lsf_surface_iterative_split_zodi.fits'\n",
    ),
    (
        f"    META_ONLY_FITS = '{OLD_PREFIX}_meta_only.fits'\n",
        "    META_ONLY_FITS = f'{DECOMP_DATA_ROOT}/{DECOMP_STEM}_meta_only.fits'\n",
    ),
]

# Setup cell (id "44b80a04"): append the two new module-level variables.
SETUP_CELL_ID = "44b80a04"
SETUP_ANCHOR_LINE = "PALACE_DIR = '../'\n"
SETUP_INSERT_AFTER_ANCHOR = (
    "\n"
    "# Input decomposition location. Switch DECOMP_DATA_ROOT to point at a\n"
    "# different corpus (e.g. 'moon_zodi_spline3' for the n_spline_knots=11\n"
    "# reduced-basis run); all downstream FITS paths derive from these two.\n"
    "DECOMP_DATA_ROOT = 'moon_zodi_spline2'\n"
    "DECOMP_STEM = 'lvmsframe_median_stack_1.2.1_p40_p70'\n"
)


def _apply_line_replacement(source_lines, needle, replacement):
    """Replace the first occurrence of `needle` (matched as a whole line
    within `source_lines`) with `replacement`. Returns the new list and a
    boolean indicating whether the substitution fired.
    """
    for i, line in enumerate(source_lines):
        if line == needle:
            new_lines = list(source_lines)
            new_lines[i] = replacement
            return new_lines, True
    return source_lines, False


def _insert_after_anchor(source_lines, anchor, block):
    """Insert `block` (a single string, may contain newlines) immediately
    after the first line that equals `anchor` (with or without a trailing
    "\n"; Jupyter omits the newline on the last element of `source`).
    Idempotent: if the block's marker already appears anywhere in source,
    no-op.
    """
    marker = "DECOMP_DATA_ROOT ="
    if any(marker in ln for ln in source_lines):
        return source_lines, False
    anchor_stripped = anchor.rstrip("\n")
    for i, line in enumerate(source_lines):
        if line == anchor or line == anchor_stripped:
            new_lines = list(source_lines)
            # Ensure the anchor line ends with a newline so the injected
            # block lands on its own lines rather than glued to the anchor.
            if not new_lines[i].endswith("\n"):
                new_lines[i] = new_lines[i] + "\n"
            insert_lines = block.splitlines(keepends=True)
            new_lines[i + 1:i + 1] = insert_lines
            # Trim the trailing "\n" on the last inserted line to preserve
            # the Jupyter convention (last source line has no newline).
            if new_lines[-1].endswith("\n"):
                new_lines[-1] = new_lines[-1].rstrip("\n")
            return new_lines, True
    return source_lines, False


def main():
    nb = json.loads(NB.read_text())

    stats_added_config = False
    stats_line_reps = {needle: 0 for needle, _ in LINE_REPLACEMENTS}

    for cell in nb["cells"]:
        if cell.get("cell_type") != "code":
            continue
        src = cell.get("source", [])
        if isinstance(src, str):
            src = src.splitlines(keepends=True)
        # Setup cell: inject the two config variables.
        if cell.get("id") == SETUP_CELL_ID:
            src, fired = _insert_after_anchor(
                src, SETUP_ANCHOR_LINE, SETUP_INSERT_AFTER_ANCHOR)
            stats_added_config = stats_added_config or fired
        # Any cell: rewrite whichever hard-coded path lines it contains.
        for needle, replacement in LINE_REPLACEMENTS:
            src, fired = _apply_line_replacement(src, needle, replacement)
            if fired:
                stats_line_reps[needle] += 1
        cell["source"] = src

    NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")

    print(f"config vars added: {stats_added_config}")
    print("line replacements applied:")
    for needle, count in stats_line_reps.items():
        preview = needle.strip().replace("\n", " ")[:80]
        print(f"  [{count}]  {preview}")
    n_unmatched = sum(1 for _, c in stats_line_reps.items() if c == 0)
    if n_unmatched:
        print(f"WARN: {n_unmatched} needle(s) never matched. Check patterns.")


if __name__ == "__main__":
    main()
