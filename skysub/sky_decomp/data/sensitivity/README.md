# Vendored absolute sensitivity curves

`sens_abs-{b,r,z}.csv` — ABSOLUTE spectrograph sensitivity in
[erg/s/cm^2/A] per [e-/s/A], on the native LVM wavelength grid per arm.

Vendored so the decomposition's per-pixel weighting has no dependency on
`$LVMCORE_DIR` being present or correct. Column 4 (0-indexed, of 8) of
`lvmcore/sensitivity/sens_percentiles-{arm}.csv`; the other seven columns are
other percentiles of the same quantity and are not used.

WHY THE PERCENTILE TABLE AND NOT `mean-sens-{arm}-v1.1.csv`: the latter cannot
carry an absolute scale by construction — it divides each arm by its own
weighted mean — so it fixes only the SHAPE within an arm and leaves both the
overall normalisation and the arm-to-arm ratios undetermined. Guessing those
ratios was wrong by 1.6x / 2.25x. See `mlp_predictor/noise.py` for the history.

Source sha256 (full 8-column files):
  sens_percentiles-b.csv  b45ea7a6a4073aff4b8c24f7ea8546d4b22a838a88687d4fe95dbf74b58e2239
  sens_percentiles-r.csv  1533afa02bdfba43659f31d0c22bb3199e3c3b0156c733e9a3e61a5ef5a50575
  sens_percentiles-z.csv  0fc356eb396dbad22574897da299eb5d79956934919adee3a1287732d4f07cbd
