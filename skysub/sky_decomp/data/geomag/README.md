# Geomagnetic grid for ionospheric pierce points

`igrf_285km_lco.npz` holds the IGRF declination and inclination (degrees) at 285 km altitude. The grid covers
lat -45..-13 deg and lon -90..-50 deg in 0.5 deg steps, for epochs 2023.0 to 2027.0 (one per year).
`skysub/mlp_predictor/ionosphere.py` interpolates it to get the magnetic declination and dip latitude
(`atan(tan I / 2)`) at each telescope's F-region pierce point. The runtime needs only numpy and scipy.

Rebuild (dev only; `ppigrf` is the `geomag` extra in `pyproject.toml`), from the `lvmsky/` repo root:

```bash
uv pip install --python .venv/bin/python -e ".[geomag]"
cd skysub/sky_decomp/data/geomag && ../../../../.venv/bin/python make_igrf_grid.py
```

Do not run `uv sync` on an existing `.venv`. The lock describes only the project extras, so a sync would
uninstall the other packages in the environment.

Values at LCO (lat -29.0, lon -70.7, 285 km): declination -1.8 deg (2023) to -2.5 deg (2027), inclination
about -30.4 deg, dip latitude about -16 deg.
