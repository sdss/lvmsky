# Moon/Zodi predictor reference fixture

`moon_zodi_predictor_reference_v1.npz` freezes three independent detector-grid
cases from the source JAX implementation: Moon below the horizon (stack row 0),
Moon high above the horizon (row 62), and strong zodiacal light (row 14469).
The reference vectors were generated with the frozen 30,000-step global
parameters and the exact high-resolution physical carriers before the
production predictor was used by the decomposition class.

The fixture SHA-256 is
`53c63d5cc9ef32a38a24a05d12a1a6fa00b02b43365d359b5274315dd16d0c4b`.
It contains the native 12,401-point wavelength grid, measured detector LSFs,
observation metadata, derived geometry, and separate Moon and Zodi vectors.
Production tests read this fixture and do not import JAX or experiment code.
