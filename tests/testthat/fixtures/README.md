# Numeric reference fixtures

## `nonparametric_granger_fieldtrip_reference.rds`

Offline WSCB-11 reference for bivariate nonparametric Granger causality.

The raw stable AR signals and multitaper cross-spectral matrix are generated
independently with NumPy/SciPy. The matrix is then passed to FieldTrip's
`sfactorization_wilson()` and `ft_connectivity_granger()` from a pinned
FieldTrip commit under GNU Octave. The committed RDS contains numeric
inputs/outputs and provenance only; no FieldTrip implementation is copied.

The fixture records:

- raw `x`/`y`, sampling rate, segment/FFT/overlap/DPSS settings;
- SciPy's one-sided 2x2 cross-spectral matrix;
- FieldTrip directional Granger spectra, transfer function, innovation
  covariance, and reconstructed spectrum; and
- FieldTrip commit, official source URLs, GNU Octave, NumPy, and SciPy
  versions, generation date, and exact commands.

Regeneration uses:

```bash
python3 generate_nonparametric_granger_fieldtrip_reference.py input \
  /tmp/wscb11-input.mat
octave --quiet generate_nonparametric_granger_fieldtrip_reference.m \
  /tmp/wscb11-input.mat /tmp/wscb11-output.mat /tmp/wscb11-fieldtrip
python3 generate_nonparametric_granger_fieldtrip_reference.py collect \
  /tmp/wscb11-input.mat /tmp/wscb11-output.mat /tmp/wscb11-reference.json
Rscript -e 'x <- jsonlite::read_json("/tmp/wscb11-reference.json", simplifyVector=TRUE); saveRDS(x, "nonparametric_granger_fieldtrip_reference.rds", version=3)'
```

The temporary FieldTrip directory must contain the two pinned official files
and minimal `ft_*` runtime shims. See the commands stored inside the fixture.
