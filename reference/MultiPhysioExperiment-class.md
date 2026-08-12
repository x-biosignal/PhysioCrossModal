# MultiPhysioExperiment class definition

The `MultiPhysioExperiment` class holds multiple `PhysioExperiment`
objects representing different signal modalities (e.g., EEG, EMG, ECG)
recorded simultaneously but potentially at different sampling rates. It
provides temporal alignment metadata and a cache for coupling analysis
results.

## Slots

- `experiments`:

  Named list of `PhysioExperiment` objects.

- `alignment`:

  [`DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  with temporal alignment metadata (one row per modality).

- `sampleMap`:

  [`DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  mapping samples across modalities.

- `couplingResults`:

  List of cached coupling analysis results.

## References

Huber, W., et al. (2015). "Orchestrating high-throughput genomic
analysis with Bioconductor." Nature Methods, 12(2), 115-121.
doi:10.1038/nmeth.3252

## See also

[`MultiPhysioExperiment()`](https://x-biosignal.github.io/PhysioCrossModal/reference/MultiPhysioExperiment.md),
[`experiments()`](https://x-biosignal.github.io/PhysioCrossModal/reference/experiments.md),
[`modalities()`](https://x-biosignal.github.io/PhysioCrossModal/reference/modalities.md),
[`alignSignals()`](https://x-biosignal.github.io/PhysioCrossModal/reference/alignSignals.md)
