# Construct a MultiPhysioExperiment object

Creates a container that holds multiple `PhysioExperiment` objects
recorded simultaneously at potentially different sampling rates.

## Usage

``` r
MultiPhysioExperiment(..., experiments = list(), alignment = NULL)
```

## Arguments

- ...:

  Named `PhysioExperiment` objects, one per modality.

- experiments:

  Alternatively, a named list of `PhysioExperiment` objects. If both
  `...` and `experiments` are provided, they are combined.

- alignment:

  Optional
  [`DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  with temporal alignment metadata. When `NULL` (the default), a default
  alignment table is built from the supplied experiments.

## Value

A
[`MultiPhysioExperiment-class`](https://x-biosignal.github.io/PhysioCrossModal/reference/MultiPhysioExperiment-class.md)
instance containing the supplied experiments, alignment metadata, and an
empty coupling results cache.

## References

Huber, W., et al. (2015). "Orchestrating high-throughput genomic
analysis with Bioconductor." Nature Methods, 12(2), 115-121.
doi:10.1038/nmeth.3252

## See also

[`experiments()`](https://x-biosignal.github.io/PhysioCrossModal/reference/experiments.md),
[`modalities()`](https://x-biosignal.github.io/PhysioCrossModal/reference/modalities.md),
[`samplingRates()`](https://x-biosignal.github.io/PhysioCrossModal/reference/samplingRates.md),
[`alignSignals()`](https://x-biosignal.github.io/PhysioCrossModal/reference/alignSignals.md),
[`couplingAnalysis()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysis.md)

## Examples

``` r
# Create two PhysioExperiment objects with different sampling rates
eeg_data <- matrix(rnorm(500 * 4), nrow = 500, ncol = 4)
emg_data <- matrix(rnorm(1000 * 2), nrow = 1000, ncol = 2)

pe_eeg <- PhysioExperiment(
  assays = list(raw = eeg_data),
  colData = S4Vectors::DataFrame(
    label = c("Fz", "Cz", "Pz", "Oz"),
    type = rep("EEG", 4)
  ),
  samplingRate = 250
)

pe_emg <- PhysioExperiment(
  assays = list(raw = emg_data),
  colData = S4Vectors::DataFrame(
    label = c("EMG1", "EMG2"),
    type = rep("EMG", 2)
  ),
  samplingRate = 1000
)

# Construct MultiPhysioExperiment
mpe <- MultiPhysioExperiment(EEG = pe_eeg, EMG = pe_emg)
mpe
#> class: MultiPhysioExperiment
#> modalities(2): EEG, EMG
#> samplingRates: EEG=250Hz, EMG=1000Hz
#>   EEG: 500 timepoints x 4 channels
#>   EMG: 1000 timepoints x 2 channels
```
