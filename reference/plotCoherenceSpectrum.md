# Plot coherence spectrum

Displays the coherence as a function of frequency from the result of
[`coherence`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md).
Optionally draws the 95\\ a horizontal dashed line.

## Usage

``` r
plotCoherenceSpectrum(
  result,
  show_threshold = TRUE,
  title = "Coherence Spectrum",
  ...
)
```

## Arguments

- result:

  A list returned by
  [`coherence`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md),
  containing at least `$coherence` and `$frequencies`. If
  `$confidence_limit` is present it is drawn as a threshold line.

- show_threshold:

  Logical; if TRUE (default) and a confidence limit is available, draw
  it on the plot.

- title:

  Character string for the plot title (default `"Coherence Spectrum"`).

- ...:

  Additional arguments (currently unused).

## Value

A `ggplot` object.

## References

Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*.
Springer-Verlag New York.
[doi:10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4)

## See also

[`coherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/coherence.md),
[`multitaperCoherence()`](https://x-biosignal.github.io/PhysioCrossModal/reference/multitaperCoherence.md),
[`plotCouplingMatrix()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotCouplingMatrix.md)

## Examples

``` r
sr <- 500
t <- seq(0, 10, length.out = sr * 10)
x <- sin(2 * pi * 20 * t) + 0.2 * rnorm(length(t))
y <- 0.8 * sin(2 * pi * 20 * t) + 0.2 * rnorm(length(t))
res <- coherence(x, y, sr = sr)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plotCoherenceSpectrum(res)
}
```
