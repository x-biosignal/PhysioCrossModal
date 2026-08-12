# Plot time-varying coupling from sliding-window analysis

Displays the peak cross-correlation over time from the result of
[`slidingCrossCorrelation`](https://x-biosignal.github.io/PhysioCrossModal/reference/slidingCrossCorrelation.md).

## Usage

``` r
plotCouplingTimecourse(result, title = "Coupling Time Course", ...)
```

## Arguments

- result:

  A list returned by
  [`slidingCrossCorrelation`](https://x-biosignal.github.io/PhysioCrossModal/reference/slidingCrossCorrelation.md),
  containing at least `$times` and `$peak_correlations`.

- title:

  Character string for the plot title (default
  `"Coupling Time Course"`).

- ...:

  Additional arguments (currently unused).

## Value

A `ggplot` object.

## References

Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*.
Springer-Verlag New York.
[doi:10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4)

## See also

[`slidingCrossCorrelation()`](https://x-biosignal.github.io/PhysioCrossModal/reference/slidingCrossCorrelation.md),
[`crossCorrelation()`](https://x-biosignal.github.io/PhysioCrossModal/reference/crossCorrelation.md),
[`plotCoherenceSpectrum()`](https://x-biosignal.github.io/PhysioCrossModal/reference/plotCoherenceSpectrum.md)

## Examples

``` r
sr <- 500
set.seed(1)
x <- rnorm(5000)
y <- c(rep(0, 10), x[1:4990])
res <- slidingCrossCorrelation(x, y, sr = sr,
                                window_sec = 1, step_sec = 0.5)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plotCouplingTimecourse(res)
}
```
