# Plot a comodulogram

Renders a comodulogram (from
[`comodulogram()`](https://x-biosignal.github.io/PhysioCrossModal/reference/comodulogram.md))
as a heatmap.

## Usage

``` r
plotComodulogram(
  result,
  title = "Comodulogram",
  mask = NULL,
  nonsignificant_alpha = 0.2
)
```

## Arguments

- result:

  A comodulogram result list.

- title:

  Plot title.

- mask:

  Optional logical significance mask or numeric p-value matrix. When
  `NULL`, uses `result$significant` if present.

- nonsignificant_alpha:

  Opacity for non-significant cells.

## Value

A `ggplot` object.

## See also

[`comodulogram()`](https://x-biosignal.github.io/PhysioCrossModal/reference/comodulogram.md),
[`modulationIndex()`](https://x-biosignal.github.io/PhysioCrossModal/reference/modulationIndex.md)

## Examples

``` r
if (FALSE) { # \dontrun{
plotComodulogram(comodulogram(sig, sr = 500,
  phase_freqs = c(4, 6, 8), amp_freqs = c(30, 50, 70)))
} # }
```
