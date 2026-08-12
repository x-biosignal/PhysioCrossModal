# Apply a computed warping to a companion signal

Carries a warping computed by
[`elasticAlign()`](https://x-biosignal.github.io/PhysioCrossModal/reference/elasticAlign.md)
onto another channel or modality (e.g. warp an EMG envelope by the
diffeomorphism that aligns the EEG).

## Usage

``` r
warpApply(signal, gamma, t = NULL)
```

## Arguments

- signal:

  Numeric signal to warp (same length as the warping grid).

- gamma:

  A warping from
  [`elasticAlign()`](https://x-biosignal.github.io/PhysioCrossModal/reference/elasticAlign.md).

- t:

  Optional grid (defaults to 0..1).

## Value

The warped signal.

## See also

[`elasticAlign()`](https://x-biosignal.github.io/PhysioCrossModal/reference/elasticAlign.md)

## Examples

``` r
t <- seq(0, 1, length.out = 50)
warpApply(sin(2 * pi * t), gamma = t^1.5, t = t)[1:3]
#> [1] 0.00000000 0.01826817 0.05167018
```
