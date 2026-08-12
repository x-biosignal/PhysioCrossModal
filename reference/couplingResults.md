# Cached coupling results on a MultiPhysioExperiment

Get or set the list of cached coupling-analysis results stored on a
[`MultiPhysioExperiment`](https://x-biosignal.github.io/PhysioCrossModal/reference/MultiPhysioExperiment.md).
The cache is populated by
[`couplingAnalysisCached()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysisCached.md),
which keys entries by a digest of the analysis arguments so repeated
identical requests are not recomputed.

## Usage

``` r
couplingResults(x)

# S4 method for class 'MultiPhysioExperiment'
couplingResults(x)

couplingResults(x) <- value

# S4 method for class 'MultiPhysioExperiment'
couplingResults(x) <- value
```

## Arguments

- x:

  A `MultiPhysioExperiment` object.

- value:

  A named list of cached results.

## Value

`couplingResults()` returns the cache list; the replacement form returns
the updated object.

## See also

[`couplingAnalysisCached()`](https://x-biosignal.github.io/PhysioCrossModal/reference/couplingAnalysisCached.md)
