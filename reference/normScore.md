# Rank Normalization Methods Using a Composite Scoring Framework

Computes a composite score to rank normalization methods based on
multiple quality metrics derived from normalized data matrices.

The function evaluates each normalization method using six criteria:
pooled coefficient of variation (PCV), within-group sample correlation,
MA-plot trend deviation, mean-SD trend deviation, relative log
expression (RLE) consistency, and total intensity consistency. These
metrics are scaled, combined into a final score, and used to rank the
normalization methods.

In addition, a correction factor derived from the coefficient of
variation of raw sample intensities is applied to the `"Log"`
normalization method. Optionally, bootstrap resampling is used to
estimate mean composite scores and percentile-based confidence
intervals.

## Usage

``` r
normScore(
  normalizedDataList,
  groupData,
  rawData,
  refGroup = NULL,
  altGroup = NULL,
  returnDetails = TRUE,
  nBoot = 1000
)
```

## Arguments

- normalizedDataList:

  Named list of normalized numeric matrices or data frames. Each element
  represents one normalization method, with rows corresponding to
  features (e.g. proteins) and columns to samples.

- groupData:

  `data.frame` or matrix describing the sample-group assignment. It must
  contain two columns, which are interpreted as sample names and group
  labels.

- rawData:

  `data.frame` or matrix of raw numeric intensities, with rows
  representing features and columns representing samples.

- refGroup:

  `character` or `NULL`. Reference group used for the MA-plot metric. If
  `NULL`, the last group in `groupData` is used.

- altGroup:

  `character` or `NULL`. Alternative group used for the metric. If
  `NULL`, the first group in `groupData` is used.

- returnDetails:

  `logical`. If `TRUE`, only the final ranking is MA-plot returned. If
  `FALSE`, detailed scores, bootstrap results, and a summaryplot are
  also returned. Default is `FALSE`.

- nBoot:

  `integer`. Number of bootstrap resamples used to estimate mean
  normScore values and confidence intervals. Default is `1000`.

## Value

If `returnDetails = TRUE`, returns a list with:

- finalRanking:

  Named numeric vector with the final normalization ranking.

If `returnDetails = FALSE`, returns a list with:

- finalRanking:

  Named numeric vector with the final normalization ranking.

- detailRanking:

  `data.frame` with scaled item-wise scores, total score, and corrected
  total score for each normalization method.

- bootstrapScore:

  `data.frame` with bootstrap mean scores and 95\\ percentile confidence
  intervals for each normalization method.

- graphic:

  A `ggplot2` object showing the bootstrap mean scores and confidence
  intervals.

## Details

The function computes the following six item scores for each
normalization method:

1.  **Item1**: mean pooled coefficient of variation across groups (PVC).

2.  **Item2**: within-group Spearman correlation summary, transformed so
    that lower values indicate better performance (Correlation).

3.  **Item3**: shape-corrected area-based deviation from the expected MA
    trend (`logFC = 0`), based on MAplot.

4.  **Item4**: area-based deviation of the mean-SD trend from
    horizontality, based on meanSD plot.

5.  **Item5**: RLE-based mean absolute percentage error relative to 1.

6.  **Item6**: quantile-based total intensity consistency metric.

Each item is scaled to the range 0 to 1 using min-max scaling and then
summed into a total score. Lower scores indicate better normalization
performance.

A correction factor based on the coefficient of variation of total raw
sample intensities is applied only to the normalization method named
`"Log"`.

When `returnDetails = FALSE`, bootstrap resampling of the six item
scores is performed.The results can be visualized using
[`plotBootstrapNormScore`](https://juliagcurras.github.io/normScore/reference/plotBootstrapNormScore.md).

## See also

[`plotNormScoreDiagnostics`](https://juliagcurras.github.io/normScore/reference/plotNormScoreDiagnostics.md),
[`plotBootstrapNormScore`](https://juliagcurras.github.io/normScore/reference/plotBootstrapNormScore.md),
[`simulateData`](https://juliagcurras.github.io/normScore/reference/simulateData.md)

## Examples

``` r

# Simulate proteomic data
simData <- simulateData(nProteins = 1000)

# Normalyze ysing NormalyzerDE package
normalizedDataList <- list(
    Norm1 = simData$logData + 0.1,
    Norm2 = simData$logData + 1,
    Norm3 = simData$logData - 1,
    Norm4 = simData$logData * 1.1,
    Norm5 = simData$logData * 0.9,
    Norm6 = simData$logData * runif(ncol(simData$logData), 0.8, 1.2)
)

# Compute ranking
result <- normScore(
    normalizedDataList,
    groupData = simData$metadata,
    rawData = simData$rawData,
    refGroup = NULL,
    altGroup = NULL,
    returnDetails = TRUE,
    nBoot = 100
)
#> Reference group set to G2.
#> Alternative group set to G1.
#> Log2-transformed raw data were added to
#>             'normalizedDataList' as 'Log'.
result
#> $finalRanking
#>        Log      Norm5      Norm2      Norm1      Norm6      Norm3      Norm4 
#> 0.07079693 1.28431816 1.91416771 2.41162482 2.63643762 3.09313202 3.65587983 
#> 
#> $detailRanking
#>           Item1 Item2     Item3     Item4     Item5      Item6    Total
#> Log   0.4713758   0.1 0.5000000 0.8142192 0.4999485 0.08450405 2.470048
#> Norm5 0.4713758   0.1 0.0000000 0.6284383 0.0000000 0.08450405 1.284318
#> Norm2 0.0000000   0.1 0.5000000 0.8142192 0.4999485 0.00000000 1.914168
#> Norm1 0.4218242   0.1 0.5000000 0.8142192 0.4999485 0.07563292 2.411625
#> Norm6 0.4713758   0.0 0.8346731 0.0000000 0.3303888 1.00000000 2.636438
#> Norm3 1.0000000   0.1 0.5000000 0.8142192 0.4999485 0.17896431 3.093132
#> Norm4 0.4713758   0.1 1.0000000 1.0000000 1.0000000 0.08450405 3.655880
#>       TotalxItem0
#> Log    0.07079693
#> Norm5  1.28431816
#> Norm2  1.91416771
#> Norm1  2.41162482
#> Norm6  2.63643762
#> Norm3  3.09313202
#> Norm4  3.65587983
#> 
#> $bootstrapScore
#>   normalization meanNormScore      ll95       ul95
#> 1           Log     0.1393238 0.0545619  0.2674648
#> 2         Norm5     3.7364414 1.1994417  7.1184431
#> 3         Norm2     5.4469503 1.3080386 10.7193999
#> 4         Norm1     6.8978096 2.2940229 13.1223678
#> 5         Norm6     7.3844707 1.7794351 14.3342643
#> 6         Norm3     8.8855229 3.1746362 16.6470613
#> 7         Norm4    10.4002546 3.3071218 19.7737763
#> 
```
