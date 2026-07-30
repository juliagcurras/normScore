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
#>        Log      Norm5      Norm6      Norm2      Norm1      Norm3      Norm4 
#> 0.06470467 1.25135716 1.89862636 1.90889209 2.38617798 3.04014697 3.63334351 
#> 
#> $detailRanking
#>           Item1 Item2     Item3     Item4     Item5      Item6    Total
#> Log   0.4713738   0.1 0.5000000 0.8090068 0.4998853 0.06196976 2.442236
#> Norm5 0.4713738   0.1 0.0000000 0.6180136 0.0000000 0.06196976 1.251357
#> Norm6 0.4713738   0.0 0.3582391 0.0000000 0.0690135 1.00000000 1.898626
#> Norm2 0.0000000   0.1 0.5000000 0.8090068 0.4998853 0.00000000 1.908892
#> Norm1 0.4218222   0.1 0.5000000 0.8090068 0.4998853 0.05546369 2.386178
#> Norm3 1.0000000   0.1 0.5000000 0.8090068 0.4998853 0.13125488 3.040147
#> Norm4 0.4713738   0.1 1.0000000 1.0000000 1.0000000 0.06196976 3.633344
#>       TotalxItem0
#> Log    0.06470467
#> Norm5  1.25135716
#> Norm6  1.89862636
#> Norm2  1.90889209
#> Norm1  2.38617798
#> Norm3  3.04014697
#> Norm4  3.63334351
#> 
#> $bootstrapScore
#>   normalization meanNormScore       ll95       ul95
#> 1           Log     0.1272803 0.04953453  0.2443309
#> 2         Norm5     3.6447895 1.19943768  6.9190588
#> 3         Norm6     5.3190620 1.00082168 10.5092757
#> 4         Norm2     5.4317108 1.30505629 10.6884795
#> 5         Norm1     6.8275021 2.26231887 12.9677234
#> 6         Norm3     8.7400366 3.11009102 16.3235006
#> 7         Norm4    10.3387300 3.25530477 19.6738705
#> 
```
