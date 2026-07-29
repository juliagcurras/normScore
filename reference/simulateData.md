# Simulate label-free proteomics-like two-group data

Simulates two-group quantitative omics data with controlled signal
structure, including abundance-dependent variance, differential
expression, sample correlation, global shifts, mean-SD dependence, and
optional MNAR missingness.

## Usage

``` r
simulateData(
  nProteins = 10000,
  nPerGroup = 20,
  muMean = 18,
  muSd = 1.2,
  muClip = c(15, 25),
  rhoWithin = 0.85,
  rhoBetween = 0.55,
  loadingSd = 0.25,
  sigmaHi = 0.05,
  sigmaLo = 0.4,
  gammaSigma = 2.5,
  propDE = 0.35,
  logFCSd = 1,
  logFCMean = 0,
  heteroLogFC = TRUE,
  fcHi = 0.55,
  fcLo = 2.5,
  gammaFC = 7,
  sampleShiftSd = 0,
  sampleShiftCap = 0.2,
  sampleSdStrength = 0,
  sampleSdRho = 0.8,
  sampleSdCap = 0.35,
  addMissing = TRUE,
  targetMissing = 0.001,
  kMnar = 1.2,
  missingBySampleSd = 0.05,
  seed = 9396
)
```

## Arguments

- nProteins:

  `integer`. Number of proteins or features to simulate. Default is
  `10000`.

- nPerGroup:

  `integer`. Number of samples per group. Default is `20`.

- muMean:

  `numeric`. Mean of the protein-wise abundance distribution on the log2
  scale. Default is `18`.

- muSd:

  `numeric`. Standard deviation of the protein-wise abundance
  distribution on the log2 scale. Default is `1.2`.

- muClip:

  `numeric` vector of length 2. Lower and upper bounds used to truncate
  simulated protein means. Default is `c(15, 25)`.

- rhoWithin:

  `numeric`. Within-group correlation target for the latent sample
  factor. Default is `0.85`.

- rhoBetween:

  `numeric`. Between-group correlation target for the latent sample
  factor. Default is `0.55`.

- loadingSd:

  `numeric`. Standard deviation of protein-specific loadings for the
  correlated latent factor. Default is `0.25`.

- sigmaHi:

  `numeric`. Residual standard deviation at high abundance. Default is
  `0.05`.

- sigmaLo:

  `numeric`. Residual standard deviation at low abundance. Default is
  `0.4`.

- gammaSigma:

  `numeric`. Controls how strongly residual variance depends on
  abundance. Default is `2.5`.

- propDE:

  `numeric`. Proportion of differentially expressed proteins. Default is
  `0.35`.

- logFCSd:

  `numeric`. Standard deviation of differential expression log-fold
  changes. Default is `1`.

- logFCMean:

  `numeric`. Mean of differential expression log-fold changes. Default
  is `0`.

- heteroLogFC:

  `logical`. If `TRUE`, the variance of log-fold changes depends on
  abundance. Default is `TRUE`.

- fcHi:

  `numeric`. Lower multiplier for abundance-dependent logFC
  heterogeneity. Default is `0.55`.

- fcLo:

  `numeric`. Upper multiplier for abundance-dependent logFC
  heterogeneity. Default is `2.5`.

- gammaFC:

  `numeric`. Controls how strongly logFC variability depends on
  abundance. Default is `7`.

- sampleShiftSd:

  `numeric`. Standard deviation of global sample shifts. Default is `0`.

- sampleShiftCap:

  `numeric`. Maximum absolute value allowed for global sample shifts.
  Default is `0.2`.

- sampleSdStrength:

  `numeric`. Strength of sample-specific mean-SD dependence. A value of
  `0` implies independence. Default is `0`.

- sampleSdRho:

  `numeric`. Correlation between sample mean rank and sample-specific SD
  effect. Must lie between `0` and `1`. Default is `0.8`.

- sampleSdCap:

  `numeric`. Maximum absolute value allowed for the log-multiplier
  controlling sample-specific SD effects. Default is `0.35`.

- addMissing:

  `logical`. Should MNAR missing values be added? Default is `TRUE`.

- targetMissing:

  `numeric`. Target overall proportion of missing Default is `0.001`.

- kMnar:

  `numeric`. Strength of abundance dependence in the MNAR missing
  values. value mechanism. Default is `1.2`.

- missingBySampleSd:

  `numeric`. Standard deviation of sample-specific missingness shifts.
  Default is `0.05`.

- seed:

  `integer`. Random seed used for reproducibility. Default is `9396`.

## Value

A list with the following elements:

- logData:

  Numeric matrix of simulated log2-scale data.

- rawData:

  Numeric matrix of simulated raw-scale data, obtained as `2^logData`.

- metadata:

  `data.frame` containing sample names and group labels.

## Details

The simulated data include several structured components:

1.  Protein-wise mean abundance on the log2 scale.

2.  A latent correlation structure inducing stronger within-group than
    between-group sample correlation.

3.  Abundance-dependent heteroscedastic residual noise to create an
    MA-like wedge pattern.

4.  Symmetric differential expression between the two groups.

5.  Optional global sample shifts affecting overall intensity.

6.  Optional sample-level mean-SD dependence.

7.  Optional MNAR missingness driven by abundance and sample-specific
    effects.

The output is intended for testing normalization methods and associated
scoring procedures under controlled simulation settings.

## Examples

``` r
simData <- simulateData(
    nProteins = 1000,
    nPerGroup = 5,
    propDE = 0.2,
    addMissing = TRUE,
    seed = 123
)

dim(simData$logData)
#> [1] 1000   10
head(simData$metadata)
#>   Samples Groups
#> 1    G1_1     G1
#> 2    G1_2     G1
#> 3    G1_3     G1
#> 4    G1_4     G1
#> 5    G1_5     G1
#> 6    G2_1     G2
```
