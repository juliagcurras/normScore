# Plot bootstrap normScore results

Generates a point-range plot showing the mean bootstrap normScore and
its 95% confidence interval for each normalization method.

## Usage

``` r
plotBootstrapNormScore(x)
```

## Arguments

- x:

  A `normScore` result object containing a `bootstrapScore` element.
  This element should be a data frame with normalization names, mean
  normScore values, and lower and upper 95% confidence limits. It can
  beobtained by setting returnDetails = `TRUE` at `normScore` main
  function.

## Value

A `ggplot` object showing bootstrap mean normScore values and 95%
confidence intervals for each normalization method.

## Examples

``` r
bootstrapScore <- data.frame(
    normalization = c("Norm1", "Norm2", "Norm3"),
    meanNormScore = c(0.8, 1.2, 1.5),
    ll95 = c(0.6, 1.0, 1.2),
    ul95 = c(1.0, 1.4, 1.8)
)

result <- list(bootstrapScore = bootstrapScore)

plotBootstrapNormScore(result)

```
