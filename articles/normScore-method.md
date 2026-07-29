# Methodological overview of normScore

## Introduction

Normalization is a critical step in quantitative omics data analysis.
Different normalization methods can alter the structure of the data in
different ways, affecting downstream results such as differential
analysis, sample clustering, and fold-change estimation.

The `normScore` package provides a composite scoring framework to
compare normalization methods using multiple complementary criteria.

Rather than relying on a single summary metric, `normScore` integrates
different aspects of normalization performance into a unified ranking.

## General rationale

A good normalization method should ideally:

- reduce unwanted technical variability
- preserve biological structure
- avoid introducing systematic distortions
- maintain sample comparability
- stabilize distributions across samples

Because these goals are not fully captured by a single number,
`normScore` uses a multi-item scoring approach.

Lower final scores indicate better normalization performance relative to
the other methods under comparison.

## Inputs required by `normScore()`

The main function of the package expects three key inputs:

1.  a named list of normalized datasets
2.  a sample annotation table with group membership
3.  the raw data matrix

Each normalized dataset should contain the same features and samples,
arranged in the same order.

The raw data are used to derive the log-transformed baseline
normalization and to compute a correction factor related to global
sample intensity variation.

## The six scoring items

The composite score is based on six item-wise metrics.

### Item 1. Pooled coefficient of variation (PCV)

This item summarizes within-group variability across samples.

For each group, the coefficient of variation is computed feature-wise,
and the mean value is summarized across groups. Lower values indicate
lower relative dispersion within biological groups.

This item captures whether a normalization method improves sample
consistency without reference to differential structure.

### Item 2. Within-group sample correlation

This item evaluates pairwise sample correlations within each group.

For each group, all unique pairwise sample-sample correlations are
extracted, and a summary statistic based on the median and interquartile
range is used. The score is transformed so that lower values indicate
better performance.

This item reflects whether samples from the same group become more
similar after normalization.

### Item 3. MA-trend deviation

This item assesses whether the relationship between fold change and
average expression shows systematic distortion.

For two selected groups, an MA-like summary is constructed using:

- `logFC`: the difference between group means
- `AveExpr`: the average of group means

A linear model is fitted, and the area between the fitted line and the
expected horizontal line `logFC = 0` is computed. An additional shape
correction factor penalizes undesirable trends in the spread of `logFC`
across the expression range.

Lower values indicate a flatter and more stable MA trend.

### Item 4. Mean-SD trend deviation

This item evaluates whether sample standard deviations depend
systematically on sample mean intensity.

For each sample, the mean and standard deviation are computed. Samples
are ordered by mean intensity, and a linear model is fitted using
standard deviation as a function of sample order.

The slope of this fitted trend is then converted into an area-based
metric measuring deviation from horizontality.

Lower values indicate a weaker dependency between sample mean and sample
variability.

### Item 5. RLE consistency

This item is inspired by the relative log expression (RLE) concept.

For each feature, the median across samples is used as a reference.
Relative expression values are obtained by dividing each observation by
the feature-wise median, and the median relative expression is then
computed sample-wise.

The final score is the mean absolute percentage error of these
sample-wise medians relative to 1.

Lower values indicate that sample-wise medians are closer to the
expected reference value.

### Item 6. Total intensity consistency

This item compares the distributional center and spread of sample
intensity profiles.

For each sample, the first quartile, median, and third quartile are
computed. The score is based on the sum of the mean absolute percentage
errors of these three summaries relative to their global sample medians.

Lower values indicate more consistent sample distributions.

## Scaling and aggregation

Because the six items are on different scales, each one is transformed
using min-max scaling before aggregation.

After scaling, the six items are summed into a total score.

This makes the final ranking easier to interpret and ensures that no
single item dominates only because of its numerical range.

## Correction applied to the log-transformed baseline

The log-transformed raw data are included in the normalization list as a
baseline reference.

A correction factor derived from the coefficient of variation of total
raw sample intensities is applied specifically to the `"Log"` method.
This is used to account for the extent of global intensity imbalance
already present in the raw data.

## Additional adjustments and weighting

During the development of the scoring framework, two additional
adjustments were introduced to improve the balance between items and
avoid dominance by specific components.

First, the contribution of the within-group correlation metric (Item 2)
was reduced. This item was observed to have relatively low
discriminative power across normalization methods, often yielding very
similar values. To prevent it from contributing disproportionately to
the final score despite its limited ability to differentiate methods,
its scaled value was down-weighted by a factor of 0.1.

Second, the correction factor applied to the log-transformed baseline
normalization (Item 0) was moderated. This factor is derived from the
coefficient of variation of total sample intensities in the raw data and
is intended to account for global intensity imbalance.

However, without adjustment, this correction could have an excessive
influence on the final score, particularly for datasets with very low or
very high values of the correction factor. To mitigate this effect, the
magnitude of the correction was reduced by scaling it down (e.g.,
dividing or adjusting its impact), ensuring that it contributes to the
ranking without dominating it.

These adjustments were introduced to maintain a balanced contribution
across items and to improve the stability and interpretability of the
final ranking.

## Bootstrap summary

When `returnDetails = FALSE`,
[`normScore()`](https://juliagcurras.github.io/normScore/reference/normScore.md)
also performs bootstrap resampling over the six item scores.

This produces:

- a mean bootstrap normScore for each normalization method
- a 95% percentile confidence interval
- a graphical summary of uncertainty

The bootstrap output provides a more stable comparative summary and
helps assess the robustness of the ranking.

## Interpretation of the final score

The final score should be interpreted as a **relative ranking
criterion** within the set of normalization methods provided by the
user.

It is not an absolute measure of normalization quality.

A method with a lower score is preferred relative to the alternatives in
the same comparison, but the interpretation depends on:

- the candidate normalization methods included
- the data structure
- the biological groups
- the intended downstream analysis

## Practical considerations

A few points should be kept in mind when using `normScore`:

- at least 100 proteins/features are required
- all normalization matrices must contain the same samples and features
- group labels should reflect the intended biological comparison
- bootstrap is optional but recommended for more stable summaries
- the final ranking is comparative, not absolute

## Conclusion

`normScore` provides a structured and extensible framework for comparing
normalization methods using multiple complementary criteria.

By combining dispersion, correlation, trend, and distributional metrics,
the package aims to offer a more robust basis for selecting a
normalization method than any single metric alone.

## Session information

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.60         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.1     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.6.1    tools_4.6.1       ragg_1.5.2        bslib_0.11.0     
#> [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
#> [25] rlang_1.3.0       fs_2.1.0
```
