#' Utility functions for internal use in normScore
#'
#' This file contains helper functions used across the package,
#' including input validation, data formatting, and auxiliary computations.
#'
#' These functions are not intended to be used directly by end users.
#' 





#' Compare Two Normalized Datasets
#'
#' @description
#' Checks whether two normalized datasets are equivalent, allowing for small
#' numerical differences.
#'
#' @param x First normalized data matrix or data frame.
#' @param y Second normalized data matrix or data frame.
#' @param tolerance `numeric`. Numerical tolerance used for comparison.
#'   Default is `1e-8`.
#'
#' @return
#' A `logical` value indicating whether `x` and `y` are considered equivalent.
#'
#' @keywords internal

isSameNormalization <- function(x, y, tolerance = 1e-8) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  if (!identical(dim(x), dim(y))) {
    return(FALSE)
  }
  
  if (!identical(rownames(x), rownames(y))) {
    return(FALSE)
  }
  
  if (!identical(colnames(x), colnames(y))) {
    return(FALSE)
  }
  
  isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE))
}






#' Ensure Presence of Log2-Transformed Data in a Normalization List
#'
#' @description
#' Checks whether a log2-transformed version of the raw data is already present
#' in a list of normalized datasets. If not found, it is added. If found under a
#' different name, it is renamed to `"Log"`.
#'
#' @param normalizedDataList Named list of normalized data matrices or data
#'   frames. Each element represents a normalization method.
#' @param rawData `data.frame` or matrix of raw numeric values from which the
#'   log2 transformation will be computed.
#' @param tolerance `numeric`. Numerical tolerance used when comparing datasets.
#'   Default is `1e-8`.
#'
#' @return
#' A named list of normalized datasets, guaranteed to contain a log2-transformed
#' version of `rawData` under the name `"Log"`. If multiple matching datasets are
#' found, only the first one is kept.
#'
#'
#' @seealso
#' \code{\link{isSameNormalization}}
#'
#' @keywords internal

addLogIfMissing <- function(normalizedDataList, rawData, tolerance = 1e-8) {
  if (any(rawData <= 0, na.rm = TRUE)) {
    stop(
      "'rawData' contains values <= 0, so log2 transformation cannot be computed safely.",
      call. = FALSE
    )
  }
  
  logData <- log(rawData, base = 2)
  
  matchingIndex <- which(vapply(
    normalizedDataList,
    function(x) isSameNormalization(x, logData, tolerance = tolerance),
    logical(1)
  ))
  
  if (length(matchingIndex) == 0) {
    normalizedDataList[["Log"]] <- logData
    message("Log2-transformed raw data were added to 'normalizedDataList' as 'Log'.")
  } else {
    names(normalizedDataList)[matchingIndex[1]] <- "Log"
    
    if (length(matchingIndex) > 1) {
      normalizedDataList <- normalizedDataList[-matchingIndex[-1]]
      message("Multiple log2-equivalent normalizations were found. The first one was kept and renamed to 'Log'.")
    }
  }
  
  return(normalizedDataList)
}






#' Coefficient of Variation
#'
#' @description
#' Computes the coefficient of variation (CV) of a numeric vector as the ratio
#' between the standard deviation and the mean. The result can be returned as
#' a proportion or as a percentage.
#'
#' @param x `numeric`. Numeric vector for which the coefficient of variation
#'   will be calculated.
#' @param proportion `logical`. If `TRUE`, the result is returned as a proportion.
#'   If `FALSE`, it is multiplied by 100 and returned as a percentage.
#'   Default is `TRUE`.
#' @param na.rm `logical`. Should missing values (`NA`) be removed before
#'   computing the mean and standard deviation? Default is `TRUE`.
#'
#' @return
#' A numeric value representing the coefficient of variation of `x`.
#'
#' @details
#' The coefficient of variation is defined as:
#'
#' \deqn{CV = \frac{sd(x)}{mean(x)}}
#'
#' If `proportion = FALSE`, the result is multiplied by 100.
#'
#' Note that if the mean of `x` is 0 or very close to 0, the result may be
#' `Inf`, `NaN`, or unstable.
#'
#'
#' @keywords internal


cv <- function(x, proportion =TRUE, na.rm = TRUE) {
  coeficiente <- (stats::sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm))
  if (!proportion){
    coeficiente <- coeficiente*100
  }
  return(coeficiente)
} 
 






#' Mean Absolute Percentage Error
#'
#' @description
#' Computes the mean absolute percentage error (MAPE) between observed
#' (`actual`) and predicted values. The result can be returned as a proportion
#' or as a percentage.
#'
#' @param actual `numeric`. Vector of observed (true) values.
#' @param predicted `numeric`. Vector of predicted values. Must have the same
#'   length as `actual`.
#' @param proportion `logical`. If `TRUE`, the result is returned as a proportion.
#'   If `FALSE`, it is multiplied by 100 and returned as a percentage.
#'   Default is `FALSE`.
#'
#' @return
#' A numeric value representing the mean absolute percentage error.
#'
#' @details
#' The mean absolute percentage error is defined as:
#'
#' \deqn{MAPE = \frac{1}{n} \sum_{i=1}^{n} \left| \frac{actual_i - predicted_i}{actual_i} \right|}
#'
#' If `proportion = FALSE`, the result is multiplied by 100.
#'
#' Note that this metric is undefined when any value in `actual` is equal to 0.
#' In such cases, the function may return `Inf` or `NaN`.
#'
#' @keywords internal


mape <- function(actual, predicted, proportion = FALSE){
  metric <- mean(abs((actual - predicted)/actual))
  metric <- ifelse(!proportion, metric*100, metric)
  return(metric)
} 






#' Compute RLE-based MAPE metric
#'
#' Computes a normScore item based on Relative Log Expression (RLE). The data
#' are centered by protein medians and transformed to linear scale. The metric
#' is then calculated as the sum of mean absolute percentage errors (MAPE)
#' comparing sample medians and interquartile values to their expected
#' references.
#'
#' @param data A matrix or data frame with proteins in rows and samples in
#'   columns, typically in log-scale.
#'
#' @return A single numeric value representing the RLE-based MAPE metric, where
#'   lower values indicate better normalization performance.
#'
#' @keywords internal

rleMAPE <- function(data) {
  rowMedians <- apply(data, 1, stats::median, na.rm = TRUE)
  
  # rleData <- as.data.frame(t(t(data) / rowMedians))
  rleData <- data - rowMedians
  rleData <- 2^rleData
  
  # sampleMedians <- apply(rleData, 2, stats::median, na.rm = TRUE)
  # sampleMedians <- 2^sampleMedians # No log for MAPE
  
  sampleQuantiles <- apply(rleData, 2, stats::quantile, na.rm = TRUE, simplify = FALSE)
  q1 <- sapply(sampleQuantiles, "[[", 2)
  sampleMedians <- sapply(sampleQuantiles, "[[", 3)
  q3 <- sapply(sampleQuantiles, "[[", 4)
  
  finalMetric <- 
    mape(actual = 1, predicted = sampleMedians, proportion = TRUE) +
    mape(actual = stats::median(q1), predicted = q1, proportion = TRUE) +
    mape(actual = stats::median(q3), predicted = q3, proportion = TRUE)
  
  # return(mape(actual = 1, predicted = sampleMedians, proportion = FALSE))
  return(finalMetric)
} 






#' Quantile-Based MAPE Consistency Metric
#'
#' @description
#' Computes a sample consistency metric based on the mean absolute percentage
#' error (MAPE) of three summary statistics across samples: first quartile (Q1),
#' median, and third quartile (Q3).
#'
#' For each sample (column), the function extracts Q1, median, and Q3.
#' It then computes the MAPE of each of these vectors with respect to their
#' global median across samples, and returns the sum of the three resulting
#' components.
#'
#' @param data `data.frame` or matrix of numeric values, where columns typically
#'   represent samples and rows represent features (e.g. genes or proteins).
#'
#' @return
#' A numeric value representing the summed quantile-based MAPE consistency
#' metric across samples.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Computes the sample-wise first quartile (Q1), median, and third
#'   quartile (Q3).
#'   \item Computes the global median of each of these three vectors.
#'   \item Computes the mean absolute percentage error between each sample-wise
#'   vector and its corresponding global median.
#'   \item Returns the sum of the three MAPE components.
#' }
#'
#' Lower values indicate greater similarity among samples in terms of their
#' central tendency and spread.
#'
#' Note that if the global median of `q1`, `sampleMedians`, or `q3` is equal to
#' 0, the corresponding MAPE value may be undefined and produce `Inf`, `NaN`,
#' or unstable results.
#'
#' @seealso
#' \code{\link{mape}}
#'
#' @keywords internal

tiMAPE <- function(data) {
  sampleQuantiles <- apply(data, 2, stats::quantile, na.rm = TRUE, simplify = FALSE)
  q1 <- sapply(sampleQuantiles, "[[", 2)
  sampleMedians <- sapply(sampleQuantiles, "[[", 3)
  q3 <- sapply(sampleQuantiles, "[[", 4)
  
  finalMetric <- 
    mape(actual = stats::median(sampleMedians), predicted = sampleMedians, proportion = TRUE) +
    mape(actual = stats::median(q1), predicted = q1, proportion = TRUE) +
    mape(actual = stats::median(q3), predicted = q3, proportion = TRUE)
  
  return(finalMetric)
} 






#' Mean Absolute Area Difference Between Predicted and Expected Lines
#'
#' @description
#' Computes an area-based metric representing the average absolute difference
#' between a predicted linear function and an expected horizontal line over
#' a given range.
#'
#' The predicted line is defined by its intercept and slope, while the expected
#' line is assumed to be horizontal at `intExpected`. The function calculates
#' the absolute area between both lines across the interval
#' from `minRange` to `maxRange` and normalizes it by the interval width.
#'
#' @param intPred `numeric`. Intercept of the predicted regression line.
#' @param coefPred `numeric`. Slope of the predicted regression line.
#' @param minRange `numeric`. Lower bound of the interval over which the area
#'   difference is computed.
#' @param maxRange `numeric`. Upper bound of the interval over which the area
#'   difference is computed.
#' @param intExpected `numeric`. Intercept of the expected horizontal line.
#'   Default is `0`.
#'
#' @return
#' A numeric value corresponding to the normalized absolute area difference
#' between the predicted line and the expected line over the specified range.
#'
#' @details
#' The function first shifts the predicted intercept so that the expected line
#' is centered at 0. It then determines whether the predicted line crosses the
#' expected line within the interval from `minRange` to `maxRange`.
#'
#' If the crossing point lies within the interval, the area is computed as the
#' sum of the two signed sub-areas on either side of the intersection.
#' If no crossing occurs within the interval, the area is computed directly over
#' the full range.
#'
#' The final area metric is normalized by dividing by:
#'
#' \deqn{maxRange - minRange}
#'
#' If `coefPred = 0`, the function returns `0`.
#'
#'
#' @keywords internal

diffAreas <- function(intPred, coefPred, minRange, maxRange, intExpected = 0){
  
  # # Moving regression lines to reach B=0, a=0
  intPred <- intPred - intExpected
  
  if (coefPred < 0){
    a = 1
    b = -1
  } else if (coefPred > 0){
    a = -1
    b = 1
  } else if (coefPred == 0){
    return(0)
  }
  
  # Cutpoint regression line with expected line
  cpX <- (-intPred)/coefPred
  cpY <- 0
  
  # Is cutpoint located inside the range?
  if (all(cpX >= minRange, cpX <= maxRange)){
    area1 <- coefPred*a*(((minRange + (intPred/coefPred))^2)/2 - ((cpX + (intPred/coefPred))^2)/2)
    area2 <- coefPred*b*(((cpX + (intPred/coefPred))^2)/2 - ((maxRange + (intPred/coefPred))^2)/2)
    areaMetric <- abs(area1+area2)
  } else if (any(cpX < minRange, cpX > maxRange)){
    areaMetric <- abs(coefPred*(((maxRange + (intPred/coefPred))^2)/2 - ((minRange + (intPred/coefPred))^2)/2))
  }
  
  areaMetric <- areaMetric/(maxRange-minRange)
  return(areaMetric)
} 






#' Mean-SD Trend Area Metric Across Samples
#'
#' @description
#' Computes an area-based metric that quantifies how far the trend between
#' sample means and sample standard deviations deviates from horizontality.
#'
#' For each sample (column), the mean and standard deviation are calculated.
#' Samples are then ordered by their mean value, and a linear model is fitted
#' using sample standard deviation as the response and sample order as the
#' predictor. The metric evaluates how far the fitted trend is from a perfectly
#' horizontal line, considering only its slope and ignoring its intercept.
#'
#' @param data `data.frame` or matrix of numeric values, where columns typically
#'   represent samples and rows represent features (e.g. genes or proteins).
#'
#' @return
#' A numeric value representing the normalized area-based deviation of the
#' fitted mean-SD trend from horizontality.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Computes the mean and standard deviation of each sample (column).
#'   \item Orders samples according to their mean values.
#'   \item Fits a linear model with sample standard deviation as the response
#'   and sample order as the predictor.
#'   \item Extracts the slope of the fitted line.
#'   \item Computes an area-based metric with \code{\link{diffAreas}} after
#'   shifting the fitted line to intercept 0, so that only deviation from
#'   horizontality is assessed.
#' }
#'
#' A value of 0 indicates a perfectly horizontal fitted trend. Larger values
#' indicate stronger deviation from horizontality.
#'
#'
#' @seealso
#' \code{\link{diffAreas}}
#'
#' @keywords internal

meanSDdiffArea <- function(data) {
  sampleMeans <- apply(data, 2, base::mean, na.rm = TRUE)
  sampleSDs <- apply(data, 2, stats::sd, na.rm = TRUE)
  
  dfAux <- data.frame(
    Mean = sampleMeans,
    SD = sampleSDs,
    Sample = colnames(data)
  )
  
  dfAux <- dfAux[order(dfAux$Mean), ]
  dfAux$Order <- seq_len(nrow(dfAux))
  
  fit <- stats::lm(SD ~ Order, data = dfAux)
  slope <- fit$coefficients["Order"]
  
  resultArea <- diffAreas(
    intPred = 0,
    coefPred = slope,
    maxRange = max(dfAux$Order),
    minRange = min(dfAux$Order),
    intExpected = 0
  )
  
  return(unname(resultArea))
} 






#' MA-Trend Area Metric with Shape Correction
#'
#' @description
#' Computes an area-based metric to quantify deviation from an expected MA trend,
#' combining both slope/intercept deviation and shape consistency across the
#' expression range.
#'
#' The function first computes an MA-like representation from two groups of
#' samples, where `logFC` is the difference between group means and `AveExpr`
#' is their average. It then fits a linear model of `logFC` as a function of
#' `AveExpr` and calculates the normalized area between the fitted trend and the
#' expected horizontal line `y = 0` using \code{\link{diffAreas}}.
#'
#' In addition, the function evaluates the shape of the MA trend by dividing
#' `AveExpr` into 10 bins, computing the interquartile range (IQR) of `logFC`
#' within each bin, and comparing the observed IQR pattern to an expected
#' decreasing trend using Spearman correlation. This produces a correction
#' factor that downweights the area metric when the expected shape is preserved.
#'
#' @param data `data.frame` or matrix of numeric values, where rows typically
#'   represent features (e.g. genes or proteins) and columns represent samples.
#' @param samplesG1 `character` vector. Names of the columns in `data`
#'   corresponding to group 1.
#' @param samplesG2 `character` vector. Names of the columns in `data`
#'   corresponding to group 2.
#'
#' @return
#' A numeric value representing the shape-corrected area-based deviation of the
#' MA trend from the expected horizontal line `logFC = 0`.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Subsets `data` to the selected samples from both groups.
#'   \item Removes rows with missing values.
#'   \item Computes `logFC` as the difference between the mean of group 2 and
#'   the mean of group 1 for each row.
#'   \item Computes `AveExpr` as the average of the two group means for each row.
#'   \item Divides `AveExpr` into 10 bins and computes the IQR of `logFC` in
#'   each bin.
#'   \item Compares the observed IQR pattern with an expected decreasing trend
#'   using Spearman correlation, and transforms the result into a correction
#'   factor ranging from 0.1 to 1.
#'   \item Fits a linear model of `logFC ~ AveExpr`.
#'   \item Computes the normalized area between the fitted line and the expected
#'   horizontal line `y = 0` using \code{\link{diffAreas}}.
#'   \item Multiplies the area metric by the shape correction factor.
#' }
#'
#' Lower values indicate a trend closer to the expected MA relationship.
#'
#' Note that the function assumes that lower variability in `logFC` is expected
#' at higher average expression values, and penalizes departures from this
#' pattern through the shape correction factor.
#'
#'
#' @seealso
#' \code{\link{diffAreas}}
#'
#' @keywords internal

maDiffArea <- function(data, samplesG1, samplesG2) {
  data <- as.data.frame(data)
  maData <- data[, c(samplesG1, samplesG2)]
  maData <- stats::na.omit(maData)
  
  maData$logFC <- apply(maData, 1, function(x) {
    mean(x[samplesG2], na.rm = TRUE) - mean(x[samplesG1], na.rm = TRUE)
  })
  
  maData$AveExpr <- apply(maData, 1, function(x) {
    (mean(x[samplesG2], na.rm = TRUE) + mean(x[samplesG1], na.rm = TRUE)) / 2
  })
  
  expressionDeciles <- stats::quantile(maData$AveExpr, probs = seq(0, 1, 0.1))
  
  iqrByBin <- sapply(seq_len(length(expressionDeciles) - 1), function(i) {
    logFCValues <- maData$logFC[
      maData$AveExpr > expressionDeciles[i] &
        maData$AveExpr <= expressionDeciles[i + 1]
    ]
    
    stats::IQR(logFCValues, na.rm = TRUE)
  })
  
  rho <- stats::cor(iqrByBin, 10:1, method = "spearman")
  rho <- (1 - rho) / 2
  shapeCorrectionFactor <- 0.1 + (1 - 0.1) * rho
  
  maData <- maData[, c("logFC", "AveExpr")]
  
  fit <- stats::lm(logFC ~ AveExpr, data = maData)
  slope <- fit$coefficients["AveExpr"]
  intercept <- fit$coefficients["(Intercept)"]
  
  resultArea <- diffAreas(
    intPred = intercept,
    coefPred = slope,
    maxRange = max(maData$AveExpr),
    minRange = min(maData$AveExpr),
    intExpected = 0
  )
  
  correctedMetric <- unname(resultArea) * shapeCorrectionFactor
  return(correctedMetric)
}






#' Coefficient of Variation per Protein Within a Group
#'
#' @description
#' Computes the coefficient of variation (CV) for each protein or feature
#' within a specified group of samples.
#'
#' The function uses a sample-group mapping table to identify the samples
#' belonging to the requested group, subsets the expression matrix accordingly,
#' and computes the coefficient of variation row-wise using \code{\link{cv}}.
#'
#' @param group A single value identifying the group for which the CV should be
#'   calculated.
#' @param groupData `data.frame` containing the mapping between samples and
#'   groups. It must contain at least the columns `Groups` and `Samples`.
#' @param data `data.frame` or matrix of numeric values, where rows typically
#'   represent proteins or features and columns represent samples.
#'
#' @return
#' A numeric vector containing the coefficient of variation for each protein
#' or feature within the selected group.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Identifies the samples belonging to the specified group using `groupData`.
#'   \item Subsets `data` to those samples.
#'   \item Computes the coefficient of variation for each row using
#'   \code{\link{cv}} with `proportion = FALSE`.
#' }
#'
#' The returned values are expressed as percentages.
#'
#'
#' @seealso
#' \code{\link{cv}}
#'
#' @keywords internal

groupProteinCV <- function(group, groupData, data) {
  groupData <- as.data.frame(groupData)
  data <- as.data.frame(data)
  
  selectedSamples <- as.vector(unlist(groupData[groupData$Groups == group, "Samples"]))
  
  subsetData <- data[, as.character(selectedSamples)]
  apply(subsetData, 1, cv, proportion = FALSE)
} 






#' Pooled Coefficient of Variation 
#'
#' @description
#' Computes the mean coefficient of variation (CV) across groups for each
#' protein or feature.
#'
#' For each group specified in `groups`, the function applies
#' \code{\link{groupProteinCV}} to calculate group-specific CV values from the
#' input data. It then computes the mean of these CV values for each protein
#' or feature across groups.
#'
#' @param data `data.frame` or matrix of numeric values, where rows typically
#'   represent proteins or features and columns represent samples.
#' @param groups A vector specifying the groups for which CV values should be
#'   calculated.
#' @param groupData `data.frame` containing the group annotation or sample-group
#'   mapping required by \code{\link{groupProteinCV}}.
#'
#' @return
#' A named numeric vector containing the pooled coefficient of
#' variation for each protein or feature.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Applies \code{\link{groupProteinCV}} to each group in `groups`.
#'   \item Combines the resulting group-wise CV values into a data frame.
#'   \item Computes the mean CV across groups for each protein or feature.
#' }
#'
#' Higher values indicate greater relative variability across samples within
#' groups.
#'
#'
#' @seealso
#' \code{\link{groupProteinCV}}
#'
#' @keywords internal

getPCV <- function(data, groups, groupData){
  groupDataProt <- as.data.frame(
    sapply(groups, groupProteinCV, 
           groupData = groupData, 
           data = data)
  )
  meanVal <- apply(groupDataProt, 2, mean, na.rm = TRUE)
  names(meanVal) <- colnames(groupDataProt)
  return(meanVal)
} 




#' Pairwise Correlations Within Groups
#'
#' @description
#' Computes pairwise correlations between samples within each group and returns
#' all unique correlation values as a single numeric vector.
#'
#' The function uses a sample-group mapping table to identify the samples
#' belonging to each group, computes the sample-sample correlation matrix within
#' each group, extracts the unique off-diagonal correlation values, and combines
#' them across groups.
#'
#' @param data `data.frame` or matrix of numeric values, where rows typically
#'   represent features (e.g. genes or proteins) and columns represent samples.
#' @param groupData `data.frame` containing the mapping between samples and
#'   groups. It must contain at least the columns `Groups` and `Samples`.
#' @param method `character`. Correlation method to be used. Passed to
#'   \code{\link[stats]{cor}}. Default is `"pearson"`.
#'
#' @return
#' A numeric vector containing all unique pairwise correlation values computed
#' within groups.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Identifies the samples belonging to each group in `groupData`.
#'   \item Subsets `data` to the samples of each group.
#'   \item Computes the sample-sample correlation matrix within each group.
#'   \item Removes diagonal values and duplicated pairwise correlations.
#'   \item Combines all unique correlation values into a single vector.
#' }
#'
#' Missing values are handled using `use = "complete.obs"` in
#' \code{\link[stats]{cor}}.
#'
#' @keywords internal

withinGroupCorrelations <- function(data, groupData, method = "spearman") {
  data <- as.data.frame(data)
  groupData <- as.data.frame(groupData)
  
  allCorrelations <- lapply(unique(groupData$Groups), function(group) {
    samplesByGroup <- groupData[groupData$Groups == group, "Samples"]
    
    groupMatrix <- data[, samplesByGroup, drop = FALSE]
    
    corMatrix <- stats::cor(groupMatrix, use = "pairwise.complete.obs", method = method)
    
    corValues <- c()
    for (i in seq_len(ncol(groupMatrix) - 1)) {
      corValues <- c(corValues, corMatrix[i, -(seq_len(i))])
    }
    
    corValues
  })
  
  unlist(allCorrelations, use.names = FALSE)
}






#' Compute within-group correlation item
#'
#' Computes the normScore correlation item for each normalization method using
#' within-group sample correlations. The item is summarized from the median and
#' interquartile range of the correlation values, and transformed so that lower
#' values indicate better performance.
#'
#' @param normalizedDataList A named list of normalized data matrices, with
#'   proteins in rows and samples in columns.
#' @param groupData A data frame containing sample-group annotation.
#' @param method Character string indicating the correlation method passed to
#'   `stats::cor()`, for example `"spearman"` or `"pearson"`.
#'
#' @return A named numeric vector with one correlation item score per
#'   normalization method.
#'
#' @keywords internal


computeCorrelation <- function(normalizedDataList, groupData, method)
  {
  dfCor <- sapply(
    normalizedDataList,
    withinGroupCorrelations,
    groupData = groupData,
    method = method, 
    simplify = FALSE, 
    USE.NAMES = TRUE
  )
  
  sapply(
    names(dfCor),
    function(normalizationName) {
      correlationValues <- dfCor[[normalizationName]]
      1 - (stats::median(correlationValues, na.rm = TRUE) - stats::IQR(correlationValues, na.rm = TRUE) / 3)
    },
    simplify = TRUE,
    USE.NAMES = TRUE
  )
  
}





#' Bootstrap Total Scores After Row Resampling
#'
#' @description
#' Computes total scores for each column after bootstrap resampling of rows.
#'
#' The function subsets the input matrix or data frame using a vector of row
#' indices, then calculates the column-wise sums of the resampled data.
#' It is intended for use as a bootstrap statistic when evaluating total scores
#' across normalizations or methods.
#'
#' @param data `data.frame` or matrix of numeric values, where rows typically
#'   represent features (e.g. proteins) and columns represent scores for
#'   different normalizations or methods.
#' @param indices `integer` vector of row indices used for bootstrap resampling.
#'
#' @return
#' A numeric vector containing one total score per column in the resampled data.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Subsets `data` using the provided row indices.
#'   \item Computes the sum of each column in the resampled data.
#'   \item Returns one total score per column.
#' }
#'
#' This function is typically used as a statistic function inside bootstrap
#' procedures.
#'
#' @keywords internal

bootstrapRowScores <- function(data, indices) {
  resampledMatrix <- data[indices, , drop = FALSE]
  totalScores <- colSums(resampledMatrix)
  return(totalScores)
}



