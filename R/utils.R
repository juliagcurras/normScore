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
#' @examples
#' x <- c(10, 12, 9, 11, 10)
#'
#' # As proportion
#' cv(x)
#'
#' # As percentage
#' cv(x, proportion = FALSE)
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
#' @examples
#' actual <- c(100, 200, 300, 400)
#' predicted <- c(110, 190, 320, 390)
#'
#' # As percentage
#' mape(actual, predicted)
#'
#' # As proportion
#' mape(actual, predicted, proportion = TRUE)
#'
#' @export
#' @keywords internal

mape <- function(actual, predicted, proportion = FALSE){
  metric <- mean(abs((actual - predicted)/actual))
  metric <- ifelse(!proportion, metric*100, metric)
  return(metric)
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
#' [`minRange`, `maxRange`] and normalizes it by the interval width.
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
#' expected line within the interval [`minRange`, `maxRange`].
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
#' @examples
#' # Predicted line crosses the expected line within the range
#' diffAreas(
#'   intPred = 2,
#'   coefPred = -0.5,
#'   minRange = 0,
#'   maxRange = 10
#' )
#'
#' # Expected line different from 0
#' diffAreas(
#'   intPred = 3,
#'   coefPred = 0.2,
#'   minRange = 1,
#'   maxRange = 8,
#'   intExpected = 1
#' )
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






#' Compute an RLE-based MAPE metric across samples
#'
#' Calculates a sample-level metric based on Relative Log Expression (RLE)
#' transformed values and mean absolute percentage error (MAPE). For each row,
#' the row median is subtracted from the original values, the result is
#' back-transformed with `2^x`, and sample-wise quartiles are computed.
#' The final metric is the sum of:
#' \itemize{
#'   \item MAPE between 1 and the sample medians
#'   \item MAPE between the median of first quartiles and each sample first quartile
#'   \item MAPE between the median of third quartiles and each sample third quartile
#' }
#'
#' Lower values indicate greater similarity across samples after RLE-based
#' centering.
#'
#' @param data A numeric matrix-like object with features in rows and samples
#'   in columns. Missing values (`NA`) are allowed.
#'
#' @return A single numeric value corresponding to the final RLE-based MAPE
#'   metric.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Computes the row-wise medians.
#'   \item Centers each row by subtracting its median.
#'   \item Applies back-transformation using `2^x`.
#'   \item Computes sample-wise first quartile, median, and third quartile.
#'   \item Combines three MAPE terms into a final summary metric.
#' }
#'
#' This function assumes that a compatible `mape()` function is available in the
#' package namespace.
#'
#' @examples
#' mat <- matrix(
#'   c(1, 2, 3,
#'     2, 3, 4,
#'     5, 6, 7),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#'
#' rleMAPE(mat)
#'
#' @seealso
#' \code{\link{mape}}
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
#' @examples
#' dataExample <- data.frame(
#'   Sample1 = c(10, 20, 30, 40, 50),
#'   Sample2 = c(11, 19, 29, 41, 49),
#'   Sample3 = c(12, 21, 28, 42, 48)
#' )
#'
#' tiMAPE(dataExample)
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
#' @examples
#' dataExample <- data.frame(
#'   Sample1 = c(10, 20, 30, 40),
#'   Sample2 = c(12, 19, 29, 41),
#'   Sample3 = c(11, 21, 31, 39)
#' )
#'
#' meanSDdiffArea(dataExample)
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
#' @examples
#' dataExample <- data.frame(
#'   G1_1 = c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
#'            20, 21, 22, 23, 24, 25, 26, 27, 28, 29),
#'   G1_2 = c(10.5, 11.5, 11.8, 13.2, 14.1, 15.3, 15.7, 17.4, 18.2, 19.1,
#'            20.1, 21.2, 21.8, 23.1, 24.0, 25.2, 26.1, 27.1, 28.2, 29.0),
#'   G2_1 = c(13.5, 14.0, 14.2, 15.8, 16.4, 17.1, 17.8, 18.9, 19.5, 20.2,
#'            21.0, 21.8, 22.5, 23.7, 24.6, 25.4, 26.5, 27.4, 28.5, 29.4),
#'   G2_2 = c(14.5, 14.3, 14.8, 16.0, 16.6, 17.3, 18.0, 19.0, 19.6, 20.3,
#'            21.1, 21.9, 22.6, 23.8, 24.7, 25.5, 26.6, 27.5, 28.6, 29.5)
#' )
#'
#' maDiffArea(
#'   data = dataExample,
#'   samplesG1 = c("G1_1", "G1_2"),
#'   samplesG2 = c("G2_1", "G2_2")
#' )
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
#' @examples
#' dataExample <- data.frame(
#'   S1 = c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
#'            20, 21, 22, 23, 24, 25, 26, 27, 28, 29),
#'   S2 = c(10.5, 11.5, 11.8, 13.2, 14.1, 15.3, 15.7, 17.4, 18.2, 19.1,
#'            20.1, 21.2, 21.8, 23.1, 24.0, 25.2, 26.1, 27.1, 28.2, 29.0),
#'   S3 = c(13.5, 14.0, 14.2, 15.8, 16.4, 17.1, 17.8, 18.9, 19.5, 20.2,
#'            21.0, 21.8, 22.5, 23.7, 24.6, 25.4, 26.5, 27.4, 28.5, 29.4),
#'   S4 = c(14.5, 14.3, 14.8, 16.0, 16.6, 17.3, 18.0, 19.0, 19.6, 20.3,
#'            21.1, 21.9, 22.6, 23.8, 24.7, 25.5, 26.6, 27.5, 28.6, 29.5)
#' )
#'
#' groupDataExample <- data.frame(
#'   Samples = c("S1", "S2", "S3", "S4"),
#'   Groups = c("A", "A", "B", "B")
#' )
#'
#' groupProteinCV(
#'   group = "A",
#'   groupData = groupDataExample,
#'   data = dataExample
#' )
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
#' @examples
#' # Example structure only
#' dataExample <- data.frame(
#'   S1 = c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
#'            20, 21, 22, 23, 24, 25, 26, 27, 28, 29),
#'   S2 = c(10.5, 11.5, 11.8, 13.2, 14.1, 15.3, 15.7, 17.4, 18.2, 19.1,
#'            20.1, 21.2, 21.8, 23.1, 24.0, 25.2, 26.1, 27.1, 28.2, 29.0),
#'   S3 = c(13.5, 14.0, 14.2, 15.8, 16.4, 17.1, 17.8, 18.9, 19.5, 20.2,
#'            21.0, 21.8, 22.5, 23.7, 24.6, 25.4, 26.5, 27.4, 28.5, 29.4),
#'   S4 = c(14.5, 14.3, 14.8, 16.0, 16.6, 17.3, 18.0, 19.0, 19.6, 20.3,
#'            21.1, 21.9, 22.6, 23.8, 24.7, 25.5, 26.6, 27.5, 28.6, 29.5)
#' )
#'
#' groupDataExample <- data.frame(
#'   Samples = c("S1", "S2", "S3", "S4"),
#'   Groups = c("A", "A", "B", "B")
#' )
#' groups <- c("A", "B")
#'
#' # Requires groupProteinCV() to be available
#' getPCV(dataExample, groups, groupDataExample)
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
#' @examples
#' dataExample <- data.frame(
#'   Sample1 = c(10, 20, 30, 40),
#'   Sample2 = c(11, 21, 29, 39),
#'   Sample3 = c(15, 18, 35, 45),
#'   Sample4 = c(14, 19, 34, 44)
#' )
#'
#' groupDataExample <- data.frame(
#'   Samples = c("Sample1", "Sample2", "Sample3", "Sample4"),
#'   Groups = c("A", "A", "B", "B")
#' )
#'
#' withinGroupCorrelations(
#'   data = dataExample,
#'   groupData = groupDataExample,
#'   method = "pearson"
#' )
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
#' @examples
#' dataExample <- data.frame(
#'   Norm1 = c(1.2, 0.8, 1.1, 0.9),
#'   Norm2 = c(1.0, 0.7, 1.3, 1.1),
#'   Norm3 = c(0.9, 0.6, 1.2, 1.0)
#' )
#'
#' indicesExample <- c(1, 2, 2, 4)
#'
#' bootstrapRowScores(dataExample, indicesExample)
#'
#' @keywords internal

bootstrapRowScores <- function(data, indices) {
  resampledMatrix <- data[indices, , drop = FALSE]
  totalScores <- colSums(resampledMatrix)
  return(totalScores)
}






#' Bootstrap Estimation of Normalization Scores
#'
#' @description
#' Computes bootstrap-based mean scores and percentile confidence intervals
#' for normalization methods from a matrix of item-wise scores.
#'
#' The function applies bootstrap resampling over rows of the input score matrix,
#' where each row typically represents a scoring component (item) and each column
#' represents a normalization method. For each resample, total scores per
#' normalization are computed using \code{\link{bootstrapRowScores}}.
#'
#' @param scoreMatrix `numeric` matrix where rows correspond to scoring items
#'   and columns correspond to normalization methods.
#' @param nBoot `integer`. Number of bootstrap resamples. Default is `1000`.
#'
#' @examples
#' scoreMatrix <- matrix(
#'   c(0.2, 0.3, 0.1,
#'     0.4, 0.2, 0.3,
#'     0.1, 0.5, 0.2),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#'
#' colnames(scoreMatrix) <- c("Norm1", "Norm2", "Norm3")
#'
#' computeBootstrapNormScore(scoreMatrix, nBoot = 100)
#'
#' @seealso
#' \code{\link{bootstrapRowScores}}
#'
#' @keywords internal

computeBootstrapNormScore <- function(scoreMatrix, nBoot = 1000) {
  bootResults <- boot::boot(
    data = scoreMatrix,
    statistic = bootstrapRowScores,
    R = nBoot
  )
  
  bootstrapMeans <- colMeans(bootResults$t)
  
  bootstrapScores <- as.data.frame(
    t(
      sapply(seq_len(ncol(scoreMatrix)), function(i) {
        ci <- boot::boot.ci(bootResults, type = "perc", index = i)
        c(
          colnames(scoreMatrix)[i],
          bootstrapMeans[i],
          ci$percent[4],
          ci$percent[5]
        )
      }, simplify = TRUE)
    )
  )
  
  colnames(bootstrapScores) <- c(
    "normalization",
    "meanNormScore",
    "ll95",
    "ul95"
  )
  bootstrapScores[-1] <- lapply(bootstrapScores[-1], as.numeric)
  
  bootstrapScores <- bootstrapScores[order(bootstrapScores$meanNormScore), ]
  
  return(bootstrapScores)
}







#' Plot Bootstrap Normalization Scores
#'
#' @description
#' Creates a summary plot of bootstrap-based normalization scores with
#' percentile confidence intervals.
#'
#' The function takes the output table generated by
#' \code{\link{computeBootstrapNormScore}} and produces a point-range style plot
#' showing the mean normScore and its 95\% confidence interval for each
#' normalization method.
#'
#' @param bootstrapScores `data.frame` containing bootstrap summary results for
#'   normalization methods. It must include the columns `"Normalization"`,
#'   `"Mean normScore"`, `"LL95%"`, and `"UL95%"`.
#'
#' @return
#' A `ggplot2` object showing the mean normScore and 95\% confidence interval
#' for each normalization method.
#'
#' @details
#' The plot displays one point per normalization method, corresponding to the
#' mean bootstrap normScore, together with a horizontal error bar representing
#' the 95\% confidence interval.
#'
#' Normalization methods are ordered by increasing mean normScore, so that
#' better-performing methods appear according to their ranking.
#'
#' @examples
#' bootstrapScores <- data.frame(
#'   normalization = c("Norm1", "Norm2", "Norm3"),
#'   meanNormScore = c(0.25, 0.40, 0.32),
#'   ll95 = c(0.20, 0.35, 0.28),
#'   ul95 = c(0.30, 0.45, 0.36)
#' )
#'
#' plotBootstrapNormScore(bootstrapScores)
#'
#' @seealso
#' \code{\link{computeBootstrapNormScore}}
#'
#' @export
#' @keywords internal 

plotBootstrapNormScore <- function(bootstrapScores) {
  
  # Checking colnames
  namesCols <- c(
    "normalization",
    "meanNormScore",
    "ll95",
    "ul95"
  )
  
  if (!all(colnames(bootstrapScores) == namesCols)){
    colnames(bootstrapScores) <- namesCols
  }
  
  # Plot!
  p <- ggplot2::ggplot(
    bootstrapScores,
    ggplot2::aes(
      x = .data$meanNormScore,
      y = stats::reorder(.data$normalization, .data$meanNormScore)
    )
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        xmin = .data$ll95,
        xmax = .data$ul95
      ),
      width = 0.2,
      orientation = "y"
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      x = "Mean normScore [95% CI]",
      y = "Normalization"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black")
    )
  
  return(p)
}


#' @importFrom rlang .data
NULL