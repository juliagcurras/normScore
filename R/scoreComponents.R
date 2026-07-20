
#' Validate and prepare inputs for normScore
#'
#' Checks the input objects used by `normScore()`, including group annotation,
#' raw data, normalized datasets, sample names, protein names, and bootstrap
#' settings. The function also aligns the order of rows and columns across
#' normalized datasets, adds the log-transformed dataset when required, and
#' identifies whether the analysis contains a single group.
#'
#' @param normalizedDataList A named list of normalized data matrices or data
#'   frames. Each element should contain proteins in rows and samples 
#'   in columns.
#' @param groupData A data frame or object coercible to a data frame containing
#'   sample-group annotation. The first column is assumed to contain sample 
#'   names and the second column group labels.
#' @param rawData A matrix or data frame containing raw intensity values, with
#'   proteins in rows and samples in columns.
#' @param refGroup Character string indicating the reference group used for
#'   group-comparison metrics. If `NULL` or invalid, it is automatically set.
#' @param altGroup Character string indicating the alternative group used for
#'   group-comparison metrics. If `NULL` or invalid, it is automatically set.
#' @param returnDetails Logical value indicating whether detailed results
#'    should be returned by `normScore()`.
#' @param nBoot Numeric value indicating the number of bootstrap resamples.
#'
#' @return A list containing validated and prepared inputs: 
#'   `normalizedDataList`, `groupData`, `rawData`, `refGroup`, `altGroup`,
#'    `returnDetails`, `nBoot`, and `singleGroup`.
#'
#' @keywords internal
#' @noRd

validateNormScoreInput <- function(
    normalizedDataList,
    groupData,
    rawData,
    refGroup,
    altGroup,
    returnDetails = NULL,
    nBoot = NULL,
    fromMainFunction = TRUE
) {
  groupInfo <- .validateGroupData(groupData, refGroup, altGroup)
  
  rawData <- .validateRawData(rawData)
  
  normalizedDataList <- .coerceNormalizedData(normalizedDataList)
  
  .validateDimnames(normalizedDataList, rawData, groupInfo$groupData)
  
  normalizedDataList <- .alignNormalizedData(normalizedDataList, rawData)
  
  normalizedDataList <- addLogIfMissing(
    normalizedDataList = normalizedDataList,
    rawData = rawData)
  
  mainArgs <- .validateMainArguments(
    returnDetails,
    nBoot,
    fromMainFunction)
  
  list(
    normalizedDataList = normalizedDataList,
    groupData = groupInfo$groupData,
    rawData = rawData,
    refGroup = groupInfo$refGroup,
    altGroup = groupInfo$altGroup,
    returnDetails = mainArgs$returnDetails,
    nBoot = mainArgs$nBoot,
    singleGroup = groupInfo$singleGroup
  )
}




#' Compute normScore item scores
#'
#' Computes the individual quality items used by `normScore()` for each
#' normalization method. These include within-group variability, sample
#' correlation, MA-plot deviation, mean-SD trend, RLE deviation, and total
#' intensity consistency. When only one group is available, the MA-plot item is
#' skipped.
#'
#' @param normalizedDataList A named list of normalized data matrices, with
#'   proteins in rows and samples in columns.
#' @param groupData A data frame containing sample-group annotation, with 
#'   sample names in the first column and group labels in the second column.
#' @param rawData A matrix or data frame containing raw intensity values. Used
#'   for consistency with the main workflow.
#' @param refGroup Character string indicating the reference group used for the
#'   MA-plot item.
#' @param altGroup Character string indicating the alternative group used for 
#'   the MA-plot item.
#' @param singleGroup Logical value indicating whether the input data contain a
#'   single group.
#'
#' @return A data frame with normalization methods in rows and normScore items
#'   in columns. If `singleGroup = TRUE`, the MA-plot item is omitted.
#'
#' @keywords internal
#' @noRd

computeNormScoreItems <- function(normalizedDataList, groupData, rawData, 
                                  refGroup, altGroup, singleGroup){
  #----- ITEM 1 - PCV
  dfPCV <- do.call(rbind, lapply(normalizedDataList, getPCV, 
                                 groups = levels(as.factor(groupData$Groups)),
                                 groupData = groupData))
  item1 <- apply(dfPCV, 1, mean, na.rm = TRUE)
  
  #----- ITEM 2 - Correlation (Spearman)
  item2 <- computeCorrelation(normalizedDataList = normalizedDataList,
                              groupData = groupData, method ="spearman") 
  
  #----- ITEM 3 - MA plot regression line vs expected 0
  if (!singleGroup){
    item3 <- vapply(normalizedDataList, maDiffArea,
      samplesG1 =  groupData[groupData$Groups == refGroup, "Samples"],
      samplesG2 = groupData[groupData$Groups == altGroup, "Samples"], 
      numeric(1))
  } else {
    item3 <- stats::setNames(rep(NA, length(normalizedDataList)), 
                             names(normalizedDataList))
  }
  
  #----- ITEM 4 - Mean-SD regression line slope deviation
  item4 <- vapply(normalizedDataList, meanSDdiffArea, numeric(1))
  
  #----- ITEM 5 - RLE: MAPE of sample medians (reference = 1)
  item5 <- vapply(normalizedDataList, rleMAPE, numeric(1))
  
  #----- ITEM 6 - Total intensity consistency
  item6 <- vapply(normalizedDataList, tiMAPE, numeric(1))
  
  #---- Score matrix
  refNames <- names(item1)
  scoreDF <- cbind(Item1 = item1[refNames], Item2 = item2[refNames],
                   Item3 = item3[refNames], Item4 = item4[refNames],
                   Item5 = item5[refNames], Item6 = item6[refNames])
  scoreDF <- as.data.frame(scoreDF)
  rownames(scoreDF) <- refNames
  
  # Selecting only items from scaled matrix
  if (singleGroup){
    itemsToSelect <- paste0("Item", c(seq_len(2), seq(4, 6)))
  } else {
    itemsToSelect <- paste0("Item", seq_len(6))
  }
  scoreDF <- scoreDF[, itemsToSelect]
  
  return(scoreDF)
}





#' Scale normScore item scores
#'
#' Applies min-max scaling to each normScore item so that item values are
#' placed on a comparable scale across normalization methods. The correlation
#' item is additionally down-weighted after scaling.
#'
#' @param scoreDF A data frame containing raw normScore item values, with
#'   normalization methods in rows and items in columns.
#'
#' @return A data frame with the same dimensions as `scoreDF`, containing 
#'  scaled item scores. Row names are preserved.
#'
#' @keywords internal
#' @noRd

scaleNormScoreItems <- function(
    scoreDF
){
  #----- Min-max scaling
  scaledScoreDF <- apply(scoreDF, 2, function(columnValues) {
    columnMin <- min(columnValues, na.rm = TRUE)
    columnMax <- max(columnValues, na.rm = TRUE)
    
    if (columnMax == columnMin) {
      return(rep(0, length(columnValues)))
    }
    
    (columnValues - columnMin) / (columnMax - columnMin)
  })
  scaledScoreDF <- as.data.frame(scaledScoreDF)
  rownames(scaledScoreDF) <- rownames(scoreDF)
  
  #---- Correcting item 2
  if ("Item2" %in% colnames(scaledScoreDF)) {
    scaledScoreDF$Item2 <- scaledScoreDF$Item2 * 0.1
  }
  
  #---- Return
  scaledScoreDF
}








#' Rank normalization methods from scaled normScore items
#'
#' Computes the total normScore for each normalization method by summing the
#' scaled item scores and ranks methods from best to worst. The log-transformed
#' reference method is additionally adjusted using the correction factor from
#' raw total intensities.
#'
#' @param scaledScoreDF A data frame containing scaled normScore item values,
#'   with normalization methods in rows and items in columns.
#' @param item0 Numeric correction factor applied to the `"Log"` method.
#'
#' @return A data frame containing the scaled item scores, total score,
#'   corrected total score, and rows ordered by increasing corrected score.
#'
#' @keywords internal
#' @noRd


rankNormScoreItems <- function(
    scaledScoreDF,
    item0
){
  #----- Final ranking 
  scaledScoreDF$Total <- rowSums(scaledScoreDF, na.rm = TRUE)
  scaledScoreDF$TotalxItem0 <- scaledScoreDF$Total
  
  scaledScoreDF[rownames(scaledScoreDF) == "Log", "TotalxItem0"] <-
    scaledScoreDF[rownames(scaledScoreDF) == "Log", "TotalxItem0"] * item0
  
  rankedScoreDF <- scaledScoreDF[order(scaledScoreDF$TotalxItem0), ]
  
  rankedScoreDF
}








#' Bootstrap Estimation of Normalization Scores
#'
#' @description
#' Computes bootstrap-based mean scores and percentile confidence intervals
#' for normalization methods from a matrix of item-wise scores.
#'
#' The function applies bootstrap resampling over rows of the input score 
#' matrix, where each row typically represents a scoring component (item) 
#' and each column represents a normalization method. For each resample, 
#' total scores per normalization are computed using 
#' \code{\link{bootstrapRowScores}}.
#'
#' @param rankedScoreDF `numeric` matrix where rows correspond to scoring items
#'   and columns correspond to normalization methods.
#' @param nBoot `integer`. Number of bootstrap resamples. Default is `1000`.
#' @param item0 Numeric correction factor applied to the `"Log"` method.
#'
#'
#' @seealso
#' \code{\link{bootstrapRowScores}}
#'
#' @keywords internal
#' @noRd

computeBootstrapNormScore <- function(rankedScoreDF, nBoot, item0){
  # Chaging format
  scoreMatrix <- t(rankedScoreDF)
  
  # Correcting Log values directly into items
  if ("Log" %in% colnames(scoreMatrix)) {
    scoreMatrix[, "Log"] <- scoreMatrix[, "Log"] * item0
  }
  
  # performing booting
  bootResults <- boot::boot(
    data = scoreMatrix,
    statistic = bootstrapRowScores,
    R = nBoot)
  bootstrapMeans <- colMeans(bootResults$t)
  
  bootstrapScores <- as.data.frame(
    t(vapply(seq_len(ncol(scoreMatrix)), function(i) {
        ci <- boot::boot.ci(bootResults, type = "perc", index = i)
        c(colnames(scoreMatrix)[i],
          bootstrapMeans[i],
          ci$percent[4],
          ci$percent[5])},
        character(4))))
  
  colnames(bootstrapScores) <- c("normalization", "meanNormScore", 
                                 "ll95", "ul95")
  bootstrapScores[-1] <- lapply(bootstrapScores[-1], as.numeric)
  bootstrapScores <- bootstrapScores[order(bootstrapScores$meanNormScore), ]
  
  return(bootstrapScores)
}


#' @importFrom rlang .data
NULL



