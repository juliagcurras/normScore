
#' Validate and prepare inputs for normScore
#'
#' Checks the input objects used by `normScore()`, including group annotation,
#' raw data, normalized datasets, sample names, protein names, and bootstrap
#' settings. The function also aligns the order of rows and columns across
#' normalized datasets, adds the log-transformed dataset when required, and
#' identifies whether the analysis contains a single group.
#'
#' @param normalizedDataList A named list of normalized data matrices or data
#'   frames. Each element should contain proteins in rows and samples in columns.
#' @param groupData A data frame or object coercible to a data frame containing
#'   sample-group annotation. The first column is assumed to contain sample names
#'   and the second column group labels.
#' @param rawData A matrix or data frame containing raw intensity values, with
#'   proteins in rows and samples in columns.
#' @param refGroup Character string indicating the reference group used for
#'   group-comparison metrics. If `NULL` or invalid, it is automatically set.
#' @param altGroup Character string indicating the alternative group used for
#'   group-comparison metrics. If `NULL` or invalid, it is automatically set.
#' @param returnDetails Logical value indicating whether detailed results should
#'   be returned by `normScore()`.
#' @param nBoot Numeric value indicating the number of bootstrap resamples.
#'
#' @return A list containing validated and prepared inputs: `normalizedDataList`,
#'   `groupData`, `rawData`, `refGroup`, `altGroup`, `returnDetails`, `nBoot`,
#'   and `singleGroup`.
#'
#' @keywords internal

validateNormScoreInput <- function(
    normalizedDataList,
    groupData,
    rawData,
    refGroup,
    altGroup,
    returnDetails,
    nBoot
    
){
  # Check 1. Group data
  groupData <- as.data.frame(groupData)
  colnames(groupData) <- c("Samples", "Groups")
  groupLevels <- levels(as.factor(groupData$Groups)) # length group levels = 1?
  singleGroup <- FALSE
  
  if (length(groupLevels) > 1){
    # Check 2. Group names
    if (is.null(refGroup) || !(refGroup %in% groupLevels)) {
      refGroup <- groupLevels[length(groupLevels)]
      message(paste0("Reference group set to ", refGroup, "."))
    }
    if (is.null(altGroup) || !(altGroup %in% groupLevels)) {
      altGroup <- groupLevels[!(groupLevels == refGroup)][1]
      message(paste0("Alternative group set to ", altGroup, "."))
    }
  } else if (length(groupLevels) == 1){
    singleGroup <- TRUE 
  } else {
    stop("Error with number of groups")
  }
  
  
  # Check 3. Number of proteins
  if(nrow(rawData) < 100){
    stop("A minimum of 100 proteins is required for normalization assessment")
  }
  
  # Check 4. Some data types
  if (!is.data.frame(rawData)){
    rawData <- as.data.frame(rawData)
  }
  isMatrix <- all(sapply(normalizedDataList, is.matrix))
  isDF <- all(sapply(normalizedDataList, is.data.frame))
  if (!all(isMatrix, isDF)){
    res <- try(lapply(normalizedDataList, as.data.frame))
    if (inherits(x = res, "try-error")){
      stop("Error related to object types in some normalizedDataList elements")
    } else {
      normalizedDataList <- res
    }
  }
  if (!is.logical(returnDetails)){
    returnDetails <- F
    warning("Wrong value for returnDetails. Set to FALSE.")
  }
  
  # Check 5. Sample names
  if (!all(colnames(rawData) %in% groupData$Samples)){
    stop("All or some sample names from rawData do not match sample names from 
         groupData.")
  }
  if (length(unique(sapply(normalizedDataList, ncol))) != 1){
    stop("Different number of samples (columns) across normalized datasets from 
         normalizedDataList.")
  }
  if (ncol(rawData) != unique(sapply(normalizedDataList, ncol))){
    stop("Different number of columns between rawData and normalized datasets
          from normalizedDataList.")
  }
  uniqueColnames <- unique(do.call(c, lapply(normalizedDataList, colnames)))
  if (length(uniqueColnames) != unique(sapply(normalizedDataList, ncol))){
    stop("Sample names do not match across normalized datasets
          from normalizedDataList.")
  }
  colnamesMatch <- all(uniqueColnames %in% colnames(rawData))
  colnamesNumber <- length(uniqueColnames) == ncol(rawData)
  if (!all(colnamesMatch, colnamesNumber)){
    stop("Sample names do not match between rawData and normalized datasets from 
         normalizedDataList.")
  }
  
  # Check 4. Protein number and order
  if (length(unique(sapply(normalizedDataList, nrow))) != 1){
    stop("Different number of proteins (rows) across normalized datasets from 
         normalizedDataList.")
  }
  if (nrow(rawData) != unique(sapply(normalizedDataList, nrow))){
    stop("Different number of rows (proteins) between rawData and normalized datasets
          from normalizedDataList.")
  }
  uniqueRownames <- unique(do.call(c, lapply(normalizedDataList, rownames)))
  if (length(uniqueRownames) != unique(sapply(normalizedDataList, nrow))){
    stop("Protein names do not match across normalized datasets
          from normalizedDataList.")
  }
  rownamesMatch <- all(uniqueRownames %in% rownames(rawData))
  rownamesNumber <- length(uniqueRownames) == nrow(rawData)
  if (!all(rownamesMatch, rownamesNumber)){
    stop("Protein names do not match between rawData and normalized datasets from 
           normalizedDataList.")
  }
  
  
  # Setting same sample and protein order. 
  refRows <- rownames(rawData)
  refCols <- colnames(rawData)
  
  normalizedDataList <- lapply(normalizedDataList, function(x) {
    as.matrix(x[refRows, refCols, drop = FALSE])
  })
  
  # Check 7. Log dataset
  normalizedDataList <- addLogIfMissing(
    normalizedDataList = normalizedDataList,
    rawData = rawData
  )
  
  # Check 8. Checking nBoot type
  if (!is.numeric(nBoot)){
    nBoot <- 500
    warning("nBoot input is not numeric. Setting nBoot to 500")
  }
  if (nBoot<100){
    nBoot <- 100
    warning("The minimum value allowed for nBoot is 100. Setting nBoot to 100.")
  }
  
  
  # Returning!
  list(
    normalizedDataList = normalizedDataList,
    groupData = groupData,
    rawData = rawData,
    refGroup = refGroup,
    altGroup = altGroup, 
    returnDetails = returnDetails,
    nBoot = nBoot, 
    singleGroup = singleGroup
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
#' @param groupData A data frame containing sample-group annotation, with sample
#'   names in the first column and group labels in the second column.
#' @param rawData A matrix or data frame containing raw intensity values. Used
#'   for consistency with the main workflow.
#' @param refGroup Character string indicating the reference group used for the
#'   MA-plot item.
#' @param altGroup Character string indicating the alternative group used for the
#'   MA-plot item.
#' @param singleGroup Logical value indicating whether the input data contain a
#'   single group.
#'
#' @return A data frame with normalization methods in rows and normScore items
#'   in columns. If `singleGroup = TRUE`, the MA-plot item is omitted.
#'
#' @keywords internal

computeNormScoreItems <- function(
    normalizedDataList,
    groupData,
    rawData,
    refGroup,
    altGroup, 
    singleGroup
){
  
  #----- ITEM 1 - PCV
  dfPCV <- do.call(
    rbind, 
      sapply(
      normalizedDataList,
      getPCV,
      groups = levels(as.factor(groupData$Groups)),
      groupData = groupData, 
      simplify = FALSE, 
      USE.NAMES =  TRUE
    )
  )
  
  item1 <- apply(dfPCV, 1, mean, na.rm = TRUE)
  
 
  #----- ITEM 2 - Correlation (Spearman)
  item2 <- computeCorrelation(
    normalizedDataList = normalizedDataList,
    groupData = groupData, 
    method ="spearman"
  ) 
  
  
  #----- ITEM 3 - MA plot regression line vs expected 0
  if (!singleGroup){
    item3 <- sapply(
      normalizedDataList,
      maDiffArea,
      samplesG1 =  groupData[groupData$Groups == refGroup, "Samples"],
      samplesG2 = groupData[groupData$Groups == altGroup, "Samples"]
    )
  } else {
    item3 <- stats::setNames(
      rep(NA, length(normalizedDataList)), 
      names(normalizedDataList)
      )
  }
  
  
  #----- ITEM 4 - Mean-SD regression line slope deviation
  item4 <- sapply(normalizedDataList, meanSDdiffArea)
  
  
  #----- ITEM 5 - RLE: MAPE of sample medians (reference = 1)
  item5 <- sapply(normalizedDataList, rleMAPE)
  
  
  #----- ITEM 6 - Total intensity consistency
  item6 <- sapply(normalizedDataList, tiMAPE)
  
  
  #---- Score matrix
  refNames <- names(item1)
  scoreDF <- cbind(
    Item1 = item1[refNames],
    Item2 = item2[refNames],
    Item3 = item3[refNames],
    Item4 = item4[refNames],
    Item5 = item5[refNames],
    Item6 = item6[refNames]
  )
  scoreDF <- as.data.frame(scoreDF)
  rownames(scoreDF) <- refNames
  
  # Selecting only items from scaled matrix
  if (singleGroup){
    itemsToSelect <- paste0("Item", c(1:2, 4:6))
  } else {
    itemsToSelect <- paste0("Item", 1:6)
  }
  scoreDF <- scoreDF[, itemsToSelect]
  
  # Return
}





#' Scale normScore item scores
#'
#' Applies min-max scaling to each normScore item so that item values are placed
#' on a comparable scale across normalization methods. The correlation item is
#' additionally down-weighted after scaling.
#'
#' @param scoreDF A data frame containing raw normScore item values, with
#'   normalization methods in rows and items in columns.
#'
#' @return A data frame with the same dimensions as `scoreDF`, containing scaled
#'   item scores. Row names are preserved.
#'
#' @keywords internal

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
#' The function applies bootstrap resampling over rows of the input score matrix,
#' where each row typically represents a scoring component (item) and each column
#' represents a normalization method. For each resample, total scores per
#' normalization are computed using \code{\link{bootstrapRowScores}}.
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

computeBootstrapNormScore <- function(
    rankedScoreDF, 
    nBoot,
    item0
  ) 
  {
  
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






#' Plot bootstrap normScore results
#'
#' Generates a point-range plot showing the mean bootstrap normScore and its
#' 95% confidence interval for each normalization method.
#'
#' @param x A `normScore` result object containing a `bootstrapScore` element.
#'   This element should be a data frame with normalization names, mean
#'   normScore values, and lower and upper 95% confidence limits. It can be 
#'   obtained by setting returnDetails = `TRUE` at `normScore` main function.
#'
#' @return A `ggplot` object showing bootstrap mean normScore values and 95%
#'   confidence intervals for each normalization method.
#'
#' @export
#'
#' @examples
#' bootstrapScore <- data.frame(
#'   normalization = c("Norm1", "Norm2", "Norm3"),
#'   meanNormScore = c(0.8, 1.2, 1.5),
#'   ll95 = c(0.6, 1.0, 1.2),
#'   ul95 = c(1.0, 1.4, 1.8)
#' )
#'
#' result <- list(bootstrapScore = bootstrapScore)
#'
#' plotBootstrapNormScore(result)
#' 
#' @export

plotBootstrapNormScore <- function(x) {
  
  if (is.null(x$bootstrapScore)) {
    stop(
      "Bootstrap scores are not available. ",
      "Run normScore(..., returnDetails = TRUE) to compute them."
    )
  }
  
  bootstrapScores <- x$bootstrapScore
  
  # Checking colnames
  namesCols <- c(
    "normalization",
    "meanNormScore",
    "ll95",
    "ul95"
  )
  
  if (!identical(colnames(bootstrapScores), namesCols)) {
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



