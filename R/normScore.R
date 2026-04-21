#' Rank Normalization Methods Using a Composite Scoring Framework
#'
#' @description
#' Computes a composite score to rank normalization methods based on multiple
#' quality metrics derived from normalized data matrices.
#'
#' The function evaluates each normalization method using six criteria:
#' pooled coefficient of variation (PCV), within-group sample correlation,
#' MA-plot trend deviation, mean-SD trend deviation, relative log expression
#' (RLE) consistency, and total intensity consistency. These metrics are scaled,
#' combined into a final score, and used to rank the normalization methods.
#'
#' In addition, a correction factor derived from the coefficient of variation
#' of raw sample intensities is applied to the `"Log"` normalization method.
#' Optionally, bootstrap resampling is used to estimate mean composite scores
#' and percentile-based confidence intervals.
#'
#' @param normalizedDataList Named list of normalized numeric matrices or data
#'   frames. Each element represents one normalization method, with rows
#'   corresponding to features (e.g. proteins) and columns to samples.
#' @param groupData `data.frame` or matrix describing the sample-group
#'   assignment. It must contain two columns, which are interpreted as sample
#'   names and group labels.
#' @param rawData `data.frame` or matrix of raw numeric intensities, with rows
#'   representing features and columns representing samples.
#' @param refGroup `character` or `NULL`. Reference group used for the MA-plot
#'   metric. If `NULL`, the last group in `groupData` is used.
#' @param altGroup `character` or `NULL`. Alternative group used for the MA-plot
#'   metric. If `NULL`, the first group in `groupData` is used.
#' @param onlyFinalRanking `logical`. If `TRUE`, only the final ranking is
#'   returned. If `FALSE`, detailed scores, bootstrap results, and a summary
#'   plot are also returned. Default is `FALSE`.
#' @param nBoot `integer`. Number of bootstrap resamples used to estimate mean
#'   normScore values and confidence intervals. Default is `1000`.
#'
#' @return
#' If `onlyFinalRanking = TRUE`, returns a list with:
#' \describe{
#'   \item{finalRanking}{Named numeric vector with the final normalization ranking.}
#' }
#'
#' If `onlyFinalRanking = FALSE`, returns a list with:
#' \describe{
#'   \item{finalRanking}{Named numeric vector with the final normalization ranking.}
#'   \item{detailRanking}{`data.frame` with scaled item-wise scores, total score,
#'   and corrected total score for each normalization method.}
#'   \item{bootstrapScore}{`data.frame` with bootstrap mean scores and 95\%
#'   percentile confidence intervals for each normalization method.}
#'   \item{graphic}{A `ggplot2` object showing the bootstrap mean scores and
#'   confidence intervals.}
#' }
#'
#' @details
#' The function computes the following six item scores for each normalization
#' method:
#' \enumerate{
#'   \item \strong{PCV}: mean pooled coefficient of variation across groups.
#'   \item \strong{Correlation}: within-group Spearman correlation summary,
#'   transformed so that lower values indicate better performance.
#'   \item \strong{MAplot}: shape-corrected area-based deviation from the
#'   expected MA trend (`logFC = 0`).
#'   \item \strong{MeanSDplot}: area-based deviation of the mean-SD trend from
#'   horizontality.
#'   \item \strong{RLEplot}: RLE-based mean absolute percentage error relative
#'   to 1.
#'   \item \strong{totalIntensity}: quantile-based total intensity consistency
#'   metric.
#' }
#'
#' Each item is scaled to the range 0 to 1 using min-max scaling and then summed
#' into a total score. Lower scores indicate better normalization performance.
#'
#' A correction factor based on the coefficient of variation of total raw sample
#' intensities is applied only to the normalization method named `"Log"`.
#'
#' In addition, the correlation item is downweighted for all methods except
#' `"CyclicLoess"` to account for its lower discrimination power in this scoring
#' framework.
#'
#' When `onlyFinalRanking = FALSE`, bootstrap resampling of the six item scores
#' is performed using \code{\link{computeBootstrapNormScore}}, and the results
#' are visualized with \code{\link{plotBootstrapNormScore}}.
#'
#' @examples
#' # Minimal example structure
#' 
#' # Simulate proteomic data
#' simData <- simulateData(nProteins = 1000)
#' 
#' # Normalyze ysing NormalyzerDE package
#' normalizedDataList <- list(
#'   Median = NormalyzerDE::medianNormalization(simData$rawData),
#'   Mean = NormalyzerDE::meanNormalization(simData$rawData),
#'   TI = NormalyzerDE::globalIntensityNormalization(simData$rawData),
#'   Quantile =  NormalyzerDE::performQuantileNormalization(simData$rawData),
#'   CyclicLoess = NormalyzerDE::performCyclicLoessNormalization(simData$rawData),
#'   VSN =  NormalyzerDE::performVSNNormalization(simData$rawData),
#'   RLR =  NormalyzerDE::performGlobalRLRNormalization(simData$rawData)
#' )
#' 
#' # Compute ranking
#' normScore(
#'   normalizedDataList,
#'   groupData = simData$metadata,
#'   rawData = simData$rawData,
#'   refGroup = NULL,
#'   altGroup = NULL,
#'   onlyFinalRanking = FALSE,
#'   nBoot = 100
#'  ) 
#' 
#'
#' @seealso
#' \code{\link{addLogIfMissing}},
#' \code{\link{getPCV}},
#' \code{\link{withinGroupCorrelations}},
#' \code{\link{maDiffArea}},
#' \code{\link{meanSDdiffArea}},
#' \code{\link{rleMAPE}},
#' \code{\link{tiMAPE}},
#' \code{\link{computeBootstrapNormScore}},
#' \code{\link{plotBootstrapNormScore}}
#'
#' @export


normScore <- function(
    normalizedDataList,
    groupData,
    rawData,
    refGroup = NULL,
    altGroup = NULL,
    onlyFinalRanking = TRUE,
    nBoot = 1000
) {
  # Input:
  # 1. List of normalized matrices (normalizedDataList)
  # 2. Sample-group annotation table (groupData)
  # 3. Raw intensities (rawData) for item 0
  # 4. Reference and alternative groups for MA metric (optional)
  # 5. Whether to return only the final ranking
  # 6. Number of bootstrap resamples
  
  #----- Initial checks ####
  groupData <- as.data.frame(groupData)
  colnames(groupData) <- c("Samples", "Groups")
  
  groupLevels <- levels(as.factor(groupData$Groups))
  scoreList <- list()
  
  if(nrow(rawData) < 100){
    stop("A minimum of 100 proteins is required for normalization assessment")
  }
  
  normalizedDataList <- addLogIfMissing(
    normalizedDataList = normalizedDataList,
    rawData = rawData
  )
  
  
  #----- ITEM 0 - correction factor ####
  totalIntensities <- colSums(rawData, na.rm = TRUE)
  item0 <- cv(totalIntensities, proportion = TRUE, na.rm = TRUE) * 3
  
  #----- ITEM 1 - PCV ####
  dfPCV <- data.frame(
    lapply(
      normalizedDataList,
      getPCV,
      groups = groupLevels,
      groupData = groupData
    )
  )
  
  dfPCV <- as.data.frame(t(dfPCV))
  dfPCV$PCV <- apply(dfPCV, 1, mean, na.rm = TRUE)
  item1 <- stats::setNames(dfPCV$PCV, rownames(dfPCV))
  
  scoreList[["PCV"]] <- item1
  
  #----- ITEM 2 - Correlation (Spearman) ####
  allCorrelationVectors <- lapply(
    normalizedDataList,
    withinGroupCorrelations,
    groupData = groupData,
    method = "spearman"
  )
  
  dfCor <- data.frame(
    sapply(
      allCorrelationVectors,
      "length<-",
      max(lengths(allCorrelationVectors))
    )
  )
  
  item2 <- sapply(
    colnames(dfCor),
    function(normalizationName) {
      correlationValues <- dfCor[, normalizationName]
      
      # Use 1 - statistic so that lower values consistently indicate better performance
      1 - (stats::median(correlationValues, na.rm = TRUE) - stats::IQR(correlationValues, na.rm = TRUE) / 3)
    },
    simplify = TRUE,
    USE.NAMES = TRUE
  )
  
  scoreList[["Correlation"]] <- item2
  
  #----- ITEM 3 - MA plot regression line vs expected 0 ####
  refGroup <- ifelse(is.null(refGroup), groupLevels[length(groupLevels)], refGroup)
  altGroup <- ifelse(is.null(altGroup), groupLevels[1], altGroup)
  
  samplesG1 <- groupData[groupData$Groups == refGroup, "Samples"]
  samplesG2 <- groupData[groupData$Groups == altGroup, "Samples"]
  
  item3 <- sapply(
    normalizedDataList,
    maDiffArea,
    samplesG1 = samplesG1,
    samplesG2 = samplesG2
  )
  
  scoreList[["MAplot"]] <- item3
  
  #----- ITEM 4 - Mean-SD regression line slope deviation ####
  item4 <- sapply(normalizedDataList, meanSDdiffArea)
  scoreList[["MeanSDplot"]] <- item4
  
  #----- ITEM 5 - RLE: MAPE of sample medians (reference = 1) ####
  item5 <- sapply(normalizedDataList, rleMAPE)
  scoreList[["RLEplot"]] <- item5
  
  #----- ITEM 6 - Total intensity consistency ####
  item6 <- sapply(normalizedDataList, tiMAPE)
  scoreList[["totalIntensity"]] <- item6
  
  #----- Join all item scores ####
  scoreDF <- dplyr::bind_cols(scoreList)
  scoreDF <- as.data.frame(scoreDF)
  rownames(scoreDF) <- names(scoreList[[1]])
  
  #----- Min-max scaling ####
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
  
  #----- Corrections ####
  # Low discrimination power for item 2 (Correlation)
  # scaledScoreDF[rownames(scaledScoreDF) != "CyclicLoess", 2] <-
  #   scaledScoreDF[rownames(scaledScoreDF) != "CyclicLoess", 2] * 0.1
  scaledScoreDF[, 2] <- scaledScoreDF[, 2] * 0.1
  
  #----- Final ranking ####
  scaledScoreDF$Total <- rowSums(scaledScoreDF)
  scaledScoreDF$TotalxItem0 <- scaledScoreDF$Total
  
  scaledScoreDF[rownames(scaledScoreDF) == "Log", "TotalxItem0"] <-
    scaledScoreDF[rownames(scaledScoreDF) == "Log", "TotalxItem0"] * item0
  
  scaledScoreDF <- scaledScoreDF[order(scaledScoreDF$TotalxItem0), ]
  
  finalRanking <- stats::setNames(scaledScoreDF$TotalxItem0, rownames(scaledScoreDF))
  names(scaledScoreDF) <- c(paste0("Item", 1:6), "Total", "TotalxItem0")
  
  if (onlyFinalRanking) {
    return(
      list(
        finalRanking = finalRanking
      )
    )
  } else {
    #----- Bootstrap over item scores ####
    scoreMatrix <- t(scaledScoreDF[, 1:6])
    
    if ("Log" %in% colnames(scoreMatrix)) {
      scoreMatrix[, "Log"] <- scoreMatrix[, "Log"] * item0
    }
    
    bootstrapScores <- computeBootstrapNormScore(
      scoreMatrix = scoreMatrix,
      nBoot = nBoot
    )
    
    bootstrapPlot <- plotBootstrapNormScore(
      bootstrapScores = bootstrapScores
    )
    
    return(
      list(
        finalRanking = finalRanking,
        detailRanking = scaledScoreDF,
        bootstrapScore = bootstrapScores,
        graphic = bootstrapPlot
      )
    )
  }
}

