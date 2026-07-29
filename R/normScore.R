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
#' @param refGroup `character` or `NULL`. Reference group used for
#'   the MA-plot metric. If `NULL`, the last group in `groupData` is used.
#' @param altGroup `character` or `NULL`. Alternative group used for the
#'   metric. If `NULL`, the first group in `groupData` is used.
#' @param returnDetails `logical`. If `TRUE`, only the final ranking is
#'   MA-plot returned. If `FALSE`, detailed scores, bootstrap results, and
#'   a summaryplot are also returned. Default is `FALSE`.
#' @param nBoot `integer`. Number of bootstrap resamples used to estimate mean
#'   normScore values and confidence intervals. Default is `1000`.
#'
#' @return
#' If `returnDetails = TRUE`, returns a list with:
#' \describe{
#'   \item{finalRanking}{Named numeric vector with the final
#'   normalization ranking.}
#' }
#'
#' If `returnDetails = FALSE`, returns a list with:
#' \describe{
#'   \item{finalRanking}{Named numeric vector with the final
#'   normalization ranking.}
#'   \item{detailRanking}{`data.frame` with scaled item-wise scores, total
#'   score, and corrected total score for each normalization method.}
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
#'   \item \strong{Item1}: mean pooled coefficient of variation across
#'   groups (PVC).
#'   \item \strong{Item2}: within-group Spearman correlation summary,
#'   transformed so that lower values indicate better performance (Correlation).
#'   \item \strong{Item3}: shape-corrected area-based deviation from the
#'   expected MA trend (`logFC = 0`), based on MAplot.
#'   \item \strong{Item4}: area-based deviation of the mean-SD trend from
#'   horizontality, based on meanSD plot.
#'   \item \strong{Item5}: RLE-based mean absolute percentage error relative
#'   to 1.
#'   \item \strong{Item6}: quantile-based total intensity consistency
#'   metric.
#' }
#'
#' Each item is scaled to the range 0 to 1 using min-max scaling and then summed
#' into a total score. Lower scores indicate better normalization performance.
#'
#' A correction factor based on the coefficient of variation of total raw sample
#' intensities is applied only to the normalization method named `"Log"`.
#'
#' When `returnDetails = FALSE`, bootstrap resampling of the six item scores
#' is performed.The results can be visualized using
#' \code{\link{plotBootstrapNormScore}}.
#'
#' @examples
#'
#' # Simulate proteomic data
#' simData <- simulateData(nProteins = 1000)
#'
#' # Normalyze ysing NormalyzerDE package
#' normalizedDataList <- list(
#'     Norm1 = simData$logData + 0.1,
#'     Norm2 = simData$logData + 1,
#'     Norm3 = simData$logData - 1,
#'     Norm4 = simData$logData * 1.1,
#'     Norm5 = simData$logData * 0.9,
#'     Norm6 = simData$logData * runif(ncol(simData$logData), 0.8, 1.2)
#' )
#'
#' # Compute ranking
#' result <- normScore(
#'     normalizedDataList,
#'     groupData = simData$metadata,
#'     rawData = simData$rawData,
#'     refGroup = NULL,
#'     altGroup = NULL,
#'     returnDetails = TRUE,
#'     nBoot = 100
#' )
#' result
#'
#' @seealso
#' \code{\link{plotNormScoreDiagnostics}},
#' \code{\link{plotBootstrapNormScore}},
#' \code{\link{simulateData}}
#'
#' @export

normScore <- function(normalizedDataList, groupData, rawData,
    refGroup = NULL, altGroup = NULL,
    returnDetails = TRUE, nBoot = 1000) {
    #----- Initial checks ####
    inputs <- validateNormScoreInput(
        normalizedDataList = normalizedDataList,
        groupData = groupData,
        rawData = rawData,
        refGroup = refGroup,
        altGroup = altGroup,
        returnDetails = returnDetails,
        nBoot = nBoot)

    item0 <- cv( #----- Item 0 as correction factor ####
        x = colSums(rawData, na.rm = TRUE),
        proportion = TRUE, na.rm = TRUE) * 4 # Correction factor

    #----- Compute items ####
    scoreDF <- computeNormScoreItems(
        normalizedDataList = inputs$normalizedDataList,
        groupData = inputs$groupData,
        rawData = inputs$rawData,
        refGroup = inputs$refGroup,
        altGroup = inputs$altGroup,
        singleGroup = inputs$singleGroup)

    #----- Scaled items ####
    scaledScoreDF <- scaleNormScoreItems(scoreDF)

    #----- Aggregate items, correct total and rank ####
    rankedScoreDF <- rankNormScoreItems(scaledScoreDF, item0 = item0)
    finalRanking <- stats::setNames(
        rankedScoreDF$TotalxItem0,
        rownames(rankedScoreDF))

    #----- Return only finalRanking
    if (!inputs$returnDetails) {
        return(list(finalRanking = finalRanking))
    } else { #----- Last output: Bootstrap over score items'
        bootstrapScores <- computeBootstrapNormScore(
            rankedScoreDF = rankedScoreDF, item0 = item0, nBoot = inputs$nBoot
        )
        #----- Return detailed results
        return(list(
            finalRanking = finalRanking, detailRanking = rankedScoreDF,
            bootstrapScore = bootstrapScores))
    }
}
