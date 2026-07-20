
#' Plot bootstrap normScore results
#'
#' Generates a point-range plot showing the mean bootstrap normScore and its
#' 95% confidence interval for each normalization method.
#'
#' @param x A `normScore` result object containing a `bootstrapScore` 
#' element. This element should be a data frame with normalization names, 
#' mean normScore values, and lower and upper 95% confidence limits. 
#' It can beobtained by setting returnDetails = `TRUE` 
#' at `normScore` main function.
#'
#' @return A `ggplot` object showing bootstrap mean normScore values 
#' and 95% confidence intervals for each normalization method.
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
    stop("Bootstrap scores are not available. ",
      "Run normScore(..., returnDetails = TRUE) to compute them.")
  }
  bootstrapScores <- x$bootstrapScore
  
  # Checking colnames
  namesCols <- c("normalization", "meanNormScore", "ll95", "ul95")
  
  if (!identical(colnames(bootstrapScores), namesCols)) {
    colnames(bootstrapScores) <- namesCols
  }
  bootstrapScores$colors <- normScorePalette(nrow(bootstrapScores))
  
  # Plot!
  p <- ggplot2::ggplot(
    bootstrapScores,
    ggplot2::aes(
      x = .data$meanNormScore,
      y = stats::reorder(.data$normalization, .data$meanNormScore),
      colour = .data$normalization)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        xmin = .data$ll95,
        xmax = .data$ul95),
      width = 0.2,
      orientation = "y"
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      x = "Mean normScore [95% CI]",
      y = "Normalization"
    ) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(
        bootstrapScores$colors,
        bootstrapScores$normalization)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"),
      legend.position = "none"
    )
  
  return(p)
}




#' Plot normScore diagnostic items
#'
#' Generates the set of diagnostic plots associated with the individual
#' normScore evaluation items. These plots provide a visual interpretation of
#' the criteria used to assess normalization performance, including raw total
#' intensity, pooled coefficient of variation, within-group correlation, MA
#' plots, mean-SD trends, RLE distributions, and intensity distributions.
#'
#' The function validates and aligns the input data before plotting. If only
#' one group is provided, the MA-plot diagnostic is skipped because it requires
#' a comparison between two groups.
#'
#' @param normalizedDataList A named list of normalized data matrices or data
#'   frames. Each element should contain proteins in rows and samples in 
#'   columns.
#' @param groupData A data frame containing sample-group annotation. The first
#'   column is assumed to contain sample names and the second column group
#'   labels.
#' @param rawData A matrix or data frame containing raw intensity values, with
#'   proteins in rows and samples in columns.
#' @param refGroup Character string indicating the reference group used for the
#'   MA-plot diagnostic. If `NULL`, it is automatically selected.
#' @param altGroup Character string indicating the alternative group used for 
#'   the MA-plot diagnostic. If `NULL`, it is automatically selected.
#'
#' @return A named list of diagnostic plots. Each element corresponds to one
#'   normScore item. If the input contains a single group, the MA-plot element 
#'   is returned as `NULL`.
#'   
#' @examples
#' # Simulate proteomic data
#' simData <- simulateData(nProteins = 1000)
#' 
#' # Normalyze ysing NormalyzerDE package
#' normalizedDataList <- list(
#'   Norm1 = simData$logData + 0.1,
#'   Norm2 = simData$logData + 1,
#'   Norm3 = simData$logData - 1,
#'   Norm4 = simData$logData * 1.1,
#'   Norm5 = simData$logData * 0.9,
#'   Norm6 = simData$logData * runif(ncol(simData$logData), 0.8, 1.2)
#' )
#' 
#' # Compute ranking
#' plotNormScoreDiagnostics(
#'   normalizedDataList = normalizedDataList,
#'   groupData = simData$metadata,
#'   rawData = simData$rawData,
#'   refGroup = NULL,
#'   altGroup = NULL
#'  ) 
#'
#' @export

plotNormScoreDiagnostics <- function(
    normalizedDataList,
    groupData,
    rawData,
    refGroup = NULL,
    altGroup = NULL
) {
  inputs <- validateNormScoreInput(
    normalizedDataList = normalizedDataList,
    groupData = groupData,
    rawData = rawData,
    refGroup = refGroup,
    altGroup = altGroup,
    fromMainFunction = FALSE)
  
  finalPlots <- list() # Just to save
  
  finalPlots[["item0"]] <- plotItem0(data = inputs$rawData)
  
  finalPlots[["item1"]] <- plotItem1(
    normalizedDataList = inputs$normalizedDataList, 
    groupData = inputs$groupData)
  
  finalPlots[["item2"]] <- plotItem2(
    normalizedDataList = inputs$normalizedDataList, 
    groupData = inputs$groupData)
  
  if (!inputs$singleGroup){
    p3 <- plotItem3(
      normalizedDataList = inputs$normalizedDataList, 
      groupData = inputs$groupData, 
      refGroup = inputs$refGroup,
      altGroup = inputs$altGroup) 
  } else {
    p3 <- NULL
  }
  finalPlots[["item3"]] <- p3
  
  p4 <- plotItem4(normalizedDataList = inputs$normalizedDataList)
  finalPlots[["item4"]] <- p4
  
  p5 <- plotItem5(normalizedDataList = inputs$normalizedDataList)
  finalPlots[["item5"]] <- p5
  
  p6 <- plotItem6(normalizedDataList = inputs$normalizedDataList)
  finalPlots[["item6"]] <- p6
  
  finalPlots
}


#' @importFrom rlang .data
NULL
