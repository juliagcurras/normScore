#' Plot functions for original criteria
#'
#' This file contains plot functions used as outputs for interpreting the
#' detailed ranking results.


#' Plot raw total intensity across samples
#'
#' Generates the diagnostic plot associated with the systematic bias penalty
#' item (item 0). Raw total intensity is computed for each sample and displayed
#' as a bar plot.
#'
#' @param data A raw intensity matrix or data frame with proteins in rows and
#'   samples in columns.
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
#' @noRd


plotItem0 <- function(
    data
) {
    # Estimating sum of intensities
    dataSum <- as.data.frame(colSums(data, na.rm = TRUE))
    colnames(dataSum) <- "Intensity"
    dataSum$Samples <- factor(
        rownames(dataSum),
        levels = rownames(dataSum),
        ordered = TRUE
    )

    # Aesthetics
    dataSum$Color <- normScorePalette(nrow(dataSum))

    # Plot
    ggplot2::ggplot(
        dataSum,
        ggplot2::aes(
            x = .data$Samples,
            y = .data$Intensity,
            fill = .data$Color
        )
    ) +
        ggplot2::geom_bar(
            stat = "identity"
            # fill = "#1786A3"
        ) +
        ggplot2::scale_fill_manual(values = dataSum$Color) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(legend.position = "none",
        axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"),
        axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"),
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)
        ) +
        ggplot2::ylab("Total intensity") +
        ggplot2::xlab("Samples") +
        ggplot2::labs(
            title = "Systematic Bias Penalty - Item 0 "
        )
}


#' Building item 1: forest plot (<5 groups)
#' 
#' @keywords internal
#' @noRd

plotItem1Forest  <- function(dfForest, groups){
    ggplot2::ggplot(dfForest, ggplot2::aes(
        x = .data$Normalization,
        y = .data$Mean, colour = .data$Group, shape = .data$Group
    )) +
    ggplot2::geom_errorbar(
        ggplot2::aes(
        ymin = .data$Lower, ymax = .data$Upper
    ),
        width = 0, linewidth = 0.7,
        position = ggplot2::position_dodge(width = 0.45)
    ) +
    ggplot2::geom_point(
        size = 2.5, na.rm = TRUE,
        position = ggplot2::position_dodge(width = 0.45)
    ) +
    ggplot2::labs(
        x = "Normalization Methods", y = "Mean of PCV (%)",
        colour = "Groups", shape = "Groups",
        title = "PCV: pooled coefficient of variation - Item 1"
    ) +
    ggplot2::scale_colour_manual(
        values = normScorePalette(length(groups))) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
        axis.line = ggplot2::element_line(linewidth = 0.5,colour = "black"),
        axis.ticks =ggplot2::element_line(linewidth = 0.5,colour = "black"),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "right"
    )
}

#' Building item 1: boxplots (>5 groups)
#' 
#' @keywords internal
#' @noRd

plotItem1Boxplot <- function(dfForest, groups){
    ggplot2::ggplot(dfForest, ggplot2::aes(
        x = .data$Normalization, y = .data$Mean
    )) +
    ggplot2::geom_boxplot(
        width = 0.65, outlier.shape = NA, na.rm = TRUE,
        alpha = 0.7, colour = "darkgrey"
    ) +
    ggplot2::geom_jitter(
        ggplot2::aes(shape = .data$Group, colour = .data$Group),
        width = 0.12, height = 0, na.rm = TRUE, size = 3, alpha = 1
    ) +
    ggplot2::scale_colour_manual(
        values = normScorePalette(length(groups))) +
    ggplot2::labs(
        x = "Normalization Methods", y = "Mean of PCV (%)",
        colour = "Groups", shape = "Groups",
        title = "PCV: pooled coefficient of variation - Item 1"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
        axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"),
        axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"),
        legend.position = "right", 
        panel.grid.minor = ggplot2::element_blank())
}


#' Plot pooled coefficient of variation
#'
#' Generates the diagnostic plot associated with the PCV item. PCV values are
#' computed for each normalization method and group, and displayed as
#' point-range or boxplot-based summaries depending on the number of groups.
#'
#' @param normalizedDataList A named list of normalized data matrices.
#' @param groupData A data frame with sample names and group labels.
#'
#' @return A `ggplot` or arranged ggplot object.
#'
#' @keywords internal
#' @noRd


plotItem1 <- function(normalizedDataList, groupData) {
    groups <- levels(as.factor(groupData$Groups))
    dfPCV <- do.call(rbind, lapply(normalizedDataList, getPCV,
        groups = groups,
        groupData = groupData, plotData = TRUE
    ))
    dfPCV <- as.data.frame(t(dfPCV))
    normNames <- colnames(dfPCV)
    groupNames <- paste0("Group ", seq_along(groups))
    dfForest <- do.call(rbind, lapply(normNames, function(normalizationName) {
        statistics <- matrix(dfPCV[[normalizationName]],
            ncol = 3, byrow = TRUE,
            dimnames = list(NULL, c("Mean", "Lower", "Upper"))
        )
        data.frame(
            Group = groupNames, Normalization = normalizationName,
            statistics, row.names = NULL, check.names = FALSE
        )
    }))
    
    # plot!
    dfForest$Normalization <- factor(dfForest$Normalization, levels = normNames)
    if (length(groups) < 5) {
        plotItem1Forest(dfForest = dfForest, groups = groups)
    } else if (length(groups) >= 5) {
        plotItem1Boxplot(dfForest = dfForest, groups = groups)
    }
}


#' Plot within-group correlations
#'
#' Generates the diagnostic plot associated with the correlation item. The plot
#' displays the distribution of within-group sample correlations for each
#' normalization method.
#'
#' @param normalizedDataList A named list of normalized data matrices.
#' @param groupData A data frame with sample names and group labels.
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
#' @noRd

plotItem2 <- function(normalizedDataList, groupData) {
    # Computing data
    dataToPlot <- computeCorrelation(
        normalizedDataList = normalizedDataList,
        groupData = groupData,
        method = "spearman",
        plotData = TRUE
    )

    # Arrange format
    dfPlot <- utils::stack(dataToPlot)
    colnames(dfPlot) <- c("Correlation", "Normalization")
    dfPlot$Normalization <- factor(dfPlot$Normalization,
        levels = colnames(dataToPlot)
    )

    ggplot2::ggplot(dfPlot, ggplot2::aes(
        x = .data$Normalization, y = .data$Correlation
    )) +
        ggplot2::geom_boxplot(
            width = 0.65, na.rm = TRUE, outlier.shape = NA,
            alpha = 0.7, colour = "black"
        ) +
        ggplot2::geom_jitter(ggplot2::aes(colour = .data$Normalization),
            width = 0.12, na.rm = TRUE, height = 0,
            size = 1, alpha = 0.3
        ) +
        ggplot2::scale_colour_manual(
            values = normScorePalette(length(levels(dfPlot$Normalization)))
        ) +
        ggplot2::labs(
            x = "Normalization Methods", y = "Intragroup correlation",
            colour = "Normalization", title = "Intragroup correlation - Item 2"
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(
        axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"),
        axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"),
        legend.position = "none"
        )
}


#' Plot MA diagnostic plots
#'
#' Generates MA plots for each normalization method using the selected
#' reference and alternative groups. These plots are used to assess
#'   intensity-dependent differences between groups.
#'
#' @param normalizedDataList A named list of normalized data matrices.
#' @param groupData A data frame with sample names and group labels.
#' @param refGroup Character string indicating the reference group.
#' @param altGroup Character string indicating the alternative group.
#'
#' @return An arranged ggplot object containing one MA plot per normalization
#'   method.
#'
#' @keywords internal
#' @noRd


plotItem3 <- function(normalizedDataList, groupData, refGroup, altGroup) {
    dataToPlot <- lapply( # Computing data
        normalizedDataList,
        maDiffArea,
        samplesG1 = groupData[groupData$Groups == refGroup, "Samples"],
        samplesG2 = groupData[groupData$Groups == altGroup, "Samples"],
        plotData = TRUE)

    minY <- floor(min(vapply(dataToPlot, function(df) {
        min(df$M, na.rm = TRUE)
    }, numeric(1)))) * 1.1
    maxY <- ceiling(max(vapply(X = dataToPlot, FUN = function(df) {
        max(df$M, na.rm = TRUE)
    }, FUN.VALUE = numeric(1)))) * 1.1
    colors <- normScorePalette(length(dataToPlot))
    names(colors) <- names(dataToPlot)

    allMAplots <- lapply(
        X = names(dataToPlot),
        FUN = function(nameNorm) {
            dfPlot <- dataToPlot[[nameNorm]]
            ggplot2::ggplot(dfPlot, ggplot2::aes(
                x = .data$A, y = .data$M)) +
                ggplot2::geom_point(
                    alpha = 0.85, na.rm = TRUE, size = 1.5,
                    fill = colors[[nameNorm]], colour = colors[[nameNorm]]) +
                ggplot2::geom_hline(
                    yintercept = 0, linetype = "dashed",
                    linewidth = 0.8, colour = "black") +
                ggplot2::ylim(c(minY, maxY)) +
                ggplot2::labs(
                    x = "A: average intensity",
                    y = paste0(
                        "M: log of mean ", altGroup,
                        " over mean ", refGroup),
                    title = nameNorm) +
                ggplot2::theme_minimal(base_size = 14) +
                ggplot2::theme(
                    axis.line = ggplot2::element_line(
                        linewidth = 0.5, colour = "black"),
                    axis.ticks = ggplot2::element_line(
                        linewidth = 0.5, colour = "black"))
        }
    )
    return(ggpubr::ggarrange(plotlist = allMAplots))
}


#' Plot mean-SD diagnostic plots
#'
#' Generates diagnostic plots showing the relationship between sample order by
#' mean intensity and standard deviation for each normalization method.
#'
#' @param normalizedDataList A named list of normalized data matrices.
#'
#' @return An arranged ggplot object containing one plot per normalization
#'   method.
#'
#' @keywords internal
#' @noRd


plotItem4 <- function(normalizedDataList) {
    dataToPlot <- lapply( # Computing data
        X = normalizedDataList,
        FUN = meanSDdiffArea,
        plotData = TRUE)

    minY <- floor(min(vapply(dataToPlot, function(df) {
        min(df$SD, na.rm = TRUE)
    }, numeric(1)))) * 1.1
    maxY <- ceiling(max(vapply(X = dataToPlot, FUN = function(df) {
        max(df$SD, na.rm = TRUE)
    }, FUN.VALUE = numeric(1)))) * 1.1

    colors <- normScorePalette(length(dataToPlot))
    names(colors) <- names(dataToPlot)

    allMeaSDplots <- lapply(X = names(dataToPlot), FUN = function(nameNorm) {
        dfPlot <- dataToPlot[[nameNorm]]
        ggplot2::ggplot(dfPlot, ggplot2::aes(
            x = .data$Order, y = .data$SD
        )) +
            ggplot2::geom_point(
                size = 2.5, alpha = 0.8, na.rm = TRUE,
                fill = colors[[nameNorm]], colour = colors[[nameNorm]]
            ) +
            ggplot2::scale_x_continuous(
                breaks = dfPlot$Order,
                labels = dfPlot$Sample
            ) +
            ggplot2::ylim(c(minY, maxY)) +
            ggplot2::labs(
                y = "Standard deviation (SD)",
                x = "Samples order by mean magnitude",
                title = nameNorm
            ) +
            ggplot2::theme_minimal(base_size = 14) +
            ggplot2::theme(
                axis.line = ggplot2::element_line(
                    linewidth = 0.5, colour = "black"
                ),
                axis.ticks = ggplot2::element_line(
                    linewidth = 0.5, colour = "black"
                ),
                panel.grid.minor = ggplot2::element_blank()
            )
    })

    return(ggpubr::ggarrange(plotlist = allMeaSDplots))
}


#' Plot RLE distributions
#'
#' Generates RLE distribution plots for each normalization method. These plots
#' are used to inspect whether sample-wise RLE distributions are centered and
#' comparable after normalization.
#'
#' @param normalizedDataList A named list of normalized data matrices.
#'
#' @return An arranged ggplot object containing one RLE plot per normalization
#'   method.
#'
#' @keywords internal
#' @noRd


plotItem5 <- function(normalizedDataList) {
    dataToPlot <- lapply( # Computing data
        X = normalizedDataList,
        FUN = rleMAPE,
        plotData = TRUE)

    # Common parameters for plotting
    minY <- floor(min(vapply(
        dataToPlot, function(df) min(df, na.rm = TRUE),
        numeric(1)
    )))
    maxY <- ceiling(max(vapply(
        dataToPlot, function(df) max(df, na.rm = TRUE),
        numeric(1)
    )))
    colors <- normScorePalette(length(dataToPlot))
    names(colors) <- names(dataToPlot)

    allRLEplots <- lapply(
        X = names(dataToPlot),
        FUN = function(nameNorm) {
            data <- as.data.frame(dataToPlot[[nameNorm]])
            dfPlot <- utils::stack(data)
            colnames(dfPlot) <- c("RLE", "Samples")
            dfPlot$Samples <- factor(dfPlot$Samples, levels = colnames(data))

            ggplot2::ggplot(dfPlot, ggplot2::aes(
                x = .data$Samples, y = .data$RLE
            )) +
                ggplot2::geom_boxplot(
                    width = 0.65, alpha = 0.7,
                    na.rm = TRUE, colour = colors[[nameNorm]]) +
                ggplot2::labs(x = "Samples", y = "RLE", title = nameNorm) +
                ggplot2::ylim(c(minY, maxY)) +
                ggplot2::theme_minimal(base_size = 14) +
                ggplot2::theme(
                    axis.line = ggplot2::element_line(
                        linewidth = 0.5, colour = "black"),
                    axis.ticks = ggplot2::element_line(
                        linewidth = 0.5, colour = "black"),
                    legend.position = "none")
        }
    )
    return(ggpubr::ggarrange(plotlist = allRLEplots))
}


#' Plot log-intensity distributions
#'
#' Generates sample-wise intensity distribution plots for each normalization
#' method. These plots are used to visually inspect global distributional
#' consistency across samples.
#'
#' @param normalizedDataList A named list of normalized data matrices.
#'
#' @return An arranged ggplot object containing one intensity distribution plot
#'   per normalization method.
#'
#' @keywords internal
#' @noRd


plotItem6 <- function(normalizedDataList) {
    # Computing data
    dataToPlot <- normalizedDataList

    # Common parameters for plotting
    minY <- floor(min(vapply(
        dataToPlot,
        function(df) min(df, na.rm = TRUE), numeric(1)
    )))
    maxY <- ceiling(max(vapply(
        dataToPlot,
        function(df) max(df, na.rm = TRUE), numeric(1)
    )))

    colors <- normScorePalette(length(dataToPlot))
    names(colors) <- names(dataToPlot)

    # Plot!
    allPlots <- lapply(
        X = names(dataToPlot),
        FUN = function(nameNorm) {
            data <- as.data.frame(dataToPlot[[nameNorm]])
            dfPlot <- utils::stack(data)
            colnames(dfPlot) <- c("Intensity", "Samples")
            dfPlot$Samples <- factor(dfPlot$Samples, levels = colnames(data))

            ggplot2::ggplot(dfPlot, ggplot2::aes(
                x = .data$Samples, y = .data$Intensity
            )) +
                ggplot2::geom_boxplot(
                    width = 0.65, alpha = 0.5, na.rm = TRUE,
                    fill = colors[[nameNorm]], colour = colors[[nameNorm]]
                ) +
                ggplot2::labs(x = "Samples", y = "Log Intensity", 
                title = nameNorm) +
                ggplot2::ylim(c(minY, maxY)) +
                ggplot2::theme_minimal(base_size = 14) +
                ggplot2::theme(
                    axis.line = ggplot2::element_line(
                        linewidth = 0.5, colour = "black"),
                    axis.ticks = ggplot2::element_line(
                        linewidth = 0.5, colour = "black"),
                    legend.position = "none")
        }
    )

    return(ggpubr::ggarrange(plotlist = allPlots))
}


#' @importFrom rlang .data
NULL
