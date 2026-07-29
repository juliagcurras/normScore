#' @noRd
.validateGroupData <- function(
    groupData,
    refGroup,
    altGroup
) {
    groupData <- as.data.frame(groupData)

    if (ncol(groupData) != 2L) {
        stop("`groupData` must contain exactly two columns.")
    }

    colnames(groupData) <- c("Samples", "Groups")
    groupLevels <- unique(as.character(groupData$Groups))

    if (length(groupLevels) == 0L) {
        stop("No valid groups were found in `groupData`.")
    }

    singleGroup <- length(groupLevels) == 1L

    if (!singleGroup) {
        groupNames <- .selectComparisonGroups(
            groupLevels,
            refGroup,
            altGroup
        )
        refGroup <- groupNames$refGroup
        altGroup <- groupNames$altGroup
    }

    list(
        groupData = groupData,
        refGroup = refGroup,
        altGroup = altGroup,
        singleGroup = singleGroup
    )
}

#' @noRd
.selectComparisonGroups <- function(
    groupLevels,
    refGroup,
    altGroup
) {
    if (is.null(refGroup) || !refGroup %in% groupLevels) {
        refGroup <- groupLevels[length(groupLevels)]
        message("Reference group set to ", refGroup, ".")
    }

    validAltGroups <- groupLevels[groupLevels != refGroup]

    if (is.null(altGroup) || !altGroup %in% validAltGroups) {
        altGroup <- validAltGroups[1]
        message("Alternative group set to ", altGroup, ".")
    }

    list(
        refGroup = refGroup,
        altGroup = altGroup
    )
}

#' @noRd
.validateRawData <- function(rawData) {
    rawData <- tryCatch(
        as.data.frame(rawData),
        error = function(e) {
            stop("`rawData` could not be converted to a data frame.")
        }
    )

    if (nrow(rawData) < 100L) {
        stop(
            "A minimum of 100 proteins is required for ",
            "normalization assessment."
        )
    }

    if (is.null(rownames(rawData)) ||
        is.null(colnames(rawData))) {
        stop("`rawData` must have row and column names.")
    }

    rawData
}

#' @noRd
.coerceNormalizedData <- function(normalizedDataList) {
    if (!is.list(normalizedDataList) ||
        length(normalizedDataList) == 0L) {
        stop("`normalizedDataList` must be a non-empty list.")
    }

    if (is.null(names(normalizedDataList)) ||
        any(names(normalizedDataList) == "")) {
        stop("All elements of `normalizedDataList` must be named.")
    }

    tryCatch(
        lapply(normalizedDataList, as.data.frame),
        error = function(e) {
            stop(
                "Some elements of `normalizedDataList` could not be ",
                "converted to data frames."
            )
        }
    )
}

#' @noRd
.validateDimnames <- function(
    normalizedDataList,
    rawData,
    groupData
) {
    if (!setequal(colnames(rawData), groupData$Samples)) {
        stop(
            "Sample names do not match between `rawData` ",
            "and `groupData`."
        )
    }

    validColumns <- vapply(
        normalizedDataList,
        function(x) {
            identical(
                sort(colnames(x)),
                sort(colnames(rawData))
            )
        },
        logical(1)
    )

    if (!all(validColumns)) {
        stop(
            "Sample names do not match between `rawData` and ",
            "the elements of `normalizedDataList`."
        )
    }

    validRows <- vapply(
        normalizedDataList,
        function(x) {
            identical(
                sort(rownames(x)),
                sort(rownames(rawData))
            )
        },
        logical(1)
    )

    if (!all(validRows)) {
        stop(
            "Protein names do not match between `rawData` and ",
            "the elements of `normalizedDataList`."
        )
    }

    invisible(TRUE)
}


#' @noRd
.alignNormalizedData <- function(
    normalizedDataList,
    rawData
) {
    refRows <- rownames(rawData)
    refCols <- colnames(rawData)

    lapply(
        normalizedDataList,
        function(x) {
            as.matrix(
                x[refRows, refCols, drop = FALSE]
            )
        }
    )
}


#' @noRd
.validateMainArguments <- function(
    returnDetails,
    nBoot,
    fromMainFunction
) {
    if (!fromMainFunction) {
        return(
            list(
                returnDetails = returnDetails,
                nBoot = nBoot
            )
        )
    }

    if (!is.logical(returnDetails) ||
        length(returnDetails) != 1L) {
        returnDetails <- FALSE
        warning(
            "`returnDetails` must be a single logical value. ",
            "Setting it to FALSE."
        )
    }

    if (!is.numeric(nBoot) || length(nBoot) != 1L) {
        nBoot <- 500
        warning(
            "`nBoot` must be a single numeric value. ",
            "Setting it to 500."
        )
    }

    if (nBoot < 100) {
        nBoot <- 100
        warning(
            "The minimum value allowed for `nBoot` is 100. ",
            "Setting it to 100."
        )
    }

    list(
        returnDetails = returnDetails,
        nBoot = nBoot
    )
}


#' @noRd
.validateMainArguments <- function(
    returnDetails,
    nBoot,
    fromMainFunction
) {
    if (!fromMainFunction) {
        return(list(
            returnDetails = returnDetails,
            nBoot = nBoot
        ))
    }

    if (!is.logical(returnDetails) ||
        length(returnDetails) != 1L ||
        is.na(returnDetails)) {
        returnDetails <- FALSE
        warning(
            "`returnDetails` must be TRUE or FALSE. ",
            "Setting it to FALSE."
        )
    }

    validNBoot <- is.numeric(nBoot) &&
        length(nBoot) == 1L &&
        is.finite(nBoot)

    if (!validNBoot) {
        nBoot <- 500L
        warning(
            "`nBoot` must be a finite numeric value. ",
            "Setting it to 500."
        )
    }

    nBoot <- as.integer(nBoot)

    if (nBoot < 100L) {
        nBoot <- 100L
        warning(
            "The minimum value allowed for `nBoot` is 100. ",
            "Setting it to 100."
        )
    }

    list(
        returnDetails = returnDetails,
        nBoot = nBoot
    )
}
