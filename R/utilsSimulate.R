#' @noRd
.simulateDataCore <- function(args) {
    groups <- rep(c("G1", "G2"), each = args$nPerGroup)
    nSamples <- length(groups)

    abundance <- .simulateProteinAbundance(args)
    deMatrix <- .simulateDifferentialExpression(
        args,
        abundance$lowWeight,
        groups
    )
    factorMatrix <- .simulateLatentStructure(
        args,
        groups
    )
    residualMatrix <- .simulateResidualNoise(
        args,
        abundance$lowWeight,
        nSamples
    )

    logData <- matrix(
        abundance$mu,
        nrow = args$nProteins,
        ncol = nSamples
    )
    logData <- logData + factorMatrix + deMatrix + residualMatrix
    logData <- .applyGlobalSampleShift(logData, args)
    logData <- .applySampleSdEffect(logData, args)
    logData <- .setSimulationDimnames(logData, groups)

    if (args$addMissing) {
        logData <- .addMnarMissingness(logData, args)
    }

    .buildSimulationOutput(logData, groups)
}


#' @noRd
.simulateProteinAbundance <- function(args) {
    mu <- stats::rnorm(
        args$nProteins,
        mean = args$muMean,
        sd = args$muSd
    )
    mu <- pmin(
        pmax(mu, args$muClip[1]),
        args$muClip[2]
    )

    denominator <- args$muClip[2] - args$muClip[1]
    lowWeight <- (args$muClip[2] - mu) / denominator
    lowWeight <- pmin(pmax(lowWeight, 0), 1)

    list(
        mu = mu,
        lowWeight = lowWeight
    )
}


#' @noRd
.simulateDifferentialExpression <- function(
    args,
    lowWeight,
    groups
) {
    nDE <- round(args$nProteins * args$propDE)
    deIndex <- integer(0)

    if (nDE > 0) {
        deIndex <- sample.int(args$nProteins, nDE)
    }

    logFC <- numeric(args$nProteins)

    if (nDE > 0) {
        logFCSd <- .computeLogFCSd(
            args,
            lowWeight,
            deIndex
        )
        logFCValues <- stats::rnorm(
            nDE,
            mean = args$logFCMean,
            sd = logFCSd
        )
        direction <- sample(
            c(-1, 1),
            nDE,
            replace = TRUE
        )
        logFC[deIndex] <- logFCValues * direction
    }

    contrast <- ifelse(groups == "G1", 0.5, -0.5)
    outer(logFC, contrast)
}


#' @noRd
.computeLogFCSd <- function(args, lowWeight, deIndex) {
    if (!args$heteroLogFC) {
        return(rep(args$logFCSd, length(deIndex)))
    }

    multiplier <- args$fcHi +
        lowWeight[deIndex]^args$gammaFC *
            (args$fcLo - args$fcHi)

    args$logFCSd * multiplier
}


#' @noRd
.simulateLatentStructure <- function(args, groups) {
    correlationMatrix <- .buildCorrelationMatrix(
        groups,
        args$rhoWithin,
        args$rhoBetween
    )
    correlationMatrix <- .stabilizeCorrelationMatrix(
        correlationMatrix
    )

    nSamples <- length(groups)
    latentFactor <- MASS::mvrnorm(
        1,
        mu = rep(0, nSamples),
        Sigma = correlationMatrix
    )
    latentFactor <- as.numeric(latentFactor)
    latentFactor <- .centerWithinGroups(
        latentFactor,
        groups
    )

    loadings <- stats::rnorm(
        args$nProteins,
        mean = 0,
        sd = args$loadingSd
    )

    outer(loadings, latentFactor)
}


#' @noRd
.buildCorrelationMatrix <- function(
    groups,
    rhoWithin,
    rhoBetween
) {
    nSamples <- length(groups)
    correlationMatrix <- matrix(
        rhoBetween,
        nrow = nSamples,
        ncol = nSamples
    )

    for (group in unique(groups)) {
        index <- which(groups == group)
        correlationMatrix[index, index] <- rhoWithin
    }

    diag(correlationMatrix) <- 1
    correlationMatrix
}


#' @noRd
.stabilizeCorrelationMatrix <- function(correlationMatrix) {
    eigenValues <- eigen(
        correlationMatrix,
        symmetric = TRUE,
        only.values = TRUE
    )$values

    minEigenValue <- min(eigenValues)

    if (minEigenValue <= 1e-8) {
        adjustment <- abs(minEigenValue) + 1e-6
        diag(correlationMatrix) <-
            diag(correlationMatrix) + adjustment
    }

    correlationMatrix
}


#' @noRd
.centerWithinGroups <- function(values, groups) {
    centeredValues <- values

    for (group in unique(groups)) {
        index <- which(groups == group)
        centeredValues[index] <-
            centeredValues[index] - mean(centeredValues[index])
    }

    centeredValues
}


#' @noRd
.simulateResidualNoise <- function(
    args,
    lowWeight,
    nSamples
) {
    sigmaVector <- args$sigmaHi +
        lowWeight^args$gammaSigma *
            (args$sigmaLo - args$sigmaHi)

    residualMatrix <- matrix(
        stats::rnorm(args$nProteins * nSamples),
        nrow = args$nProteins,
        ncol = nSamples
    )

    residualMatrix * sigmaVector
}


#' @noRd
.applyGlobalSampleShift <- function(logData, args) {
    if (args$sampleShiftSd <= 0) {
        return(logData)
    }

    nSamples <- ncol(logData)
    sampleShift <- stats::rnorm(
        nSamples,
        mean = 0,
        sd = args$sampleShiftSd
    )
    sampleShift <- sampleShift - mean(sampleShift)

    if (.isFiniteScalar(args$sampleShiftCap)) {
        sampleShift <- pmin(
            pmax(sampleShift, -args$sampleShiftCap),
            args$sampleShiftCap
        )
        sampleShift <- sampleShift - mean(sampleShift)
    }

    sweep(logData, 2, sampleShift, "+")
}


#' @noRd
.applySampleSdEffect <- function(logData, args) {
    if (args$sampleSdStrength == 0) {
        return(logData)
    }

    sampleMeans <- colMeans(logData, na.rm = TRUE)
    zMean <- .safeScale(
        rank(sampleMeans, ties.method = "average")
    )
    randomNoise <- .safeScale(
        stats::rnorm(ncol(logData))
    )

    rho <- max(0, min(1, args$sampleSdRho))
    zSd <- rho * zMean +
        sqrt(1 - rho^2) * randomNoise

    logMultiplier <- args$sampleSdStrength * zSd

    if (.isFiniteScalar(args$sampleSdCap)) {
        logMultiplier <- pmin(
            pmax(logMultiplier, -args$sampleSdCap),
            args$sampleSdCap
        )
    }

    sdScale <- exp(logMultiplier)
    centeredData <- sweep(logData, 2, sampleMeans, "-")
    scaledData <- sweep(centeredData, 2, sdScale, "*")
    sweep(scaledData, 2, sampleMeans, "+")
}


#' @noRd
.safeScale <- function(x) {
    scaledValues <- as.numeric(scale(x))

    if (anyNA(scaledValues)) {
        return(rep(0, length(x)))
    }

    scaledValues
}


#' @noRd
.isFiniteScalar <- function(x) {
    !is.null(x) &&
        length(x) == 1L &&
        is.finite(x)
}


#' @noRd
.setSimulationDimnames <- function(logData, groups) {
    rownames(logData) <- paste0(
        "P",
        sprintf("%05d", seq_len(nrow(logData)))
    )

    sampleNumber <- stats::ave(
        seq_along(groups),
        groups,
        FUN = seq_along
    )
    colnames(logData) <- paste0(
        groups,
        "_",
        sampleNumber
    )

    logData
}


#' @noRd
.addMnarMissingness <- function(logData, args) {
    nProteins <- nrow(logData)
    nSamples <- ncol(logData)

    sampleShift <- stats::rnorm(
        nSamples,
        mean = 0,
        sd = args$missingBySampleSd
    )
    sampleShift <- sampleShift - mean(sampleShift)

    shiftMatrix <- matrix(
        sampleShift,
        nrow = nProteins,
        ncol = nSamples,
        byrow = TRUE
    )

    objective <- function(intercept) {
        probability <- stats::plogis(
            intercept - args$kMnar * logData +
                shiftMatrix
        )
        mean(probability, na.rm = TRUE) -
            args$targetMissing
    }

    intercept <- stats::uniroot(
        objective,
        interval = c(-50, 50)
    )$root

    probability <- stats::plogis(
        intercept - args$kMnar * logData +
            shiftMatrix
    )
    randomValues <- matrix(
        stats::runif(nProteins * nSamples),
        nrow = nProteins,
        ncol = nSamples
    )

    logData[randomValues < probability] <- NA_real_
    logData
}


#' @noRd
.buildSimulationOutput <- function(logData, groups) {
    metadata <- data.frame(
        Samples = colnames(logData),
        Groups = factor(
            groups,
            levels = c("G1", "G2")
        )
    )

    list(
        logData = logData,
        rawData = 2^logData,
        metadata = metadata
    )
}
