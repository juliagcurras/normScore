#' Simulate Label-Free Proteomics-Like Two-Group Data
#'
#' @description
#' Simulates two-group quantitative omics data with controlled signal structure,
#' including protein-specific mean abundance, within- and between-group sample
#' correlation, abundance-dependent residual variance, differential expression,
#' sample-level global shifts, mean-SD dependence, and optional MNAR missingness.
#'
#' The function is designed to generate synthetic datasets suitable for testing
#' normalization and scoring workflows without including large example datasets
#' in the package.
#'
#' @param nProteins `integer`. Number of proteins or features to simulate.
#'   Default is `10000`.
#' @param nPerGroup `integer`. Number of samples per group. Default is `20`.
#' @param muMean `numeric`. Mean of the protein-wise abundance distribution on
#'   the log2 scale. Default is `18`.
#' @param muSd `numeric`. Standard deviation of the protein-wise abundance
#'   distribution on the log2 scale. Default is `1.2`.
#' @param muClip `numeric` vector of length 2. Lower and upper bounds used to
#'   truncate simulated protein means. Default is `c(15, 25)`.
#' @param rhoWithin `numeric`. Within-group correlation target for the latent
#'   sample factor. Default is `0.85`.
#' @param rhoBetween `numeric`. Between-group correlation target for the latent
#'   sample factor. Default is `0.55`.
#' @param loadingSd `numeric`. Standard deviation of protein-specific loadings
#'   for the correlated latent factor. Default is `0.25`.
#' @param sigmaHi `numeric`. Residual standard deviation at high abundance.
#'   Default is `0.05`.
#' @param sigmaLo `numeric`. Residual standard deviation at low abundance.
#'   Default is `0.4`.
#' @param gammaSigma `numeric`. Controls how strongly residual variance depends
#'   on abundance. Default is `2.5`.
#' @param propDE `numeric`. Proportion of differentially expressed proteins.
#'   Default is `0.35`.
#' @param logFCSd `numeric`. Standard deviation of differential expression
#'   log-fold changes. Default is `1`.
#' @param logFCMean `numeric`. Mean of differential expression log-fold changes.
#'   Default is `0`.
#' @param heteroLogFC `logical`. If `TRUE`, the variance of log-fold changes
#'   depends on abundance. Default is `TRUE`.
#' @param fcHi `numeric`. Lower multiplier for abundance-dependent logFC
#'   heterogeneity. Default is `0.55`.
#' @param fcLo `numeric`. Upper multiplier for abundance-dependent logFC
#'   heterogeneity. Default is `2.5`.
#' @param gammaFC `numeric`. Controls how strongly logFC variability depends on
#'   abundance. Default is `7`.
#' @param sampleShiftSd `numeric`. Standard deviation of global sample shifts.
#'   Default is `0`.
#' @param sampleShiftCap `numeric`. Maximum absolute value allowed for global
#'   sample shifts. Default is `0.2`.
#' @param sampleSdStrength `numeric`. Strength of sample-specific mean-SD
#'   dependence. A value of `0` implies independence. Default is `0`.
#' @param sampleSdRho `numeric`. Correlation between sample mean rank and
#'   sample-specific SD effect. Must lie between `0` and `1`. Default is `0.8`.
#' @param sampleSdCap `numeric`. Maximum absolute value allowed for the
#'   log-multiplier controlling sample-specific SD effects. Default is `0.35`.
#' @param addMissing `logical`. Should MNAR missing values be added?
#'   Default is `TRUE`.
#' @param targetMissing `numeric`. Target overall proportion of missing values.
#'   Default is `0.001`.
#' @param kMnar `numeric`. Strength of abundance dependence in the MNAR missing
#'   value mechanism. Default is `1.2`.
#' @param missingBySampleSd `numeric`. Standard deviation of sample-specific
#'   missingness shifts. Default is `0.05`.
#' @param seed `integer`. Random seed used for reproducibility. Default is `9396`.
#'
#' @return
#' A list with the following elements:
#' \describe{
#'   \item{logData}{Numeric matrix of simulated log2-scale data.}
#'   \item{rawData}{Numeric matrix of simulated raw-scale data, obtained as `2^logData`.}
#'   \item{metadata}{`data.frame` containing sample names and group labels.}
#' }
#'
#' @details
#' The simulated data include several structured components:
#' \enumerate{
#'   \item Protein-wise mean abundance on the log2 scale.
#'   \item A latent correlation structure inducing stronger within-group than
#'   between-group sample correlation.
#'   \item Abundance-dependent heteroscedastic residual noise to create an
#'   MA-like wedge pattern.
#'   \item Symmetric differential expression between the two groups.
#'   \item Optional global sample shifts affecting overall intensity.
#'   \item Optional sample-level mean-SD dependence.
#'   \item Optional MNAR missingness driven by abundance and sample-specific effects.
#' }
#'
#' The output is intended for testing normalization methods and associated
#' scoring procedures under controlled simulation settings.
#'
#' @examples
#' simData <- simulateData(
#'   nProteins = 1000,
#'   nPerGroup = 5,
#'   propDE = 0.2,
#'   addMissing = TRUE,
#'   seed = 123
#' )
#' 
#'
#' dim(simData$logData)
#' head(simData$metadata)
#'
#' @export

simulateData <- function(
    nProteins = 10000,
    nPerGroup = 20,
    
    # Protein means (log2 scale)
    muMean = 18,
    muSd = 1.2,
    muClip = c(15, 25),
    
    # Correlation structure (within-group > between-group)
    rhoWithin = 0.85,
    rhoBetween = 0.55,
    loadingSd = 0.25,
    
    # MA wedge: residual variance depends on abundance
    sigmaHi = 0.05,
    sigmaLo = 0.4,
    gammaSigma = 2.5,
    
    # Differential expression (symmetric between groups)
    propDE = 0.35,
    logFCSd = 1,
    logFCMean = 0,
    heteroLogFC = TRUE,
    fcHi = 0.55,
    fcLo = 2.5,
    gammaFC = 7,
    
    # Global sample shift
    sampleShiftSd = 0,
    sampleShiftCap = 0.2,
    
    # Mean-SD dependence by sample
    sampleSdStrength = 0,
    sampleSdRho = 0.8,
    sampleSdCap = 0.35,
    
    # Missing values
    addMissing = TRUE,
    targetMissing = 0.001,
    kMnar = 1.2,
    missingBySampleSd = 0.05,
    
    seed = 9396
) {
  set.seed(seed)
  
  #----- Group labels and sample names ####
  groups <- rep(c("G1", "G2"), each = nPerGroup)
  nSamples <- length(groups)
  
  #----- Protein mean abundances ####
  mu <- stats::rnorm(nProteins, muMean, muSd)
  mu <- pmin(pmax(mu, muClip[1]), muClip[2])
  
  lowAbundanceWeight <- (muClip[2] - mu) / (muClip[2] - muClip[1])
  lowAbundanceWeight <- pmin(pmax(lowAbundanceWeight, 0), 1)
  
  #----- Symmetric differential expression (+/- logFC / 2) ####
  nDE <- round(nProteins * propDE)
  deIndex <- if (nDE > 0) sample.int(nProteins, nDE) else integer(0)
  
  logFC <- rep(0, nProteins)
  
  if (nDE > 0) {
    logFCSdVector <- rep(logFCSd, nDE)
    
    if (heteroLogFC) {
      fcMultiplier <- fcHi + (lowAbundanceWeight[deIndex]^gammaFC) * (fcLo - fcHi)
      logFCSdVector <- logFCSd * fcMultiplier
    }
    
    logFC[deIndex] <- stats::rnorm(nDE, logFCMean, logFCSdVector) * sample(c(-1, 1), nDE, TRUE)
  }
  
  groupContrast <- ifelse(groups == "G1", 0.5, -0.5)
  deMatrix <- outer(logFC, groupContrast)
  
  #----- Correlated latent factor for sample correlation ####
  correlationMatrix <- matrix(rhoBetween, nSamples, nSamples)
  diag(correlationMatrix) <- 1
  
  group1Index <- which(groups == "G1")
  group2Index <- which(groups == "G2")
  
  correlationMatrix[group1Index, group1Index] <- rhoWithin
  diag(correlationMatrix[group1Index, group1Index]) <- 1
  
  correlationMatrix[group2Index, group2Index] <- rhoWithin
  diag(correlationMatrix[group2Index, group2Index]) <- 1
  
  eigenValues <- eigen(correlationMatrix, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenValues) <= 1e-8) {
    correlationMatrix <- correlationMatrix + diag(abs(min(eigenValues)) + 1e-6, nSamples)
  }
  
  sampleLatentFactor <- as.numeric(
    MASS::mvrnorm(1, mu = rep(0, nSamples), Sigma = correlationMatrix)
  )
  
  sampleLatentFactor[group1Index] <- sampleLatentFactor[group1Index] - mean(sampleLatentFactor[group1Index])
  sampleLatentFactor[group2Index] <- sampleLatentFactor[group2Index] - mean(sampleLatentFactor[group2Index])
  
  proteinLoadings <- stats::rnorm(nProteins, 0, loadingSd)
  factorMatrix <- outer(proteinLoadings, sampleLatentFactor)
  
  #----- Heteroscedastic residual noise (MA wedge) ####
  sigmaVector <- sigmaHi + (lowAbundanceWeight^gammaSigma) * (sigmaLo - sigmaHi)
  residualMatrix <- matrix(stats::rnorm(nProteins * nSamples), nrow = nProteins, ncol = nSamples) * sigmaVector
  
  #----- Global sample shifts ####
  sampleShift <- rep(0, nSamples)
  
  if (sampleShiftSd > 0) {
    sampleShift <- stats::rnorm(nSamples, mean = 0, sd = sampleShiftSd)
    sampleShift <- sampleShift - mean(sampleShift)
    
    if (!is.null(sampleShiftCap) && is.finite(sampleShiftCap)) {
      sampleShift <- pmin(pmax(sampleShift, -sampleShiftCap), sampleShiftCap)
      sampleShift <- sampleShift - mean(sampleShift)
    }
  }
  
  #----- Base log2 matrix ####
  logData <- matrix(mu, nrow = nProteins, ncol = nSamples) + factorMatrix + deMatrix + residualMatrix
  
  if (any(sampleShift != 0)) {
    logData <- sweep(logData, 2, sampleShift, "+")
  }
  
  #----- Mean-SD dependence by sample ####
  if (sampleSdStrength != 0) {
    sampleMeans <- colMeans(logData, na.rm = TRUE)
    zMean <- as.numeric(scale(rank(sampleMeans, ties.method = "average")))
    if (anyNA(zMean)) {
      zMean <- rep(0, nSamples)
    }
    
    randomNoise <- as.numeric(scale(stats::rnorm(nSamples)))
    if (anyNA(randomNoise)) {
      randomNoise <- rep(0, nSamples)
    }
    
    rho <- max(0, min(1, sampleSdRho))
    zSd <- rho * zMean + sqrt(1 - rho^2) * randomNoise
    
    logMultiplier <- sampleSdStrength * zSd
    if (!is.null(sampleSdCap) && is.finite(sampleSdCap)) {
      logMultiplier <- pmin(pmax(logMultiplier, -sampleSdCap), sampleSdCap)
    }
    
    sdScale <- exp(logMultiplier)
    centeredLogData <- sweep(logData, 2, sampleMeans, "-")
    logData <- sweep(centeredLogData, 2, sdScale, "*")
    logData <- sweep(logData, 2, sampleMeans, "+")
  }
  
  #----- Row and column names ####
  rownames(logData) <- paste0("P", sprintf("%05d", seq_len(nProteins)))
  colnames(logData) <- paste0(groups, "_", stats::ave(seq_along(groups), groups, FUN = seq_along))
  
  #----- Optional MNAR missingness ####
  missingInfo <- NULL
  
  if (addMissing) {
    sampleMissingShift <- stats::rnorm(nSamples, 0, missingBySampleSd)
    sampleMissingShift <- sampleMissingShift - mean(sampleMissingShift)
    
    missingFunction <- function(intercept) {
      missingProb <- stats::plogis(
        intercept - kMnar * logData +
          matrix(sampleMissingShift, nrow = nProteins, ncol = nSamples, byrow = TRUE)
      )
      mean(missingProb, na.rm = TRUE) - targetMissing
    }
    
    interceptHat <- stats::uniroot(missingFunction, interval = c(-50, 50))$root
    
    missingProb <- stats::plogis(
      interceptHat - kMnar * logData +
        matrix(sampleMissingShift, nrow = nProteins, ncol = nSamples, byrow = TRUE)
    )
    
    missingMask <- matrix(stats::runif(nProteins * nSamples), nrow = nProteins, ncol = nSamples) < missingProb
    logData[missingMask] <- NA_real_
    
    missingInfo <- list(
      probability = missingProb,
      mask = missingMask,
      sampleShift = sampleMissingShift
    )
  }
  
  metadata <- data.frame(
    Samples = colnames(logData),
    Groups = factor(groups, levels = c("G1", "G2"))
  )
  
  return(
    list(
      logData = logData,
      rawData = 2^logData,
      metadata = metadata
      # missingInfo = missingInfo
    )
  )
}
