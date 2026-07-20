#' Simulate label-free proteomics-like two-group data
#'
#' Simulates two-group quantitative omics data with controlled signal
#' structure, including abundance-dependent variance, differential
#' expression, sample correlation, global shifts, mean-SD dependence,
#' and optional MNAR missingness.
#'
#' @param nProteins `integer`. Number of proteins or features to simulate.
#'   Default is `10000`.
#' @param nPerGroup `integer`. Number of samples per group. 
#' Default is `20`.
#' @param muMean `numeric`. Mean of the protein-wise abundance distribution on
#'   the log2 scale. Default is `18`.
#' @param muSd `numeric`. Standard deviation of the protein-wise abundance
#'   distribution on the log2 scale. Default is `1.2`.
#' @param muClip `numeric` vector of length 2. Lower and upper bounds used to
#'   truncate simulated protein means. Default is `c(15, 25)`.
#' @param rhoWithin `numeric`. Within-group correlation target for the 
#'   latent sample factor. Default is `0.85`.
#' @param rhoBetween `numeric`. Between-group correlation target for the 
#'   latent sample factor. Default is `0.55`.
#' @param loadingSd `numeric`. Standard deviation of protein-specific loadings
#'   for the correlated latent factor. Default is `0.25`.
#' @param sigmaHi `numeric`. Residual standard deviation at high abundance.
#'   Default is `0.05`.
#' @param sigmaLo `numeric`. Residual standard deviation at low abundance.
#'   Default is `0.4`.
#' @param gammaSigma `numeric`. Controls how strongly residual variance 
#'   depends on abundance. Default is `2.5`.
#' @param propDE `numeric`. Proportion of differentially expressed proteins.
#'   Default is `0.35`.
#' @param logFCSd `numeric`. Standard deviation of differential expression
#'   log-fold changes. Default is `1`.
#' @param logFCMean `numeric`. Mean of differential expression log-fold
#'    changes. Default is `0`.
#' @param heteroLogFC `logical`. If `TRUE`, the variance of log-fold
#'    changes depends on abundance. Default is `TRUE`.
#' @param fcHi `numeric`. Lower multiplier for abundance-dependent logFC
#'   heterogeneity. Default is `0.55`.
#' @param fcLo `numeric`. Upper multiplier for abundance-dependent logFC
#'   heterogeneity. Default is `2.5`.
#' @param gammaFC `numeric`. Controls how strongly logFC variability depends
#'    on abundance. Default is `7`.
#' @param sampleShiftSd `numeric`. Standard deviation of global sample shifts.
#'   Default is `0`.
#' @param sampleShiftCap `numeric`. Maximum absolute value allowed for global
#'   sample shifts. Default is `0.2`.
#' @param sampleSdStrength `numeric`. Strength of sample-specific mean-SD
#'   dependence. A value of `0` implies independence. Default is `0`.
#' @param sampleSdRho `numeric`. Correlation between sample mean rank and
#'   sample-specific SD effect. Must lie between `0` and `1`. 
#'   Default is `0.8`.
#' @param sampleSdCap `numeric`. Maximum absolute value allowed for the
#'   log-multiplier controlling sample-specific SD effects. Default is `0.35`.
#' @param addMissing `logical`. Should MNAR missing values be added?
#'   Default is `TRUE`.
#' @param targetMissing `numeric`. Target overall proportion of missing
#'   Default is `0.001`.
#' @param kMnar `numeric`. Strength of abundance dependence in the MNAR 
#'    missing values. value mechanism. Default is `1.2`.
#' @param missingBySampleSd `numeric`. Standard deviation of sample-specific
#'   missingness shifts. Default is `0.05`.
#' @param seed `integer`. Random seed used for reproducibility. 
#' Default is `9396`.
#'
#' @return
#' A list with the following elements:
#' \describe{
#'   \item{logData}{Numeric matrix of simulated log2-scale data.}
#'   \item{rawData}{Numeric matrix of simulated raw-scale data, obtained 
#'   as `2^logData`.}
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
#'   \item Optional MNAR missingness driven by abundance and sample-specific 
#'   effects.
#' }
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
#' dim(simData$logData)
#' head(simData$metadata)
#'
#' @export

simulateData <- function(
    nProteins = 10000,
    nPerGroup = 20,
    muMean = 18,
    muSd = 1.2,
    muClip = c(15, 25),
    rhoWithin = 0.85,
    rhoBetween = 0.55,
    loadingSd = 0.25,
    sigmaHi = 0.05,
    sigmaLo = 0.4,
    gammaSigma = 2.5,
    propDE = 0.35,
    logFCSd = 1,
    logFCMean = 0,
    heteroLogFC = TRUE,
    fcHi = 0.55,
    fcLo = 2.5,
    gammaFC = 7,
    sampleShiftSd = 0,
    sampleShiftCap = 0.2,
    sampleSdStrength = 0,
    sampleSdRho = 0.8,
    sampleSdCap = 0.35,
    addMissing = TRUE,
    targetMissing = 0.001,
    kMnar = 1.2,
    missingBySampleSd = 0.05,
    seed = 9396
) {
  args <- as.list(environment())
  args$seed <- NULL

  withr::with_seed(
    seed,
    .simulateDataCore(args)
  )
}
