
# isSameNormalization ####

test_that("isSameNormalization returns TRUE for identical datasets", {
  x <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("P1", "P2"), c("S1", "S2"))
  )
  
  y <- x
  
  expect_true(isSameNormalization(x, y))
})

test_that("isSameNormalization returns FALSE for numerical differences above tolerance", {
  x <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("P1", "P2"), c("S1", "S2"))
  )
  
  y <- matrix(
    c(1, 2, 3, 4.01),
    nrow = 2,
    dimnames = list(c("P1", "P2"), c("S1", "S2"))
  )
  
  expect_false(isSameNormalization(x, y, tolerance = 1e-8))
})

test_that("isSameNormalization returns FALSE for different dimensions", {
  x <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("P1", "P2"), c("S1", "S2"))
  )
  
  y <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3
  )
  
  expect_false(isSameNormalization(x, y))
})

test_that("isSameNormalization returns FALSE for different row names", {
  x <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("P1", "P2"), c("S1", "S2"))
  )
  
  y <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("P1", "P3"), c("S1", "S2"))
  )
  
  expect_false(isSameNormalization(x, y))
})


# addLogIfMissing function ####
test_that("addLogIfMissing adds Log when missing", {
  rawData <- data.frame(
    S1 = c(100, 200),
    S2 = c(120, 240)
  )
  
  normalizedList <- list(
    VSN = rawData / 10
  )
  
  result <- addLogIfMissing(normalizedList, rawData)
  
  expect_true("Log" %in% names(result))
  expect_equal(result[["Log"]], log(rawData, base = 2))
})

test_that("addLogIfMissing renames existing matching normalization to Log", {
  rawData <- data.frame(
    S1 = c(100, 200),
    S2 = c(120, 240)
  )
  
  normalizedList <- list(
    customName = log(rawData, base = 2)
  )
  
  result <- addLogIfMissing(normalizedList, rawData)
  
  expect_length(result, 1)
  expect_true("Log" %in% names(result))
  expect_false("customName" %in% names(result))
})




# CV function ####

test_that("cv returns the coefficient of variation as proportion", {
  x <- c(10, 12, 8, 10)
  
  expected <- stats::sd(x) / mean(x)
  
  expect_equal(cv(x), expected)
})

test_that("cv returns percentage when proportion is FALSE", {
  x <- c(10, 12, 8, 10)
  
  expected <- (stats::sd(x) / mean(x)) * 100
  
  expect_equal(cv(x, proportion = FALSE), expected)
})

test_that("cv handles missing values correctly", {
  x <- c(10, 12, NA, 10)
  
  expected <- stats::sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  
  expect_equal(cv(x, na.rm = TRUE), expected)
})






# MAPE function ####
test_that("mape computes mean absolute percentage error as percentage", {
  actual <- c(100, 200, 400)
  predicted <- c(110, 180, 420)
  
  expected <- mean(abs((actual - predicted) / actual)) * 100
  
  expect_equal(mape(actual, predicted, proportion = FALSE), expected)
})

test_that("mape computes mean absolute percentage error as proportion", {
  actual <- c(100, 200, 400)
  predicted <- c(110, 180, 420)
  
  expected <- mean(abs((actual - predicted) / actual))
  
  expect_equal(mape(actual, predicted, proportion = TRUE), expected)
})




# rleMAPE function ####
test_that("rleMAPE returns a single non-negative numeric value", {
  mat <- matrix(
    c(1, 2, 3,
      2, 3, 4,
      5, 6, 7),
    nrow = 3,
    byrow = TRUE
  )
  
  result <- rleMAPE(mat)
  
  expect_type(result, "double")
  expect_length(result, 1)
  expect_false(is.na(result))
  expect_gte(result, 0)
})

test_that("rleMAPE is zero when all columns are identical", {
  mat <- matrix(
    c(1, 1, 1,
      2, 2, 2,
      3, 3, 3,
      4, 4, 4),
    nrow = 4,
    byrow = TRUE
  )
  
  result <- rleMAPE(mat)
  
  expect_equal(result, 0)
})
test_that("rleMAPE handles missing values when enough data are available", {
  mat <- matrix(
    c(1,  1, NA,
      2,  2,  2,
      3, NA,  3,
      4,  4,  4),
    nrow = 4,
    byrow = TRUE
  )
  
  result <- rleMAPE(mat)
  
  expect_type(result, "double")
  expect_length(result, 1)
  expect_false(is.na(result))
  expect_gte(result, 0)
})





# tiMAPE function ####
test_that("tiMAPE returns 0 when all samples have identical distributions", {
  dataExample <- data.frame(
    S1 = c(10, 20, 30, 40),
    S2 = c(10, 20, 30, 40),
    S3 = c(10, 20, 30, 40)
  )
  
  expect_equal(tiMAPE(dataExample), 0)
})

test_that("tiMAPE returns a positive value when sample quantiles differ across samples", {
  dataExample <- data.frame(
    S1 = c(10, 20, 30, 40),
    S2 = c(12, 22, 32, 42),
    S3 = c(8, 18, 28, 38)
  )
  
  expect_true(tiMAPE(dataExample) > 0)
})

test_that("tiMAPE matches the manually computed expected value", {
  dataExample <- data.frame(
    S1 = c(10, 20, 30, 40),
    S2 = c(12, 22, 32, 42),
    S3 = c(8, 18, 28, 38)
  )
  
  sampleQuantiles <- apply(dataExample, 2, stats::quantile, na.rm = TRUE, simplify = FALSE)
  q1 <- sapply(sampleQuantiles, "[[", 2)
  sampleMedians <- sapply(sampleQuantiles, "[[", 3)
  q3 <- sapply(sampleQuantiles, "[[", 4)
  
  expected <- 
    mean(abs((stats::median(sampleMedians) - sampleMedians) / stats::median(sampleMedians))) +
    mean(abs((stats::median(q1) - q1) / stats::median(q1))) +
    mean(abs((stats::median(q3) - q3) / stats::median(q3)))
  
  expect_equal(tiMAPE(dataExample), expected)
})





# diffAreas function ####
test_that("diffAreas returns zero when slope is zero", {
  expect_equal(
    diffAreas(intPred = 5, coefPred = 0, minRange = 0, maxRange = 10),
    0
  )
})

test_that("diffAreas returns a non-negative value", {
  result <- diffAreas(intPred = 2, coefPred = -0.5, minRange = 0, maxRange = 10)
  
  expect_true(result >= 0)
})



# meanSDdiffArea function ####

test_that("meanSDdiffArea returns 0 when sample standard deviations are identical", {
  dataExample <- data.frame(
    S1 = c(10, 20, 30, 40),
    S2 = c(11, 21, 31, 41),
    S3 = c(12, 22, 32, 42)
  )
  
  expect_equal(meanSDdiffArea(dataExample), 0)
})

test_that("meanSDdiffArea returns a non-negative value", {
  dataExample <- data.frame(
    S1 = c(10, 20, 30, 40),
    S2 = c(10, 20, 30, 50),
    S3 = c(10, 20, 30, 60)
  )
  
  expect_true(meanSDdiffArea(dataExample) >= 0)
})

test_that("meanSDdiffArea matches the manually computed expected value", {
  dataExample <- data.frame(
    S1 = c(10, 20, 30, 40),
    S2 = c(10, 20, 30, 50),
    S3 = c(10, 20, 30, 60)
  )
  
  sampleMeans <- apply(dataExample, 2, base::mean, na.rm = TRUE)
  sampleSDs <- apply(dataExample, 2, stats::sd, na.rm = TRUE)
  
  dfAux <- data.frame(
    Mean = sampleMeans,
    SD = sampleSDs,
    Sample = colnames(dataExample)
  )
  
  dfAux <- dfAux[order(dfAux$Mean), ]
  dfAux$Order <- seq_len(nrow(dfAux))
  
  fit <- stats::lm(SD ~ Order, data = dfAux)
  slope <- fit$coefficients["Order"]
  
  expected <- diffAreas(
    intPred = 0,
    coefPred = slope,
    maxRange = max(dfAux$Order),
    minRange = min(dfAux$Order),
    intExpected = 0
  )
  
  expect_equal(meanSDdiffArea(dataExample), unname(expected))
})





# getPCV function ####
test_that("getPCV returns the mean CV for each group", {
  dataExample <- data.frame(
    S1 = c(10, 20, 30),
    S2 = c(12, 22, 28),
    S3 = c(100, 200, 300),
    S4 = c(110, 210, 290)
  )
  
  groupDataExample <- data.frame(
    Samples = c("S1", "S2", "S3", "S4"),
    Groups = c("A", "A", "B", "B")
  )
  
  groupsExample <- c("A", "B")
  
  expected <- c(
    A = mean(groupProteinCV("A", groupDataExample, dataExample), na.rm = TRUE),
    B = mean(groupProteinCV("B", groupDataExample, dataExample), na.rm = TRUE)
  )
  
  expect_equal(
    getPCV(data = dataExample, groups = groupsExample, groupData = groupDataExample),
    expected
  )
})

test_that("getPCV returns a named numeric vector", {
  dataExample <- data.frame(
    S1 = c(10, 20, 30),
    S2 = c(12, 22, 28),
    S3 = c(100, 200, 300),
    S4 = c(110, 210, 290)
  )
  
  groupDataExample <- data.frame(
    Samples = c("S1", "S2", "S3", "S4"),
    Groups = c("A", "A", "B", "B")
  )
  
  result <- getPCV(
    data = dataExample,
    groups = c("A", "B"),
    groupData = groupDataExample
  )
  
  expect_type(result, "double")
  expect_equal(names(result), c("A", "B"))
})

test_that("getPCV returns 0 when all within-group CV values are 0", {
  dataExample <- data.frame(
    S1 = c(10, 20, 30),
    S2 = c(10, 20, 30),
    S3 = c(100, 200, 300),
    S4 = c(100, 200, 300)
  )
  
  groupDataExample <- data.frame(
    Samples = c("S1", "S2", "S3", "S4"),
    Groups = c("A", "A", "B", "B")
  )
  
  expect_equal(
    getPCV(
      data = dataExample,
      groups = c("A", "B"),
      groupData = groupDataExample
    ),
    c(A = 0, B = 0)
  )
})



# groupProteinCV function ####
test_that("groupProteinCV computes protein-wise CV within the selected group", {
  dataExample <- data.frame(
    S1 = c(10, 20, 30),
    S2 = c(12, 22, 28),
    S3 = c(100, 200, 300),
    S4 = c(110, 210, 290)
  )
  
  groupDataExample <- data.frame(
    Samples = c("S1", "S2", "S3", "S4"),
    Groups = c("A", "A", "B", "B")
  )
  
  expected <- apply(dataExample[, c("S1", "S2")], 1, cv, proportion = FALSE)
  
  expect_equal(
    groupProteinCV(group = "A", groupData = groupDataExample, data = dataExample),
    expected
  )
})

test_that("groupProteinCV returns one value per protein", {
  dataExample <- data.frame(
    S1 = c(10, 20, 30),
    S2 = c(12, 22, 28),
    S3 = c(100, 200, 300),
    S4 = c(110, 210, 290)
  )
  
  groupDataExample <- data.frame(
    Samples = c("S1", "S2", "S3", "S4"),
    Groups = c("A", "A", "B", "B")
  )
  
  result <- groupProteinCV(group = "A", groupData = groupDataExample, data = dataExample)
  
  expect_type(result, "double")
  expect_length(result, nrow(dataExample))
})

test_that("groupProteinCV uses only samples from the requested group", {
  dataExample <- data.frame(
    S1 = c(10, 20, 30),
    S2 = c(10, 20, 30),
    S3 = c(100, 200, 300),
    S4 = c(100, 200, 300)
  )
  
  groupDataExample <- data.frame(
    Samples = c("S1", "S2", "S3", "S4"),
    Groups = c("A", "A", "B", "B")
  )
  
  expect_equal(
    groupProteinCV(group = "A", groupData = groupDataExample, data = dataExample),
    c(0, 0, 0)
  )
  
  expect_equal(
    groupProteinCV(group = "B", groupData = groupDataExample, data = dataExample),
    c(0, 0, 0)
  )
})





# withinGroupCorrelations function ####
test_that("withinGroupCorrelations returns numeric correlations", {
  dataExample <- data.frame(
    S1 = 1:4,
    S2 = c(1:3, 5),
    S3 = c(2, 4:6),
    S4 = 10:13,
    S5 = 4:1,
    S6 = c(1, 7, 2, 3)
  )
  
  groupData <- data.frame(
    Samples = c("S1", "S2", "S3", "S4", "S5", "S6"),
    Groups = c("A", "A", "A", "B", "B", "B")
  )
  
  result <- withinGroupCorrelations(dataExample, groupData)
  
  expect_type(result, "double")
  expect_length(result, 6)
})


test_that("withinGroupCorrelations output change for different methods", {
  dataExample <- data.frame(
    S1 = c(1, 2, 3, 4),
    S2 = c(1, 2, 3, 5),
    S3 = c(21, 15, 10, 6),
    S4 = c(1, 3, 7, 9)
  )
  
  groupData <- data.frame(
    Samples = c("S1", "S2", "S3", "S4"),
    Groups = c("A", "A", "B", "B")
  )
  
  resultS <- withinGroupCorrelations(dataExample, groupData, method = "spearman")
  resultP <- withinGroupCorrelations(dataExample, groupData, method = "pearson")
  
  expect_all_false(resultS == resultP)
  expect_all_true(resultS[1]>0 & resultS[2]<0)
  expect_all_true(resultP[1]>0 & resultP[2]<0)
})




# bootstrapRowScores ####
test_that("bootstrapRowScores returns column sums for the resampled rows", {
  dataExample <- data.frame(
    Norm1 = c(1, 2, 3),
    Norm2 = c(4, 5, 6)
  )
  
  indicesExample <- c(1, 3)
  
  expected <- colSums(dataExample[indicesExample, , drop = FALSE])
  
  expect_equal(
    bootstrapRowScores(data = dataExample, indices = indicesExample),
    expected
  )
})

test_that("bootstrapRowScores allows repeated row indices", {
  dataExample <- data.frame(
    Norm1 = c(1, 2, 3),
    Norm2 = c(4, 5, 6)
  )
  
  indicesExample <- c(2, 2, 3)
  
  expected <- colSums(dataExample[indicesExample, , drop = FALSE])
  
  expect_equal(
    bootstrapRowScores(data = dataExample, indices = indicesExample),
    expected
  )
})

test_that("bootstrapRowScores returns one value per column", {
  dataExample <- data.frame(
    Norm1 = c(1, 2, 3),
    Norm2 = c(4, 5, 6),
    Norm3 = c(7, 8, 9)
  )
  
  result <- bootstrapRowScores(data = dataExample, indices = c(1, 2))
  
  expect_type(result, "double")
  expect_length(result, ncol(dataExample))
})





# computeBootstrapNormScore  ####
test_that("computeBootstrapNormScore returns a data frame with the expected columns", {
  set.seed(123)
  
  scoreMatrix <- matrix(
    c(0.2, 0.3, 0.1,
      0.4, 0.2, 0.3,
      0.1, 0.5, 0.2),
    nrow = 3,
    byrow = TRUE
  )
  
  colnames(scoreMatrix) <- c("Norm1", "Norm2", "Norm3")
  
  result <- computeBootstrapNormScore(scoreMatrix, nBoot = 20)
  
  expect_s3_class(result, "data.frame")
  expect_equal(
    colnames(result),
    c("normalization", "meanNormScore", "ll95", "ul95")
  )
})

test_that("computeBootstrapNormScore returns one row per normalization method", {
  set.seed(123)
  
  scoreMatrix <- matrix(
    c(0.2, 0.3, 0.1,
      0.4, 0.2, 0.3,
      0.1, 0.5, 0.2),
    nrow = 3,
    byrow = TRUE
  )
  
  colnames(scoreMatrix) <- c("Norm1", "Norm2", "Norm3")
  
  result <- computeBootstrapNormScore(scoreMatrix, nBoot = 20)
  
  expect_equal(nrow(result), ncol(scoreMatrix))
  expect_equal(sort(result$normalization), sort(colnames(scoreMatrix)))
})

test_that("computeBootstrapNormScore returns results sorted by meanNormScore", {
  set.seed(123)
  
  scoreMatrix <- matrix(
    c(0.2, 0.3, 0.1,
      0.4, 0.2, 0.3,
      0.1, 0.5, 0.2),
    nrow = 3,
    byrow = TRUE
  )
  
  colnames(scoreMatrix) <- c("Norm1", "Norm2", "Norm3")
  
  result <- computeBootstrapNormScore(scoreMatrix, nBoot = 20)
  
  expect_true(all(diff(result$meanNormScore) >= 0))
})






# plotBootstrapNormScore  ####

test_that("plotBootstrapNormScore returns a ggplot object", {
  bootstrapScores <- data.frame(
    normalization = c("Norm1", "Norm2", "Norm3"),
    meanNormScore = c(0.25, 0.40, 0.32),
    ll95 = c(0.20, 0.35, 0.28),
    ul95 = c(0.30, 0.45, 0.36)
  )
  
  result <- plotBootstrapNormScore(bootstrapScores)
  
  expect_s3_class(result, "ggplot")
})

test_that("plotBootstrapNormScore works when input column names differ", {
  bootstrapScores <- data.frame(
    A = c("Norm1", "Norm2", "Norm3"),
    B = c(0.25, 0.40, 0.32),
    C = c(0.20, 0.35, 0.28),
    D = c(0.30, 0.45, 0.36)
  )
  
  result <- plotBootstrapNormScore(bootstrapScores)
  
  expect_s3_class(result, "ggplot")
})

test_that("plotBootstrapNormScore keeps the expected axis labels", {
  bootstrapScores <- data.frame(
    normalization = c("Norm1", "Norm2", "Norm3"),
    meanNormScore = c(0.25, 0.40, 0.32),
    ll95 = c(0.20, 0.35, 0.28),
    ul95 = c(0.30, 0.45, 0.36)
  )
  
  result <- plotBootstrapNormScore(bootstrapScores)
  
  expect_equal(result$labels$x, "Mean normScore [95% CI]")
  expect_equal(result$labels$y, "Normalization")
})










