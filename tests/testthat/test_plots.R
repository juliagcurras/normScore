

test_that("plotNormScoreDiagnostics returns all diagnostic plots for two groups", {
  groupData <- data.frame(
    Samples = paste0("S", 1:6),
    Groups = rep(c("A", "B"), each = 3)
  )
  
  rawData <- matrix(
    runif(120 * 6, 100, 1000),
    nrow = 120,
    ncol = 6,
    dimnames = list(paste0("P", 1:120), groupData$Samples)
  )
  
  normalizedDataList <- list(
    Norm1 = log2(rawData),
    Norm2 = log2(rawData + 1)
  )
  
  plots <- plotNormScoreDiagnostics(
    normalizedDataList = normalizedDataList,
    groupData = groupData,
    rawData = rawData,
    refGroup = "B",
    altGroup = "A"
  )
  
  expect_type(plots, "list")
  expect_named(plots, paste0("item", 0:6))
  expect_false(is.null(plots$item3))
})


test_that("plotNormScoreDiagnostics skips item3 for single-group input", {
  groupData <- data.frame(
    Samples = paste0("S", 1:6),
    Groups = rep("A", 6)
  )
  
  rawData <- matrix(
    runif(120 * 6, 100, 1000),
    nrow = 120,
    ncol = 6,
    dimnames = list(paste0("P", 1:120), groupData$Samples)
  )
  
  normalizedDataList <- list(
    Norm1 = log2(rawData),
    Norm2 = log2(rawData + 1)
  )
  
  plots <- plotNormScoreDiagnostics(
    normalizedDataList = normalizedDataList,
    groupData = groupData,
    rawData = rawData,
    refGroup = NULL,
    altGroup = NULL
  )
  
  expect_type(plots, "list")
  expect_named(plots, paste0("item", c(0:2, 4:6)))
  expect_null(plots$item3)
})



# plotBootstrapNormScore  ####

test_that("plotBootstrapNormScore returns a ggplot object", {
  bootstrapScore <- data.frame(
    normalization = c("Norm1", "Norm2", "Norm3"),
    meanNormScore = c(0.8, 1.2, 1.5),
    ll95 = c(0.6, 1.0, 1.2),
    ul95 = c(1.0, 1.4, 1.8)
  )
  
  result <- list(bootstrapScore = bootstrapScore)
  
  p <- plotBootstrapNormScore(result)
  
  expect_s3_class(p, "ggplot")
})


test_that("plotBootstrapNormScore stops if bootstrapScore is missing", {
  result <- list(finalRanking = c(Norm1 = 0.8, Norm2 = 1.2))
  
  expect_error(
    plotBootstrapNormScore(result),
    "Bootstrap scores are not available"
  )
})


test_that("plotBootstrapNormScore renames bootstrapScore columns when needed", {
  bootstrapScore <- data.frame(
    method = c("Norm1", "Norm2", "Norm3"),
    mean = c(0.8, 1.2, 1.5),
    lower = c(0.6, 1.0, 1.2),
    upper = c(1.0, 1.4, 1.8)
  )
  
  result <- list(bootstrapScore = bootstrapScore)
  
  p <- plotBootstrapNormScore(result)
  
  expect_s3_class(p, "ggplot")
})


test_that("plotBootstrapNormScore maps expected variables", {
  bootstrapScore <- data.frame(
    normalization = c("Norm1", "Norm2", "Norm3"),
    meanNormScore = c(0.8, 1.2, 1.5),
    ll95 = c(0.6, 1.0, 1.2),
    ul95 = c(1.0, 1.4, 1.8)
  )
  
  result <- list(bootstrapScore = bootstrapScore)
  
  p <- plotBootstrapNormScore(result)
  
  expect_equal(p$labels$x, "Mean normScore [95% CI]")
  expect_equal(p$labels$y, "Normalization")
})



