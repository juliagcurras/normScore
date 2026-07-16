test_that("plotItem0 returns a ggplot object", {
  rawData <- matrix(
    runif(120 * 6, 100, 1000),
    nrow = 120,
    ncol = 6,
    dimnames = list(paste0("P", 1:120), paste0("S", 1:6))
  )
  
  p <- plotItem0(rawData)
  
  expect_s3_class(p, "ggplot")
})


test_that("plotItem1 returns a ggplot object", {
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
  
  p <- plotItem1(
    normalizedDataList = normalizedDataList,
    groupData = groupData
  )
  
  expect_s3_class(p, "ggplot")
})


test_that("plotItem2 returns a ggplot object", {
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
  
  p <- plotItem2(
    normalizedDataList = normalizedDataList,
    groupData = groupData
  )
  
  expect_s3_class(p, "ggplot")
})


test_that("plotItem3 returns an arranged ggplot object", {
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
  
  p <- plotItem3(
    normalizedDataList = normalizedDataList,
    groupData = groupData,
    refGroup = "B",
    altGroup = "A"
  )
  
  expect_s3_class(p, "gg")
})


test_that("plotItem4 returns an arranged ggplot object", {
  rawData <- matrix(
    runif(120 * 6, 100, 1000),
    nrow = 120,
    ncol = 6,
    dimnames = list(paste0("P", 1:120), paste0("S", 1:6))
  )
  
  normalizedDataList <- list(
    Norm1 = log2(rawData),
    Norm2 = log2(rawData + 1)
  )
  
  p <- plotItem4(normalizedDataList)
  
  expect_s3_class(p, "gg")
})


test_that("plotItem5 returns an arranged ggplot object", {
  rawData <- matrix(
    runif(120 * 6, 100, 1000),
    nrow = 120,
    ncol = 6,
    dimnames = list(paste0("P", 1:120), paste0("S", 1:6))
  )
  
  normalizedDataList <- list(
    Norm1 = log2(rawData),
    Norm2 = log2(rawData + 1)
  )
  
  p <- plotItem5(normalizedDataList)
  
  expect_s3_class(p, "gg")
})


test_that("plotItem6 returns an arranged ggplot object", {
  rawData <- matrix(
    runif(120 * 6, 100, 1000),
    nrow = 120,
    ncol = 6,
    dimnames = list(paste0("P", 1:120), paste0("S", 1:6))
  )
  
  normalizedDataList <- list(
    Norm1 = log2(rawData),
    Norm2 = log2(rawData + 1)
  )
  
  p <- plotItem6(normalizedDataList)
  
  expect_s3_class(p, "gg")
})

