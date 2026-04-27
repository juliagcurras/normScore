test_that("normScore returns only finalRanking when returnDetails is TRUE", {
  simData <- simulateData(
    nProteins = 200,
    nPerGroup = 3,
    addMissing = FALSE,
    seed = 123
  )
  
  normalizedDataList <- list(
    Median = simData$logData,
    Shifted = simData$logData + 0.1
  )
  
  result <- normScore(
    normalizedDataList = normalizedDataList,
    groupData = simData$metadata,
    rawData = simData$rawData,
    returnDetails = TRUE
  )
  
  expect_true(is.list(result))
  expect_equal(names(result), "finalRanking")
  expect_type(result$finalRanking, "double")
  expect_true("Log" %in% names(result$finalRanking))
  expect_true(names(result[[1]][1]) == "Log")
})

test_that("normScore returns detailed output when returnDetails is FALSE", {
  
  simData <- simulateData(
    nProteins = 200,
    nPerGroup = 3,
    addMissing = FALSE,
    seed = 123
  )
  
  normalizedDataList <- list(
    Median = simData$logData,
    Shifted = simData$logData + 0.1
  )
  
  result <- normScore(
    normalizedDataList = normalizedDataList,
    groupData = simData$metadata,
    rawData = simData$rawData,
    returnDetails = FALSE,
    nBoot = 20
  )
  
  expect_true(is.list(result))
  expect_equal(
    names(result),
    c("finalRanking", "detailRanking", "bootstrapScore")
  )
  expect_s3_class(result$detailRanking, "data.frame")
  expect_s3_class(result$bootstrapScore, "data.frame")
})


test_that("normScore includes Log in the ranking even if it is not provided in normalizedDataList", {
  simData <- simulateData(
    nProteins = 200,
    nPerGroup = 3,
    addMissing = FALSE,
    seed = 123
  )
  
  normalizedDataList <- list(
    Median = simData$logData + 0.2,
    Shifted = simData$logData + 0.1
  )
  
  result <- normScore(
    normalizedDataList = normalizedDataList,
    groupData = simData$metadata,
    rawData = simData$rawData,
    returnDetails = TRUE
  )
  
  expect_true("Log" %in% names(result$finalRanking))
})



test_that("normScore fails when fewer than 100 proteins are provided", {
  simData <- simulateData(
    nProteins = 50,
    nPerGroup = 3,
    addMissing = FALSE,
    seed = 123
  )
  
  normalizedDataList <- list(
    Log = simData$logData
  )
  
  expect_error(
    normScore(
      normalizedDataList = normalizedDataList,
      groupData = simData$metadata,
      rawData = simData$rawData,
      returnDetails = TRUE
    ),
    "minimum of 100 proteins"
  )
})



test_that("normScore returns finalRanking sorted in increasing order", {
  simData <- simulateData(
    nProteins = 200,
    nPerGroup = 3,
    addMissing = FALSE,
    seed = 123
  )
  
  normalizedDataList <- list(
    Log = simData$logData,
    Shifted = simData$logData + 0.1,
    Shifted2 = simData$logData + 0.2
    )
  
  result <- normScore(
    normalizedDataList = normalizedDataList,
    groupData = simData$metadata,
    rawData = simData$rawData,
    returnDetails = FALSE
  )
  
  expect_true(all(diff(result$finalRanking) >= 0))
})



