

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



# MSE function ####
test_that("mse computes mean squared error correctly", {
  actual <- c(1, 2, 3)
  predicted <- c(1, 2, 4)
  
  expected <- mean((actual - predicted)^2)
  
  expect_equal(mse(actual, predicted), expected)
})

test_that("mse is zero for identical vectors", {
  x <- c(2, 4, 6)
  
  expect_equal(mse(x, x), 0)
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












