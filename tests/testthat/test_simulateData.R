

# simulateData ####

test_that("simulateData returns expected output structure", {
  sim <- simulateData(
    nProteins = 100,
    nPerGroup = 3,
    addMissing = FALSE,
    seed = 123
  )
  
  expect_true(is.list(sim))
  expect_true(all(c("logData", "rawData", "metadata") %in% names(sim)))
  expect_equal(nrow(sim$logData), 100)
  expect_equal(ncol(sim$logData), 6)
  expect_equal(nrow(sim$metadata), 6)
})

test_that("simulateData returns same output using same seed", {
  sim1 <- simulateData(nProteins = 100, seed = 123)
  sim2 <- simulateData(nProteins = 100, seed = 123)
  
  expect_equal(sim1, sim2)
})