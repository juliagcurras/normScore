# Tests for score component functions


# validateNormScoreInput  ####
test_that("validateNormScoreInput prepares valid two-group input", {
    groupData <- data.frame(
        Samples = paste0("S", 1:6),
        Groups = rep(c("A", "B"), each = 3)
    )

    set.seed(123)
    rawData <- matrix(
        runif(120 * 6, min = 100, max = 1000),
        nrow = 120,
        ncol = 6,
        dimnames = list(paste0("P", 1:120), groupData$Samples)
    )

    normalizedDataList <- list(
        Norm1 = as.data.frame(log2(rawData)),
        Norm2 = as.data.frame(log2(rawData + 1))
    )

    result <- validateNormScoreInput(
        normalizedDataList = normalizedDataList,
        groupData = groupData,
        rawData = rawData,
        refGroup = "B",
        altGroup = "A",
        returnDetails = TRUE,
        nBoot = 100
    )

    expect_type(result, "list")
    expect_named(
        result,
        c(
            "normalizedDataList",
            "groupData",
            "rawData",
            "refGroup",
            "altGroup",
            "returnDetails",
            "nBoot",
            "singleGroup"
        )
    )

    expect_false(result$singleGroup)
    expect_equal(result$refGroup, "B")
    expect_equal(result$altGroup, "A")
    expect_true(result$returnDetails)
    expect_equal(result$nBoot, 100)

    expect_true(all(vapply(result$normalizedDataList, is.matrix, logical(1))))
    expect_true("Log" %in% names(result$normalizedDataList))
})


test_that("validateNormScoreInput detects single-group input", {
    groupData <- data.frame(
        Samples = paste0("S", 1:6),
        Groups = rep("A", 6)
    )

    set.seed(123)
    rawData <- matrix(
        runif(120 * 6, min = 100, max = 1000),
        nrow = 120,
        ncol = 6,
        dimnames = list(paste0("P", 1:120), groupData$Samples)
    )

    normalizedDataList <- list(
        Norm1 = log2(rawData),
        Norm2 = log2(rawData + 1)
    )

    result <- validateNormScoreInput(
        normalizedDataList = normalizedDataList,
        groupData = groupData,
        rawData = rawData,
        refGroup = NULL,
        altGroup = NULL,
        returnDetails = TRUE,
        nBoot = 100
    )

    expect_true(result$singleGroup)
    expect_null(result$refGroup)
    expect_null(result$altGroup)
})


test_that("validateNormScoreInput stops when rawData has fewer than 100
        proteins", {
    groupData <- data.frame(
        Samples = paste0("S", 1:6),
        Groups = rep(c("A", "B"), each = 3)
    )

    rawData <- matrix(
        runif(50 * 6, min = 100, max = 1000),
        nrow = 50,
        ncol = 6,
        dimnames = list(paste0("P", 1:50), groupData$Samples)
    )

    normalizedDataList <- list(
        Norm1 = log2(rawData),
        Norm2 = log2(rawData + 1)
    )

    expect_error(
        validateNormScoreInput(
            normalizedDataList = normalizedDataList,
            groupData = groupData,
            rawData = rawData,
            refGroup = "B",
            altGroup = "A",
            returnDetails = TRUE,
            nBoot = 100
        ),
        "minimum of 100 proteins"
    )
})


test_that("validateNormScoreInput reorders normalized datasets to
        match rawData", {
    groupData <- data.frame(
        Samples = paste0("S", 1:6),
        Groups = rep(c("A", "B"), each = 3)
    )

    set.seed(123)
    rawData <- matrix(
        runif(120 * 6, min = 100, max = 1000),
        nrow = 120,
        ncol = 6,
        dimnames = list(paste0("P", 1:120), groupData$Samples)
    )

    shuffledData <- log2(rawData)[rev(rownames(rawData)), rev(colnames(rawData))]

    normalizedDataList <- list(
        Log = shuffledData,
        Norm1 = shuffledData + 0.5,
        Norm2 = log2(rawData) + 1
    )

    result <- validateNormScoreInput(
        normalizedDataList = normalizedDataList,
        groupData = groupData,
        rawData = rawData,
        refGroup = "B",
        altGroup = "A",
        returnDetails = TRUE,
        nBoot = 100
    )

    expect_equal(rownames(result$normalizedDataList$Norm1), rownames(rawData))
    expect_equal(colnames(result$normalizedDataList$Norm1), colnames(rawData))
})


# computeNormScoreItems  ####
test_that("computeNormScoreItems returns expected item columns
        for two groups", {
    groupData <- data.frame(
        Samples = paste0("S", 1:6),
        Groups = rep(c("A", "B"), each = 3)
    )

    set.seed(123)
    rawData <- matrix(
        runif(120 * 6, min = 100, max = 1000),
        nrow = 120,
        ncol = 6,
        dimnames = list(paste0("P", 1:120), groupData$Samples)
    )

    normalizedDataList <- list(
        Norm1 = log2(rawData),
        Norm2 = log2(rawData + 1),
        Norm3 = log2(rawData * 1.1)
    )

    result <- computeNormScoreItems(
        normalizedDataList = normalizedDataList,
        groupData = groupData,
        rawData = rawData,
        refGroup = "B",
        altGroup = "A",
        singleGroup = FALSE
    )

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), length(normalizedDataList))
    expect_equal(rownames(result), names(normalizedDataList))
    expect_equal(colnames(result), paste0("Item", 1:6))
})


test_that("computeNormScoreItems skips Item3 for single-group input", {
    groupData <- data.frame(
        Samples = paste0("S", 1:6),
        Groups = rep("A", 6)
    )

    set.seed(123)
    rawData <- matrix(
        runif(120 * 6, min = 100, max = 1000),
        nrow = 120,
        ncol = 6,
        dimnames = list(paste0("P", 1:120), groupData$Samples)
    )

    normalizedDataList <- list(
        Norm1 = log2(rawData),
        Norm2 = log2(rawData + 1),
        Norm3 = log2(rawData * 1.1)
    )

    result <- computeNormScoreItems(
        normalizedDataList = normalizedDataList,
        groupData = groupData,
        rawData = rawData,
        refGroup = NULL,
        altGroup = NULL,
        singleGroup = TRUE
    )

    expect_s3_class(result, "data.frame")
    expect_equal(rownames(result), names(normalizedDataList))
    expect_equal(colnames(result), paste0("Item", c(1:2, 4:6)))
    expect_false("Item3" %in% colnames(result))
})


# scaleNormScoreItems  ####
test_that("scaleNormScoreItems applies min-max scaling and down-weights Item2", {
    scoreDF <- data.frame(
        Item1 = c(10, 20, 30),
        Item2 = c(1, 2, 3),
        Item3 = c(5, 5, 5),
        row.names = c("Norm1", "Norm2", "Norm3")
    )

    result <- scaleNormScoreItems(scoreDF)

    expect_s3_class(result, "data.frame")
    expect_equal(dim(result), dim(scoreDF))
    expect_equal(rownames(result), rownames(scoreDF))

    expect_equal(result$Item1, c(0, 0.5, 1))
    expect_equal(result$Item2, c(0, 0.05, 0.1))
    expect_equal(result$Item3, c(0, 0, 0))
})


# rankNormScoreItems  ####
test_that("rankNormScoreItems computes totals, applies item0 to Log,
        and ranks rows", {
    scaledScoreDF <- data.frame(
        Item1 = c(0.1, 0.2, 0.3),
        Item2 = c(0.1, 0.2, 0.3),
        row.names = c("Norm1", "Log", "Norm2")
    )

    result <- rankNormScoreItems(
        scaledScoreDF = scaledScoreDF,
        item0 = 3
    )

    expect_s3_class(result, "data.frame")
    expect_true(all(c("Total", "TotalxItem0") %in% colnames(result)))

    expect_equal(result["Norm1", "Total"], 0.2)
    expect_equal(result["Log", "Total"], 0.4)
    expect_equal(result["Log", "TotalxItem0"], 1.2)

    expect_equal(rownames(result), c("Norm1", "Norm2", "Log"))
})


# computeBootstrapNormScore  ####
test_that("computeBootstrapNormScore returns a data frame with the
        expected columns", {
    set.seed(123)

    scoreMatrix <- matrix(
        c(
            0.2, 0.3, 0.1,
            0.4, 0.2, 0.3,
            0.1, 0.5, 0.2
        ),
        nrow = 3,
        byrow = TRUE
    )

    colnames(scoreMatrix) <- c("Item1", "Item2", "Item3")
    rownames(scoreMatrix) <- c("Log", "Norm1", "Norm2")
    item0 <- c(0.9, 0.95, 0.92)

    result <- computeBootstrapNormScore(
        rankedScoreDF = scoreMatrix,
        nBoot = 100, item0 = item0
    )

    expect_s3_class(result, "data.frame")
    expect_equal(
        colnames(result),
        c("normalization", "meanNormScore", "ll95", "ul95")
    )
    expect_equal(nrow(result), ncol(scoreMatrix))
    expect_equal(sort(result$normalization), sort(rownames(scoreMatrix)))
})


test_that("computeBootstrapNormScore returns results sorted by meanNormScore", {
    set.seed(123)

    scoreMatrix <- matrix(
        c(
            0.2, 0.3, 0.1,
            0.4, 0.2, 0.3,
            0.1, 0.5, 0.2
        ),
        nrow = 3,
        byrow = TRUE
    )

    colnames(scoreMatrix) <- c("Item1", "Item2", "Item3")
    rownames(scoreMatrix) <- c("Log", "Norm1", "Norm2")
    item0 <- c(0.9, 0.95, 0.92)

    result <- computeBootstrapNormScore(
        rankedScoreDF = scoreMatrix,
        nBoot = 100, item0 = item0
    )

    expect_true(all(diff(result$meanNormScore) >= 0))
})
