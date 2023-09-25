# Define character vectors used for unit tests, which may be changed
keyFields <- c('testID', 'wmmVersion')
vectorFields <- c('x', 'y', 'z')
vectorDotFields <- paste0(vectorFields, 'Dot')
testFields <- c(
  vectorFields,
  vectorDotFields
)
calculatedFields <- paste0(testFields, 'Calculated')
testthatFields <- c(keyFields, testFields)

# Calculate magnetic field values, warning about 'is in'
testthat::test_that("warns: is in the blackout zone", {
  for (rn in which(testData$h < 2000)) {
    testthat::expect_warning(
      testData[
        rn
        , (calculatedFields) := GetMagneticFieldWMM(
          lon = lon,
          lat = lat,
          height = height * 1e3,
          time = year,
          wmmVersion = wmmVersion
        )
        , by = testID
      ]
     , "Location is in the blackout zone")
  }
})

# Calculate magnetic field values, warning about 'is approaching'
testthat::test_that("warns: is approaching the blackout zone", {
  for (rn in which(testData$h >= 2000 & testData$h < 6000)) {
    testthat::expect_warning(
      testData[
        rn
        , (calculatedFields) := GetMagneticFieldWMM(
          lon = lon,
          lat = lat,
          height = height * 1e3,
          time = year,
          wmmVersion = wmmVersion
        )
        , by = testID
      ]
     , "Location is approaching the blackout zone")
  }
})

# Calculate magnetic field values, no warning
testthat::test_that("no warnings", {
  testthat::expect_silent(
    testData[
      h >= 6000
      , (calculatedFields) := GetMagneticFieldWMM(
        lon = lon,
        lat = lat,
        height = height * 1e3,
        time = year,
        wmmVersion = wmmVersion
      )
      , by = testID
    ])
})

testthat::test_that("nothing missed", {
  expect_false(anyNA(testData[, c(calculatedFields), with = FALSE]))
})

# Copy table to help with set of test fields
calculatedData <- data.table::copy(testData)[
  , mget(c(keyFields, calculatedFields))
]
data.table::setnames(
  calculatedData,
  calculatedFields,
  testFields
)

# Perform unit tests
testthat::test_that('WMM Test Values, Not 2005', {
  expect_true(
    all.equal(
      testData[
        wmmVersion != 'WMM2005'
        , mget(vectorFields)
      ],
      calculatedData[
        wmmVersion != 'WMM2005'
        , mget(vectorFields)
      ],
      tolerance = 5e-6
    )
  )
})

testthat::test_that('WMM Test Values, 2005 only', {
  expect_true(
    all.equal(
      testData[
        wmmVersion == 'WMM2005'
        , mget(vectorFields)
      ],
      calculatedData[
        wmmVersion == 'WMM2005'
        , mget(vectorFields)
      ],
      # Setting tolerance to 5e-5 due to small precision provided in 2005
      # test values.
      tolerance = 5e-5
    )
  )
})

testthat::test_that('WMM Test Values (secular), not 2005', {
  expect_true(
    all.equal(
      testData[
        wmmVersion != 'WMM2005'
        , mget(vectorDotFields)
      ],
      calculatedData[
        wmmVersion != 'WMM2005'
        , mget(vectorDotFields)
      ],
      # Setting tolerance to 0.001 due to small precision provided for secular
      # variation test values.
      tolerance = 1e-3
    )
  )
})

testthat::test_that('WMM Test Values (secular), 2005 only', {
  expect_true(
    all.equal(
      testData[
        wmmVersion == 'WMM2005'
        , mget(vectorDotFields)
      ],
      calculatedData[
        wmmVersion == 'WMM2005'
        , mget(vectorDotFields)
      ],
      # Setting tolerance to 0.01 due to small precision provided for secular
      # variation test values.
      tolerance = 1e-2
    )
  )
})
