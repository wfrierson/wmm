# Define path to WMM test data
folderExtdata <- file.path(
  system.file(package = 'wmm'),
  'extdata'
)

pathTestData <- file.path(
  folderExtdata,
  'WMMTestValues.csv'
)

# Import WMM test data
testData <- data.table::fread(
  pathTestData,
  sep = '|',
  header = TRUE,
  stringsAsFactors = FALSE
)[
  , testID := .I
]

# Define character vectors used for unit tests, which may be changed
keyField <- c('testID', 'wmmVersion')
testFields <- c('x', 'y', 'z')
calculatedFields <- paste0(testFields, 'Calculated')
testthatFields <- c(keyField, testFields)

# Calculate magnetic field values
testData[
  , (calculatedFields) := GetMagneticFieldWMM(
    lon = lon,
    lat = lat,
    height = height * 1e3,
    time = year,
    wmmVersion = wmmVersion
  )
  , by = testID
]

# Copy table to help with set of test fields
calculatedData <- data.table::copy(testData)[
  , mget(c(keyField, calculatedFields))
]
data.table::setnames(
  calculatedData,
  calculatedFields,
  testFields
)

# Perform unit tests
testthat::context('Testing WMM benchmarks...')
testthat::test_that('WMM Test Values, Not 2005', {
  expect_true(
    all.equal(
      testData[
        wmmVersion != ' WMM2005'
        , mget(testthatFields)
      ],
      calculatedData[wmmVersion != ' WMM2005'],
      tolerance = 5e-4
    )
  )
})

testthat::test_that('WMM Test Values, 2005 only', {
  expect_true(
    all.equal(
      testData[
        wmmVersion == ' WMM2005'
        , mget(testthatFields)
      ],
      calculatedData[wmmVersion == ' WMM2005'],
      # Setting tolerance to 0.5 due to small precision provided in 2005
      # test values.
      tolerance = 5e-1
    )
  )
})

