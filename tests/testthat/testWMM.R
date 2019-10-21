testthat::context('Testing WMM...')

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
keyFields <- c('testID', 'wmmVersion')
vectorFields <- c('x', 'y', 'z')
vectorDotFields <- paste0(vectorFields, 'Dot')
testFields <- c(
  vectorFields,
  vectorDotFields
)
calculatedFields <- paste0(testFields, 'Calculated')
testthatFields <- c(keyFields, testFields)

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
  , mget(c(keyFields, calculatedFields))
]
data.table::setnames(
  calculatedData,
  calculatedFields,
  testFields
)

# Perform unit tests
testthat::context('Testing WMM main field benchmarks...')
testthat::test_that('WMM Test Values, Not 2005', {
  expect_true(
    all.equal(
      testData[
        wmmVersion != ' WMM2005'
        , mget(vectorFields)
      ],
      calculatedData[
        wmmVersion != ' WMM2005'
        , mget(vectorFields)
      ],
      tolerance = 2.5e-5
    )
  )
})

testthat::test_that('WMM Test Values, 2005 only', {
  expect_true(
    all.equal(
      testData[
        wmmVersion == ' WMM2005'
        , mget(vectorFields)
      ],
      calculatedData[
        wmmVersion == ' WMM2005'
        , mget(vectorFields)
      ],
      # Setting tolerance to 0.5 due to small precision provided in 2005
      # test values.
      tolerance = 5e-1
    )
  )
})

testthat::context('Testing WMM secular variation field benchmarks...')
testthat::test_that('WMM Test Values (secular)', {
  expect_true(
    all.equal(
      testData[
        , mget(vectorDotFields)
        ],
      calculatedData[
        , mget(vectorDotFields)
        ],
      # Setting tolerance to 0.007 due to small precision provided secular
      # variation test values.
      tolerance = 7e-3
    )
  )
})
