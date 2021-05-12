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
