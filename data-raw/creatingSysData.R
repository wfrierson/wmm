###############################################################################
## The following objects are built and then saved in R/sysdata.rda

###############################################################################
## Define generic constants

# Constant to convert from radians to degrees
.kRadToDegree <- 180 / pi

###############################################################################
## Define special constants

# Flattening constant, f
.kFlatteningConstant <- 1 / 298.257223563

# Eccentricity, e
.kEccentricity <- sqrt(.kFlatteningConstant * (2 - .kFlatteningConstant))

# Earth's semi-major axis, A
.kEarthSemimajorAxis <- 6378137

# Geomagnetic reference radius, a
.kGeomagneticRadius <- 6371200

###############################################################################
## Build lookup table used to speed up calculation of P_{n,m}(mu)

# Import WMM Gauss coefficients
.folderExtdata <- file.path(
  'inst',
  'extdata'
)

# The underlying data is found here:
# https://www.ngdc.noaa.gov/geomag/WMM/data/WMMReports/wmm2000.pdf
# https://www.ngdc.noaa.gov/geomag/WMM/data/WMM2005/TRWMM_2005.pdf
# https://www.ngdc.noaa.gov/geomag/WMM/data/WMM2010/WMM2010COF.zip
# https://www.ngdc.noaa.gov/geomag/WMM/data/WMM2015/WMM2015COF.zip
# https://www.ngdc.noaa.gov/geomag/WMM/data/WMM2015/WMM2015v2COF.zip
#
# The following files are in the extdata folder and are formatted versions of
# the WMM coefficients.
.filenamesWMM <- c(
  'WMM2000.csv',
  'WMM2005.csv',
  'WMM2010.csv',
  'WMM2015.csv',
  'WMM2015v2.csv'
)

.pathsWMM <- file.path(
  .folderExtdata,
  .filenamesWMM
)

.kCoefficientsWMM <- data.table::rbindlist(
  lapply(
    .pathsWMM,
    function(path) {
      versionWMM <- tools::file_path_sans_ext(basename(path))

      output <- data.table::fread(
        path,
        sep = '|',
        header = TRUE,
        colClasses = rep('numeric', 6)
      )[
        , version := (versionWMM)
      ]

      return(output)
    }
  ),
  use.names = TRUE
)

# For programming ease, replace NA values for h and h_dot with 0
.kCoefficientsWMM[
  is.na(h)
  , h := 0
][
  is.na(h_dot)
  , h_dot := 0
]


data.table::setkey(
  .kCoefficientsWMM,
  n, m, version
)

# Partition data.table into list of data.tables split by WMM version.
.kCoefficientsWMM <- split(.kCoefficientsWMM, by = 'version')

# Restate each data.table as a matrix with indices (n, m) for each Gauss
# coefficient
.RestateGaussCoefficient <- function(coefName, coefTable = .kCoefficientsWMM) {
  output <- lapply(
    coefTable,
    function(coefWMM)
      as.matrix(
        data.table::dcast(coefWMM, n ~ m, value.var = coefName)[
          , -c('n')
        ]
      )
  )

  return(output)
}

.kCoefficientsWMMg <- .RestateGaussCoefficient('g')
.kCoefficientsWMMh <- .RestateGaussCoefficient('h')
.kCoefficientsWMMgDot <- .RestateGaussCoefficient('g_dot')
.kCoefficientsWMMhDot <- .RestateGaussCoefficient('h_dot')

# n is degree & m is order.
# Note: nDegree = 13 needed to calculate P_Schmidt_muDeriv, even though only the
# first 12 degrees are summed.
.kLegendreTemplate <- data.table::data.table(n = 1:13)
.kLegendreTemplate <- .kLegendreTemplate[
  , .(m = 0:n)
  , by = n
]

# Reshape .kLegendreTemplate and cast as matrix with indices (n, m) and values
# equal to m. This will be used as a template to store computed Legendre values,
# and so the values of order m are not important.
#
# Note: The column names are intentionally off from the column number in order
# to be consistent with the Legendre indices. The column names will be used for
# the order m.
.kLegendreTemplate <- as.matrix(
  data.table::dcast(.kLegendreTemplate, n ~ m, value.var = 'm')[
    , -c('n')
  ]
)

###############################################################################
# Index sequence to compute all needed associated legendre functions

.kLegendreSequence <- list(
  # n values
  'n' = c(
    1:13, 1:13, 2:13, 3:13, 4:13, 5:13, 6:13, 7:13, 8:13, 9:13, 10:13, 11:13,
    12:13
  )
  , #m values
  'm' = c(
    rep(0, 13), rep(1, 13), rep(2, 12), rep(3, 11), rep(4, 10), rep(5, 9),
    rep(6, 8), rep(7, 7), rep(8, 6), rep(9, 5), rep(10, 4), rep(11, 3),
    rep(12, 2)
  )
)

# Define index used for recursion:
#    Indices <= 5 mean the first 5 associated Legendre functions are used.
#    Index = 6 means the constant m recursion relation is used.
#    Index == 7 means the constant n recursion relation is used.
.kLegendreSequence <- data.table::as.data.table(.kLegendreSequence)[
  , index := ifelse(
    n == 1 & m == 0, 1, ifelse(
    n == 1 & m == 1, 2, ifelse(
    n == 2 & m == 0, 3, ifelse(
    n == 2 & m == 1, 4, ifelse(
    n == 2 & m == 2, 5, ifelse(
    n > 2 & m <= 1, 6, 7
  ))))))
]

# Calculate Schmidt normalization factor
.kLegendreSequence[
  , normalizationFactor := sqrt(2 * factorial(n - m) / factorial(n + m))
]

# Restate .kLegendreSequence back as list of vectors to use in downstream
# loop, i.e., don't lookup values in data.table, just pull the needed values
# in order of .kLegendreSequence.
.kLegendreSequence <- as.list(.kLegendreSequence)

###############################################################################
## Create grid of constants representing Legendre degree and order indices

.kDegreeIndexMatrix <- outer(
  1:13,
  0:13,
  FUN = function(n, m) n
)

.kOrderIndexMatrix <- outer(
  1:13,
  0:13,
  FUN = function(n, m) m
)

# Prevent calculating complex numbers for a few lines of code by removing
# values for unneeded Legendre indices
filterUnneededIndices <- as.vector(.kDegreeIndexMatrix < .kOrderIndexMatrix)

.kDegreeIndexMatrix[filterUnneededIndices] <- NA
.kOrderIndexMatrix[filterUnneededIndices] <- NA

###############################################################################
## Save Objects

# Store objects not directly accessible to user
usethis::use_data(
  .kLegendreTemplate,
  .kCoefficientsWMMg,
  .kCoefficientsWMMh,
  .kCoefficientsWMMgDot,
  .kCoefficientsWMMhDot,
  .kEccentricity,
  .kEarthSemimajorAxis,
  .kGeomagneticRadius,
  .kRadToDegree,
  .kLegendreSequence,
  .kDegreeIndexMatrix,
  .kOrderIndexMatrix
,internal = TRUE, overwrite = TRUE)
