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

# n is degree & m is order.
# Note: nDegree = 13 needed to calculate P_Schmidt_muDeriv, even though only the
# first 12 degrees are summed.
.kLegendreIndices <- data.table::data.table(n = 1:13, key = "n")
.kLegendreIndices <- .kLegendreIndices[
  , .(m = 0:n)
  , by = eval(data.table::key(.kLegendreIndices))
]
data.table::setkey(.kLegendreIndices, n, m)

# Define index used for recursion, see P_recursive for details.
.kLegendreIndices[, index := .I]

# index == 6 means the constant m recursion relation is used.
.kLegendreIndices[
  .kLegendreIndices[!J(0:2)][data.table::CJ(unique(n), 0:1)]
  , index := 6
]

# index == 7 means the constant n recursion relation is used.
.kLegendreIndices[index > 6, index := 7]

# Create explicit row number to improve calculation speed.
.kLegendreIndices[, rowNumb := 1:104]

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

###############################################################################
## Save Objects

# Store objects not directly accessible to user
usethis::use_data(
  .kLegendreIndices,
  .kCoefficientsWMM,
  .kEccentricity,
  .kEarthSemimajorAxis,
  .kGeomagneticRadius,
  .kRadToDegree,
  .kLegendreSequence
,internal = TRUE, overwrite = TRUE)
