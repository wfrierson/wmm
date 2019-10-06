###############################################################################
## The following objects are built and then saved in R/sysdata.rda

###############################################################################
## Define generic constants

# Constant to convert from radians to degrees
.radToDegree <- 180 / pi

###############################################################################
## Define specialconstants

# Flattening constant
.f <- 1 / 298.257223563

# Eccentricity
.e <- sqrt(.f * (2 - .f))

# Earth's semi-major axis
.A <- 6378137

# Geomagnetic reference radius
.a <- 6371200

###############################################################################
## Build lookup table used to speed up calculation of P_{n,m}(mu)

# Import WMM Gauss coefficients
.folder.extdata <- file.path(
  'inst',
  'extdata'
)

# The underlying data is found here:
# https://www.ngdc.noaa.gov/geomag/WMM/data/WMM2005/TRWMM_2005.pdf
# https://www.ngdc.noaa.gov/geomag/WMM/data/WMM2010/WMM2010COF.zip
# https://www.ngdc.noaa.gov/geomag/WMM/data/WMM2015/WMM2015COF.zip
# https://www.ngdc.noaa.gov/geomag/WMM/data/WMM2015/WMM2015v2COF.zip
#
# The following files are in the extdata folder and are formatted versions of
# the WMM coefficients.
.filenames.wmm <- c(
  'WMM2005.csv',
  'WMM2010.csv',
  'WMM2015.csv',
  'WMM2015v2.csv'
)

.paths.wmm <- file.path(
  .folder.extdata,
  .filenames.wmm
)

.coefficientsWMM <- data.table::rbindlist(
  lapply(
    .paths.wmm,
    function(path) {
      versionWMM <- tools::file_path_sans_ext(basename(path))

      output <- data.table::fread(
        path,
        sep = '|',
        header = TRUE
      )[
        , version := (versionWMM)
      ]

      return(output)
    }
  ),
  use.names = TRUE
)
data.table::setkey(
  .coefficientsWMM,
  n, m, version
)

# n is degree & m is order.
# Note: nDegree = 13 needed to calculate P_Schmidt_muDeriv, even though only the
# first 12 degrees are summed.
.P <- data.table::data.table(n = 1:13, key = "n")
.P <- .P[, .(m = 0:n), by = eval(data.table::key(.P))]
data.table::setkey(.P, n, m)

# Define index used for recursion, see P_recursive for details.
.P[, index := .I]

# index == 6 means the constant m recursion relation is used.
.P[.P[!J(0:2)][data.table::CJ(unique(n), 0:1)], index := 6]

# index == 7 means the constant n recursion relation is used.
.P[index > 6, index := 7]

# Create explicit row number to improve calculation speed.
.P[, rowNumb := 1:104]

###############################################################################
## Save Objects

# Store objects not directly accessible to user
usethis::use_data(
  .P,
  .coefficientsWMM,
  .f,
  .e,
  .A,
  .a,
  .radToDegree
,internal = TRUE, overwrite = TRUE)
