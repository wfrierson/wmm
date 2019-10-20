#' Calculate Expected Magnetic Field from WMM2015
#'
#' Calculate the magnetic field for a given location and time using the fitted spherical harmonic model from the 2015 World Magnetic Model.
#'
#' @param lon GPS longitude
#' @param latGD GPS latitude, geodetic
#' @param latGC GPS latitude, geocentric
#' @param radius Radius of curvature of prime vertical at given geodetic latitude
#' @param time Annualized date time. E.g., 2015-02-01 = (2015 + 32/365) = 2015.088
#' @param highestDegree Highest degree used to compute the WMM magnetic field. Note: This is a diagnostic. \code{highestDegree} should always be 12.
#' @param wmmVersion String representing WMM version to use. Must be consistent with \code{time} and one of the following: 'derived', 'WMM2000', 'WMM2005', 'WMM2010', 'WMM2015', 'WMM2015v2'. Default 'derived' value will infer the latest WMM version consistent with \code{time}.
#'
#' @return Expected magnetic field from WMM2015 expressed as a vector, \eqn{m_{\lambda_t,\varphi_t,h_t,t}^{WMM}}{m_wmm(lambda_t, phi_t, h_t, t)}
#'
#' @import data.table
.CalculateMagneticField <- function(
  lon,
  latGD,
  latGC,
  radius,
  time,
  highestDegree = 12,
  wmmVersion = 'derived'
) {
  # NULLing out data.table-related names before using them to make
  # devtools::check() & CRAN happy
  J <- NULL
  . <- NULL
  n <- NULL
  m <- NULL
  P <- NULL
  P_Schmidt <- NULL
  P_Schmidt_muDeriv <- NULL
  xGeocentric <- NULL
  yGeocentric <- NULL
  zGeocentric <- NULL
  xDotGeocentric <- NULL
  yDotGeocentric <- NULL
  zDotGeocentric <- NULL
  g <- NULL
  h <- NULL

  deltaLatitude <- latGC - latGD   # Leave in degrees for deltaLatitude
  lon <- lon / .kRadToDegree
  latGC <- latGC / .kRadToDegree
  latGD <- latGD / .kRadToDegree

  # Sequential calculations
  mu <- sin(latGC)
  legendreTable <- .RunLegendreProcedure(mu)

  # Derive other legendre functions
  .CalculateSchmidtLegendre(legendreTable)
  data.table::setorder(legendreTable, m, n)

  .CalculateSchmidtLegendreDerivative(legendreTable, mu)
  data.table::setkey(legendreTable, n, m)

  legendreTable <- legendreTable[1:90]

  # Lookup Gauss coefficients
  legendreTable[
    , c('g', 'h', 'gDot0', 'hDot0') :=
      .CalculateGaussCoef(n, m, time, wmmVersion)
    , by = .(n,m)
  ]

  # Compute geocentric WMM summand values (i.e., to be summed)
  .CalculateGeocentricFieldSummand(
    legendreTable,
    radius,
    lon,
    latGC
  )

  # Sum of geocentric values for relevant degrees
  legendreTable <- legendreTable[
    ,.(
      xGeocentric = sum(xGeocentric),
      yGeocentric = sum(yGeocentric),
      zGeocentric = sum(zGeocentric),
      xDotGeocentric = sum(xDotGeocentric),
      yDotGeocentric = sum(yDotGeocentric),
      zDotGeocentric = sum(zDotGeocentric)
    )
  ]

  # Package the geocentric sums with deltaLatitude for downstream inputs
  geocentricField <- list(
    legendreTable$xGeocentric,
    legendreTable$yGeocentric,
    legendreTable$zGeocentric,
    deltaLatitude
  )

  geocentricDotField <- list(
    legendreTable$xDotGeocentric,
    legendreTable$yDotGeocentric,
    legendreTable$zDotGeocentric,
    deltaLatitude
  )

  # Convert predicted magnetic field from geocentric to geodetic coordinates
  geodentricField <- do.call(
    .ConvertGeocentricToGeodeticFieldComponents,
    geocentricField
  )

  geodentricDotField <- do.call(
    .ConvertGeocentricToGeodeticFieldComponents,
    geocentricDotField
  )

  output <- union(
    geodentricField,
    geodentricDotField
  )

  names(output) <- c(
    'x', 'y', 'z',
    'xDot', 'yDot', 'zDot'
  )

  return(output)
}

#' Calculate Expected Magnetic Field from WMM
#'
#' Function that takes in geodetic GPS location and annualized time, and returns the expected magnetic field from WMM.
#'
#' @param lon GPS longitude
#' @param lat GPS latitude, geodetic
#' @param height GPS height in meters above ellipsoid
#' @param time Annualized date time. E.g., 2015-02-01 = (2015 + 32/365) = 2015.088
#' @param wmmVersion String representing WMM version to use. Must be consistent with \code{time} and one of the following: 'derived', 'WMM2000', 'WMM2005', 'WMM2010', 'WMM2015', 'WMM2015v2'. Default 'derived' value will infer the latest WMM version consistent with \code{time}.
#'
#' @return list of calculated main field and secular variation vector components in nT and nT/yr, resp.: \code{x}, \code{y}, \code{z}, \code{xDot}, \code{yDot}, \code{zDot}
#' @export
#'
#' @examples
#' GetMagneticFieldWMM(
#'    lon = 240,
#'    lat = -80,
#'    height = 1e5,
#'    time = 2017.5,
#'    wmmVersion = 'WMM2015'
#' )
#'
#' ## Expected output
#' # X = 5683.51754 95763 nT
#' # Y = 14808.84920 23104 nT
#' # Z = -50163.01336 54779 nT
#'
#' ## Calculated Output
#' # X = 5683.518 nT
#' # Y = 14808.85 nT
#' # Z = -50163.01 nT
GetMagneticFieldWMM <- function(
  lon,
  lat,
  height,
  time,
  wmmVersion = 'derived'
) {
  geocentric <- .ConvertGeodeticToGeocentricGPS(lat, height)

  output <- .CalculateMagneticField(
    lon = lon,
    latGD = lat,
    latGC = geocentric[['latitude_GC']][1],
    radius = geocentric[['radius_GC']][1],
    time = time,
    wmmVersion = wmmVersion
  )
  return(output)
}
