#' Calculate Expected Magnetic Field from WMM2015
#'
#' Calculate the magnetic field for a given location and time using the fitted spherical harmonic model from the 2015 World Magnetic Model.
#'
#' @param lon GPS longitude
#' @param latGD GPS latitude, geodetic
#' @param latGC GPS latitude, geocentric
#' @param radius Radius of curvature of prime vertical at given geodetic latitude
#' @param time Annualized date time. E.g., 2015-02-01 = (2015 + 32/365) = 2015.088
#' @param wmmVersion String representing WMM version to use. Must be consistent with \code{time} and one of the following: 'derived', 'WMM2000', 'WMM2005', 'WMM2010', 'WMM2015', 'WMM2015v2'. Default 'derived' value will infer the latest WMM version consistent with \code{time}.
#'
#' @return Expected magnetic field from WMM2015 expressed as a vector, \eqn{m_{\lambda_t,\varphi_t,h_t,t}^{WMM}}{m_wmm(lambda_t, phi_t, h_t, t)}
.CalculateMagneticField <- function(
  lon,
  latGD,
  latGC,
  radius,
  time,
  wmmVersion = 'derived'
) {
  # Check consistency of given time and wmmVersion
  derivedVersionInfo <- .CheckVersionWMM(t = time, wmmVersion = wmmVersion)

  if(wmmVersion == 'derived') {
    # Assume first value of 'version' output from .DeriveVersionInfo is the most
    # appropriate.
    # E.g., if the derived reference year is 2015, use the out-of-cycle version.
    wmmVersion <- derivedVersionInfo[['version']][1]
  }

  deltaLatitude <- latGC - latGD   # Leave in degrees for deltaLatitude
  lon <- lon / .kRadToDegree
  latGC <- latGC / .kRadToDegree
  latGD <- latGD / .kRadToDegree

  # Run workhorse of algorithm, i.e., use recursion relations to sequentially
  # calculate different Legendre values to later be summed
  mu <- sin(latGC)
  legendreTable <- .RunLegendreProcedure(mu)

  # Get Gauss coefficients given input time and reference year
  gaussCoef <- .CalculateGaussCoef(
    t = time,
    t0 = derivedVersionInfo[['year']],
    wmmVersion = wmmVersion
  )

  geocentricFields <- .CalculateGeocentricFieldSum(
    legendreTable,
    gaussCoef,
    radius,
    lon,
    latGC,
    deltaLatitude
  )

  # Convert predicted magnetic field from geocentric to geodetic coordinates
  geodentricField <- do.call(
    .ConvertGeocentricToGeodeticFieldComponents,
    geocentricFields[['Main Field']]
  )

  geodentricDotField <- do.call(
    .ConvertGeocentricToGeodeticFieldComponents,
    geocentricFields[['Secular Variation']]
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
