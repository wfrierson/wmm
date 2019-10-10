#' Calculate Expected Magnetic Field from WMM2015
#'
#' Calculate the magnetic field for a given location and time using the fitted spherical harmonic model from the 2015 World Magnetic Model.
#'
#' @param lon GPS longitude
#' @param latGD GPS latitude, geodetic
#' @param latGC GPS latitude, geocentric
#' @param radius Radius of curvature of prime vertical at given geodetic latitude
#' @param timeInYears Annualized date time. E.g., 2015-02-01 = (2015 + 32/365) = 2015.088
#' @param highestDegree Highest degree used to compute the WMM magnetic field. Note: This is a diagnostic. \code{highestDegree} should always be 12.
#' @return Expected magnetic field from WMM2015 expressed as a vector, \eqn{m_{\lambda_t,\varphi_t,h_t,t}^{WMM}}{m_wmm(lambda_t, phi_t, h_t, t)}
#' @examples
#' # This code reproduces the output of the numerical example described in WMM 2015
#' lon <- 240
#' latGD <- -80
#' height <- 1e5
#'
#' geocentric <- .ConvertGeodeticToGeocentricGPS(latGD, height)
#'
#' magField_WMM(
#'     lon = lon
#'     , latGD = latGD
#'     , latGC = geocentric[['latitude_GC']][1]
#'     , radius = geocentric[['radius_GC']][1]
#'     , timeInYears = 2017.5
#' )
#'
#' ## Expected output
#' # X = 5683.51754 95763 nT
#' # Y = 14808.84920 23104 nT
#' # Z = -50163.01336 54779 nT
#'
#' ## Output a/o 2017-10-03
#' # X = 5683.518 nT
#' # Y = 14808.85 nT
#' # Z = -50163.01 nT
#' @seealso \code{\link{getMagneticFieldWMM}}
#'
#' @import data.table
.CalculateMagneticField <- function(
  lon,
  latGD,
  latGC,
  radius,
  timeInYears,
  highestDegree = 12,
  wmmVersion = 'derived'
) {
  deltaLatitude <- latGC - latGD   # Leave in degrees for deltaLatitude
  lon <- lon / .kRadToDegree
  latGC <- latGC / .kRadToDegree
  latGD <- latGD / .kRadToDegree

  legendreTable <- data.table::copy(.kLegendreIndices)

  # Index sequence to compute all needed associated legendre functions
  legendreSequence <- data.table::copy(.kLegendreSequence)

  # Sequential calculations
  mu <- sin(latGC)
  .RunLegendreProcedure(legendreTable, legendreSequence, mu)

  # Derive other legendre functions
  legendreTable[
    , P_Schmidt := ifelse(
      m == 0,
      P,
      sqrt(2 * factorial(n - m) / factorial(n + m)) * P
    )
    ]
  data.table::setorder(legendreTable, m, n)
  legendreTable[
    , P_Schmidt_muDeriv := (
      (n + 1) * mu * P_Schmidt -
        sqrt((n + 1)^2 - m^2) * data.table::shift(P_Schmidt, type = 'lead')
    ) / (1 - mu^2)
    ]
  data.table::setkey(legendreTable, n, m)

  # Lookup Gauss coefficients
  legendreTable[
    , c('g', 'h') := .CalculateGaussCoef(n, m, timeInYears, wmmVersion)
    ]

  # Compute geocentric WMM summand values (i.e., to be summed)
  legendreTable[
    , xGeocentric := -((.kGeomagneticRadius / radius) ^ (n + 2)) *
      (g * cos(m * lon) + h * sin(m * lon)) * P_Schmidt_muDeriv * cos(latGC)
    ]
  legendreTable[
    , yGeocentric := ((.kGeomagneticRadius / radius) ^ (n + 2)) * m *
      (g * sin(m * lon) - h * cos(m * lon)) * P_Schmidt / cos(latGC)
    ]
  legendreTable[
    , zGeocentric := -(n + 1) * ((.kGeomagneticRadius / radius) ^ (n + 2)) *
      (g * cos(m * lon) + h * sin(m * lon)) * P_Schmidt
    ]

  legendreTable <- legendreTable[!J(13)][
    J(1:highestDegree)
    ,.(
      xGeocentric = sum(xGeocentric)
      , yGeocentric = sum(yGeocentric)
      , zGeocentric = sum(zGeocentric)
    )
    ]
  geocentricField <- list(
    legendreTable$xGeocentric,
    legendreTable$yGeocentric,
    legendreTable$zGeocentric,
    deltaLatitude
  )

  output <- do.call(
    .ConvertGeocentricToGeodeticFieldComponents,
    geocentricField
  )

  return(output)
}

#' Calculate Expected Magnetic Field from WMM2015 (user-friendly)
#'
#' Convenience function that takes in geodetic GPS location and annualized time, and returns the expected magnetic field from WMM2015.
#' This is a user-friendly version of \code{magField_WMM}.
#'
#' @param lon GPS longitude
#' @param lat GPS latitude, geodetic
#' @param height GPS height in meters above ellipsoid
#' @param timeInYears Annualized date time. E.g., 2015-02-01 = (2015 + 32/365) = 2015.088
#' @return Expected magnetic field from WMM2015 expressed as a vector, \eqn{m_{\lambda_t,\varphi_t,h_t,t}^{WMM}}{m_wmm(lambda_t, phi_t, h_t, t)}
#' @seealso \code{\link{magField_WMM}}, \code{\link{getDeclinationWMM}}
#' @export
GetMagneticFieldWMM <- function(
  lon,
  lat,
  height,
  timeInYears,
  wmmVersion = 'derived'
) {
  geocentric <- .ConvertGeodeticToGeocentricGPS(lat, height)

  output <- magField_WMM(
    lon = lon,
    latGD = lat,
    latGC = geocentric[['latitude_GC']][1],
    radius = geocentric[['radius_GC']][1],
    timeInYears = timeInYears,
    wmmVersion
  )
  return(output)
}
