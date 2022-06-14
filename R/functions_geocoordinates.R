#' Radius of curvature of prime vertical
#'
#' Calculate radius of curvature of prime vertical at given geodetic latitude.
#'
#' @noRd
#'
#' @param latitudeGD Geodetic latitude in decimal degrees
#' @return Radius of curvature of prime vertical at given geodetic latitude
.CalculateRadiusCurvature <- function(latitudeGD) {
  output <- .kEarthSemimajorAxis / sqrt(
    1 - .kEccentricity^2 * sin(latitudeGD / .kRadToDegree)^2
  )
  return(output)
}

#' Convert from Geodetic to Geocentric Coordinates
#'
#' Convert geodetic coordinates to geocentric coordinates
#'
#' @noRd
#'
#' @param latitudeGD Geodetic latitude in decimal degrees
#' @param height Height in meters above ellipsoid (not mean sea level)
#' @return List with first element as geocentric latitude in decimal degrees and second element as geocentric radius
.ConvertGeodeticToGeocentricGPS <- function(latitudeGD, height) {
  radiusCurvature <- .CalculateRadiusCurvature(latitudeGD)
  latitudeGD <- latitudeGD / .kRadToDegree

  p <- (radiusCurvature + height) * cos(latitudeGD)
  z <- (radiusCurvature * (1 - .kEccentricity^2) + height) * sin(latitudeGD)
  r <- sqrt(p^2 + z^2)
  latitudeGC <- asin(z / r) * .kRadToDegree

  output <- list(latitudeGC, r)
  names(output) <- c('latitude_GC', 'radius_GC')
  return(output)
}

#' Geocentric Coordinates to Geodetic Coordinates
#'
#' Convert Geocentric Coordinates to Geodetic Coordinates.
#'
#' @noRd
#'
#' @param xGeocentric X-coordinate in geocentric system
#' @param yGeocentric Y-coordinate in geocentric system
#' @param zGeocentric Z-coordinate in geocentric system
#' @param deltaLat (Geocentric Latitude - Geodetic Latitude) in decimal degrees
#' @return Vector of length 3 representing geodetic coordinates consistent with given geocentric data
.ConvertGeocentricToGeodeticFieldComponents <- function(
  xGeocentric,
  yGeocentric,
  zGeocentric,
  deltaLat
) {
  # Get difference in latitude angle
  deltaLat <- deltaLat / .kRadToDegree
  cosDeltLat <- cos(deltaLat)
  sinDeltLat <- sin(deltaLat)

  xGeodetic <- xGeocentric * cosDeltLat - zGeocentric * sinDeltLat
  yGeodetic <- yGeocentric
  zGeodetic <- xGeocentric * sinDeltLat + zGeocentric * cosDeltLat

  output <- list(
    'x' = xGeodetic,
    'y' = yGeodetic,
    'z' = zGeodetic
  )
  return(output)
}

#' Calculate sum of geocentric field components
#'
#' @noRd
#'
#' @param legendreTable \code{data.table} modified by \code{.CalculateSchmidtLegendreDerivative}
#' @param gaussCoef Gauss coefficients as calculated by \code{.CalculateGaussCoef}
#' @param radius Radius of curvature of prime vertical at given geodetic latitude
#' @param lon GPS longitude
#' @param latGC GPS latitude, geocentric
#' @param deltaLatitude (Geocentric Latitude - Geodetic Latitude) in decimal degrees
.CalculateGeocentricFieldSum <- function(
  legendreTable,
  gaussCoef,
  radius,
  lon,
  latGC,
  deltaLatitude
) {
  cosLonM <- outer(
    1:13,
    0:13,
    FUN = function(n, m) cos(lon * m)
  )

  sinLonM <- outer(
    1:13,
    0:13,
    FUN = function(n, m) sin(lon * m)
  )

  radiusRatioPower <- (.kGeomagneticRadius / radius) ^ (.kDegreeIndexMatrix + 2)

  gaussCoefGCosLonM <- gaussCoef[['g']] * cosLonM
  gaussCoefGSinLonM <- gaussCoef[['g']] * sinLonM
  gaussCoefHCosLonM <- gaussCoef[['h']] * cosLonM
  gaussCoefHSinLonM <- gaussCoef[['h']] * sinLonM

  gaussCoefGDot0CosLonM <- gaussCoef[['gDot0']] * cosLonM
  gaussCoefGDot0SinLonM <- gaussCoef[['gDot0']] * sinLonM
  gaussCoefHDot0CosLonM <- gaussCoef[['hDot0']] * cosLonM
  gaussCoefHDot0SinLonM <- gaussCoef[['hDot0']] * sinLonM

  xGeocentric <- -radiusRatioPower * (
    gaussCoefGCosLonM + gaussCoefHSinLonM
  ) * legendreTable[['Derivative Schmidt P']] * cos(latGC)

  yGeocentric <- radiusRatioPower * .kOrderIndexMatrix * (
    gaussCoefGSinLonM - gaussCoefHCosLonM
  ) * legendreTable[['Schmidt P']] / cos(latGC)

  zGeocentric <- -(.kDegreeIndexMatrix + 1) * radiusRatioPower * (
    gaussCoefGCosLonM + gaussCoefHSinLonM
  ) * legendreTable[['Schmidt P']]

  xDotGeocentric <- -radiusRatioPower * (
    gaussCoefGDot0CosLonM + gaussCoefHDot0SinLonM
  ) * legendreTable[['Derivative Schmidt P']] * cos(latGC)

  yDotGeocentric <- radiusRatioPower * .kOrderIndexMatrix * (
    gaussCoefGDot0SinLonM - gaussCoefHDot0CosLonM
  ) * legendreTable[['Schmidt P']] / cos(latGC)

  zDotGeocentric <- -(.kDegreeIndexMatrix + 1) * radiusRatioPower * (
    gaussCoefGDot0CosLonM + gaussCoefHDot0SinLonM
  ) * legendreTable[['Schmidt P']]

  # Package the sum of geocentric values with deltaLatitude to later calculate
  # geodentric values
  geocentricField <- list(
    sum(xGeocentric[-13, -14], na.rm = TRUE),
    sum(yGeocentric[-13, -14], na.rm = TRUE),
    sum(zGeocentric[-13, -14], na.rm = TRUE),
    deltaLatitude
  )

  geocentricDotField <- list(
    sum(xDotGeocentric[-13, -14], na.rm = TRUE),
    sum(yDotGeocentric[-13, -14], na.rm = TRUE),
    sum(zDotGeocentric[-13, -14], na.rm = TRUE),
    deltaLatitude
  )

  output <- list(
    'Main Field' = geocentricField,
    'Secular Variation' = geocentricDotField
  )

  return(output)
}
