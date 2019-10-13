#' Radius of curvature of prime vertical
#'
#' Calculate radius of curvature of prime vertical at given geodetic latitude.
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
    'X' = xGeodetic,
    'Y' = yGeodetic,
    'Z' = zGeodetic
  )
  return(output)
}

#' Calculate summand of geocentric field component
#'
#' @param legendreTable \code{data.table} modified by \code{.CalculateSchmidtLegendreDerivative}
#' @param radius Radius of curvature of prime vertical at given geodetic latitude
#' @param lon GPS longitude
#' @param latGC GPS latitude, geocentric
.CalculateGeocentricFieldSummand <- function(
  legendreTable,
  radius,
  lon,
  latGC
) {
  # NULLing out data.table-related names before using them to make
  # devtools::check() & CRAN happy
  n <- NULL
  m <- NULL
  P <- NULL
  P_Schmidt <- NULL
  P_Schmidt_muDeriv <- NULL
  g <- NULL
  h <- NULL
  gDot0 <- NULL
  hDot0 <- NULL

  legendreTable[
    , `:=` (
      # Summands of main field vector components
      xGeocentric = -((.kGeomagneticRadius / radius) ^ (n + 2)) *
        (g * cos(m * lon) + h * sin(m * lon)) * P_Schmidt_muDeriv * cos(latGC),
      yGeocentric = ((.kGeomagneticRadius / radius) ^ (n + 2)) * m *
        (g * sin(m * lon) - h * cos(m * lon)) * P_Schmidt / cos(latGC),
      zGeocentric = -(n + 1) * ((.kGeomagneticRadius / radius) ^ (n + 2)) *
        (g * cos(m * lon) + h * sin(m * lon)) * P_Schmidt,

      # Summands of secular variation vector components
      xDotGeocentric = -((.kGeomagneticRadius / radius) ^ (n + 2)) *
        (gDot0 * cos(m * lon) + hDot0 * sin(m * lon)) * P_Schmidt_muDeriv *
        cos(latGC),
      yDotGeocentric = ((.kGeomagneticRadius / radius) ^ (n + 2)) * m *
        (gDot0 * sin(m * lon) - hDot0 * cos(m * lon)) * P_Schmidt / cos(latGC),
      zDotGeocentric = -(n + 1) * ((.kGeomagneticRadius / radius) ^ (n + 2)) *
        (gDot0 * cos(m * lon) + hDot0 * sin(m * lon)) * P_Schmidt
    )
  ]
}
