#' Calculate Expected Magnetic Elements from WMM2025
#'
#' Calculate the magnetic elements (i.e., horizontal intensity, total intensity, inclination, declination, and their secular variation) for given magnetic orthogonal components
#'
#' @noRd
#'
#' @param orthComps named \code{list} containing magnetic orthogonal components
#'
#' @return Expected magnetic components from WMM2025. \code{list}
#'
.CalculateMagneticElements <- function(
  orthComps
) {
  h <- sqrt(orthComps[['x']]^2 + orthComps[['y']]^2)
  f <- sqrt(h^2 + orthComps[['z']]^2)
  i <- atan2(orthComps[['z']], h)
  d <- atan2(orthComps[['y']], orthComps[['x']])

  hDot <- (
    orthComps[['x']] * orthComps[['xDot']] +
      orthComps[['y']] * orthComps[['yDot']]
  ) / h

  fDot <- (
    orthComps[['x']] * orthComps[['xDot']] +
      orthComps[['y']] * orthComps[['yDot']] +
      orthComps[['z']] * orthComps[['zDot']]
  ) / f

  iDot <- (
    h * orthComps[['zDot']] -
      orthComps[['z']] * hDot
  ) / f^2

  dDot <- (
    orthComps[['x']] * orthComps[['yDot']] -
      orthComps[['y']] * orthComps[['xDot']]
  ) / h^2

  output <- list(
    'h' = h,
    'f' = f,
    'i' = i * .kRadToDegree,
    'd' = d * .kRadToDegree,
    'hDot' = hDot,
    'fDot' = fDot,
    'iDot' = iDot * .kRadToDegree,
    'dDot' = dDot * .kRadToDegree
  )

  return(output)
}

#' Calculate Expected Magnetic Field from WMM2025
#'
#' Calculate the magnetic field for a given location and time using the fitted spherical harmonic model from the 2025 World Magnetic Model.
#'
#' @noRd
#'
#' @param lon GPS longitude
#' @param latGD GPS latitude, geodetic
#' @param latGC GPS latitude, geocentric
#' @param radius Radius of curvature of prime vertical at given geodetic latitude
#' @param time Annualized date time. E.g., 2015-02-01 = (2015 + 32/365) = 2015.088
#' @param wmmVersion String representing WMM version to use. Must be consistent with \code{time} and one of the following: 'derived', 'WMM2000', 'WMM2005', 'WMM2010', 'WMM2015', 'WMM2015v2', 'WMM2020', 'WMM2025'. Default 'derived' value will infer the latest WMM version consistent with \code{time}.
#'
#' @return Expected magnetic field from WMM2025, \eqn{m_{\lambda_t,\varphi_t,h_t,t}^{WMM}}{m_wmm(lambda_t, phi_t, h_t, t)}. \code{list}
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
  legendreTable <- .CalcLegendre(mu)

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

  magElements <- .CalculateMagneticElements(output)

  output <- c(
    output,
    magElements
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
#' @param time Annualized date time. E.g., 2015-02-01 = (2015 + 32/365) = 2015.088; optionally an object (length 1) of class 'POSIXt' or 'Date'
#' @param wmmVersion String representing WMM version to use. Must be consistent with \code{time} and one of the following: 'derived', 'WMM2000', 'WMM2005', 'WMM2010', 'WMM2015', 'WMM2015v2', 'WMM2020', 'WMM2025'. Default 'derived' value will infer the latest WMM version consistent with \code{time}.
#'
#' @return \code{list} of calculated main field and secular variation vector components in nT and nT/yr, resp. The magnetic element intensities (i.e., horizontal and total intensities, h & f) are in nT and the magnetic element angles (i.e., inclination and declination, i & d) are in degrees, with their secular variation in nT/yr and deg/yr, resp.: \code{x}, \code{y}, \code{z}, \code{xDot}, \code{yDot}, \code{zDot}, \code{h}, \code{f}, \code{i}, \code{d}, \code{hDot}, \code{fDot}, \code{iDot}, \code{dDot}
#' @export
#'
#' @examples
#' GetMagneticFieldWMM(
#'    lon = 240,
#'    lat = -80,
#'    height = 0,
#'    time = 2025,
#'    wmmVersion = 'WMM2025'
#' )
#'
#' ## Expected output
#' # x = 6117.5 nT
#' # y = 15751.9 nT
#' # z = -52022.5 nT
#' # xDot = 33.3 nT/yr
#' # yDot = -8.6 nT/yr
#' # zDot = 95.5 nT/yr
#' # h = 16898.1 nT
#' # f = 54698.2 nT
#' # i = -72 deg
#' # d = 68.78 deg
#' # hDot = 4.0 nT/yr
#' # fDot = -89.6 nT/yr
#' # iDot = 0.03 deg/yr
#' # dDot = -0.12 deg/yr
#'
#' ## Calculated output
#' #$x
#' #[1] 6117.548
#'
#' #$y
#' #[1] 15751.91
#'
#' #$z
#' #[1] -52022.52
#'
#' #$xDot
#' #[1] 33.27253
#'
#' #$yDot
#' #[1] -8.646338
#'
#' #$zDot
#' #[1] 95.54722
#'
#' #$h
#' #[1] 16898.13
#'
#' #$f
#' #[1] 54698.17
#'
#' #$i
#' #[1] -72.00499
#'
#' #$d
#' #[1] 68.77539
#'
#' #$hDot
#' #[1] 3.985647
#'
#' #$fDot
#' #[1] -89.64208
#'
#' #$iDot
#' #[1] 0.03489031
#'
#' #$dDot
#' #[1] -0.1157767
#'
GetMagneticFieldWMM <- function(
  lon,
  lat,
  height,
  time,
  wmmVersion = 'derived'
) {
  # Validating inputs
  if (
    !(
      !is.na(lat) & !is.na(lon) & !is.na(height) &
      !is.infinite(lat) & !is.infinite(lon) & !is.infinite(height) &
      is.numeric(lat) & is.numeric(lon) & is.numeric(height)
    )
  ) {
    stop(
      "Please check that 'lat', 'lon', and 'height' inputs are numeric."
    )
  }

  geocentric <- .ConvertGeodeticToGeocentricGPS(lat, height)

  if (!is.numeric(time)) {
    if (inherits(time, c("POSIXt", "Date"))) {
      YrJul <- with(
        as.POSIXlt(time),
        c(1900 + year, yday + hour/24 + min/1440 + sec/86400)
      )
      # https://www.timeanddate.com/date/leapyear.html#rules
      leapyear <- +((YrJul[1] %% 4 == 0) & ((!YrJul[1] %% 100 == 0) | (YrJul[1] %% 400 == 0)))
      ydays <- 365 + leapyear
      time <- YrJul[1] + YrJul[2]/ydays
    } else {
      stop(
        "Please check that 'time' input is numeric or in POSIXt or Date formats. Unrecognized 'time' input: ",
        paste(sQuote(class(time)), collapse = ", ")
      )
    }
  }

  output <- .CalculateMagneticField(
    lon = lon,
    latGD = lat,
    latGC = geocentric[['latitude_GC']][1],
    radius = geocentric[['radius_GC']][1],
    time = time,
    wmmVersion = wmmVersion
  )

  h <- output[['h']]
  .CheckBlackoutZone(h)

  return(output)
}

#' @rdname GetMagneticFieldWMM
#' @export
wmm <- GetMagneticFieldWMM
