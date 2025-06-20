.datatable.aware = TRUE

#' wmm: R Implementation of World Magnetic Model
#'
#' The \code{wmm} package calculates magnetic field at a given location and time according to the World Magnetic Model.
#'
#' @section WMM functions:
#' This package has 1 exported function, \code{\link{GetMagneticFieldWMM}}, which returns a list of:
#' \itemize{
#'   \item Main field and secular variation vector components in nT and nT/yr, resp.
#'   \item Magnetic element intensities (i.e., horizontal and total intensities, \code{h} & \code{f}) in nT with their secular variation in nT/yr
#'   \item Magnetic element angles (i.e., inclination and declination, \code{i} & \code{d}) in degrees with their secular variation in deg/yr
#' }
#'
#' \code{GetMagneticFieldWMM(lambda_t, phi_t, h_t, t)} = (\code{x}, \code{y}, \code{z}, \code{xDot}, \code{yDot}, \code{zDot}, \code{h}, \code{f}, \code{i}, \code{d}, \code{hDot}, \code{fDot}, \code{iDot}, \code{dDot})
#'
#' @section Acknowledgments:
#' Thanks to:
#' \itemize{
#'     \item The WMM team past, present, and future for making the Gauss coefficients public domain
#'     \item Alex Breeze for tech reviewing the original version of this code, years ago
#' }
#'
#' @noRd
#'
#' @name wmm

NULL
