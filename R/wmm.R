.datatable.aware = TRUE

#' wmm: R Implementation of World Magnetic Model
#'
#' The \code{wmm} package calculates magnetic field at a given location and time according to the World Magnetic Model.
#'
#' @section WMM functions:
#' This package has 1 exported function, \code{\link{GetMagneticFieldWMM}}, which returns a list of calculated main field and secular variation vector components in nT and nT/yr, resp.: \code{x}, \code{y}, \code{z}, \code{xDot}, \code{yDot}, \code{zDot}.
#'
#' @section Acknowledgements:
#' Thanks to:
#' \itemize{
#'     \item The WMM team past, present, and future for making the Gauss coefficients public domain
#'     \item Alex Breeze for tech reviewing the original version of this code, years ago
#' }
#'
#' @docType package
#' @name wmm

NULL
