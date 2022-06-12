#' Lookup Table for Gauss coefficients g & h
#'
#' Find Gauss coefficient \eqn{g_{n,m}(t)} consistent with the World Magnetic Model.
#'
#' @noRd
#'
#' @param t Annualized date time. E.g., 2015-02-01 = (2015 + 32/365) = 2015.088
#' @param t0 Annualized reference time associated with \code{t}
#' @param wmmVersion String representing WMM version to use. Must be consistent with \code{time} and one of the following: 'derived', 'WMM2000', 'WMM2005', 'WMM2010', 'WMM2015', 'WMM2015v2', 'WMM2020'. Default 'derived' value will infer the latest WMM version consistent with \code{time}.
#'
#' @return vector of Gauss coefficients, \eqn{g_{n,m}(t)} and \eqn{h_{n,m}(t)}
.CalculateGaussCoef <- function(t, t0, wmmVersion = 'derived') {
  gaussG <- .kLegendreTemplate
  gaussH <- .kLegendreTemplate
  gaussGDot <- .kLegendreTemplate
  gaussHDot <- .kLegendreTemplate

  gaussG[1:12, 1:13] <- .kCoefficientsWMMg[[wmmVersion]] +
    (t - t0) * .kCoefficientsWMMgDot[[wmmVersion]]
  gaussH[1:12, 1:13] <- .kCoefficientsWMMh[[wmmVersion]] +
    (t - t0) * .kCoefficientsWMMhDot[[wmmVersion]]
  gaussGDot[1:12, 1:13] <- .kCoefficientsWMMgDot[[wmmVersion]]
  gaussHDot[1:12, 1:13] <- .kCoefficientsWMMhDot[[wmmVersion]]

  output <- list(
    'g' = gaussG,
    'h' = gaussH,
    'gDot0' = gaussGDot,
    'hDot0' = gaussHDot
  )

  return(output)
}
