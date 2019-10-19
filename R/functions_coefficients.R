#' Lookup Table for Gauss coefficients g & h
#'
#' Find Gauss coefficient \eqn{g_{n,m}(t)} consistent with the World Magnetic Model.
#'
#' @param n Degree of associated Legendre function to compute
#' @param m Order of associated Legendre function to compute
#' @param t Annualized date time. E.g., 2015-02-01 = (2015 + 32/365) = 2015.088
#' @param wmmVersion String representing WMM version to use. Must be consistent with \code{time} and one of the following: 'derived', 'WMM2000', 'WMM2005', 'WMM2010', 'WMM2015', 'WMM2015v2'. Default 'derived' value will infer the latest WMM version consistent with \code{time}.
#'
#' @return vector of Gauss coefficients, \eqn{g_{n,m}(t)} and \eqn{h_{n,m}(t)}
#'
#' @import data.table
#' @importFrom utils tail
.CalculateGaussCoef <- function(n, m, t, wmmVersion = 'derived') {
  # Check consistency of given time t and wmmVersion
  .CheckVersionWMM(t = t, wmmVersion = wmmVersion)

  # Get WMM versions and reference year associated with given time t
  derivedVersionInfo <- .DeriveVersionInfo(t)

  # Derive reference year
  t0 <- derivedVersionInfo[['year']]

  if(wmmVersion == 'derived') {
    # Assume first value of 'version' output from .DeriveVersionInfo is the most
    # appropriate.
    # E.g., if the derived reference year is 2015, use the out-of-cycle version.
    wmmVersion <- derivedVersionInfo[['version']][1]
  }

  # Get g coefficient
  g0 <- .kCoefficientsWMMg[[wmmVersion]][n, as.character(m)]
  gDot0 <- .kCoefficientsWMMgDot[[wmmVersion]][n, as.character(m)]
  coefG <- g0 + (t - t0) * gDot0

  h0 <- .kCoefficientsWMMh[[wmmVersion]][n, as.character(m)]
  hDot0 <- .kCoefficientsWMMhDot[[wmmVersion]][n, as.character(m)]
  coefH <- h0 + (t - t0) * hDot0

  output <- list(
    'g' = coefG,
    'h' = coefH,
    'gDot0' = gDot0,
    'hDot0' = hDot0
  )

  return(output)
}
