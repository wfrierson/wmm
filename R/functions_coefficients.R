#' Lookup Table for Gauss coefficients g & h
#'
#' Find Gauss coefficient \eqn{g_{n,m}(t)} consistent with the World Magnetic Model.
#'
#' @param n Degree of associated Legendre function to compute
#' @param m Order of associated Legendre function to compute
#' @param t Annualized date time. E.g., 2015-02-01 = (2015 + 32/365) = 2015.088
#' @param wmmVersion
#'
#' @return vector of Gauss coefficients, \eqn{g_{n,m}(t)} and \eqn{h_{n,m}(t)}
#'
#' @import data.table
.CalculateGaussCoef <- function(n, m, t, wmmVersion = 'derived') {
  .CheckVersionWMM(t = t, wmmVersion = wmmVersion)
  derivedVersionInfo <- .DeriveVersionInfo(t)
  t0 <- derivedVersionInfo[['year']]

  if(wmmVersion == 'derived') {
    # Assume last value of 'version' output from .DeriveVersionInfo is the most
    # appropriate.
    # E.g., if the derived reference year is 2015, use the out-of-cycle version.
    wmmVersion <- tail(derivedVersionInfo[['version']], 1)
  }

  # Rename degree and order to avoid using the same name fields in
  # .kCoefficientsWMM.
  nDegree <- n
  mOrder <- m
  coefficientsWMM <- data.table::copy(
    .kCoefficientsWMM[J(nDegree, mOrder, wmmVersion)]
  )

  # Get g coefficient
  g0 <- coefficientsWMM$g
  gDot0 <- coefficientsWMM$g_dot
  coefG <- g0 + (t - t0) * gDot0

  h0 <- coefficientsWMM$h
  hDot0 <- coefficientsWMM$h_dot
  coefH <- h0 + (t - t0) * hDot0

  output <- list(
    'g' = coefG,
    'h' = coefH
  )

  return(output)
}
