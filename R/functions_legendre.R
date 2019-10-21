#' Run recursion to compute associated Legendre functions
#'
#' Use recursion relations to compute the associated Legendre function, \eqn{P_{n,m}(\mu)}{P_{n,m}(mu)}. User supplies degree \eqn{n} and order \eqn{m} as well as associated Legendre functions with smaller degree and order indices for recursion.
#' When \eqn{n \leq 2}{n <= 2}, \eqn{P_{n,m}(\mu)}{P_{n,m}(mu)} is directly calculated with known functions (i.e., no recursion). When \eqn{n > 2}, the following recursion relations are used based on the order \eqn{m}:
#' \deqn{P_{n > 2, m \leq 1}(\mu) = \frac{(2n - 1) \cdot \mu \cdot P_{n-1,m}(\mu) - (n - 1 + m) \cdot P_{n-2,m}(\mu)}{(n - m)}}{P_{n > 2, m <= 1}(mu) = ((2 * n - 1) * mu * Pn_1 - (n - 1 + m) * Pn_2) / (n - m)}
#' \deqn{P_{n > 2, m > 1}(\mu) = \frac{2\mu(m - 1) \cdot P_{n,m-1}(\mu)}{\sqrt{1 - \mu^2}} - (n + m - 1) \cdot (n - m + 2) \cdot P_{n,m-2}(\mu)}{P_{n > 2, m > 1}(mu) = 2 * mu * (m - 1) * Pm_1 / sqrt(1 - mu^2) - (n + m - 1) * (n - m + 2) * Pm_2}
#'
#' @param n Degree of associated Legendre function to compute
#' @param m Order of associated Legendre function to compute
#' @param mu Argument of associated Legendre function to compute
#' @param index Index from \code{.kLegendreIndices} associated with \code{n} and \code{m}
#' @param Pn_1 \eqn{P_{n-1,m}(\mu)}{P_{n-1,m}(mu)}
#' @param Pn_2 \eqn{P_{n-2,m}(\mu)}{P_{n-2,m}(mu)}
#' @param Pm_1 \eqn{P_{n,m-1}(\mu)}{P_{n,m-1}(mu)}
#' @param Pm_2 \eqn{P_{n,m-2}(\mu)}{P_{n,m-2}(mu)}
#'
#' @return \eqn{P_{n,m}(\mu)}{P_{n,m}(mu)}, scalar
.CalculateRecursiveLegendre <- function(
  n,
  m,
  mu,
  index,
  Pn_1 = NULL,
  Pn_2 = NULL,
  Pm_1 = NULL,
  Pm_2 = NULL
) {
  output <- switch(
    index,
    # n == 1 & m == 0, index == 1
    mu,
    # n == 1 & m == 1, index == 2
    sqrt(1 - mu^2),
    # n == 2 & m == 0, index == 3
    (3 * mu^2 - 1) / 2,
    # n == 2 & m == 1, index == 4
    3 * mu * sqrt(1 - mu^2),
    # n == 2 & m == 2, index == 5
    3 * (1 - mu^2),
    # n > 2 & m <= 1, index == 6
    ((2 * n - 1) * mu * Pn_1 - (n - 1 + m) * Pn_2) / (n - m),
    # n > 2 & m > 1, index == 7
    2 * mu * (m - 1) * Pm_1 / sqrt(1 - mu^2) - (n + m - 1) * (n - m + 2) * Pm_2
  )

  return(output)
}

#' Compute Associated Legendre Functions Given Sequence of (degree, order) Indices
#'
#' Procedure that computes the associated Legendre function, \eqn{P_{n,m}(\mu)}{P_{n,m}(mu)}, given a sequence of (degree, order) indices and function argument \eqn{\mu}{mu}. This is computed via recursive relationships for Legendre functions.
#'
#' @param mu Function argument to \eqn{P_{n,m}(\mu)}{P_{n,m}(mu)}
.RunLegendreProcedure <- function(mu) {
  legendreP <- .kLegendreTemplate
  legendreSchmidtP <- .kLegendreTemplate
  legendreDerivSchmidtP <- .kLegendreTemplate

  for(
    x in seq_along(.kLegendreSequence[['n']])
  ) {
    nDegree <- .kLegendreSequence[['n']][x]
    mOrder <- .kLegendreSequence[['m']][x]
    index <- .kLegendreSequence[['index']][x]

    legendreValue <- .CalculateRecursiveLegendre(
      nDegree, mOrder, mu, index,
      Pn_1 = legendreP[(nDegree - 1), as.character(mOrder)],
      Pn_2 = legendreP[(nDegree - 2), as.character(mOrder)],
      Pm_1 = legendreP[nDegree, as.character(mOrder - 1)],
      Pm_2 = legendreP[nDegree, as.character(mOrder - 2)]
    )

    legendreP[nDegree, as.character(mOrder)] <- legendreValue
  }

  # Compute semi-normalized Schmidt Legendre values
  legendreSchmidtP <- .kNormalizationFactors * legendreP

  # Get the next n value of the Schmidt semi-normalized P calculation
  legendreSchmidtPNext <- rbind(
    legendreSchmidtP[-1, ],
    matrix(
      rep(0, 14),
      ncol = ncol(legendreSchmidtP),
      dimnames = dimnames(legendreSchmidtP)
    )
  )

  legendreDerivSchmidtP <- (
    (.kDegreeIndexMatrix + 1) * mu * legendreSchmidtP -
      sqrt((.kDegreeIndexMatrix + 1)^2 - .kOrderIndexMatrix^2) *
      legendreSchmidtPNext
  ) / (1 - mu^2)

  output <- list(
    'P' = legendreP,
    'Schmidt P' = legendreSchmidtP,
    'Derivative Schmidt P' = legendreDerivSchmidtP
  )

  return(output)
}
