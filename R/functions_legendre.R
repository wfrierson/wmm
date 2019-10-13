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
#' @param Pn_1 \eqn{P_{n-1,m}(\mu)}{P_{n-1,m}(mu)}
#' @param Pn_2 \eqn{P_{n-2,m}(\mu)}{P_{n-2,m}(mu)}
#' @param Pm_1 \eqn{P_{n,m-1}(\mu)}{P_{n,m-1}(mu)}
#' @param Pm_2 \eqn{P_{n,m-2}(\mu)}{P_{n,m-2}(mu)}
#' @return \eqn{P_{n,m}(\mu)}{P_{n,m}(mu)}, scalar
#'
#' @import data.table
.CalculateRecursiveLegendre <- function(
  n,
  m,
  mu,
  Pn_1 = NULL,
  Pn_2 = NULL,
  Pm_1 = NULL,
  Pm_2 = NULL
) {
  # NULLing out data.table-related names before using them to make
  # devtools::check() & CRAN happy
  J <- NULL

  # Rename degree and order to avoid using the same name fields in
  # .kLegendreIndices.
  nDegree <- n
  mOrder <- m

  index <- .kLegendreIndices[J(nDegree, mOrder)]$index

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
#' Procedure that recursively computes associated Legendre function, \eqn{P_{n,m}(\mu)}{P_{n,m}(mu)}, given a sequence of (degree, order) indices and function argument \eqn{\mu}{mu}.
#'
#' @param legendreTable copy of internal data.table \code{.kLegendreIndices} to store intermediate function values
#' @param legendreSequence Sequence of (degree, order) indices contained in a list
#' @param mu Function argument to \eqn{P_{n,m}(\mu)}{P_{n,m}(mu)}
#'
#' @import data.table
.RunLegendreProcedure <- function(legendreTable, legendreSequence, mu) {
  # NULLing out data.table-related names before using them to make
  # devtools::check() & CRAN happy
  J <- NULL

  invisible(lapply(
    # The following vector is equivalent to code below, which is constant:
    # seq_along(legendreSequence[['n']])
    1:103,
    function(x) {
      # Rename degree and order to avoid using the same name fields in
      # .kCoefficientsWMM.
      nDegree <- legendreSequence[[1]][x]
      mOrder <- legendreSequence[[2]][x]
      rowNumb <- legendreTable[J(nDegree, mOrder)]$rowNumb

      legendreValue <- if(nDegree <= 2) {
        .CalculateRecursiveLegendre(
          nDegree, mOrder, mu
        )
      } else if(mOrder <= 1) {
        .CalculateRecursiveLegendre(
          nDegree, mOrder, mu,
          Pn_1 = legendreTable[J(nDegree - 1, mOrder)]$P,
          Pn_2 = legendreTable[J(nDegree - 2, mOrder)]$P
        )
      } else {
        .CalculateRecursiveLegendre(
          nDegree, mOrder, mu,
          Pm_1 = legendreTable[J(nDegree, mOrder - 1)]$P,
          Pm_2 = legendreTable[J(nDegree, mOrder - 2)]$P
        )
      }

      data.table::set(legendreTable, rowNumb, 'P', legendreValue)
    }
  ))
}

#' Calculate Schmidt semi-normalized Legendre function
#'
#' @param legendreTable \code{data.table} modified by \code{.RunLegendreProcedure}
.CalculateSchmidtLegendre <- function(legendreTable) {
  legendreTable[
    , P_Schmidt := ifelse(
      m == 0,
      P,
      sqrt(2 * factorial(n - m) / factorial(n + m)) * P
    )
  ]
}

#' Calculate mu-derivative of Schmidt semi-normalized Legendre function
#'
#' @param legendreTable \code{data.table} modified by \code{.CalculateSchmidtLegendre}
#' @param mu Function argument to \eqn{P_{n,m}(\mu)}{P_{n,m}(mu)}
.CalculateSchmidtLegendreDerivative <- function(legendreTable, mu) {
  legendreTable[
    , P_Schmidt_muDeriv := (
      (n + 1) * mu * P_Schmidt -
        sqrt((n + 1)^2 - m^2) * data.table::shift(P_Schmidt, type = 'lead')
    ) / (1 - mu^2)
  ]
}
