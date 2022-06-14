#' Compute Legendre Components
#'
#' Function that computes the components of the associated Legendre function, \eqn{P_{n,m}(\mu)}{P(mu, n, m)}, only dependent on (degree, order) indices. This function is only used to precompute values.
#'
#' The underlying equation used is: \deqn{P(x, n, m)=(-1)^m * 2^n * (1-x^2)^(m/2) * sum(for m <= k <= n: k!/(k-m)! * x^(k-m) * choose(n, k) * choose((n+k-1)/2, n))}
#'
#' @noRd
#'
#' @param n degree of associated Legendre function
#' @param m order of associated Legendre function
.CalcLegendreComponents <- function(n, m) {
  # Since these values will be in a product, just zero out the non-sense indices
  if (m > n) {
    return(NA)
  }

  mSequence <- seq(from = m, to = n)

  # Keep only the indices where the binomial coefficient,
  # choose((n + k - 1)/2, n), is non-zero.
  mSequence <- mSequence[which((n + mSequence - 1) %% 2 == 1)]

  mRange <- length(mSequence)
  output <- vector(mode = 'numeric', length = mRange)

  for (
    index in seq_along(mSequence)
  ) {
    k <- mSequence[index]

    # Note: When calcuating the magnetic field, the associated values of mu
    # will be multiplied and summed: (1 - mu^2)^(m/2) * mu^mSequence
    output[index] <- 2^n * factorial(k) / factorial(k - m) * choose(n, k) *
      choose((n + k - 1)/2, n)
  }

  return(output)
}

#' Calculate Polynomial Components for Associated Legendre Function
#'
#' Function that computes the polynomial components of \code{mu} that are paired with the output of \code{.CalcLegendreComponents} to create the indiviudal components of the associated Legendre function, \eqn{P_{n,m}(\mu)}{P_{n,m}(mu)}.
#'
#' The underlying equation used is: \deqn{P(x, n, m)=(-1)^m * 2^n * (1-x^2)^(m/2) * sum(for m <= k <= n: k!/(k-m)! * x^(k-m) * choose(n, k) * choose((n+k-1)/2, n))}
#'
#' @noRd
#'
#' @param mu Function argument to \eqn{P_{n,m}(\mu)}{P_{n,m}(mu)}
.CalcPolynomialComponents <- function(mu) {
  output <- (1 - mu^2)^.kSelectedExponentsM * mu^.kSelectedIndicesM
  output <- replace(output, which(is.na(output)), 0)

  return(output)
}

#' Compute Associated Legendre Functions Given Sequence of (degree, order) Indices
#'
#' Procedure that computes the associated Legendre function, \eqn{P_{n,m}(\mu)}{P_{n,m}(mu)}, given a sequence of (degree, order) indices and function argument \eqn{\mu}{mu}. This is computed via a closed-form equation.
#'
#' The underlying equation used is: \deqn{P(x, n, m)=(-1)^m * 2^n * (1-x^2)^(m/2) * sum(for m <= k <= n: k!/(k-m)! * x^(k-m) * choose(n, k) * choose((n+k-1)/2, n))}
#'
#' @noRd
#'
#' @param mu Function argument to \eqn{P_{n,m}(\mu)}{P_{n,m}(mu)}
.CalcLegendre <- function(mu) {
  polynomialComponents <- .CalcPolynomialComponents(mu)

  legendreP <- rowSums(.kLegendreComponents * polynomialComponents, dims = 2)

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
  ) / ((1 - mu) * (1 + mu))

  output <- list(
    'P' = legendreP,
    'Schmidt P' = legendreSchmidtP,
    'Derivative Schmidt P' = legendreDerivSchmidtP
  )

  return(output)
}
