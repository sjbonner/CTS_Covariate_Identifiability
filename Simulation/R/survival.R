survival <- function(z, s.params) {
  ## Returns survival probabilities for given weights
  ## Logistic function of weight

  ## Set default bounds if not provided
  if (is.null(s.params$lower)) {
    s.params$lower <- 1
  }

  if (is.null(s.params$upper)) {
    s.params$upper <- 1
  }

  ## Compute linear predictor
  eta <- s.params$beta[1] + s.params$beta[2] * z

  ## Compute survival probability
  phi <- s.params$lower + (s.params$upper - s.params$lower) / (1 + exp(-eta))

  return(phi)
}
