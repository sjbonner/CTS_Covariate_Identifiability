generate_z <- function(a, Ncap, d_params) {
  ## Create vector of means
  if (a == 1) {
    drift <- d_params$drift
  } else {
    drift <- head(d_params$drift, -(a - 1))
  }

  mu <- cumsum(c(d_params$mu_init, drift))

  ## Create var-cov matrix
  v <- lower.tri(matrix(nrow = Ncap - a + 1, ncol = Ncap - a + 1), diag = T)
  Sigma <- d_params$sigmasq * v %*% t(v)

  ## Generate z's
  z <- c(rep(NA, a - 1), mvrnorm(1, mu, Sigma))

  return(as.vector(z))
}

generate_d <- function(a, z, Ncap, s_params) {
  ## Compute survival probabilities
  phi <- survival(z[a:(Ncap - 1)], s_params)

  ## Compute probability of dying on remaining occasions
  probs <- cumprod(c(1, phi)) * c(1 - phi, 1)

  ## Sample last occasion individual is alive
  d <- sample(a:Ncap, size = 1, prob = probs)

  return(d)
}

generate_c <- function(a, z, d, Ncap, c_params) {

  ## If only alive for one occasion
  if (a == d) {
    ch <- c(rep(0, a - 1), 1, rep(0, Ncap - a))
  }

  ## Generate captures
  ch <- c(
    rep(0, a - 1), # No captures before first capture
    1, # First capture
    (runif(d - a, 0, 1) <= c_params$p[(a + 1):d]), # Until last alive
    rep(0, Ncap - d)
  ) # After death

  ## Simulate recovery
  if (d < Ncap) {
    ch[d + 1] <- 2 * rbinom(1, 1, c_params$lambda[d])
  }

  return(ch)
}

generate_ind <- function(id, Ncap, d_params, s_params, c_params) {
  # Generates a full history for a single individual
  # Individual must be captured at least once to appear in dataset
  # id is a number here but may be more complex

  ## Generate time of first capture
  a <- sample(1:(Ncap - 1), 1)

  ## Simulate covariate values
  z <- generate_z(a, Ncap, d_params)

  ## Simulte survival
  d <- generate_d(a, z, Ncap, s_params)

  ## Simulate captures
  ch <- generate_c(a, z, d, Ncap, c_params)

  ## Time of last capture
  b <- max((1:Ncap)[ch == 1])

  ## Censor covariate values
  z <- ifelse(ch == 1, z, NA)

  return(list(id = id, z = z, ch = ch, a = a, b = b, d = d))
}

generate_data <- function(Nind, Ncap, d_params, s_params, c_params) {
  # Generates a dataset using the following information
  # Nind= number of individuals
  # Ncap= number of capture occasions
  # d_params= parameters of the drift process
  # s_params= parameters of the survival function
  # c_params= parameters of the capture function

  # Create id numbers
  id <- 1:Nind

  # Create capture-recapture process for each individual using function generate_ind()
  capture_dat <- lapply(id, generate_ind, Ncap, d_params, s_params, c_params)
  return(capture_dat)
}
