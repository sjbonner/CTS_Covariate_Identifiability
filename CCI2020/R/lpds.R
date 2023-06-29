dlpd_p <- function(p, ptrue, lambdatrue, phitrue, normalize = TRUE, scale = FALSE){
  # Density of the limiting posterior distribution for p
  
  # Initialized output vector
  lpd <- rep(0, length(p))
  
  # Compute transparent reparametrization
  rho1 <- phitrue * ptrue
  rho2 <- (1- phitrue) * lambdatrue
  
  # Compute lpd
  phi <- rho1/p
  valid <- intersect(which(phi > rho1), which(phi < 1-rho2))
  
  alpha <- log(phi[valid]) - log(1 - phi[valid])
  
  lpd[valid] <- dt(alpha/10, 1, 0) / (1-phi[valid])^(1 + (lambdatrue > 0))
           
  if(scale){
    # Set maximum to 1
    lpd <- lpd/max(lpd, na.rm = TRUE)
  }
  else if(normalize){
    # Set integral to 1
    c <- integrate(dlpd_p,rho1/(1-rho2), 1, ptrue = ptrue, lambdatrue = lambdatrue, phitrue = phitrue, normalize = FALSE)$value
    lpd <- lpd/c
  }
  
  return(lpd)
}

Elpd_p <- function(ptrue, lambdatrue, phitrue){
  integrand <- function(p, ptrue, phitrue, lambdatrue){
    p * dlpd_p(p, ptrue, lambdatrue, phitrue, scale = FALSE, normalize = FALSE)
  }
  
  integrate(integrand, 0, 1, ptrue = ptrue, lambdatrue = lambdatrue, phitrue = phitrue)$value /
    integrate(dlpd_p, 0, 1, ptrue = ptrue, lambdatrue = lambdatrue, phitrue = phitrue, normalize = FALSE)$value
}

Vlpd_p <- function(ptrue, lambdatrue, phitrue){
  integrand <- function(p, ptrue, phitrue, lambdatrue){
      p^2 * dlpd_p(p, ptrue, lambdatrue, phitrue, scale = FALSE, normalize = FALSE)
  }
  
  EP <- Elpd_p(ptrue, lambdatrue, phitrue)
  
  E2P <- integrate(integrand, 0, 1, ptrue = ptrue, lambdatrue = lambdatrue, phitrue = phitrue)$value /
    integrate(dlpd_p, 0, 1, ptrue = ptrue, lambdatrue = lambdatrue, phitrue = phitrue, normalize = FALSE)$value

  return(E2P - EP^2)
}

## Lambda
dlpd_lambda <- function(lambda, ptrue, lambdatrue, phitrue, normalize = TRUE, scale = FALSE){
  # Density of the limiting posterior distribution for lambda
  
  # Initialized output vector
  lpd <- rep(0, length(lambda))
  
  # Compute transparent reparametrization
  rho1 <- phitrue * ptrue
  rho2 <- (1- phitrue) * lambdatrue
  
  # Compute lpd
  phi <- 1 - rho2/lambda
  valid <- intersect(which(phi > rho1), which(phi < 1-rho2))
  
  alpha <- log(phi[valid]) - log(1 - phi[valid])
  
  lpd[valid] <- dt(alpha/10, 1, 0) / phi[valid]^2
  
  if(scale){
    # Set maximum to 1
    lpd <- lpd/max(lpd, na.rm = TRUE)
  }
  else if(normalize){
    # Set integral to 1
    c <- integrate(dlpd_lambda,0,1, ptrue = ptrue, lambdatrue = lambdatrue, phitrue = phitrue, normalize = FALSE)
    lpd <- lpd/c
  }
  
  return(lpd)
}

Elpd_lambda <- function(ptrue, lambdatrue, phitrue){
  integrand <- function(lambda, ptrue, phitrue, lambdatrue){
    lambda * dlpd_lambda(lambda, ptrue, lambdatrue, phitrue, scale = FALSE, normalize = FALSE)
  }
  
  integrate(integrand, 0, 1, ptrue = ptrue, lambdatrue = lambdatrue, phitrue = phitrue)$value /
    integrate(dlpd_lambda, 0, 1, ptrue = ptrue, lambdatrue = lambdatrue, phitrue = phitrue, normalize = FALSE)$value
}

Vlpd_lambda <- function(ptrue, lambdatrue, phitrue){
  integrand <- function(lambda, ptrue, phitrue, lambdatrue){
    lambda^2 * dlpd_lambda(lambda, ptrue, lambdatrue, phitrue, scale = FALSE, normalize = FALSE)
  }
  
  ELambda <- Elpd_lambda(ptrue, lambdatrue, phitrue)
  
  E2Lambda <- integrate(integrand, 0, 1, ptrue = ptrue, lambdatrue = lambdatrue, phitrue = phitrue)$value /
    integrate(dlpd_lambda, 0, 1, ptrue = ptrue, lambdatrue = lambdatrue, phitrue = phitrue, normalize = FALSE)$value
  
  return(E2Lambda - ELambda^2)
}

## Alpha
dlpd_alpha <- function(alpha, ptrue, lambdatrue, phitrue, normalize = TRUE, scale = FALSE){
  # Density of the limiting posterior distribution for alpha
  
  # Initialized output vector
  lpd <- rep(0, length(alpha))
  
  # Compute transparent reparametrization
  rho1 <- phitrue * ptrue
  rho2 <- (1- phitrue) * lambdatrue
  
  # Compute lpd
  phi <- (1 + exp(-alpha)) ^ -1
  valid <- intersect(which(phi > rho1), which(phi < 1-rho2))

  lpd[valid] <- dt(alpha[valid]/10, 1, 0) / (phi[valid] * (1 - phi[valid]))
  
  if(scale){
    # Set maximum to 1
    lpd <- lpd/max(lpd, na.rm = TRUE)
  }
  else if(normalize){
    # Set integral to 1
    c <- integrate(dlpd_lambda,0,1, ptrue = ptrue, lambdatrue = lambdatrue, phitrue = phitrue, normalize = FALSE)
    lpd <- lpd/c
  }
  
  return(lpd)
}

Elpd_alpha <- function(ptrue, lambdatrue, phitrue){
  integrand <- function(alpha, ptrue, phitrue, lambdatrue){
    alpha * dlpd_alpha(alpha, ptrue, lambdatrue, phitrue, scale = FALSE, normalize = FALSE)
  }
  
  integrate(integrand, 0, 1, ptrue = ptrue, lambdatrue = lambdatrue, phitrue = phitrue)$value /
    integrate(dlpd_alpha, 0, 1, ptrue = ptrue, lambdatrue = lambdatrue, phitrue = phitrue, normalize = FALSE)$value
}

Vlpd_alpha <- function(ptrue, lambdatrue, phitrue){
  integrand <- function(alpha, ptrue, phitrue, lambdatrue){
    alpha^2 * dlpd_alpha(alpha, ptrue, lambdatrue, phitrue, scale = FALSE, normalize = FALSE)
  }
  
  Ealpha <- Elpd_alpha(ptrue, lambdatrue, phitrue)
  
  E2alpha <- integrate(integrand, 0, 1, ptrue = ptrue, lambdatrue = lambdatrue, phitrue = phitrue)$value /
    integrate(dlpd_alpha, 0, 1, ptrue = ptrue, lambdatrue = lambdatrue, phitrue = phitrue, normalize = FALSE)$value
  
  return(E2alpha - Ealpha^2)
}