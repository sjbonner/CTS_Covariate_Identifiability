##' Run truncated complete data model with recaptures and recoveries
##'
##' Run model with truncated likelihood
##' @title Run model with truncated likelihood
##' @param k Truncation factor
##' @param indata Input data
##' @param model Model ("logit","scaled_logit", or "generalized_logit")
##' @param inits Initial values
##' @param priors Specification of hyperparameters
##' @param pars Parameters to monitor 
##' @param chains Number of chains
##' @param burnin Length of burnin period
##' @param sampling Length of sampling period
##' @export
##' @import spark
##' @import nimble
##' @import coda
##' @return 
##' @author Simon Bonner
run_complete_mrr_model <- function(k,
                                   indata,
                                   model = c("logit", "scaled", "generalized"),
                                   inits = NULL,
                                   priors = NULL,
                                   pars = NULL,
                                   chains = 3,
                                   burnin = 1000,
                                   sampling = 1000) {

  ## 0) Format data for spark
  Nind <- length(indata)
  Ncap <- length(indata[[1]]$ch)

  ch <- t(sapply(indata, "[[", "ch"))
  z <- t(sapply(indata, "[[", "z"))

  spark_data <- list(
    chmat = ch,
    other = z,
    freq = rep(1, length(indata))
  )

  ## 1) Truncate data
  trunc_data <- spark(indata = spark_data, k = k, datatype = "livedead")

  ## 2) Construct delta vector (occasions between release and recapture/recovery)
  delta <- ifelse(
    trunc_data$recapture + trunc_data$recovery > 0,
    pmax(trunc_data$recapture,trunc_data$recovery ) - trunc_data$release,
    pmin(Ncap, Ncap - trunc_data$release)
  )

  ## 3) Format covariate matrix
  Z <- matrix(NA, nrow = trunc_data$nrelease, ncol = k + 1)

  for (i in 1:trunc_data$nrelease) {
    Z[i, 1] <- trunc_data$other[trunc_data$ind[i], trunc_data$release[i]]
    if (trunc_data$recapture[i] > 0) {
      Z[i, delta[i] + 1] <- trunc_data$other[trunc_data$ind[i], trunc_data$recapture[i]]
    }
  }

  ## 4) Order data so that recaptures come before recoveries come before neither
  tmp_recap <- which(trunc_data$recapture > 0)
  tmp_recov <- which(trunc_data$recovery > 0)
           
  index <- c(tmp_recap, tmp_recov, (1:trunc_data$nrelease)[-c(tmp_recap,tmp_recov)])
  
  nrecapture <- length(tmp_recap)
  nrecovery <- length(tmp_recov)

  ## 5) Construct data list
  trunc_jags_data <- list(
    nocc = Ncap,
    nrecapture = nrecapture,
    nrecovery = nrecovery,
    release = trunc_data$release[index],
    nrelease = trunc_data$nrelease,
    delta = delta[index],
    Z = Z[index, ],
    dummy = rep(1, nrow(trunc_data$chmat))
  )

  ## Set initial values
  if (is.null(inits)) {
    inits <- list(
      beta.phi = c(0, 0),
      p = .5,
      lambda = .5,
      mu.z = rep(0, Ncap - 1),
      tau.z = 1
    )

    if (model == "scaled_logit") {
      inits$upper <- 1
    }
    if (model == "generalized_logit") {
      inits$lower <- 0
      inits$upper <- 1
    }
  }
  
  ## Set hyperparameters

  ## Covariate model
  if(is.null(priors[["mu"]])){
    trunc_jags_data$mu.z.hyper <- c(0,1)
  }
  else{
    trunc_jags_data$mu.z.hyper <- priors$mu
  }

  if(is.null(priors[["sigma"]])){
    trunc_jags_data$sigma.z.hyper <- c(.74^2, 4)
  }
  else{
    trunc_jags_data$sigma.z.hyper <- priors$sigma
  }
    
  ## Survival
  if(is.null(priors[["phi"]])){
    trunc_jags_data$beta.phi.hyper <- rbind(c(0, .01, 1), c(0, .16, 1))
  }
  else{
    trunc_jags_data$beta.phi.hyper <- priors$phi
  }

  ## Capture
  if(is.null(priors[["p"]])){
    trunc_jags_data$p.hyper <- c(0,1)
  }
  else{
    trunc_jags_data$p.hyper <- priors$p
  }

  ## Recovery
  if(is.null(priors[["lambda"]])){
    trunc_jags_data$lambda.hyper <- c(0,1)
  }
  else{
    trunc_jags_data$lambda.hyper <- priors$lambda
  }
  
  ## Identify parameters to monitor
  if(is.null(pars)){
    pars <- c("beta.phi", "p", "lambda")
    
    if (model == "scaled_logit") {
      pars <- c(pars, "upper")
    }
    if (model == "generalized_logit") {
      pars <- c(pars, "lower", "upper")
    }
  }
  
  ## Run model
  model_file <- system.file("Nimble",paste0("complete_mrr_", model, ".R"),
                            package="CCI2020")
  
  ## Run model
  nimble_model <- readBUGSmodel(model_file,
    data = trunc_jags_data,
    inits = inits
  )

  time <- system.time(
    trunc_jags_mcmc <- nimbleMCMC(
      model = nimble_model,
      monitors = pars,
      nchains = chains,
      nburnin = 0,
      niter = burnin + sampling,
      samplesAsCodaMCMC = TRUE
    )
  )

  ## Return summaries
  output <- list(
    summary = summary(window(trunc_jags_mcmc, start = burnin + 1)),
    mcmc = trunc_jags_mcmc,
    time = time
  )

  if (chains > 1) {
    output$ess <- effectiveSize(window(trunc_jags_mcmc, start = burnin + 1))
    output$diag <- gelman.diag(window(trunc_jags_mcmc, start = burnin + 1))
  }

  return(output)
}

#' @export
drecap <- nimbleFunction(
  run = function(d = double(0),
                 phi = double(1),
                 p = double(0)) {
    
    # Compute probability of being recaptured on the d-th
    # occasion following release.
    
    if (d == 1) {
      # Recaptured on subsequent occasion
      P <- phi[1] * p
    }
    else{
      # Recaptured later than subsequent occasion
      P <- prod(phi[1:d]) *
        (1 - p)^(d-1) *
        p
    }  
    
    return(P)
    returnType(double(0))
  })


#' @export
drecov <- nimbleFunction(
  run = function(d = double(0),
                 phi = double(1),
                 p = double(0),
                 lambda = double(0)) {
    
    # Compute probability of being recovered on the d-th
    # occasion following release.
    
    if (d == 1) {
      # Recovered on subsequent occasion
      P <- (1 - phi[1]) * lambda
    }
    else{
      # Recovered later than subsequent occasion
      P <- prod(phi[1:(d - 1)]) * (1 - phi[d]) *
        (1 - p)^(d - 1) *
        lambda
    }  
    
    return(P)
    returnType(double(0))
  })

#' @export
devade <- nimbleFunction(
  run = function(d = double(0),
                 phi = double(1),
                 p = double(0),
                 lambda = double(0)) {
    
    # Compute probability of being neither recaptured nor recovered by the 
    # d-th occasion following release (inclusive). 
    
    P <- phi[d] * (1 - p) +
      (1 - phi[d]) * (1 - lambda)
    
    if(d > 1){
      for(k in 1:(d-1)){
        P <- (1 - phi[d - k]) * (1 - lambda) +
          phi[d - k] * (1 - p) * P
      }
    }
    
    return(P)
    returnType(double(0))
  })

