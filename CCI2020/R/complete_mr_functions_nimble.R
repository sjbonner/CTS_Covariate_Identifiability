##' Run truncated complete data model with recaptures only
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
run_complete_mr_model <- function(k,
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
    chmat = ifelse(ch == 1, 1, 0),
    other = z,
    freq = rep(1, length(indata))
  )

  ## 1) Truncate data
  trunc_data <- spark(indata = spark_data, k = k)

  ## 2) Construct delta vector (occasions between release and recapture/recovery)
  delta <- ifelse(
    trunc_data$recapture > 0,
    trunc_data$recapture - trunc_data$release,
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

  ## 4) Order data so that recaptures come first
  tmp <- which(trunc_data$recapture > 0)
  index <- c(tmp, (1:trunc_data$nrelease)[-tmp])
  nrecapture <- length(tmp)

  ## 5) Construct data list
  trunc_jags_data <- list(
    nocc = Ncap,
    nrecapture = nrecapture,
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

  ## Identify parameters to monitor
  if(is.null(pars)){
    pars <- c("beta.phi", "p")
    
    if (model == "scaled_logit") {
      pars <- c(pars, "upper")
    }
    if (model == "generalized_logit") {
      pars <- c(pars, "lower", "upper")
    }
  }
  
  ## Run model
  model_file <- system.file("Nimble",paste0("complete_mr_", model, ".R"),
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
