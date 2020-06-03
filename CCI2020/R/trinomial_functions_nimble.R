##' Run model with trinomial likelihood
##'
##' Run model with trinomial likelihood
##' @title 
##' @param indata Input data
##' @param model Model ("logit","scaled_logit", or "generalized_logit")
##' @param inits Initial values
##' @param pars Parameters to monitor 
##' @param coda_dir Directory for saving samples
##' @param chains Number of chains
##' @param burnin Length of burnin period
##' @param sampling Length of sampling period
##' @export
##' @return 
##' @author Simon Bonner
run_trinomial <- function(indata,
                          model = c("logit", "scaled", "generalized"),
                          inits = NULL,
                          pars = NULL,
                          coda_dir,
                          chains = 3,
                          burnin = 1000,
                          sampling = 1000) {

  ## 1) Format data

  ## Extract basic parameters
  Nind <- length(indata)
  Ncap <- length(indata[[1]]$ch)

  ## Extract matrices of capture histories and covariates
  ch <- t(sapply(indata, "[[", "ch"))
  z <- t(sapply(indata, "[[", "z"))

  ## Identify releases
  tmp <- which(ch[, -Ncap] == 1, arr.ind = TRUE)
  trinomial_data <- list(
    nrelease = nrow(tmp),
    z = sapply(1:nrow(tmp), function(i) z[tmp[i, 1], tmp[i, 2]]),
    outcome = sapply(1:nrow(tmp), function(i) ch[tmp[i, 1], tmp[i, 2] + 1]) + 1
  )

  ## Set initial values
  if(is.null(inits)){
    inits <- list(
      p = .5,
      lambda = .5,
      beta.phi = c(0, 0)
    )
    
    if (model == "scaled_logit") {
      inits$upper <- 1
    }
    
    if (model == "generalized_logit") {
      inits$lower <- 0
      inits$upper <- 1
    }
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
  model_file <- file.path("JAGS", paste0("trinomial_", model, ".R"))
  coda_file <- file.path(coda_dir, paste0("trinomial_", model, ".rds"))

  nimble_model <- readBUGSmodel(model_file,
    data = trinomial_data,
    inits = inits)

  ## Burn-in
  time <- system.time(
    trinomial_mcmc <- nimbleMCMC(
      model = nimble_model,
      monitors = pars,
      nchains = chains,
      nburnin = 0,
      niter = burnin + sampling,
      samplesAsCodaMCMC = TRUE
    )
  )

  ## Save output
  saveRDS(trinomial_mcmc, file = coda_file)

  ## Return summaries
  list(
    summary = summary(trinomial_mcmc),
    ess = effectiveSize(trinomial_mcmc),
    diag = gelman.diag(trinomial_mcmc),
    mcmc = trinomial_mcmc,
    time = time
  )
}
