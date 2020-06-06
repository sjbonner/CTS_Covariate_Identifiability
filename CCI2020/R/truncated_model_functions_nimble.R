##' Run model with truncated likelihood
##'
##' Run model with truncated likelihood
##' @title Run model with truncated likelihood
##' @param k Truncation factor
##' @param indata Input data
##' @param model Model ("logit","scaled_logit", or "generalized_logit")
##' @param inits Initial values
##' @param pars Parameters to monitor 
##' @param coda_dir Directory for saving samples
##' @param chains Number of chains
##' @param burnin Length of burnin period
##' @param sampling Length of sampling period
##' @export
##' @import spark
##' @import nimble
##' @import coda
##' @return 
##' @author Simon Bonner
run_trunc_model <- function(k,
                            indata,
                            model = c("logit", "scaled", "generalized"),
                            inits = NULL,
                            pars = NULL,
                            coda_dir,
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

  ## 2) Construct drecap
  drecap <- ifelse(
    trunc_data$recapture > 0,
    trunc_data$recapture - trunc_data$release,
    pmin(k, Ncap - trunc_data$release)
  )

  ## 3) Format covariate matrix
  Z <- matrix(NA, nrow = trunc_data$nrelease, ncol = k + 1)

  for (i in 1:trunc_data$nrelease) {
    Z[i, 1] <- trunc_data$other[trunc_data$ind[i], trunc_data$release[i]]
    if (trunc_data$recapture[i] > 0) {
      Z[i, drecap[i] + 1] <- trunc_data$other[trunc_data$ind[i], trunc_data$recapture[i]]
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
    drecap = drecap[index],
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

  ## Identify parameters to monitor
  if(is.null(pars)){
    pars <- c("beta.phi", "p", "mu.z", "sigma.z")
    
    if (model == "scaled_logit") {
      pars <- c(pars, "upper")
    }
    if (model == "generalized_logit") {
      pars <- c(pars, "lower", "upper")
    }
  }

  ## Run model
  model_file <- system.file("JAGS",paste0("cts_cov_jags_", model, ".R"),
                            package="CCI2020")

  coda_file <- file.path(coda_dir, paste0("trunc_coda_", k, "_", model, ".rds"))

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

  ## Save output
  saveRDS(trunc_jags_mcmc, file = coda_file)

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
