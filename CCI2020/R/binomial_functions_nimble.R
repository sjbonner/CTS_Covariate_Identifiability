##' Run model with binomial likelihood
##'
##' Run model with binomial likelihood
##' @title Run model with binomial likelihood
##' @param indata Input data
##' @param model Model ("logit","scaled_logit", or "generalized_logit")
##' @param inits Initial values
##' @param pars Parameters to monitor 
##' @param coda_dir Directory for saving samples
##' @param chains Number of chains
##' @param burnin Length of burnin period
##' @param sampling Length of sampling period
##' @export
##' @import nimble
##' @import coda
##' @return 
##' @author Simon Bonner
run_binomial <- function(indata,
                         model=c("logit","scaled","generalized"),
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
  ch <- ifelse(t(sapply(indata,"[[","ch")) == 1, 1, 0)
  z <- t(sapply(indata,"[[","z"))

  ## Identify releases
  tmp <- which(ch[,-Ncap] == 1, arr.ind = TRUE)

  binomial_data <- list(
    nrelease = nrow(tmp),
    z = sapply(1:nrow(tmp), function(i) z[tmp[i,1],tmp[i,2]]),
    outcome = sapply(1:nrow(tmp), function(i) ch[tmp[i,1],tmp[i,2] + 1]) + 1)

  ## Set initial values
  if(is.null(inits)){
    inits <- list(
      p = .5,
      beta.phi = c(0,0))

    if(model == "scaled_logit"){
      inits$upper <- 1
    }
    
    if(model == "generalized_logit"){
      inits$lower <- 0
      inits$upper <- 1
    }
  }
  
  ## Identify parameters to monitor
  if(is.null(pars)){
    pars <- c("beta.phi","p")
    
    if (model == "scaled_logit") {
      pars <- c(pars, "upper")
    }
    if (model == "generalized_logit") {
      pars <- c(pars, "lower", "upper")
    }
  }

  ## Run model
  model_file <- system.file("JAGS",paste0("binomial_",model,".R"),
                            package="CCI2020")
  coda_file <- file.path(coda_dir, paste0("binomial_",model,".rds"))

  ## Run model
  nimble_model <- readBUGSmodel(model_file,
                                data = binomial_data,
                                inits = inits)
  
  time <- system.time(
    binomial_mcmc <- nimbleMCMC(
      model = nimble_model,
      monitors = pars,
      nchains = chains,
      nburnin = 0,
      niter = burnin + sampling,
      samplesAsCodaMCMC = TRUE)
  )
  
  ## Save output
  saveRDS(binomial_mcmc, file = coda_file)

  ## Return summaries
  output <- list(summary = summary(window(binomial_mcmc, start = burnin + 1)),
       mcmc = binomial_mcmc,
       time = time)

  if (chains > 1) {
    output$ess <- effectiveSize(window(binomial_mcmc, start = burnin + 1))
    output$diag <- gelman.diag(window(binomial_mcmc, start = burnin + 1))
  }

  return(output)
}
