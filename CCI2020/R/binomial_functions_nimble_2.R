##' Run model with binomial likelihood (alternate input format)
##'
##' Run model with binomial likelihood (alternate input format)
##' @title Run model with binomial likelihood (alternate input format)
##' @param ch Matrix of capture histories
##' @param Zphi Array of (partially) observed covariate values for survival
##' @param Zp Array of (partially) observed covariate values for capture
##' @param intercepts If TRUE then add intercepts (matrices of 1s) to both Zphi and Zp
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
##' @import abind
##' @return 
##' @author Simon Bonner
run_binomial_2 <- function(ch = NULL,
                         Zphi = NULL,
                         Zp = NULL,
                         intercepts = TRUE,
                         model=c("logit","scaled","generalized"),
                         inits = NULL,
                         pars = NULL,
                         coda_dir,
                         chains = 3,
                         burnin = 1000,
                         sampling = 1000) {

  ## Extract basic values
  Nind <- nrow(ch)
  Ncap <- ncol(ch)

  ## Identify releases
  tmp <- which(ch[,-Ncap] == 1, arr.ind = TRUE)

  ## Add intercepts to design matrices
  if(intercepts){
    Zphi <- abind(matrix(1, nrow = nrow(Zphi), ncol = ncol(Zphi)),
                  Zphi,
                  along = 3)
    
    Zp <- abind(matrix(1, nrow = nrow(Zphi), ncol = ncol(Zphi)),
                Zp,
                along = 3)
  }

  ## Format data for nimble
  binomial_data <- list(
    nrelease = nrow(tmp),
    occasion = tmp[,2],
    Zphi = t(sapply(1:nrow(tmp), function(i) Zphi[tmp[i,1],tmp[i,2],])),
    nZphi = dim(Zp)[3],
    Zp = t(sapply(1:nrow(tmp), function(i) Zp[tmp[i,1],tmp[i,2]+1,])),
    nZp = dim(Zp)[3],
    outcome = sapply(1:nrow(tmp), function(i) ch[tmp[i,1],tmp[i,2] + 1]) + 1)

  ## Set initial values
  if(is.null(inits)){
    inits <- list(beta.p = rep(0,binomial_data$nZp),
                  beta.phi = rep(0, binomial_data$nZphi))
    
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
    pars <- c("beta.p","beta.phi")
    
    if (model == "scaled_logit") {
      pars <- c(pars, "upper")
    }
    if (model == "generalized_logit") {
      pars <- c(pars, "lower", "upper")
    }
  }

  ## Run model
  model_file <- system.file("JAGS",paste0("binomial_",model,"_2.R"),
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
