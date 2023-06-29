##' Run one replicate of the simulation
##'
##' Run one replicate of the simulation
##' @title Run one replicate of the simulation
##'
##' @param pars Vector of simulation parameters including sim (scenario number), rep (replicate number wihin scenario), nind (number of individuals), ncap (number of individuals), model (either logit, scaled_logit, or generalized_logit), chains (number of chains), burnin (number of burnin samples), samples (number of samples to retain for computing posterior summary statistics)
##' @param models Character vector of models to run including binomial, trinomial, complete_mr, and/or complete_mrr
##' @param gendata If TRUE then generate data. Otherwise, reuse existing data.
##' @param verbose If true then print logging information to console
##'
##' @export
##' @importFrom MASS mvrnorm
##' @return
##' @author Simon Bonner
run_replicate <- function(pars, 
                          models = NULL, 
                          gendata = TRUE,
                          verbose = TRUE) {
  
  ## Setup exit condition
  on.exit(gc())

  if (verbose) {
    cat(date(), ": Setting parameters values...\n")
  }
  
  ## Default to running all models
  if(is.null(models)){
    models <- c("binomial","trinomial","complete_mr","complete_mrr")
  }

  ## 1) Covariate model
  d_params <- list(
    mu_init = .5, # Mean initial covariate value
    drift = rnorm(pars$ncap - 1, 0, .25), # Drifts,
    sigmasq = .5
  ) # Individual variance

  ## 2) Survival
  
  # Default parameters
  if(!("alpha" %in% names(pars)))
    pars$alpha <- -.5
  
  if(!("beta" %in% names(pars)))
    pars$beta <- 4

  if(!("lower" %in% names(pars)))
    pars$lower <- .3
  
  if(!("upper" %in% names(pars)))
    pars$upper <- .7
  
  # Construct parameter list
  s_params <- list(
    lower = ifelse(pars$model == "generalized_logit", pars$lower, 0),
    upper = ifelse(pars$model == "logit", 1, pars$upper),
    beta = c(pars$alpha, pars$beta)
  )

  ## 3) Capture and recovery
  if(!("p" %in% names(pars)))
    pars$p <- .6
  
  if(!("lambda" %in% names(pars)))
    pars$lambda <- .3
  
  c_params <- list(
    p = c(NA, rep(pars$p, pars$ncap - 1)),
    lambda = rep(pars$lambda, pars$ncap - 1)
  )

  ## Simulate Data for nind Individuals ##
  if (verbose) {
    cat(date(), ": Simulating data...\n")
  }
  
  datafile <- file.path("Data", paste0(pars$model, "_data_", pars$sim, "_", pars$rep, ".rds"))
  
  if(gendata){
    sim_data <- generate_data(pars$nind, pars$ncap, d_params, s_params, c_params)
    saveRDS(sim_data, datafile)
  }
  else{
    if(file.exists(datafile)){
      sim_data <- readRDS(datafile)  
    }
    else{
      stop("Cannot find input data file:",datafile,".\n")
    }
  } 
  
  on.exit(rm("sim_data"), add = TRUE, after = FALSE)
  
  ## Fit trinomial model ##
  if("trinomial" %in% models){
    if (verbose) {
      cat(date(), ": Running trinomial model...\n")
    }
    
    trinomial_time <- system.time(
      trinomial_out <- run_trinomial(indata = sim_data,
                                     model = pars$model,
                                     chains = pars$chains,
                                     burnin = pars$burnin,
                                     sampling = pars$samples))
    
    saveRDS(
      trinomial_out,
      file.path("Output", paste0("trinomial_", pars$model, "_out_", pars$sim, "_", pars$rep, ".rds"))
    )
    on.exit(rm("trinomial_out"), add = TRUE, after = FALSE)
  }
    
  ## Fit binomial model ##
  if("binomial" %in% models){
    if (verbose) {
      cat(date(), ": Running binomial model...\n")
    }
    
    binomial_time <- system.time(
      binomial_out <- run_binomial(indata = sim_data,
                                   model = pars$model,
                                   chains = pars$chains,
                                   burnin = pars$burnin,
                                   sampling = pars$samples))
    
    saveRDS(
      binomial_out,
      file.path("Output", paste0("binomial_", pars$model, "_out_", pars$sim, "_", pars$rep, ".rds"))
    )
    on.exit(rm("binomial_out"), add = TRUE, after = FALSE)
  }
  
  ## Fit complete MR data model ##
  if("complete_mr" %in% models){
    if (verbose) {
      cat(date(), ": Running complete MR data model...\n")
    }
    
    complete_mr_time <- system.time(
      complete_mr_out <- run_complete_mr_model(k = pars$ncap,
                                               indata = sim_data,
                                               model = pars$model,
                                               chains = pars$chains,
                                               burnin = pars$burnin,
                                               sampling = pars$samples))
    
    saveRDS(
      complete_mr_out,
      file.path("Output", paste0("complete_mr_", pars$model, "_out_", pars$sim, "_", pars$rep, ".rds"))
    )
    on.exit(rm("complete_mr_out"), add = TRUE, after = FALSE)
  }
      
  ## Fit complete MRR data model ##
  if("complete_mrr" %in% models){
    if (verbose) {
      cat(date(), ": Running complete MRR data model...\n")
    }
    
    complete_mrr_time <- system.time(
      complete_mrr_out <- run_complete_mrr_model(k = pars$ncap,
                                                 indata = sim_data,
                                                 model = pars$model,
                                                 chains = pars$chains,
                                                 burnin = pars$burnin,
                                                 sampling = pars$samples))
    
    saveRDS(
      complete_mrr_out,
      file.path("Output", paste0("complete_mrr_", pars$model, "_out_", pars$sim, "_", pars$rep, ".rds"))
    )
    on.exit(rm("complete_mrr_out"), add = TRUE, after = FALSE)
  }
  
  ## Extract summary ##
  if (verbose) {
    cat(date(), ": Computing summary...\n")
  }

  as_tibble.summary.mcmc <- function(summary, model, y) {
    cbind(summary[[1]], summary[[2]]) %>%
      as_tibble(rownames = "Parameter") %>%
      add_column(Model = model, y = y, .before = 1)
  }
  
  summary_list <- NULL
  
  if("trinomial" %in% models)
    summary_list$trinomial <- as_tibble(trinomial_out$summary, "Trinomial", 3)
  
  if("binomial" %in% models)
    summary_list$binomial <- as_tibble(binomial_out$summary, "Binomial", 4)
  
  if("complete_mr" %in% models)
    summary_list$complete_mr <- as_tibble(complete_mr_out$summary, "Complete MR", 1)
  
  if("complete_mrr" %in% models)
    summary_list$complete_mrr <- as_tibble(complete_mrr_out$summary, "Complete MRR", 2)
  
  
  summary_all <- summary_list %>%
    bind_rows()
  
  
  ## Extract run times ##
  time_list <- NULL
  
  if("binomial" %in% models)
    time_list$binomial <- c(binomial_time, Model = "Binomial")
  
  if("trinomial" %in% models)
    time_list$trinomial <- c(trinomial_time, Model = "Trinomial")
  
  if("complete_mr" %in% models)
    time_list$complete_mr <- c(complete_mr_time, Model = "Complete MR")
  
  if("complete_mrr" %in% models)
    time_list$complete_mrr <- c(complete_mrr_time, Model = "Complete MRR")
  
  time <- time_list %>%
    bind_rows()
  
  saveRDS(
    list(summary = summary_all,
         time = time),
    file.path("Output", paste0("summary_", pars$sim, "_", pars$rep, ".rds"))
  )
  on.exit(rm("summary_all"), add = TRUE, after = FALSE)
  
  ## Clean-up
  if (verbose) {
    cat(date(), ": Cleaning up...\n")
  }

  return(1)
}
