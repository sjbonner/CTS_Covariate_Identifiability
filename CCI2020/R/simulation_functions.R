##' Run one replicate of the simulation
##'
##' Run one replicate of the simulation
##' @title Run one replicate of the simulation
##' @param pars Vector of simulation parameters including sim (scenario number), rep (replicate number wihin scenario), nind (number of individuals), ncap (number of individuals), model (either logit, scaled_logit, or generalized_logit), chains (number of chains), burnin (number of burnin samples), samples (number of samples to retain for computing posterior summary statistics)
##' @param verbose If true then print logging information to console
##' @export
##' @import MASS
##' @return
##' @author Simon Bonner
run_replicate <- function(pars, verbose = TRUE) {
  
  ## Setup exit condition
  on.exit(gc())

  if (verbose) {
    cat(date(), ": Setting parameters values...\n")
  }

  ## 1) Covariate model
  d_params <- list(
    mu_init = .5, # Mean initial covariate value
    drift = rnorm(pars$ncap - 1, 0, .25), # Drifts,
    sigmasq = .5
  ) # Individual variance

  ## 2) Survival
  s_params <- list(
    lower = ifelse(pars$model == "generalized_logit", .3, 0),
    upper = ifelse(pars$model == "logit", 1, .7),
    beta = c(-.5, 4)
  )

  ## 3) Capture and recovery
  c_params <- list(
    p = c(NA, rep(.6, pars$ncap - 1)),
    lambda = rep(.3, pars$ncap - 1)
  )

  ## Simulate Data for nind Individuals ##
  if (verbose) {
    cat(date(), ": Simulating data...\n")
  }
  
  sim_data <- generate_data(pars$nind, pars$ncap, d_params, s_params, c_params)
  saveRDS(sim_data, file.path("Data", paste0(pars$model, "_data_", pars$sim, "_", pars$rep, ".rds")))
  on.exit(rm("sim_data"), add = TRUE, after = FALSE)
  
  ## Fit trinomial model ##
  if (verbose) {
    cat(date(), ": Running trinomial model...\n")
  }

  trinomial_out <- run_trinomial(indata = sim_data,
                                 model = pars$model,
                                 coda_dir = "Trinomial",
                                 chains = pars$chains,
                                 burnin = pars$burnin,
                                 sampling = pars$samples)
  saveRDS(
    trinomial_out,
    file.path("Output", paste0("trinomial_", pars$model, "_out_", pars$sim, "_", pars$rep, ".rds"))
  )
  on.exit(rm("trinomial_out"), add = TRUE, after = FALSE)
  

  ## Fit binomial model ##
  if (verbose) {
    cat(date(), ": Running binomial model...\n")
  }

  binomial_out <- run_binomial(indata = sim_data,
                               model = pars$model,
                               coda_dir = "Binomial",
                               chains = pars$chains,
                               burnin = pars$burnin,
                               sampling = pars$samples)
  saveRDS(
    binomial_out,
    file.path("Output", paste0("binomial_", pars$model, "_out_", pars$sim, "_", pars$rep, ".rds"))
  )
  on.exit(rm("binomial_out"), add = TRUE, after = FALSE)
  
  ## Fit complete data model ##
  if (verbose) {
    cat(date(), ": Running complete data model...\n")
  }

  complete_out <- run_trunc_model(k = pars$ncap,
                                  indata = sim_data,
                                  model = pars$model,
                                  coda_dir = "Complete_Data",
                                  chains = pars$chains,
                                  burnin = pars$burnin,
                                  sampling = pars$samples)
  saveRDS(
    complete_out,
    file.path("Output", paste0("complete_", pars$model, "_out_", pars$sim, "_", pars$rep, ".rds"))
  )
  on.exit(rm("complete_out"), add = TRUE, after = FALSE)
  ## Extract summary ##
  if (verbose) {
    cat(date(), ": Computing summary...\n")
  }

  as_tibble.summary.mcmc <- function(summary, model, y) {
    cbind(summary[[1]], summary[[2]]) %>%
      as_tibble(rownames = "Parameter") %>%
      add_column(Model = model, y = y, .before = 1)
  }
  
  summary_all <- bind_rows(
    as_tibble(trinomial_out$summary, "Trinomial", 3),
    as_tibble(binomial_out$summary, "Binomial", 4),
    as_tibble(complete_out$summary, "Complete", 1)
  )

  saveRDS(
    summary_all,
    file.path("Output", paste0("summary_", pars$sim, "_", pars$rep, ".rds"))
  )
  on.exit(rm("summary_all"), add = TRUE, after = FALSE)
  
  ## Clean-up
  if (verbose) {
    cat(date(), ": Cleaning up...\n")
  }

  return(1)
}
