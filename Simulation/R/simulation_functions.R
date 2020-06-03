run_replicate <- function(pars, verbose = TRUE) {

  ## Extract Parameters ##
  attach(pars)

  if (verbose) {
    cat(date(), ": Setting parameters values...\n")
  }

  ## 1) Covariate model
  d_params <- list(
    mu_init = .5, # Mean initial covariate value
    drift = rnorm(ncap - 1, 0, .25), # Drifts,
    sigmasq = .5
  ) # Individual variance

  ## 2) Survival
  s_params <- list(
    lower = ifelse(model == "generalized_logit", .3, 0),
    upper = ifelse(model == "logit", 1, .7),
    beta = c(-.5, 4)
  )

  ## 3) Capture and recovery
  c_params <- list(
    p = c(NA, rep(.6, ncap - 1)),
    lambda = rep(.3, ncap - 1)
  )

  ## Simulate Data for nind Individuals ##
  if (verbose) {
    cat(date(), ": Simulating data...\n")
  }

  sim_data <- generate_data(nind, ncap, d_params, s_params, c_params)
  saveRDS(sim_data, file.path("Data", paste0(model, "_data_", sim, "_", rep, ".rds")))

  ## Fit trinomial model ##
  if (verbose) {
    cat(date(), ": Running trinomial model...\n")
  }

  ## trinomial_out <- run_trinomial(indata = sim_data,
  ##                                model = model,
  ##                                coda_dir = "Trinomial",
  ##                                chains = chains,
  ##                                burnin = burnin,
  ##                                sampling = samples)
  ## saveRDS(
  ##   trinomial_out,
  ##   file.path("Output", paste0("trinomial_", model, "_out_", sim, "_", rep, ".rds"))
  ## )

  ## Fit binomial model ##
  if (verbose) {
    cat(date(), ": Running binomial model...\n")
  }

  binomial_out <- run_binomial(indata = sim_data,
                               model = model,
                               coda_dir = "Binomial",
                               chains = chains,
                               burnin = burnin,
                               sampling = samples)
  saveRDS(
    binomial_out,
    file.path("Output", paste0("binomial_", model, "_out_", sim, "_", rep, ".rds"))
  )

  ## Fit alternative trinomial model ##
  if (verbose) {
    cat(date(), ": Running alternative trinomial model...\n")
  }

  alt_trinomial_out <- run_trunc_model(k = 2,
                                       indata = sim_data,
                                       model = model,
                                       coda_dir = "Alt_Trinomial",
                                       chain = chains,
                                       burnin = burnin,
                                       sampling = samples)
  saveRDS(
    alt_trinomial_out,
    file.path("Output", paste0("alt_trinomial_", model, "_out_", sim, "_", rep, ".rds"))
  )

  ## Fit complete data model ##
  if (verbose) {
    cat(date(), ": Running complete data model...\n")
  }

  complete_out <- run_trunc_model(k = ncap,
                                  indata = sim_data,
                                  model = model,
                                  coda_dir = "Complete_Data",
                                  chains = chains,
                                  burnin = burnin,
                                  sampling = samples)
  saveRDS(
    complete_out,
    file.path("Output", paste0("complete_", model, "_out_", sim, "_", rep, ".rds"))
  )

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
    as_tibble(alt_trinomial_out$summary, "Alt. Trinomial", 2),
    as_tibble(complete_out$summary, "Complete", 1)
  )

  saveRDS(
    summary_all,
    file.path("Output", paste0("summary_", sim, "_", rep, ".rds"))
  )

  ## Clean-up
  if (verbose) {
    cat(date(), ": Cleaning up...\n")
  }

  rm(list = c(
    "sim_data",
    "trinomial_out",
    "binomial_out",
    "complete_out",
    "summary_all"
  ))

  gc()

  return(summary_all)
}
