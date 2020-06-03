read_mcmc <- function(file) {
  cat(file, "\n")
  tmp <- readRDS(file)
  tmp$mcmc %>%
    ggs()
}

load_all <- function(mymodel, mynind, myrep, burnin = 2000, thin = 10) {
  sim <- filter(pars_mat, model == mymodel, nind == mynind, rep == myrep) %>%
    pull(sim)

  tibble(Model = c("binomial", "trinomial", "alt_trinomial", "complete")) %>%
    group_by(Model) %>%
    do(read_mcmc(file.path("Output", paste0(.$Model, "_", mymodel, "_out_", sim, "_", myrep, ".rds")))) %>%
    filter(Chain == 1) %>% # See notes of February 18
    filter(Iteration > burnin) %>% # Remove burnin
    filter(Iteration %% thin == 1) %>% # Thin by factor of 10
    add_column(model = mymodel, .before = 1)
}

posterior_density <- function(data, par, truth, xlim = NULL) {
  myplot <- data %>%
    filter(Parameter == par) %>%
    ggplot(aes(x = value)) +
    geom_density() +
    # geom_histogram(aes(y = ..density..),bins = 100) +
    facet_wrap(vars(Model)) +
    geom_vline(
      data = filter(truth, Parameter == par),
      aes(xintercept = Value), lty = 2
    )

  if (!is.null(xlim)) {
    myplot <- myplot + xlim(xlim)
  }

  print(myplot)
}

posterior_density_grid <- function(data, pars, truth) {
  data %>%
    filter(Parameter %in% pars) %>%
    ggplot(aes(x = value, colour = Model)) +
    geom_density() +
    facet_grid(model ~ Parameter, scales = "free") +
    geom_vline(
      data = filter(truth, Parameter %in% pars),
      aes(xintercept = Value), lty = 2
    )
}

scatter_plot <- function(data, par1, par2, truth, xlim = NULL, ylim = NULL) {
  myplot <- data %>%
    filter(Parameter %in% c(par1, par2)) %>%
    spread(key = Parameter, value = value) %>%
    ggplot(aes_string(x = par1, y = par2)) +
    geom_point() +
    geom_point(
      data = spread(truth, key = Parameter, value = Value),
      colour = "red"
    ) +
    facet_grid(model ~ Model)

  if (is.null(xlim)) {
    myplot <- myplot +
      xlim(c(0, 1))
  }

  if (is.null(ylim)) {
    myplot <- myplot +
      ylim(c(0, 1))
  }

  print(myplot)
}
