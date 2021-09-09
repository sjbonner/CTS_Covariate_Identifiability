read_mcmc <- function(file) {
  cat(file, "\n")
  tmp <- readRDS(file)
  tmp$mcmc %>%
    ggs()
}

load_all <- function(mymodel, mynind, myrep, dir = "Output", burnin = 2000, thin = 1) {
  sim <- filter(pars_mat, model == mymodel, nind == mynind, rep == myrep) %>%
    pull(sim)

  tibble(Model = c("binomial", "trinomial", "alt_trinomial", "complete")) %>%
    group_by(Model) %>%
    do(read_mcmc(file.path(dir, paste0(.$Model, "_", mymodel, "_out_", sim, "_", myrep, ".rds")))) %>%
    #filter(Chain == 1) %>% # See notes of February 18
    filter(Iteration > burnin) %>% # Remove burnin
    filter(Iteration %% thin == 0) %>% # Thin
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

posterior_density_grid <- function(data,
                                   truth,
                                   trans = "identity",
                                   scaled = TRUE,
                                   x.breaks = NULL,
                                   y.breaks = NULL,
                                   ylab = "Density") {
  
    myplot <- data %>%
        ggplot(aes(x = value, colour = Model))

    if(scaled)
        myplot <- myplot + stat_density(aes(y = ..scaled..),geom = "line", position = "identity")
    else
        myplot <- myplot + stat_density(geom = "line", position = "identity")

    myplot <- myplot +
        scale_x_continuous(name = "", n.breaks = x.breaks) + 
        scale_y_continuous(name = ylab, trans = trans, n.breaks = y.breaks) + 
        facet_grid(link_name_tex ~ Parameter_name, scales = "free", labeller = label_parsed) +
        geom_vline(data = truth, aes(xintercept = Value), lty = 2)

    myplot
}

posterior_caterpillar_grid <- function(data, truth){
  ## Compute summary statstics
  summ <- data %>%
    group_by(link_name_tex,Model,Parameter_name) %>%
    summarize(Mean = mean(value),
              Q2.5 = quantile(value,.025),
              Q25 = quantile(value,.25),
              Q75 = quantile(value,.75),
              Q97.5 = quantile(value,.975))

  ## Draw plot
  summ %>%
    ggplot(aes(x = Mean, y = Model)) +
    geom_point() + 
    facet_grid(link_name_tex ~ Parameter_name,
               scale = "free_x",labeller = label_parsed) +
    geom_errorbarh(aes(xmin = Q25, xmax = Q75), lwd = 2, height = 0) +
    geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, height = 0)) +
    geom_vline(data=truth,aes(xintercept = Value), lty = 2) +
    xlab("") +
    ylab("Model") +
    scale_y_discrete(limits = rev(levels(data$Model)))
}

scatter_plot <- function(data,
                         par1,
                         par2,
                         truth = NULL,
                         xlim = NULL,
                         ylim = NULL,
                         xlab = NULL,
                         ylab = NULL) {

  ## Create data 
  plot_df <- data %>%
    filter(Parameter %in% c(par1, par2)) %>%
    select(-Parameter_name) %>%
    spread(key = Parameter, value = value)

  myplot <- plot_df %>% 
    ggplot(aes_string(x = par1, y = par2)) +
    geom_point() +
    facet_grid(link_name ~ Model)

  if (!is.null(truth)) {
    myplot <- myplot +
      geom_point(
        data = spread(select(truth,-Parameter_name), key = Parameter, value = Value),
        colour = "red"
      )
  }
  
  if (!is.null(xlim)) {
    myplot <- myplot +
      xlim(c(0, 1))
  }

  if (!is.null(ylim)) {
    myplot <- myplot +
      ylim(c(0, 1))
  }

  if (!is.null(xlab)){
    myplot <- myplot +
      xlab(xlab)
  }

  if(!is.null(ylab)){
    myplot <- myplot +
      ylab(ylab)
  }

  print(myplot)
}

density_plot_2d <- function(data,
                         par1,
                         par2,
                         truth = NULL,
                         xlim = NULL,
                         ylim = NULL,
                         xlab = NULL,
                         ylab = NULL) {

  ## Create data 
  plot_df <- data %>%
    filter(Parameter %in% c(par1, par2)) %>%
    select(-Parameter_name) %>%
    spread(key = Parameter, value = value)

  myplot <- plot_df %>% 
    ggplot(aes_string(x = par1, y = par2)) +
    geom_density_2d(contour_var = "ndensity") +
    facet_grid(link_name ~ Model)

  if (!is.null(truth)) {
    myplot <- myplot + 
      geom_point(
        data = spread(select(truth,-Parameter_name), key = Parameter, value = Value),
        colour = "red"
      ) 
  }
  
  if (!is.null(xlim)) {
    myplot <- myplot +
      xlim(c(0, 1))
  }

  if (!is.null(ylim)) {
    myplot <- myplot +
      ylim(c(0, 1))
  }

  if (!is.null(xlab)){
    myplot <- myplot +
      xlab(xlab)
  }

  if(!is.null(ylab)){
    myplot <- myplot +
      ylab(ylab)
  }

  print(myplot)
}
