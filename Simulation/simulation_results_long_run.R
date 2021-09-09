## Load packages
library(tidyverse)
library(ggmcmc)
library(latex2exp)
library(coda)

## Source code
source("R/simulation_results_functions.R")

## Load parameter matrix
pars_mat <- read_csv("parameter_matrix_long.csv")

## Parameters common to simulation scenarios
pars <- c("p","lambda","beta.phi[1]","beta.phi[2]","upper","lower")

## Blank plotting window
x11(width=6,height=4)

## Pretty names
par_names <- TeX(paste0("$",
                        c("p",
                          "\\lambda",
                          "\\alpha",
                          "\\beta",
                          "\\gamma_1",
                          "\\gamma_0"),
                        "$"))

link_names <- c("Logit",
                "Scaled Logit",
                "Generalized Logit")

model_names <- c("Alt. Trinomial",
                 "Binomial",
                 "Complete",
                 "Trinomial")

## Define data set with true parameter values
truth <- tibble(Parameter =pars,
                Parameter_name =
                  factor(pars,
                         levels = pars,
                         labels = par_names),
                Value = c(.6, .3, -.5, 4, .7, .3))

## Parameters defining simulation
nind <- 10000
rep <- "long"

## Load complete output
all_mcmc_logit <- load_all("logit", nind, rep, dir = "Output", burnin = 10000, thin = 1)
all_mcmc_generalized <- load_all("generalized_logit", nind, rep, dir = "Output", burnin = 10000, thin = 1)
all_mcmc_scaled <- load_all("scaled_logit", nind, rep, dir = "Output", burnin = 10000, thin = 1)

all_mcmc <- bind_rows(all_mcmc_logit,
                      all_mcmc_scaled,
                      all_mcmc_generalized)

rm(all_mcmc_logit,all_mcmc_generalized,all_mcmc_scaled)
gc()

## Rename model to link (it was a silly choice)
all_mcmc <- all_mcmc %>%
  rename(Link = model)

## Convert chain to integer
all_mcmc <- all_mcmc %>%
  mutate(Chain = as.factor(Chain))

## Keep only parameters of interest and remove alternative trinomial
all_mcmc <- all_mcmc %>%
  filter(Parameter %in% pars)

## Add pretty names
all_mcmc$Parameter_name <- factor(all_mcmc$Parameter,
                                  levels = pars,
                                  labels = par_names)

all_mcmc$link_name <- factor(all_mcmc$Link,
                         levels = c("logit",
                                    "scaled_logit",
                                    "generalized_logit"),
                         labels = link_names)

all_mcmc$link_name_tex <- factor(all_mcmc$Link,
                         levels = c("logit",
                                    "scaled_logit",
                                    "generalized_logit"),
                         labels = TeX(link_names))

all_mcmc$Model <- factor(all_mcmc$Model,
                         levels = c("alt_trinomial",
                                    "binomial",
                                    "complete",
                                    "trinomial"),
                         labels = model_names)

## Individual Density Plots
## 1) Logit
filter(all_mcmc,link_name == "Logit",
       Model != "Alt. Trinomial") %>%
  posterior_density_grid(truth)

## 2) Scaled Logit

## Construct limiting posterior density for upper bound on survival
rho <- truth %>%
  select(Parameter,Value) %>%
  spread(key = Parameter, value = Value) %>%
  mutate(rho0 = lower * p,
         rho1 = upper *p) %>%
  select(rho0,rho1) %>%
  as.numeric()

limiting_posterior_1 <- tibble(Model = "Binomial",
                               x = seq(0,1,length = 100),
                               pi = (x > rho[2]) * (-log(rho[2])/(2 * x))) %>%
  crossing(tibble(Link = c("scaled_logit", "generalized_logit"),
                  link_name = c("Scaled Logit", "Generalized Logit"),
                  link_name_tex = factor(Link,
                                         levels = c("logit",
                                                    "scaled_logit",
                                                    "generalized_logit"),
                                         labels = TeX(link_names)))) %>%
  crossing(tibble(Parameter = c("upper","p"),
                  Parameter_name = pull(filter(truth, Parameter %in% c("upper","p")),"Parameter_name")))

limiting_posterior_2 <- tibble(Link = "generalized_logit",
                               link_name = "Generalized Logit",
                               link_name_tex = factor(Link,
                                                      levels = c("logit",
                                                        "scaled_logit",
                                                        "generalized_logit"),
                                                      labels = TeX(link_names)),
                               Model = "Binomial",
                               x = seq(0,1,length= 100),
                               pi = (x > rho[1]) * (x < rho[1]/rho[2])/(rho[1]/rho[2]-rho[1])/4,
                               Parameter = "lower",
                               Parameter_name = pull(filter(truth, Parameter %in% c("lower")),"Parameter_name"))

limiting_posterior <- bind_rows(limiting_posterior_1,
                                limiting_posterior_2)
  
filter(all_mcmc,link_name == "Scaled Logit") %>%
  posterior_density_grid(truth) +
  geom_line(data = filter(limiting_posterior, link_name == "Scaled Logit"),
            mapping = aes(x = x, y = pi),
            lty = 2)

## 3) Generalized Logit
filter(all_mcmc,link_name == "Generalized Logit") %>%
  posterior_density_grid(truth) +
  geom_line(data = filter(limiting_posterior, link_name == "Generalized Logit"),
            mapping = aes(x = x, y = pi),
            lty = 2)

## Combined density plot
filter(all_mcmc,
       Model != "Alt. Trinomial") %>%
posterior_density_grid(truth,
                                        #trans = "sqrt",
                       scaled = TRUE,
                       x.breaks = 3,
                       y.breaks = 3,
                       ylab = "Density (Scaled)") +
    geom_line(data = limiting_posterior,
            mapping = aes(x = x, y = pi),
            lty = 2) + 
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90))

dev.copy2pdf(file = "Figures/combined_density_plot.pdf")

## Caterpillar
all_mcmc %>%
    filter(Model != "Alt. Trinomial") %>%
    posterior_caterpillar_grid(truth)

dev.copy2pdf(file = "Figures/combined_caterpillar_plot.pdf")

## Scatterplots and density plots of sampled values

## 1) p vs upper
id_region <- tibble(Model = "Binomial",
                    p = seq(rho[2],1,length = 100),
                    lower = rho[1]/p,
                    upper = rho[2]/p)
all_mcmc %>%
    filter(Link != "logit",
           Model != "Alt. Trinomial",
           Iteration %% 100 == 0) %>%
    scatter_plot("p", "upper", truth, c(0,1), c(0,1),TeX("$p$"), TeX("$\\gamma_1$")) +
    geom_line(aes(x=p, y = upper), data = id_region, colour = "red")

dev.copy2pdf(file = "Figures/p_vs_upper_scatterplot.pdf")

all_mcmc %>%
    filter(Link != "logit",
           Model != "Alt. Trinomial") %>%
    density_plot_2d("p", "upper", truth, c(0,1), c(0,1),TeX("$p$"), TeX("$\\gamma_1$")) +
    geom_line(aes(x=p, y = upper), data = id_region, colour = "red")

dev.copy2pdf(file = "Figures/p_vs_upper_contourplot.pdf")

## 2) p vs lower

x11(width = 6, height = 2.25)

all_mcmc %>%
    filter(Link == "generalized_logit",
           Model != "Alt. Trinomial",
           Iteration %% 100 == 0) %>%
    scatter_plot("p", "lower", truth, c(0,1), c(0,1),TeX("$p$"), TeX("$\\gamma_0$")) +
    geom_line(aes(x=p, y = lower), data = id_region, colour = "red")

dev.copy2pdf(file = "Figures/p_vs_lower_scatterplot.pdf")

all_mcmc %>%
    filter(Link == "generalized_logit",
           Model != "Alt. Trinomial") %>%
    density_plot_2d("p", "lower", truth, c(0,1), c(0,1),TeX("$p$"), TeX("$\\gamma_0$")) +
    geom_line(aes(x=p, y = lower), data = id_region, colour = "red")

dev.copy2pdf(file = "Figures/p_vs_lower_contourplot.pdf")

## ## 3) beta_0 vs beta_
## filter(all_mcmc, Link == "generalized_logit") %>%
##   density_plot_2d("beta.phi[1]", "beta.phi[2]", truth, NULL, NULL,TeX("$\\beta_1$"), TeX("$\\beta_2$"))

## if(beamer)
##   dev.copy2pdf(file = "Figures/beta_1_vs_beta_2_contourplot_beamer.pdf")
## if(!beamer)
##   dev.copy2pdf(file = "Figures/beta_1_vs_beta_2_contourplot.pdf")

## tmp <- all_mcmc %>%
##   select(-Parameter_name) %>%
##   filter(Parameter %in% c("beta.phi[1]","beta.phi[2]")) %>%
##   spread(key=Parameter, value = value)

## tmp %>%
##   ggplot(aes(x = `beta.phi[1]`, y = `beta.phi[2]`)) +
##   geom_point() +
##     geom_point(data = spread(select(truth,-Parameter_name),
##                              key = Parameter, value = Value),
##                colour = "red") + 
##   facet_wrap(vars(Model))

## Compute Gelman diagnostics and effective sample sizes

link <- "generalized_logit"
model <- "Binomial"

coda <- lapply(1:3,function(i){
  all_mcmc %>%
  filter(Link == link,
         Model == model,
         Chain == i) %>%
  ungroup() %>%
  select(Iteration, Parameter, value) %>%
  spread(key = Parameter, value = value) %>%
  select(-Iteration) %>%
    as.mcmc()}) %>%
  as.mcmc.list()
  
(ess <- effectiveSize(coda))
100 * ess/300000

gelman.diag(coda)
