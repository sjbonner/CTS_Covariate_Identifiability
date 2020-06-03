## Load packages
library(tidyverse)
library(ggmcmc)

## Source code
source("R/simulation_results_functions.R")

## Load parameter matrix
pars_mat <- read_csv("parameter_matrix.csv")

## Define data set with true parameter values
truth <- tibble(Parameter = c("p","lambda","beta.phi[1]","beta.phi[2]","lower","upper"),
                Value = c(.6, .3, -.5, 4, .3, .7))

## Load all summaries
## load_summary <- function(sim, rep){
##   file <- file.path("Output",paste0("summary_",sim,"_",rep,".rds"))
##   readRDS(file)           
## }

tmp <- pars_mat %>%
  rowid_to_column(var = "ID") %>%
  select(ID,sim,rep,nind,model) %>%
  mutate(file = file.path("Output",paste0("summary_",sim,"_",rep,".rds")),
         exists = file.exists(file)) 

filter(tmp, !exists)

summaries <- tmp %>%
  filter(exists)%>%
  group_by(sim,rep,nind,model) %>%
  do(
    readRDS(.$file)           
  )

(tmp <- summaries %>%
   filter(Parameter == "beta.phi[1]") %>%
   group_by(model,Model, nind) %>%
   summarize(Count = n()))

## Summarize across simulations
results <- summaries %>%
  group_by(sim, model, Model, Parameter, nind) %>%
  summarize(MeanSD = mean(SD))


## Plot mean standard deviation versus sample size
pars <- c("p","beta.phi[1]","beta.phi[2]","lower","upper")

filter(results, Parameter %in% pars) %>%
  ggplot(aes(x = nind, y=MeanSD, colour=Model)) +
  geom_point() +
  geom_line() + 
  geom_hline(yintercept = 0) +
  facet_grid(Parameter ~ model, scale = "free") +
  xlab("Sample Size")

## Examine results for one simulated data set per model
nind <- 10000
rep <- 2

## 1) Logit
all_mcmc_logit <- load_all("logit", nind, rep)
pars <- c("p","beta.phi[1]","beta.phi[2]","lower","upper")
posterior_density_grid(all_mcmc_logit, pars, truth)

## 2) Scaled Logit
all_mcmc_scaled <- load_all("scaled_logit", nind, rep)
pars <- c("p","beta.phi[1]","beta.phi[2]","lower","upper")
posterior_density_grid(all_mcmc_scaled, pars, truth)

## 3) Generalized Logit
all_mcmc_generalized <- load_all("generalized_logit", nind, rep)
pars <- c("p","beta.phi[1]","beta.phi[2]","lower","upper")
posterior_density_grid(all_mcmc_generalized, pars, truth)

pars <- c("p","beta.phi[1]","beta.phi[2]","lower","upper")
all_mcmc <- bind_rows(all_mcmc_logit,
                      all_mcmc_scaled,
                      all_mcmc_generalized) %>%
    filter(Parameter %in% pars)
  
all_mcmc$Model <- factor(all_mcmc$Model,
                         levels = c("binomial","trinomial","alt_trinomial","complete"))

all_mcmc$model <- factor(all_mcmc$model,
                         levels = c("logit","scaled_logit","generalized_logit"))

## Combined density plot
posterior_density_grid(all_mcmc, pars, truth)

## Caterpillar

all_mcmc %>%
  group_by(model,Model,Parameter) %>%
  summarize(Mean = mean(value),
            Q2.5 = quantile(value,.025),
            Q25 = quantile(value,.25),
            Q75 = quantile(value,.75),
            Q97.5 = quantile(value,.975)) %>%
  ggplot(aes(x = Mean, y = Model)) +
  geom_point() + 
  facet_grid(model ~ Parameter, scale = "free_x") +
  geom_errorbarh(aes(xmin = Q25, xmax = Q75), lwd = 2, height = 0) +
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, height = 0)) +
  geom_vline(data = filter(truth, Parameter %in% pars),
             aes(xintercept = Value), lty = 2) +
  xlab("")


## Scatterplots of sampled values

## 1) p vs upper
filter(all_mcmc, model != "logit") %>%
  scatter_plot("p", "upper", truth, c(0,1), c(0,1))

## 2) p vs lower
filter(all_mcmc, model == "generalized_logit") %>%
  scatter_plot("p", "lower", truth, c(0,1), c(0,1))

## 3) beta_0 vs beta_1
all_mcmc %>%
  filter(Chain == 1) %>% # See notes of February 18
  filter(Iteration > 2000 ) %>%    # Remove burn-in
  filter(Iteration %% 10 == 1) %>% # Thin
  filter(Parameter %in% c("beta.phi[1]","beta.phi[2]")) %>%
  spread(key=Parameter, value = value) %>%
  ggplot(aes(x = `beta.phi[1]`, y = `beta.phi[2]`)) +
  geom_point() +
  geom_point(data = spread(truth,key = Parameter, value = Value),
             colour = "red") + 
  facet_wrap(vars(Model))

