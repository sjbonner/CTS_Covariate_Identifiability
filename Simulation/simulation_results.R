## Load packages
library(tidyverse)
library(ggmcmc)
library(latex2exp)

## If true generate 16x9 aspect ratio figure for beamer
beamer <- FALSE

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
   summarize(Count = n()) %>%
   arrange(Count))

## Summarize across simulations
results <- summaries %>%
  group_by(sim, model, Model, Parameter, nind) %>%
  summarize(MeanSD = mean(SD))

## Keep only parameters common to simulation scenarios
pars <- c("p","lambda","beta.phi[1]","beta.phi[2]","lower","upper")

results <- results %>%
  filter(Parameter %in% pars)

## Add nicer model and parameter names
results$model_name <- factor(results$model,
                             levels = c("logit",
                                        "scaled_logit",
                                        "generalized_logit"),
                             labels = TeX(c("Logit",
                                        "Scaled Logit",
                                        "Generalized Logit")))

results$Parameter_name = factor(results$Parameter,
                                levels = c("p",
                                           "lambda",
                                           "beta.phi[1]",
                                           "beta.phi[2]",
                                           "upper",
                                           "lower"),
                                 labels = TeX(
                                   paste0("$",
                                          c("p",
                                            "\\lambda",
                                            "\\alpha",
                                            "\\beta",
                                            "\\gamma_1",
                                            "\\gamma_0"),
                                          "$")))

## Compute SD of limiting posterior distribution
rho <- truth %>%
  filter(Parameter %in% c("p","upper")) %>%
  summarize(Value = prod(Value)) %>%
  pull("Value")

limiting_sd_scaled <- results %>%
  filter(model == "scaled_logit",
         Parameter %in% c("p","upper"),
         Model == "Binomial") %>%
  mutate(LVar = -(1-rho^2)/(2 * log(rho)) - ((1-rho)/log(rho))^2,
         LSD = sqrt(LVar))

limiting_sd_generalized <- results %>%
  filter(model == "generalized_logit",
         Parameter %in% c("p"),
         Model == "Binomial") %>%
  mutate(LVar = rho - (rho * log(rho)/(1 - rho))^2,
         LSD = sqrt(LVar))


## Plot mean standard deviation versus sample size
if(beamer)
    x11(width=7,height=3)
if(!beamer)
    x11(width = 7, height = 6)

results %>%
    filter(Model != "Alt. Trinomial") %>%
  ggplot(aes(x = nind, y=MeanSD, colour=Model)) +
  geom_point() +
  geom_line() + 
  geom_hline(yintercept = 0) +
  facet_grid(Parameter_name ~ model_name,
             scale = "free",
             labeller = label_parsed) +
  xlab("Sample Size") +
  ylab("Mean Posterior Standard Deviation") +
  scale_y_continuous(n.breaks = 3) + 
  theme(strip.text.y.right = element_text(angle = 0)) # + 
  ## geom_line(data = bind_rows(limiting_sd_scaled,
  ##                            limiting_sd_generalized),
  ##           mapping = aes(x = nind, y = LSD),
  ##           lty = 2)

if(beamer)
  dev.copy2pdf(file = "Figures/simulation_summary_beamer.pdf")
if(!beamer)
  dev.copy2pdf(file = "Figures/simulation_summary.pdf")
  
