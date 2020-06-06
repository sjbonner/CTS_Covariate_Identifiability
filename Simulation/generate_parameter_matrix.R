## Load packages
library(tidyverse)

## Set simulation parameters
debug <- FALSE

if(debug){
  nrep <- 2
  chains <- 3
  burnin <- 100
  samples <- 100
}
if(!debug){
  nrep <- 50
  chains <- 3
  burnin <- 2000
  samples <- 20000
}

## Define simulation parameters
pars_mat <- crossing(nind = c(500,1000,5000,10000),
                     ncap = 10,
                     model = c("logit","scaled_logit","generalized_logit"),
                     chains = chains,
                     burnin = burnin,
                     samples = samples) %>%
  arrange(nind,ncap,model) %>%
  rowid_to_column("sim") %>%
  crossing(rep = 1:nrep) %>%
  rowid_to_column("task")
  
write_csv(pars_mat,"parameter_matrix.csv")
