## Options
options(echo = TRUE)

## Load packages
library(tidyverse)
library(MASS)

## Load package containing common files
## library(CCI2020)

library(devtools)
devtools::load_all("../CCI2020")

## Source simulation specific code
source("R/survival.R")
source("R/generate_data.R")
source("R/simulation_functions.R")

## Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
k <- args[1]
pars_mat_name <- args[2]

cat(k)
cat(pars_mat_name)

## Load parameter matrix
pars_mat <- read_csv(pars_mat_name)

## Run replicate
test <- run_replicate(pars_mat[k,])
