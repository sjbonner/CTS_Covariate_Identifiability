## Load packages
library(tidyverse)
library(parallel)
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
args = commandArgs(trailingOnly = TRUE)
k = args[1]

## Load parameter matrix
pars_mat <- read_csv("parameter_matrix.csv")

## Run replicate
test <- run_replicate(pars_mat[k,])
