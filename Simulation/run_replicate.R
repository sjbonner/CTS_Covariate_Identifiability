## Load packages
library(tidyverse)
library(parallel)
library(MASS)
library(spark)
library(nimble)
library(coda)

## Source code
source("R/survival.R")
source("R/generate_data.R")
source("R/truncated_model_functions_nimble.R")
source("R/trinomial_functions_nimble.R")
source("R/binomial_functions_nimble.R")
source("R/simulation_functions.R")

## Read command line arguments
args = commandArgs(trailingOnly = TRUE)
k = args[1]

## Load parameter matrix
pars_mat <- read_csv("parameter_matrix.csv")

## Run replicate
test <- run_replicate(pars_mat[k,])
