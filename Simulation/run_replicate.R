## Options
options(echo = TRUE)

## Load packages
library(tidyverse)
library(MASS)

## Load installed manuscript package 
library(CCI2020)

## Load manuscript package from working directory
## library(devtools)
## devtools::install("../CCI2020")

## Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
k <- args[1]
pars_mat_name <- args[2]

cat(k,"\n")
cat(pars_mat_name,"\n")

## Load parameter matrix
pars_mat <- read_csv(pars_mat_name)

## Run replicate
test <- run_replicate(pars_mat[k,])
