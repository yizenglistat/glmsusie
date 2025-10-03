rm(list=ls())
setwd("~/Documents/glmcs")


library(glmcs)
library(glue)
library(data.table)
source("files/prepare_data.R")

temp <- sapply(list.files("R", pattern="*.R$", full.names=TRUE, ignore.case=TRUE), function(x) source(x))
rm(temp)

standardize <- TRUE
decompose <- TRUE
shrinkage <- TRUE

susie_only <- FALSE
regression <- "gaussian"

L <- 25

ukb_metab2_full_path <- glue("data/ukb_metab2_full.csv")
ukb_metab2_full <- fread(ukb_metab2_full_path)

data_list <- prepare_data(ukb_metab2_full, phase="phase2")
phase_name <- "phase2full"

n <- 1000#length(data_list$y)
p <- ncol(data_list$X)

fit <- glmcs(
	X=data_list$X[1:n,],
	y=data_list$y[1:n],
	L=L,
	standardize=standardize,
	decompose=decompose,
	shrinkage=shrinkage,
	family=gaussian()
)