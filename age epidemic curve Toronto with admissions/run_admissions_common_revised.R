#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(cmdstanr); library(dplyr); library(magrittr)})
a <- commandArgs(trailingOnly=TRUE)
if(length(a)!=1L) stop("Pass the 1-based task index")
idx <- as.integer(a[1]); if(is.na(idx)) stop("Invalid index")
project <- "/home/somerse5/wastewater_admissions_common_prior_sensitivity"
set.seed(2917); ids <- sample(1:3000, 500)
if(length(ids)!=500L || anyDuplicated(ids)) stop("Expected 500 unique IDs from previous admissions workflow")
if(idx<1L || idx>length(ids)) stop("Task index outside IDs file")
waste_id <- ids[idx]
out <- "/scratch/somerse5/wastewater_admissions_common_prior_sensitivity/fits"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
set_cmdstan_path(Sys.getenv("CMDSTAN"))
e <- new.env(parent=globalenv())
load(file.path(project,"epidemic_model_nonlag_data.RData"),envir=e)
init <- with(e, case_age_merged$case[case_age_merged$earliest_week_end_date == (date_range[1]-7)])
ratio <- with(e, data_foranalysis_full$analysis_d[[waste_id]] %>%
  arrange(earliest_week_end_date) %>%
  filter(earliest_week_end_date >= date_range[1], earliest_week_end_date <= date_range[2], !is.na(ratio_cases)) %$%
  ratio_v_u_fixed)
Y <- t(as.matrix(e$case_age_matrix[,-1,drop=FALSE]))
A <- t(as.matrix(e$hosp_age_matrix[,-1,drop=FALSE]))
I <- nrow(Y); J <- ncol(Y)
stopifnot(length(init)==I, length(ratio)==J, identical(dim(A),dim(Y)), all(A<=Y), all(is.finite(ratio)), all(ratio>0))
data <- list(I=I,J=J,Y=Y,A=A,ratio=as.numeric(ratio),mean_X0=init*2,
 lower_X0=init,sd_X0=c(786,233,124,30),phi_p=-log(.5)/.4,
 phi_tau=-log(.5)/.4,tau_mean_prior=c(-4,-4,-2,-2))
stopifnot(length(data$sd_X0)==I, length(data$tau_mean_prior)==I)
message("Task ",idx," / ",length(ids),"; wastewater ID ",waste_id,"; I=",I," J=",J)
mod <- cmdstan_model(file.path(project,"model_admissions_common_revised_priors.stan"),
 exe_file=file.path(project,"model_admissions_common_revised_priors"),compile=FALSE)
fit <- mod$sample(data=data,seed=12345L+waste_id,chains=4,parallel_chains=4,
 iter_warmup=1000,iter_sampling=2000,thin=1,save_warmup=FALSE,
 adapt_delta=.95,max_treedepth=12,output_dir=out)
# Do not silently publish a partial fit if a chain failed.
if(length(fit$output_files())!=4L || any(!file.exists(fit$output_files()))) stop("Missing chain output")
draws <- fit$draws()
if(dim(draws)[2]!=4L || dim(draws)[1]!=2000L) stop("Incomplete Stan sampling")
wb_rhat <- fit$summary()
sampler_diagnostics <- fit$sampler_diagnostics()
filename <- file.path(out,sprintf("admissions_common_revised_wastewatersample_%d.RData",waste_id))
tmp <- tempfile(pattern=".saving_",tmpdir=out,fileext=".RData")
save(draws,wb_rhat,sampler_diagnostics,file=tmp,compress=TRUE)
if(!file.rename(tmp,filename)) stop("Could not save ",filename)
unlink(fit$output_files())
message("Saved ",filename)
