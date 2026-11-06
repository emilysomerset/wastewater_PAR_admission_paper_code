#!/usr/bin/env Rscript
source_dir <- '/scratch/somerse5/wastewater_model/results_noadm_common_rw_long'
dest_dir <- '/scratch/somerse5/wastewater_model/results_noadm_common_rw_long_reduced'
dir.create(dest_dir, recursive=TRUE, showWarnings=FALSE)
files <- sort(list.files(source_dir, pattern='^stan_samples_common_rw_long_wastewatersample_[0-9]+\\.RData$', full.names=TRUE))
if(length(files)!=500L) stop('Expected 500 full fits; found ', length(files))
for(f in seq_along(files)) {
  id <- sub('^stan_samples_common_rw_long_wastewatersample_([0-9]+)\\.RData$', '\\1', basename(files[f]))
  target <- file.path(dest_dir, sprintf('stan_samples_common_rw_wastewatersample_%s.RData',id))
  if(file.exists(target)) { if(f %% 50L == 0L) message('Checked ',f,'/500'); next }
  e <- new.env(parent=emptyenv()); load(files[f], envir=e)
  d <- e$draws
  if(length(dim(d))!=3L || dim(d)[1]!=2000L || dim(d)[2]!=4L) stop('Unexpected draw dimensions: ',basename(files[f]))
  set.seed(2917L + as.integer(id))
  idx <- sort(sample.int(dim(d)[1],250L,replace=FALSE))
  draws <- d[idx,,,drop=FALSE]
  wb_rhat <- e$wb_rhat
  if(is.null(wb_rhat)) stop('Missing wb_rhat: ',basename(files[f]))
  # Retain original full-fit diagnostics; these are NOT computed on the subset.
  tmp <- tempfile(pattern='.partial_',tmpdir=dest_dir,fileext='.RData')
  save(draws,wb_rhat,file=tmp,compress=TRUE)
  if(!file.rename(tmp,target)) stop('Failed saving ',target)
  rm(e,d,draws,wb_rhat); gc(verbose=FALSE)
  if(f %% 25L == 0L) message('Reduced ',f,'/500 fits')
}
message('Finished. Reduced files: ',length(list.files(dest_dir,pattern='^stan_samples_common_rw_wastewatersample_[0-9]+\\.RData$')))
