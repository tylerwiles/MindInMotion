
# Description: Before this script is run, the matlab script with the same name
# should have already be run. From the time series of stride intervals from the
# mind in motion dataset, the Hurst exponent is calculated for the final 50,
# 100, 150, 200 strides (when the number of strides allows) and then the time
# series is dropped from 2-10% randomly or contiguously. Each of these
# combinations of time series lengths and % drops are run multiple times (i.e.,
# 100 repetitions) but the average of them are exported. This code is written to
# be run on Hipergator but minor changes will run locally.

options(scipen = 999)
rm(list = ls())

library(data.table)
library(doParallel)
library(dplyr)
library(here)
library(foreach)
library(parallel)
library(Rcpp)
library(RcppArmadillo)

my.directory = "/blue/dferris/twiles/HSIM_DATA/H_Simulation/MIM"
# my.directory = "R:/Ferris-Lab/twiles/H_Simulation/MIM" # For Testing
my.files = list.files(my.directory, full.names = FALSE)

drops = seq(0.02, 0.10, by = 0.02)
min.vec = c(50, 100, 150, 200)

n.files = length(my.files)
n.drops = length(drops)
n.min = length(min.vec)
n.rep = 100 # number of repetitions for each i / min.contacts / drop

# Parallel Setup --------------------------------------------------------------------

# Rcpp::sourceCpp(here::here('FUNCTIONS', 'CPP', 'bayesH.cpp')) # For Testing
cpp.path = here::here('FUNCTIONS', 'CPP', 'bayesH.cpp')

slurm.cpus = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", ""))
if (is.na(slurm.cpus) || slurm.cpus < 1L) {
  # physical cores preferred (avoid hyperthreads)
  n.cores = max(1L, parallel::detectCores(logical = FALSE))
} else {
  n.cores = slurm.cpus
}
print('done 1')
cl = makeCluster(n.cores)
print('done 2')
registerDoParallel(cl)
print('done 3')

clusterExport(cl, varlist = c("cpp.path"), envir = environment()) # export the C++ file path to workers
print('done 4')

# sync library paths
lib.paths = .libPaths()
clusterCall(cl, function(p) .libPaths(p), lib.paths)
print('done 5')

clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(Rcpp)
    library(RcppArmadillo)
  })
  Rcpp::sourceCpp(cpp.path, rebuild = FALSE, showOutput = FALSE, verbose = FALSE)
  TRUE
})
print('done 6')

# Run MIM Data ------------------------------------------------------

result.lists = foreach(i = seq_len(n.files),
                       .export = c("drops", "min.vec", "n.drops", "n.min", "n.rep",
                                   "my.directory", "my.files"),
                       .packages = c("data.table"),
                       .combine = "c") %dopar% {
                         
                         base.filename = my.files[i]
                         full.filename = file.path(my.directory, base.filename)
                         
                         png.filename = strsplit(base.filename, ".", fixed = TRUE)[[1]][1]
                         
                         cat("Analyzing:", png.filename)
                         
                         dat = data.table::fread(full.filename, header = FALSE)
                         dat = as.numeric(dat[[1L]])
                         
                         # Max possible rows for this file: one row per (min.contacts, drop)
                         max.rows.i = n.min * n.drops
                         rows.list.i = vector("list", length = max.rows.i)
                         row.idx.i = 0L
                         
                         for (k in seq_len(n.min)) { # Loop over minimum number of stride intervals
                           
                           min.contacts = min.vec[k]
                           
                           if (length(dat) >= min.contacts) { # only run for this threshold if we have enough contacts
                             
                             dat.cut = tail(dat, min.contacts) # Cut to last min.contacts
                             
                             h.orig.val = median(bayesH(dat.cut, 200)) # H for original series for this (file, min.contacts) pair
                             
                             n.total = length(dat.cut) # total strides for this file
                             
                             for (j in seq_len(n.drops)) { # Loop over % drops
                               
                               this.drop = drops[j]
                               
                               n.keep = n.total - floor(n.total * this.drop)
                               n.drop = n.total - n.keep
                               
                               strides.val = n.total
                               strides.cut.val = n.drop
                               strides.new.length.val = n.keep
                               
                               # collect H values over n.rep repetitions
                               h.random.vec = numeric(n.rep)
                               h.contig.vec = numeric(n.rep)
                               
                               for (r in seq_len(n.rep)) { # Loop over number of repetitions
                                 
                                 # random removal
                                 keep.idx.random = sort(sample.int(n.total, n.keep))
                                 dat.cut.random = dat.cut[keep.idx.random]
                                 h.random.vec[r] = median(bayesH(dat.cut.random, 200))
                                 
                                 # contiguous removal
                                 start.idx = sample.int(n.total - n.drop + 1L, 1L)
                                 
                                 keep.mask = rep(TRUE, n.total)
                                 if (n.drop > 0L) {
                                   keep.mask[start.idx:(start.idx + n.drop - 1L)] = FALSE
                                 }
                                 
                                 dat.cut.contiguous  = dat.cut[keep.mask]
                                 h.contig.vec[r] = median(bayesH(dat.cut.contiguous, 200))
                               }
                               
                               # Average across reps
                               h.random.mean = mean(h.random.vec)
                               h.contig.mean = mean(h.contig.vec)
                               
                               row.idx.i = row.idx.i + 1L # Store ONE row per (file, min.contacts, drop)
                               
                               rows.list.i[[row.idx.i]] = list(id = png.filename,
                                                               rep = n.rep, # Number of repetitions that were averaged
                                                               stride.intervals = strides.val, # Number of original strides
                                                               drops = this.drop, # % of trial to be dropped
                                                               stride.intervals.cut = strides.cut.val, # Number of strides cut
                                                               stride.intervals.new.length = strides.new.length.val, # New number of strides
                                                               h.original = h.orig.val, # H from the original time series
                                                               h.random = h.random.mean, # H from the randomly cut timeseries
                                                               h.contiguous = h.contig.mean # H from the contiguously cut timeseries
                               )
                               
                             }
                             
                           }
                           
                         }
                         
                         # Trim unused entries for this file
                         if (row.idx.i < length(rows.list.i)) {
                           rows.list.i = rows.list.i[seq_len(row.idx.i)]
                         }
                         
                         rows.list.i
                       }

stopCluster(cl)

# Flatten and make table
rows.list.all = result.lists
results = do.call(rbind.data.frame, c(rows.list.all, list(stringsAsFactors = FALSE)))
results = results[, c("id", "rep", "stride.intervals", "drops",
                      "stride.intervals.cut", "stride.intervals.new.length",
                      "h.original", "h.random", "h.contiguous")]

write.csv(results, here::here("DATA", "H_Simulation", "H_Simulation_MIM.csv"), row.names = FALSE)
