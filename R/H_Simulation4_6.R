


# Words

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

# Sweep Setup -----------------------------------------------------------------------

# Stride interval mean and standard deviation from Gaitprint
si.mean = 1.05
si.sd = 0.09

# Create a list of parameter sweeps
strides = seq(50, 200, 50)
drops = seq(0.02, .1, 0.02)
h.target = c(seq(0.4, 0.6, 0.1))

params = expand.grid(strides = strides,
                     drops = drops,
                     h.target = h.target,
                     KEEP.OUT.ATTRS = FALSE)

# Create a specified number of times we want to repeat our list of parameter sweeps and assign an id value
repetitions = 3000
params = params[rep(seq_len(nrow(params)), each = repetitions), ] # replicate the rows
rownames(params) = NULL # reset row names
params$id = rep(seq_len(repetitions), times = nrow(params) / repetitions) # assign an id column 1:repetitions repeating across all combinations
print('params done')

# Parallel Setup --------------------------------------------------------------------

# n_cores = parallel::detectCores() - 1 # Use if not on cluster for testing
slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", ""))
if (is.na(slurm_cpus) || slurm_cpus < 1L) {
  # physical cores preferred (avoid hyperthreads)
  n_cores <- max(1L, parallel::detectCores(logical = FALSE))
} else {
  n_cores <- slurm_cpus
}
print('done 1')
cl = makeCluster(n_cores)
print('done 2')
registerDoParallel(cl)
print('done 3')

cpp_path = here::here('FUNCTIONS', 'CPP', 'bayesH.cpp')
r_path = here::here('FUNCTIONS', 'R', 'fgn_sim.R')
sampen_path = here::here('FUNCTIONS', 'CPP', 'Ent_Samp.cpp')

# Testing versions
# Rcpp::sourceCpp(here::here('FUNCTIONS', 'CPP', 'bayesH.cpp'))
# source(here::here('FUNCTIONS', 'R', 'fgn_sim.R'))
# Rcpp::sourceCpp(here::here('FUNCTIONS', 'CPP', 'Ent_Samp.cpp'))
print('done 4')

## workers need: paths + params + scalars
clusterExport(cl, c("cpp_path","r_path","sampen_path","params","si.sd","si.mean"), envir = environment())
print('done 5')

## compile/load Rcpp code and source R on each worker
lib_paths = .libPaths()
clusterCall(cl, function(p) .libPaths(p), lib_paths)

clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(Rcpp)
    library(RcppArmadillo)
  })
  source(r_path)  # fgn_sim.R
  Rcpp::sourceCpp(cpp_path,    rebuild = TRUE, showOutput = FALSE, verbose = FALSE) # bayesH.cpp
  Rcpp::sourceCpp(sampen_path, rebuild = TRUE, showOutput = FALSE, verbose = FALSE) # Ent_Samp.cpp
  TRUE
})
print('done 6')

# Sweep it --------------------------------------------------------------------------

results = foreach(i = 1:nrow(params),
                  .combine = rbind,
                  .export = character(0),  # don't re-export; already clusterExport'ed
                  .noexport = c('bayesH','Ent_Samp')) %dopar% {
                    
                    strides.temp = params$strides[i]
                    drops.temp = params$drops[i]
                    h.target.temp = params$h.target[i]
                    id.temp = params$id[i]
                    
                    # Create a time series that mimics the behavior of typical stride intervals
                    dat = fgn_sim(n = strides.temp, H = h.target.temp) * si.sd + si.mean
                    h.original = median(bayesH(dat, 200)) # Take the median of 200 samples for H
                    # tol = 0.25 * sd(dat)
                    # sampen.original = Ent_Samp(dat, 2, tol) # Run Sample Entropy on 25% the SD of the time series
                    
                    # Randomly drop percentage number of datapoints
                    dat.rand = dat[-sample(seq_along(dat), size = floor(length(dat) * drops.temp))]
                    h.test.random = median(bayesH(dat.rand, 200)) # Take the median of 200 samples for H
                    # tol = 0.25 * sd(dat.rand)
                    # sampen.test.random = Ent_Samp(dat.rand, 2, tol) # Run Sample Entropy on 25% the SD of the time series
                    
                    # Drop contiguous chunks of the percentage of datapoints
                    temp.length = length(dat)
                    chunk.length = max(1, min(floor(temp.length * drops.temp), temp.length - 1)) # Obtain number of datapoints to remove
                    chunk.start = sample(temp.length - chunk.length + 1, 1) # Randomly determine the starting point of that chunk
                    chunk.index = chunk.start:(chunk.start + chunk.length - 1) # Obtain the indices to remove
                    dat.contiguous = dat[-chunk.index] # Remove indices
                    h.test.contiguous = median(bayesH(dat.contiguous, 200)) # Take the median of 200 samples for H
                    # tol = 0.25 * sd(dat.contiguous)
                    # sampen.test.contiguous = Ent_Samp(dat.contiguous, 2, tol) # Run Sample Entropy on 25% the SD of the time series
                    
                    data.frame(id = id.temp,
                               strides = strides.temp,
                               drops = drops.temp,
                               strides.cut = length(dat) - length(dat.rand),
                               strides.new.length = length(dat.rand),
                               h.target = h.target.temp,
                               h.original = h.original,
                               h.test.random = h.test.random,
                               h.test.contiguous = h.test.contiguous,
                               # sampen.original = sampen.original,
                               # sampen.test.random = sampen.test.random,
                               # sampen.test.contiguous = sampen.test.contiguous,
                               stringsAsFactors = FALSE)
                    
                  }

print('sweep done')

stopCluster(cl)

write.csv(results, here::here("DATA", "H_Simulation", "H_Simulation.csv"), row.names = FALSE)
