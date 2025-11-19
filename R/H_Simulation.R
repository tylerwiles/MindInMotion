
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

# Sweep Setup -----------------------------------------------------------------------

# Commonly observed stride interval mean and standard deviation
si.mean = 1
si.sd   = 0.2

# Create a list of parameter sweeps
strides = seq(32, 512, 32)
drops = seq(0.01, .1, 0.01)
h.target = c(seq(0.5, 0.9, 0.1), 0.999)
# For Testing
# strides = seq(32, 97, 32)
# drops = seq(0.05, .15, 0.05)
# h.target = c(seq(0.5, 0.6, 0.1), 0.999)

params = expand.grid(strides = strides,
                     drops = drops,
                     h.target = h.target,
                     KEEP.OUT.ATTRS = FALSE)

# Create a specified number of times we want to repeat our list of parameter sweeps and assign an id value
repetitions = 1000
params = params[rep(seq_len(nrow(params)), each = repetitions), ] # replicate the rows
rownames(params) = NULL # reset row names
params$id = rep(seq_len(repetitions), times = nrow(params) / repetitions) # assign an id column 1:repetitions repeating across all combinations
print('params done')
# Parallel Setup --------------------------------------------------------------------

# n_cores = max(1, parallel::detectCores() - 1)
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
r_path = here::here('FUNCTIONS', 'R',   'fgn_sim.R')
print('done 4')

## workers need: paths + params + scalars
clusterExport(cl, c("cpp_path","r_path","params","si.sd","si.mean"), envir = environment())
print('done 5')

## compile/load Rcpp code and source R on each worker
clusterEvalQ(cl, {
  suppressPackageStartupMessages(library(Rcpp))
  source(r_path)
  Rcpp::sourceCpp(cpp_path, rebuild = TRUE, showOutput = FALSE, verbose = FALSE)
  TRUE
})
print('done 6')

set.seed(9232025)

# Sweep it --------------------------------------------------------------------------

results = foreach(i = 1:nrow(params),
                  .combine = rbind,
                  .export = character(0),  # don't re-export; already clusterExport'ed
                  .noexport = 'bayesH') %dopar% {
                    
                    strides.temp = params$strides[i]
                    drops.temp = params$drops[i]
                    h.target.temp = params$h.target[i]
                    id.temp = params$id[i]
                    
                    # Create a time series that mimics the behavior of typical stride intervals
                    dat = fgn_sim(n = strides.temp, H = h.target.temp) * si.sd + si.mean
                    h.original = median(bayesH(dat, 200)) # Take the median of 200 samples for H
                    
                    # Randomly drop percentage number of datapoints
                    dat.rand = dat[-sample(seq_along(dat), size = floor(length(dat) * drops.temp))]
                    h.test.random = median(bayesH(dat.rand, 200)) # Take the median of 200 samples for H
                    
                    # Drop contiguous chunks of the percentage of datapoints
                    temp.length = length(dat)
                    chunk.length = max(1, min(floor(temp.length * drops.temp), temp.length - 1)) # Obtain number of datapoints to remove
                    chunk.start = sample(temp.length - chunk.length + 1, 1) # Randomly determine the starting point of that chunk
                    chunk.index = chunk.start:(chunk.start + chunk.length - 1) # Obtain the indices to remove
                    dat.contig = dat[-chunk.index] # Remove indices
                    h.test.contiguous = median(bayesH(dat.contig, 200)) # Take the median of 200 samples for H
                    
                    data.frame(id = id.temp,
                               strides = strides.temp,
                               drops = drops.temp,
                               strides.cut = length(dat.rand),
                               strides.new.length = length(dat.rand),
                               h.target = h.target.temp,
                               h.original = h.original,
                               h.test.random = h.test.random,
                               h.test.contiguous = h.test.contiguous,
                               stringsAsFactors = FALSE)

                  }

print('sweep done')

stopCluster(cl)

# Split results into quarters and export
split.points = ceiling(seq(0, nrow(results), length.out = 5))
for (i in 1:4) {
  part = results[(split.points[i] + 1):split.points[i + 1], ]
  write.csv(part, here::here("DATA", paste0("H_Simulation_Q", i, ".csv")), row.names = FALSE)
}
