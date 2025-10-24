

options(scipen = 999)
rm(list=ls())

library(data.table)
library(dplyr)
library(here)
library(ggplot2)

dat = data.table::fread(here('DATA', 'min_strides.csv'), header = TRUE)

dat$group = ifelse(grepl("^H", dat$id), "healthy", "unhealthy") # Specify clinical(?) group

dat$speed = gsub("p", ".", dat$speed)

dat$id = factor(dat$id)
dat$treadmill = factor(dat$treadmill)
dat$speed = factor(dat$speed)
dat$terrain = factor(dat$terrain,
                     levels = c('flat', 'low', 'med', 'high', 'NaN'))
dat$trial = factor(dat$trial)
dat$group = factor(dat$group)

hist(dat$n_stride_intervals, breaks = 50)
