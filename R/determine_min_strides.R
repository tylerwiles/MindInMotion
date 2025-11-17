

options(scipen = 999)
rm(list=ls())

library(data.table)
library(dplyr)
library(here)
library(ggplot2)

dat = data.table::fread(here('DATA', 'min_contacts.csv'), header = TRUE)

dat$speed = gsub("p", ".", dat$speed)

dat$id = factor(dat$id)
dat$treadmill = factor(dat$treadmill)
dat$speed = factor(dat$speed)
dat$terrain = factor(dat$terrain,
                     levels = c('flat', 'low', 'med', 'high', 'NaN'))
dat$trial = factor(dat$trial)

hist(dat$n.steps, breaks = 50)
hist(dat$n.contacts.left, breaks = 50)
hist(dat$n.contacts.right, breaks = 50)
hist(dat$n.contacts.diff, breaks = 50)

table(dat$id[dat$n.steps < 50])
table(dat$id[dat$n.steps < 50 & dat$n.contacts.diff > 2 | dat$n.contacts.diff < -2])

hist(dat$n.steps[dat$terrain == 'flat'], breaks = 50)
hist(dat$n.steps[dat$terrain == 'low'], breaks = 50)
hist(dat$n.steps[dat$terrain == 'med'], breaks = 50)
hist(dat$n.steps[dat$terrain == 'high'], breaks = 50)

hist(dat$n.contacts.diff[dat$terrain == 'flat'], breaks = 50)
hist(dat$n.contacts.diff[dat$terrain == 'low'], breaks = 50)
hist(dat$n.contacts.diff[dat$terrain == 'med'], breaks = 50)
hist(dat$n.contacts.diff[dat$terrain == 'high'], breaks = 50)

