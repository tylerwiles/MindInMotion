

options(scipen = 999)
rm(list=ls())

library(data.table)
library(dplyr)
library(here)
library(ggplot2)
library(tidyr)
library(ggridges)

dat = data.table::fread(here('DATA', 'min_contacts.csv'), header = TRUE)

dat = dat %>%
  pivot_longer(cols = matches("^n\\.contacts(\\.(left|right))?\\.seg\\d+$"),
               names_to = c("measure", "segment"),
               names_pattern = "^(n\\.contacts(?:\\.(?:left|right))?)\\.seg(\\d+)$",
               values_to = "value") %>%
  mutate(segment = paste0("seg", segment)) %>%
  pivot_wider(names_from  = measure,
              values_from = value)

dat = dat[dat$treadmill == 'TM',]
dat = dat[!duplicated(dat[, c("id", "treadmill", "speed", "terrain", "trial")]), ] # Ignore segments

dat$speed = gsub("p", ".", dat$speed)

dat$id = factor(dat$id)
dat$treadmill = factor(dat$treadmill)
dat$speed = factor(dat$speed)
dat$terrain = factor(dat$terrain,
                     levels = c('flat', 'low', 'med', 'high', 'NaN'))
dat$trial = factor(dat$trial)
dat$segment = factor(dat$segment)

# Get difference between number of left and right contacts
dat$n.contacts.diff = dat$n.contacts.left - dat$n.contacts.right

hist(dat$n.contacts.total, breaks = 50)
hist(dat$n.contacts.left, breaks = 50)
hist(dat$n.contacts.right, breaks = 50)
hist(dat$n.contacts.diff, breaks = 50)

nrow(dat[dat$n.contacts.total < 50, ])
nrow(dat[dat$n.contacts < 50, ])


table(dat$id[dat$n.contacts.total < 50])
table(dat$id[dat$n.contacts.total < 50 & dat$n.contacts.diff > 2 | dat$n.contacts.diff < -2])

hist(dat$n.contacts.total[dat$terrain == 'flat'], breaks = 50)
hist(dat$n.contacts.total[dat$terrain == 'low'], breaks = 50)
hist(dat$n.contacts.total[dat$terrain == 'med'], breaks = 50)
hist(dat$n.contacts.total[dat$terrain == 'high'], breaks = 50)

ggplot(dat, aes(x = n.contacts, y = segment, fill = segment)) +
  geom_density_ridges() +
  labs(title = "Ridgeline Plot Example",
       x = "Sepal Length",
       y = "Species")

# Are there a signfiicantly different number of contacts per segment?
dat.temp = dat[dat$n.contacts >= 32, ]
m0 = aov(n.contacts ~ segment, data = dat.temp)
summary(m0)
t.test(x = dat.temp$n.contacts[dat.temp$segment == 'seg1'],
       y = dat.temp$n.contacts[dat.temp$segment == 'seg2'])
t.test(x = dat.temp$n.contacts[dat.temp$segment == 'seg1'],
       y = dat.temp$n.contacts[dat.temp$segment == 'seg3'])
t.test(x = dat.temp$n.contacts[dat.temp$segment == 'seg2'],
       y = dat.temp$n.contacts[dat.temp$segment == 'seg3'])

mean(dat.temp$n.contacts[dat.temp$segment == 'seg1'])
mean(dat.temp$n.contacts[dat.temp$segment == 'seg2'])
mean(dat.temp$n.contacts[dat.temp$segment == 'seg3'])
sd(dat.temp$n.contacts[dat.temp$segment == 'seg1'])
sd(dat.temp$n.contacts[dat.temp$segment == 'seg2'])
sd(dat.temp$n.contacts[dat.temp$segment == 'seg3'])

dat[dat$n.contacts >= 32, ] %>%
  group_by(segment, treadmill, id) %>%
  summarise(min_contacts = min(n.contacts),
            max_contacts = max(n.contacts),
            .groups = "drop") %>%
  group_by(segment, treadmill) %>%
  summarise(mean_min_contacts = mean(min_contacts),
            mean_max_contacts = mean(max_contacts),
            .groups = "drop")



dat[dat$n.contacts.total >= 32, ] %>%
  group_by(segment, treadmill, id) %>%
  summarise(min_contacts = min(n.contacts.total),
            max_contacts = max(n.contacts.total),
            .groups = "drop") %>%
  group_by(segment, treadmill) %>%
  summarise(mean_min_contacts = mean(min_contacts),
            mean_max_contacts = mean(max_contacts),
            .groups = "drop")





ggplot(dat, aes(x = n.contacts.total, y = terrain, fill = terrain)) +
  geom_density_ridges() +
  labs(title = "Ridgeline Plot Example",
       x = "Sepal Length",
       y = "Species")









