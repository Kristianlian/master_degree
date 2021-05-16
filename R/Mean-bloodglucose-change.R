#### Blood glucose change analysis

## Author: Kristian Lian
## Project: Ribose

## This script analyses log-fold change in plasma glucose per supplement. 


# Packages
library(tidyverse);library(readxl);library(nlme);library(lme4);library(knitr);library(broom);library(emmeans);library(dplyr)

## Data

gluc.dat <- read_excel("./data/glucose/fingerstick.xlsx", na = "NA")
gluc.dat$sample_time <- as.character(gluc.dat$sample_time)

## Handling the data by creating a new factor called time from sample_time. This factor combines any observation at 0 to baseline, etc. 
# The code also turns glucose and placebo values from character to numeric for later analysis

glu.dat <- gluc.dat %>%
        mutate(time = if_else(sample_time == "0",
                              "baseline",
                              if_else(sample_time == "45",
                                      "min45",
                                      if_else(sample_time == "90",
                                              "min90",
                                              if_else(sample_time == "120",
                                                      "min120",
                                                      if_else(sample_time == "135",
                                                              "min135",
                                                              if_else(sample_time == "150",
                                                                      "min150",
                                                                      if_else(sample_time == "270",
                                                                              "min270", sample_time))))))),
               time = factor(time, levels = c("baseline", "min45", "min90", "min120", "min135", "min150", "min270"))) %>%
        mutate(glu = as.numeric(glu),
               lak = as.numeric(lak)) %>%
        print()

## Change data
# The code beneath summarizes the mean glucose values at each time, grouped by subject, time and supplement, creating a wider data set with observations of 
# participants glucose measurements per time point.
# Then, mutate() is used to calculate change scores, where each timepoint is log-transformed and compared to baseline. baseline = baseline - mean(baseline,
# na.rm = TRUE) mean centers the baseline values. Subject, supplement, baseline and change scores are then selected and pivoted for modeling. 

change_dat.glu <- glu.dat %>%
        dplyr::select(subject, time, glu, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(glu = mean(glu, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = glu) %>%
        #print()
        
        ungroup() %>%
        mutate(change.45 = log(min45)-log(baseline),
               change.90 = log(min90)-log(baseline),
               change.120 = log(min120)-log(baseline),
               change.135 = log(min135)-log(baseline),
               change.150 = log(min150)-log(baseline),
               change.270 = log(min270)-log(baseline),
               baseline = baseline - mean(baseline, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.45, change.90, change.120, change.135, 
               change.150, change.270) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.45:change.270)) %>%
        print()

## Linear mixed effects model
# This model tries to explain the change by time and supplement, accounting for potential differences in baseline values and that the same participants
# are measured at multiple time points. 
# It produces results on both the time effect and the difference between the groups at any timepoint. We are interested in the difference between groups.

m1 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat.glu)

plot(m1)

summary(m1)

### Get estimated means from the model, these are average increase at 
# pre = 0 (the average pre value)
# These are log-fold change values (changeble in the mutate function)


confint.m1 <- confint(emmeans(m1, specs = ~"supplement|time")) %>%
        data.frame()



