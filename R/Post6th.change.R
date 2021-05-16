#### Humac analysis pre vs. post sixth RT session



# This script provides log-fold change score calculation of humac test prior to and after RT session 6.


# Packages
library(readxl);library(tidyverse);library(nlme);library(lme4);library(broom);library(knitr);library(emmeans)

## Handling the data by creating a new factor called time from timepoint. This factor combines any observation at T1 and T2 to baseline, etc. 
# The code also sorts the order of the factor time, from baseline to session 6, using time = factor(time, levels c()), and sets placebo to be compared to 
# glucose via supplement = factor(supplement, levels = c()). Acute code is called to set a new factor named acute, so that its possible to divid post 5th
# session data from post 6th session data

humac <- read_excel("./data/tests/ribose.humac.xlsx", na = "NA") %>%
        
        
        mutate(time = if_else(timepoint == "D-1", 
                              "baseline", 
                              if_else(timepoint %in% c("D4", "D5"), 
                                      "test1", 
                                      if_else(timepoint %in% c("D8", "D9"), 
                                              "test2", 
                                              if_else(timepoint %in% c("T3", "T4") & acute %in% c("rest", "post30min", "post2h"),
                                                      "test3",
                                                      if_else(acute == "post23h", "test4", timepoint)))))) %>%
        mutate(time = factor(time, levels = c("baseline", "test1", "test2", "test3", "test4")), 
               acute = factor(acute, levels = c("rest", "post30min", "post2h", "post23h")), 
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        print()


## Change-data
# The code beneath summarizes the mean values at each time, grouped by subject, time and supplement, creating a wider data set with observations of 
# participants glucose measurements per time point.
# Then, mutate() is used to calculate change scores, where each timepoint is log-transformed and compared to baseline. baseline = baseline - mean(baseline,
# na.rm = TRUE) mean centers the baseline values. Subject, supplement, baseline and change scores are then selected and pivoted for modeling. The data set is
# filtered according to test exercise (isometric, isokinetic 60 or isokinetic 240)

# Isometric
isom.dat <- humac %>%
        filter(test == "isom",
               time %in% c("test3", "test4")) %>%
        print()

change_dat <- isom.dat %>%
        dplyr::select(subject, acute, supplement, peak.torque) %>%
        group_by(subject, acute, supplement) %>%
        summarise(peak.torque = mean(peak.torque, na.rm = TRUE)) %>%
        pivot_wider(names_from = acute, 
                    values_from = peak.torque) %>%
        #print()
        
        ungroup() %>%
        mutate(change.2 = log(post30min)-log(rest),
               change.3 = log(post2h)-log(rest),
               change.4 = log(post23h)-log(rest),
               baseline = rest - mean(rest, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.4)) %>%
        print()

# Isok.60

isok60.dat <- humac %>%
        filter(test == "isok.60",
               time %in% c("test3", "test4")) %>%
        print()

change_dat2 <- isok60.dat %>%
        dplyr::select(subject, acute, supplement, peak.torque) %>%
        group_by(subject, acute, supplement) %>%
        summarise(peak.torque = mean(peak.torque, na.rm = TRUE)) %>%
        pivot_wider(names_from = acute, 
                    values_from = peak.torque) %>%
        #print()
        
        ungroup() %>%
        mutate(change.2 = log(post30min)-log(rest),
               change.3 = log(post2h)-log(rest),
               change.4 = log(post23h)-log(rest),
               baseline = rest - mean(rest, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.4)) %>%
        print()

## Isok.240

isok240.dat <- humac %>%
        filter(test == "isok.240",
               time %in% c("test3", "test4")) %>%
        print()

change_dat3 <- isok240.dat %>%
        dplyr::select(subject, acute, supplement, peak.torque) %>%
        group_by(subject, acute, supplement) %>%
        summarise(peak.torque = mean(peak.torque, na.rm = TRUE)) %>%
        pivot_wider(names_from = acute, 
                    values_from = peak.torque) %>%
        #print()
        
        ungroup() %>%
        mutate(change.2 = log(post30min)-log(rest),
               change.3 = log(post2h)-log(rest),
               change.4 = log(post23h)-log(rest),
               baseline = rest - mean(rest, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.4)) %>%
        print()


## Linear mixed effects model
# This model tries to explain the change by time and supplement, accounting for potential differences in baseline values and that the same participants
# are measured at multiple time points. 
# It produces results on both the time effect and the difference between the groups at any timepoint. We are interested in the difference between groups.

# Isometric
m1 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat)
plot(m1)

summary(m1)

# Isok.60
m2 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat2)
plot(m2)

summary(m2)

# Isok.240
m3 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat3)
plot(m3)

summary(m3)

## Fold-change estimated means
# Gets estimated means from the model, these are average increase at pre = 0 (the average pre value).
# These are log-fold change values (changeble with the mutate function)

# Isometric

confint.m1 <- confint(emmeans(m1, specs = ~"supplement|time")) %>%
        data.frame() %>%
        print()

# Isok.60

confint.m2 <- confint(emmeans(m2, specs = ~"supplement|time")) %>%
        data.frame() %>%
        print()

# Isok.240

confint.m3 <- confint(emmeans(m3, specs = ~"supplement|time")) %>%
        data.frame()


