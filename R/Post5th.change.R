#### Humac analysis pre vs. post fifth RT session


## Author: Kristian Lian

# Purpose: This script plots mean torque per supplement (both through intervention and pre vs. post) results from the ribose project, 
# and analyses the data per test (isometric, isokinetic 60, isokinetic 240) in a linear model.

## Time-points
# D-1: Baseline, before any supplementation or training
# D4, D5, D8 and D9: Day 4, 5, 8 and 9 of the intervention, humac testing of the leg that performed
# RT the preceding day
# T3: Post testing leg #1 (leg that started the intervention). Leg #1 is tested four times at T3/T4:
# Test 1 leg 1: 1.5hrs after protein ingestion, 45min before RT (T3)
# Test 2 leg 1: 30min after RT (T3)
# Test 3 leg 1: 2hrs after RT (T3)
# Test 4 leg 1: ~23hrs after RT (T4)
# Test 1 serve as a post test for the 5 RT sessions and pre test before the sixth session, test 2,
# 3, and 4 serve as post test following sixth session
# T4 and 13 follow the same design for leg #2

## Data
# Date of testing
# Subject
# Test type: isok.60 (isokinetic 60), isok.240 (isokinetic 240), isom (isometric)
# Peak.torque: Highest peak torque from each test
# Leg: left or right leg
# Supplement: glucose or placebo

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

rest.dat <- humac %>%
        filter(acute == "rest" ) %>%
        print()

## Baseline analysis - comparison of the two legs
# A baseline analysis comparing peak torque for each exercise at baseline between the two legs via a paired t.test, and providing a summary of mean peak
# torque and sd

# Isometric
base.isom <- humac %>%
        filter(time == "baseline",
               test == "isom") %>%
        select(subject, time, test, supplement, peak.torque) %>%
        group_by(supplement) %>%
        pivot_wider(names_from = supplement,
                    values_from = peak.torque) %>%
        print()

isom.ttest <- t.test(base.isom$glucose, base.isom$placebo, paired = TRUE)

isom.summary <- humac %>%
        filter(time == "baseline",
               test == "isom") %>%
        select(subject, time, test, supplement, peak.torque) %>%
        group_by(supplement) %>%
        mutate(m = mean(peak.torque),
               s = sd(peak.torque)) %>%
        print()

# Isok 60

base.60 <- humac %>%
        filter(time == "baseline",
               test == "isok.60") %>%
        select(subject, time, test, supplement, peak.torque) %>%
        group_by(supplement) %>%
        pivot_wider(names_from = supplement,
                    values_from = peak.torque) %>%
        print()

isok60.ttest <- t.test(base.60$glucose, base.60$placebo, paired = TRUE)

isok60.summary <- humac %>%
        filter(time == "baseline",
               test == "isok.60") %>%
        select(subject, time, test, supplement, peak.torque) %>%
        group_by(supplement) %>%
        mutate(m = mean(peak.torque),
               s = sd(peak.torque)) %>%
        print()

# Isok 240

base.240 <- humac %>%
        filter(time == "baseline",
               test == "isok.240") %>%
        select(subject, time, test, supplement, peak.torque) %>%
        group_by(supplement) %>%
        pivot_wider(names_from = supplement,
                    values_from = peak.torque) %>%
        print()

isok240.ttest <- t.test(base.240$glucose, base.240$placebo, paired = TRUE)

isok240.summary <- humac %>%
        filter(time == "baseline",
               test == "isok.240") %>%
        select(subject, time, test, supplement, peak.torque) %>%
        group_by(supplement) %>%
        mutate(m = mean(peak.torque),
               s = sd(peak.torque)) %>%
        print()


## Change-data
# The code beneath summarizes the mean values at each time, grouped by subject, time and supplement, creating a wider data set with observations of 
# participants glucose measurements per time point.
# Then, mutate() is used to calculate change scores, where each timepoint is log-transformed and compared to baseline. baseline = baseline - mean(baseline,
# na.rm = TRUE) mean centers the baseline values. Subject, supplement, baseline and change scores are then selected and pivoted for modeling. The data set is
# filtered according to test exercise (isometric, isokinetic 60 or isokinetic 240)

# Isometric
isom.dat <- rest.dat %>%
        filter(test == "isom") %>%
        print()

change_dat <- isom.dat %>%
        dplyr::select(subject, time, supplement, peak.torque) %>%
        group_by(subject, time, supplement) %>%
        summarise(peak.torque = mean(peak.torque, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = peak.torque) %>%
        
        ungroup() %>%
        mutate(change.2 = log(test1)-log(baseline),
               change.3 = log(test2)-log(baseline),
               change.4 = log(test3)-log(baseline),
               baseline = baseline - mean(baseline, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.4)) %>%
        print()

# Isok.60

isok60.dat <- rest.dat %>%
        filter(test == "isok.60") %>%
        print()

change_dat2 <- isok60.dat %>%
        dplyr::select(subject, time, supplement, peak.torque) %>%
        group_by(subject, time, supplement) %>%
        summarise(peak.torque = mean(peak.torque, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = peak.torque) %>%
        
        ungroup() %>%
        mutate(change.2 = log(test1)-log(baseline),
               change.3 = log(test2)-log(baseline),
               change.4 = log(test3)-log(baseline),
               baseline = baseline - mean(baseline, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.4)) %>%
        print()

## Isok.240

isok240.dat <- rest.dat %>%
        filter(test == "isok.240") %>%
        print()

change_dat3 <- isok240.dat %>%
        dplyr::select(subject, time, supplement, peak.torque) %>%
        group_by(subject, time, supplement) %>%
        summarise(peak.torque = mean(peak.torque, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = peak.torque) %>%
        
        ungroup() %>%
        mutate(change.2 = log(test1)-log(baseline),
               change.3 = log(test2)-log(baseline),
               change.4 = log(test3)-log(baseline),
               baseline = baseline - mean(baseline, na.rm = TRUE),
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

# Mean of all subjects 

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
        data.frame() 

# Isok.60

confint.m2 <- confint(emmeans(m2, specs = ~"supplement|time")) %>%
        data.frame() %>%
        print()

# Isok.240

confint.m3 <- confint(emmeans(m3, specs = ~"supplement|time")) %>%
        data.frame()

## Emmeans figures

# Isom
confint.m1 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "Test 1", "change.3" = "Test 2",
                                  "change.4" = "Test 3")) +
        labs(x = "", y = "Isometric \n(nm change)\n", fill = "Supplement") +
        theme_classic() +
        theme(axis.text.x = element_text(size=8))

# Isok 60

confint.m2 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "Test 1", "change.3" = "Test 2",
                                  "change.4" = "Test 3")) +
        labs(x = "", y = "Isokinetic 60 \n(nm change)\n", fill = "Supplement") +
        theme_classic() +
        theme(axis.text.x = element_text(size=8))


# Isok 240

confint.m3 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "Test 1", "change.3" = "Test 2",
                                  "change.4" = "Test 3")) +
        labs(x = "Time-Point", y = "Isokinetic 240 \n(nm change)\n", fill = "Supplement") +
        theme_classic() +
        theme(axis.text.x = element_text(size=8))




