#### Humac figure post sixth RT session


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



## Data handling

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


## Change data

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
        mutate(change.2 = post30min-rest,
               change.3 = post2h-rest,
               change.4 = post23h-rest,
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
        mutate(change.2 = post30min-rest,
               change.3 = post2h-rest,
               change.4 = post23h-rest,
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
        mutate(change.2 = post30min-rest,
               change.3 = post2h-rest,
               change.4 = post23h-rest,
               baseline = rest - mean(rest, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.4)) %>%
        print()


## Create a model
# This model tries to explain the change by time and supplement, accounting for baseline.
# It produces results on both the time effect and the difference between the groups at any timepoint

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

### Get estimated means from the model, these are average increase at 
# pre = 0 (the average pre value)
# remember that these are absolute (kg) values (changeble in the mutate function)

# Isometric

confint.m1 <- confint(emmeans(m1, specs = ~"supplement|time")) %>%
        data.frame()

# Isok.60

confint.m2 <- confint(emmeans(m2, specs = ~"supplement|time")) %>%
        data.frame()

# Isok.240

confint.m3 <- confint(emmeans(m3, specs = ~"supplement|time")) %>%
        data.frame()

## Emmeans figure 

pos <- position_dodge(width = 0.2)

# Isom
isom.plot <- confint.m1 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "30min post", "change.3" = "2hr post",
                                  "change.4" = "23hr post")) +
        labs(x = "", y = "Isometric \n(nm)\n", title = "Peak torque change post 6th RT session", fill = "Supplement") +
        theme_classic() +
        #theme(plot.background = element_rect(fill = "gray80")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(legend.position = "none")

# Isok 60
isok60.plot <- confint.m2 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "30min post", "change.3" = "2hr post",
                                  "change.4" = "23hr post")) +
        labs(x = "", y = "Isokinetic 60 Peak Torque \n(nm)\n", fill = "Supplement") +
        theme_classic() +
        #theme(plot.background = element_rect(fill = "gray80")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(legend.position = "none")

# Isok 240
isok240.plot <- confint.m3 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "30min post", "change.3" = "2hr post",
                                  "change.4" = "23hr post")) +
        labs(x = "", y = "Isokinetic 240 Peak Torque \n(nm)\n", fill = "Supplement") +
        theme_classic() +
        #theme(plot.background = element_rect(fill = "gray80")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Joined
library(cowplot)

plot_grid(isom.plot, isok60.plot, isok240.plot, nrow = 3, rel_widths = c(0.5,0.5,1))
