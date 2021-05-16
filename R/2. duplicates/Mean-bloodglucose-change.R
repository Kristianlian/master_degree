#### Blood glucose change analysis

## Author: Kristian Lian
## Project: Ribose

## This script analyses and plots the difference in blood glucose between supplement groups, with x axis
# scaling (using sample_time instead of the mutated factor "time" as in the original analysis)

## Data
# Subject
# Timepoints: T1 (baseline, before supplements and RT), T3 (posttest leg #1), and T4 (posttest leg #2)
# Sample_time: Minutes following protein ingestion (0-270)
# Time: time of day
# Glu: blood glucose, mmol/L
# Lak: blood lactat, mmol/L
# Supplement: glucose or placebo

# Packages
library(tidyverse);library(readxl);library(nlme);library(lme4);library(knitr);library(broom);library(emmeans);library(dplyr)

## Data handling

gluc.dat <- read_excel("./data/glucose/fingerstick.xlsx", na = "NA")
gluc.dat$sample_time <- as.character(gluc.dat$sample_time)

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
# Filter out participant 101, 102, 103? No pre value from T3 and T4, only from T1 - or use the pre value from T1?
# Filter out participant 107, missing value

change_dat.glu <- glu.dat %>%
       # filter(subject != "107") %>%
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

## Create a model
# This model tries to explain the change by time and supplement, accounting for baseline.
# It produces results on both the time effect and the difference between the groups at any timepoint

m1 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat.glu)
plot(m1)

summary(m1)


### Get estimated means from the model, these are average increase at 
# pre = 0 (the average pre value)
# remember that these are fold change values (changeble in the mutate function)


confint.m1 <- confint(emmeans(m1, specs = ~"supplement|time")) %>%
        data.frame()


## Emmeans figure
pos <- position_dodge(width = 0.2)

confint.m1 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before = 2) %>%
        mutate(time.c = gsub("change.", "", time), 
               time.c = as.numeric(time.c)) %>%
        #mutate(time = factor(time, levels = c("change.45", "change.90", "change.120", "change.135", "change.150", "change.270"))) %>%
        #mutate(time = as.numeric(time)) %>%
        ggplot(aes(time.c, exp(emmean), group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)),
                      size = 0.5,
                      position = pos) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 2) +
        #scale_x_discrete(labels=c("change.45" = "45", "change.90" = "90",
         #                         "change.120" = "120", "change.135" = "135", 
          #                        "change.150" = "150", "change.270" = "270")) +
        labs(x = "Time-point", y = "Blood glucose levels \n(fold change mmol/L)\n", fill = "Supplement") +
        theme_classic() +
       # theme(plot.background = element_rect(fill = "gray80")) +
        scale_x_continuous(limits = c(0, 300), breaks = c(0, 45, 90, 120, 135, 150, 270),
                           expand = expand_scale(0)) +
        theme(axis.text.x = element_text(angle = 45))
       


## Table
cbind(coef(summary(m1)), data.frame(confint(m1))[1:13, ])

tab_model(m1,
          string.pred = "Coefficient",
          string.ci = " Conf.Int (95%)",
          string.p = "P-value",
          dv.labels = "Change",
          show.re.var = FALSE,
          show.icc = FALSE)


#

