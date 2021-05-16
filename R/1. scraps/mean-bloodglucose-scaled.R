#### Glucose #### Glucose analysis 3

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

## Fitting a model

m1 <- lmer(glu ~ sample_time * supplement + (1|subject), data = glu.dat)
m1.log <- lmer(log(glu) ~ sample_time * supplement + (1|subject), data = glu.dat)

plot(m1)
plot(m1.log)
summary(m1)
confint(m1)

## Figure

# With caption 
pos <- position_dodge(width = 0.2)

emmeans(m1, specs = ~ "sample_time|supplement") %>%
        data.frame() %>%
        mutate(sample_time = as.numeric(as.character(sample_time))) %>%
        ggplot(aes(sample_time, emmean, group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin= lower.CL, ymax = upper.CL), 
                      position = pos,
                      width = 0.1) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos,
                   size = 3) +
        labs(x = "Time-point", y = "Blood glucose \n(mmol/L)\n", color = "Supplement",
             caption = "Values are mean and CL") +
        theme_bw() +
        theme(plot.background = element_rect(fill = "gray80")) +
       # scale_y_continuous(limits = c(1, 8), breaks = c(1, 2, 3, 4, 5, 6, 7, 8),
                     #      expand = expand_scale(0)) +
        scale_x_continuous(limits = c(0, 300), breaks = c(0, 45, 90, 120, 135, 150, 270),
                           expand = expand_scale(0))

# Without caption and background color - probably better for adding arrows and such in power point

emmeans(m1, specs = ~ "sample_time|supplement") %>%
        data.frame() %>%
        mutate(sample_time = as.numeric(as.character(sample_time))) %>%
        ggplot(aes(sample_time, emmean, group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin= lower.CL, ymax = upper.CL),
                      position = position_dodge(width = 0.2), 
                      width = 0.1) +
        geom_line(position = position_dodge(width = 0.2)) +
        geom_point(shape = 21, position = position_dodge(width = 0.2),
                   size = 3) +
        labs(x = NULL, y = "Blood glucose \n(mmol/L)\n", fill = "Supplement") +
        theme_classic() +
        #theme(plot.background = element_rect(fill = "gray80")) +
        # scale_y_continuous(limits = c(1, 8), breaks = c(1, 2, 3, 4, 5, 6, 7, 8),
        #      expand = expand_scale(0)) +
        scale_x_continuous(limits = c(0, 300), breaks = c(0, 45, 90, 120, 135, 150, 270),
                           expand = expand_scale(0))


