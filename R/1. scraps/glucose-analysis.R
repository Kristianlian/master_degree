#### Glucose analysis 2

## Author: Kristian Lian
## Project: Ribose

## This script analyses and plots the difference in blood glucose between supplement groups 

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
                                      "45min",
                                      if_else(sample_time == "90",
                                              "90min",
                                              if_else(sample_time == "120",
                                                      "120min",
                                                      if_else(sample_time == "135",
                                                              "135min",
                                                              if_else(sample_time == "150",
                                                                      "150min",
                                                                      if_else(sample_time == "270",
                                                                              "270min", sample_time))))))),
               time = factor(time, levels = c("baseline", "45min", "90min", "120min", "135min", "150min", "270min"))) %>%
        mutate(glu = as.numeric(glu),
               lak = as.numeric(lak))


        
## Fitting model

# Without control

m1 <- lmer(glu ~ time * supplement + (1|subject), data = glu.dat)
#summary(m1)
plot(m1)
#confint(m1)

# Try log transforming
m2 <- lmer(log(glu) ~ time * supplement + (1|subject), data = glu.dat)
plot(m2)
summary(m2)
confint(m2)


# Looks better? A bit more spread, but still a couple data points far away

# Figure: mean blood glucose levels per supplement, presented as change in mmol/L
emmeans(m1, specs = ~ "time|supplement") %>%
        data.frame() %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin= lower.CL, ymax = upper.CL), position = position_dodge(width = 0.2), 
                      width = 0.05) +
        geom_line(position = position_dodge(width = 0.2)) +
        geom_point(shape = 21, position = position_dodge(width = 0.2),
                   size = 3) +
        scale_x_discrete(breaks = c("baseline", "45min", "90min", "120min", "135min", "150min", "270min"),
                           labels = c("Baseline", "45min", "90min", "120min", "135min", "150min", "270min")) +
        labs(x = "Time-point", y = "Blood glucose \n(mmol/L)\n", color = "Supplement",
             caption = "Values are mean and CL") +
        theme_bw() +
        theme(plot.background = element_rect(fill = "gray80"))



pairs(emmeans(m2, specs = ~ "time|supplement"), reverse = TRUE) %>%
        confint()



cbind(coef(summary(m2)), data.frame(confint(m2))[3:16, ])


# To do:
# Add arrows for protein ingestion, RT (as e.g. a shaded area), and statistics to figure (power.point)
# Scale timeline to actual distance between measuring



