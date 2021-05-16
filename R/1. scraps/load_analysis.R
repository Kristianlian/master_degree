#### Training load analysis

## Author: Kristian Lian
## Project: Ribose

# Purpose: This script plots and analyses mean load pre vs. post and from session 1-5 per supplement.

## Time-points
# T1 and T2: First session of either leg #1 (T1) and leg #2 (T2), aka baseline
# D3&4: Session 2
# D5&6: Session 3
# D7&8: Session 4
# D9&10: Session 5
# T3&T4: Post testing + session 6

## Data
# Subject
# Exercise type: unilateral leg press and -knee extension
# 3x10 per exercise
# Load: kg resistance per set
# Leg: left or right leg
# Supplement: glucose or placebo

## Packages
library(readxl);library(tidyverse);library(knitr);library(lme4);library(broom);library(emmeans)

## Data handling

tot.load <- read_excel("./data/training/ribose.training.xlsx", na = "NA") %>%
        mutate(time = if_else(timepoint %in% c("T1", "T2"),
                              "baseline",
                              if_else(timepoint %in% c("D3", "D4"),
                                      "session2",
                                      if_else(timepoint %in% c("D5", "D6"),
                                              "session3",
                                              if_else(timepoint %in% c("D7", "D8"),
                                                      "session4",
                                                      if_else(timepoint %in% c("D9", "D10"),
                                                              "session5",
                                                              if_else(timepoint %in% c("T3", "T4"),
                                                                      "post", timepoint))))))) %>%
        mutate(time = factor(time, levels = c("baseline", "session1", "session2", "session3", "session4", "session5", "post")),
               supplement = factor(supplement, levels = c("placebo", "glucose")))

## Figures

# Total load pre vs. post
pos <- position_dodge(width = 0.2)


tot.load %>%
        filter(time %in% c("baseline", "post")) %>%
        ggplot(aes(time, load, color = supplement)) +
        geom_boxplot() +
        labs(x = "Time-point", y = "Load (kg)", color = "Supplement",
             title = "Mean session load developement per supplement",
             subtitle = "Different patterns per supplement")

# Total load session 1-5

tot.load %>%
        filter(time != "post") %>%
        ggplot(aes(time, load, color = supplement)) +
        geom_boxplot() +
        labs(x = "Time-point", y = "Volume (kg)", color = "Supplement",
             title = "Mean session load developement per supplement",
             subtitle = "Different patterns per supplement")

## Models

# Total load pre vs. post
m1 <- lmer(load ~ time * supplement + (1|subject),
           data = filter(tot.load, time %in% c("baseline", "post")))
plot(m1)
# Try log transforming, data seems uneven
m1l <- lmer(log(load) ~ time * supplement + (1|subject),
            data = filter(tot.load, time %in% c("baseline", "post")))
plot(m1l)

summary(m1l)

confint(m1l)

# Get mean estimates

emmeans(m1l, specs = ~ "time|supplement") %>%
        data.frame() %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
        
        geom_errorbar(aes(ymin= lower.CL, ymax = upper.CL), position = position_dodge(width = 0.2),
                      width = 0.05) +
        #geom_line(position = position_dodge(width = 0.2)) +
        geom_point(shape = 21, position = position_dodge(width = 0.2),
                   size = 3) +
        labs(x = "Time-point", y = "Total session load \n(log-kg per session)\n", fill = "Supplement") +
        theme_classic() +
        theme(plot.background = element_rect(fill = "gray80"))

pairs(emmeans(m1l, specs = ~ "time|supplement"), reverse = TRUE) %>%
        confint()


cbind(coef(summary(m1l)), data.frame(confint(m1l))[3:6, ])



# Total load session 1-5
m2 <- lmer(load ~ time * supplement + (1|subject),
           data = filter(tot.load, time != "post"))
plot(m2)
# Try log transforming, data seems uneven
m2l <- lmer(log(load) ~ time * supplement + (1|subject),
            data = filter(tot.load, time != "post"))
plot(m2l)

summary(m2l)

confint(m2l)

# Get mean estimates

emmeans(m2l, specs = ~ "time|supplement") %>%
        data.frame() %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
        
        geom_errorbar(aes(ymin= lower.CL, ymax = upper.CL), position = position_dodge(width = 0.2),
                      width = 0.05) +
        #geom_line(position = position_dodge(width = 0.2)) +
        geom_point(shape = 21, position = position_dodge(width = 0.2),
                   size = 3) +
        labs(x = "Time-point", y = "Total session load \n(log-kg per session)\n", fill = "Supplement") +
        theme_classic() +
        theme(plot.background = element_rect(fill = "gray80"))

pairs(emmeans(m2l, specs = ~ "time|supplement"), reverse = TRUE) %>%
        confint()


cbind(coef(summary(m2l)), data.frame(confint(m2l))[3:12, ])



