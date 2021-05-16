#### Total training session volume pre vs. post

## Author: Kristian Lian
## Project: Ribose

## Purpose: This script plots and analyses mean development of total training session volume per supplement, both pre vs. post and
# from session 1-5

## Time-points
# T1 and T2: First session of either leg #1 (T1) and leg #2 (T2), aka baseline
# D3&4: Session 2
# D5&6: Session 3
# D7&8: Session 4
# D9&10: Session 5
# T3&T4: Post testing + session 6

## Data contains
# Subject
# Leg: left or right
# lp.volume: Session volume for leg press
# ke.volume: Session volume for knee extension
# tot.volume: Total volume for both leg press and knee extension
# Supplement: placebo or glucose
# Weight: Subjects weight
# Height: Subjects height


## Packages
library(readxl);library(tidyverse);library(knitr);library(lme4);library(broom);library(emmeans)


## Data handling
tot.vol <- read_excel("./data/training/ribose.volume.xlsx", na = "NA") %>%
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

# Total volume pre vs. post
pos <- position_dodge(width = 0.2)


tot.vol %>%
        filter(time %in% c("baseline", "post"),
               subject != "115") %>%
        ggplot(aes(time, tot.volume, color = supplement)) +
        geom_boxplot() +
        labs(x = "Time-point", y = "Volume (kg)", color = "Supplement",
             title = "Mean session volume developement per supplement",
             subtitle = "Different patterns per supplement")

# Total volume session 1-5

tot.vol %>%
        filter(time != "post") %>%
        ggplot(aes(time, tot.volume, color = supplement)) +
        geom_boxplot() +
        labs(x = "Time-point", y = "Volume (kg)", color = "Supplement",
             title = "Mean session volume developement per supplement",
             subtitle = "Different patterns per supplement")
        
        
## Models

# Total volume pre vs. post
m1 <- lmer(tot.volume ~ time * supplement + (1|subject),
           data = filter(tot.vol, time %in% c("baseline", "post")))
plot(m1)
# Try log transforming, data seems uneven
m1l <- lmer(log(tot.volume) ~ time * supplement + (1|subject),
           data = filter(tot.vol, time %in% c("baseline", "post")))
plot(m1l)

# Didn´t see a whole lot of difference totally, the spread changes a bit but still seems uneven

summary(m1)

confint(m1)

# Get mean estimates

emmeans(m1, specs = ~ "time|supplement") %>%
        data.frame() %>%
        ggplot(aes(time, emmean, color = supplement, group = supplement)) +
        
        geom_errorbar(aes(ymin= lower.CL, ymax = upper.CL), width = 0.2) +
        geom_point() +      
        geom_line() 

pairs(emmeans(m1, specs = ~ "time|supplement"), reverse = TRUE) %>%
        confint()


cbind(coef(summary(m1)), data.frame(confint(m1))[3:6, ])



# Total volume session 1-5
m2 <- lmer(tot.volume ~ time * supplement + (1|subject),
           data = filter(tot.vol, time != "post"))
plot(m2)
# Try log transforming, data seems uneven
m2l <- lmer(log(tot.volume) ~ time * supplement + (1|subject),
            data = filter(tot.vol, time != "post"))
plot(m2l)

# Looks better

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
        labs(x = "Time-point", y = "Total session volum \n(log-kg per session)\n", fill = "Supplement") +
        theme_classic() +
        theme(plot.background = element_rect(fill = "gray80"))


pairs(emmeans(m2l, specs = ~ "time|supplement"), reverse = TRUE) %>%
        confint()


cbind(coef(summary(m2l)), data.frame(confint(m2l))[3:12, ])
