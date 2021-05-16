#### Volume session 1-5 analysis

## Packages
library(readxl);library(tidyverse);library(knitr);library(lme4);library(broom);library(emmeans)


## Data handling
# Volume
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



## Change data

change_dat.vol <- tot.vol %>%
        dplyr::select(subject, time, tot.volume, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(tot.volume = mean(tot.volume, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = tot.volume) %>%
        #print()
        
        ungroup() %>%
        mutate(change.2 = log(session2)-log(baseline),
               change.3 = log(session3)-log(baseline),
               change.4 = log(session4)-log(baseline),
               change.5 = log(session5)-log(baseline),
               change.post = log(post)-log(baseline),
               baseline = baseline - mean(baseline, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4, change.5, change.post) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.post)) %>%
        print()


## Create a model
# This model tries to explain the change by time and supplement, accounting for baseline.
# It produces results on both the time effect and the difference between the groups at any timepoint

# Mean of all subjects 
m1 <- lmerTest::lmer(change ~ baseline + time + supplement:time + (1|subject),
                     data = change_dat.vol)
plot(m1)

summary(m1)


### Get estimated means from the model, these are average increase at 
# pre = 0 (the average pre value)
# remember that these are absolute (kg) values (changeble in the mutate function)

confint.m1 <- confint(emmeans(m1, specs = ~"supplement|time")) %>%
        data.frame()


## Emmeans figure mean all subjects
# exp() funksjonen tilbaketransformerer dataene, altså viser tallene nå prosentvis forandring (1.1 = 10%)
pos <- position_dodge(width = 0.2)

 confint.m1 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, exp(emmean), group = supplement, fill = supplement)) +
        #geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)),
         #             position = pos,
          #            width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "Session 2", "change.3" = "Session 3",
                                  "change.4" = "Session 4", "change.5" = "Session 5", 
                                  "change.post" = "Session 6")) +
        labs(x = "Time-point", y = "Absolute mean change in total session volume \n(log-kg)\n", fill = "Supplement") +
        theme_classic() +
        theme(plot.background = element_rect(fill = "gray80"))


## Table
cbind(coef(summary(m1)), data.frame(confint(m1))[1:11, ])

# Comments: 
# 1. Include baseline in the figure, otherwise done


