#### Load change score analysis

# Packages
library(readxl);library(tidyverse);library(knitr);library(lme4);library(broom);library(emmeans)


# Data handling

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

saveRDS(tot.load, "./data/derivedData/tot.load.RDS")

## Change data

change_dat.load <- tot.load %>%
        dplyr::select(subject, time, load, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(load = mean(load, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = load) %>%
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

saveRDS(change_dat.load, "./data/derivedData/change_dat.load.RDS")
## Create a model
# This model tries to explain the change by time and supplement, accounting for baseline.
# It produces results on both the time effect and the difference between the groups at any timepoint

m1 <- lmerTest::lmer(change ~ baseline + time + supplement:time + (1|subject),
                     data = change_dat.load)
plot(m1)

summary(m1)


saveRDS(m1, "./data/derivedData/m1.load.RDS")

### Get estimated means from the model, these are average increase at 
# pre = 0 (the average pre value)
# remember that these are absolute (kg) values (changeble in the mutate function)

confint.m1 <- confint(emmeans(m1, specs = ~"supplement|time")) %>%
        data.frame()

saveRDS(confint.m1, "./data/derivedData/confint.load.RDS")

## Emmeans figure
pos <- position_dodge(width = 0.2)


load.fig <- confint.m1 %>%
        data.frame() %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.2" = "Session 2", "change.3" = "Session 3",
                                  "change.4" = "Session 4", "change.5" = "Session 5", 
                                  "change.post" = "Session 6")) +
        labs(x = "Time-point", y = "Absolute mean change in total session load \n(log-kg)\n", fill = "Supplement") +
        theme_classic() +
        theme(plot.background = element_rect(fill = "gray80"))

saveRDS(load.fig, "./data/derivedData/load.fig.RDS")

cbind(coef(summary(m1)), data.frame(confint(m1))[1:11, ])

# Comments: 
# 1. Include baseline in figure, or is it good as is?
        

