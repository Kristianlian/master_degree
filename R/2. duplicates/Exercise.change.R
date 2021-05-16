#### Volume session 1-5 analysis, per exercise

## Packages
library(readxl);library(tidyverse);library(knitr);library(lme4);library(broom);library(emmeans);library(sjPlot)

# Data
lp.vol <- read_excel("./data/training/ribose.volume.xlsx", na = "NA") %>%
        select(subject, timepoint, lp.volume, supplement) 

ke.vol <- read_excel("./data/training/ribose.volume.xlsx", na = "NA") %>%
        select(subject, timepoint, ke.volume, supplement) 

## Data handling
# Leg press
lp.volh <- lp.vol %>%
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
                                                                      "session6", timepoint))))))) %>%
        mutate(time = factor(time, levels = c("baseline", "session1", "session2", "session3", "session4", "session5", "session6")),
               supplement = factor(supplement, levels = c("placebo", "glucose")))

# Knee extension
ke.volh <- ke.vol %>%
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
                                                                      "session6", timepoint))))))) %>%
        mutate(time = factor(time, levels = c("baseline", "session1", "session2", "session3", "session4", "session5", "session6")),
               supplement = factor(supplement, levels = c("placebo", "glucose")))

# Leg press
# log() = log transforms the data, gives more precise analysis

change_dat.vol <- lp.volh %>%
        dplyr::select(subject, time, lp.volume, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(lp.volume = mean(lp.volume, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = lp.volume) %>%
        #print()
        
        ungroup() %>%
        mutate(change.2 = log(session2)-log(baseline),
               change.3 = log(session3)-log(baseline),
               change.4 = log(session4)-log(baseline),
               change.5 = log(session5)-log(baseline),
               change.6 = log(session6)-log(baseline),
               baseline = baseline - mean(baseline, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4, change.5, change.6) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.6)) %>%
        print()

# Knee extension

change_dat.vol2 <- ke.volh %>%
        dplyr::select(subject, time, ke.volume, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(ke.volume = mean(ke.volume, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = ke.volume) %>%
       # print()
        
        ungroup() %>%
        mutate(change.2 = log(session2)-log(baseline),
               change.3 = log(session3)-log(baseline),
               change.4 = log(session4)-log(baseline),
               change.5 = log(session5)-log(baseline),
               change.6 = log(session6)-log(baseline),
               baseline = baseline - mean(baseline, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4, change.5, change.6) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.6)) %>%
        print()

## Create a model
# This model tries to explain the change by time and supplement, accounting for baseline.
# It produces results on both the time effect and the difference between the groups at any timepoint

# Leg press
m1 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat.vol)
plot(m1)

summary(m1)

# Knee extension
m2 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat.vol2)
plot(m2)

summary(m2)


### Get estimated means from the model, these are average increase at 
# pre = 0 (the average pre value)
# remember that these are absolute (kg) values (changeble in the mutate function)

# Leg press
confint.m1 <- confint(emmeans(m1, specs = ~"supplement|time")) %>%
        data.frame()

# Knee extension
confint.m2 <- confint(emmeans(m2, specs = ~"supplement|time")) %>%
        data.frame()

## Emmeans figure mean all subjects
# exp() funksjonen tilbaketransformerer dataene, altså viser tallene nå prosentvis forandring (1.1 = 10%)
pos <- position_dodge(width = 0.2)

# Leg press
confint.m1 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, exp(emmean), group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)),
                     position = pos,
                     width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "Session 2", "change.3" = "Session 3",
                                  "change.4" = "Session 4", "change.5" = "Session 5", 
                                  "change.6" = "Session 6")) +
        labs(x = "Time-point", y = "Leg press training volume \n(fold change kg x reps)\n", fill = "Supplement") +
        theme_classic() +
        #theme(plot.background = element_rect(fill = "gray80")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Knee extension
confint.m2 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, exp(emmean), group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "Session 2", "change.3" = "Session 3",
                                  "change.4" = "Session 4", "change.5" = "Session 5", 
                                  "change.6" = "Session 6")) +
        labs(x = "Time-point", y = "Knee extension training volume \n(fold change kg x reps)\n", fill = "Supplement") +
        theme_classic() +
        #theme(plot.background = element_rect(fill = "gray80")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

# If you want % change, add code below above ggplot
#mutate(change = (exp(emmean)-1)*100,
#      lower.ci = (exp(lower.CL)-1)*100,
#     upper.ci = (exp(upper.CL)-1)*100) %>%

## Table
cbind(coef(summary(m1)), data.frame(confint(m1))[1:11, ])

tab_model(m1, m2,
          pred.labels = c("Baseline", "Baseline-S2", "Baseline-S3", "Baseline-S4", 
                          "Baseline-S5", "Baseline-S6", "Group diff S2", "Group diff S3",
                          "Group diff S4", "Group diff S5", "Group diff S6"),
          string.pred = "Coefficient",
          string.ci = " Conf.Int (95%)",
          string.p = "P-value",
          dv.labels = c("LP change", "KE change"),
          show.re.var = FALSE,
          show.icc = FALSE)

# Fiks labels tabell, er verdiene i absolutte nå? ellers goooooooooooooood

