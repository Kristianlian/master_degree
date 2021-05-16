### Total volume change, both absolute values and fold-change

# This script provides log-fold change score calculation of training volume, and illustration of absolute and fold change in training volume. 
# Absolute changes are presented as mean change ± SD per supplement, fold changes are represented as mean change ± CI. 

# Packages
library(readxl);library(tidyverse);library(knitr);library(lme4);library(broom);library(emmeans)

# Data
tot.vol <- read_excel("./data/training/ribose.volume.xlsx", na = "NA") %>%
        select(subject, timepoint, tot.volume, supplement) 

## Handling the data by creating a new factor called time from timepoint. This factor combines any observation at T1 and T2 to baseline, etc. 
# The code also sorts the order of the factor time, from baseline to session 6, using time = factor(time, levels c()), and sets placebo to be compared to 
# glucose via supplement = factor(supplement, levels = c()).

tot.volh <- tot.vol %>%
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

## Baseline analysis - comparison of the two legs
# A baseline analysis comparing training volume at baseline sessions between the two legs via a paired t.test, and providing a summary of mean training
# volume and sd

base.vol <- tot.volh %>%
        filter(time == "baseline") %>%
        select(subject, time, supplement, tot.volume) %>%
        group_by(supplement) %>%
        pivot_wider(names_from = supplement,
                    values_from = tot.volume) %>%
        print()

vol.ttest <- t.test(base.vol$glucose, base.vol$placebo, paired = TRUE)

vol.summary <- tot.volh %>%
        filter(time == "baseline") %>%
        select(subject, time, supplement, tot.volume) %>%
        group_by(supplement) %>%
        mutate(m = mean(tot.volume),
               s = sd(tot.volume)) %>%
        print()

## Absolute data - summarising mean and SD for barplot
# Mean and SD are calculated from absolute (kg) total volume, creating a basis for a barplot. 

vol.exp <- tot.volh %>%
        select(subject, supplement, time, tot.volume) %>% 
        group_by(supplement, time) %>%
        summarise(mean.vol = mean(tot.volume),
                  sd.vol = sd(tot.volume)) %>%
        print()

totvol.barplot <- ggplot(vol.exp, aes(fill = supplement, y = mean.vol, x = time)) +
        annotate("text", x = c("session4", "session5", "session6"),
                 y = c(8500, 8500, 8700), label = "†") +
        geom_bar(position = "dodge", stat = "identity") +
        geom_errorbar(aes( ymin = mean.vol - sd.vol, ymax = mean.vol + sd.vol),
                     width = 0.2,
                     position = position_dodge(width = 1)) +
        labs(x = "", y = "Training volume \n(kg)\n", fill = "Supplement") +
        scale_x_discrete(labels=c("baseline" = "Baseline", "session2" = "Session 2", "session3" = "Session 3",
                                  "session4" = "Session 4", "session5" = "Session 5", 
                                  "session6" = "Session")) +
        #scale_y_continuous(limits = c(0, 10), breaks = c(0, 1, 2, 3, 4, 5, 6,
         #                                                7, 8, 9, 10),
          #                 expand = expand_scale(0)) +
        theme_classic() +
        theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1))


saveRDS(totvol.barplot, "./data/derivedData/totvol.barplot.RDS")

## Change-data
# The code beneath summarizes the mean values at each time, grouped by subject, time and supplement, creating a wider data set with observations of 
# participants glucose measurements per time point.
# Then, mutate() is used to calculate change scores, where each timepoint is log-transformed and compared to baseline. baseline = baseline - mean(baseline,
# na.rm = TRUE) mean centers the baseline values. Subject, supplement, baseline and change scores are then selected and pivoted for modeling.

change_dat <- tot.volh %>%
        dplyr::select(subject, time, tot.volume, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(tot.volume = mean(tot.volume, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = tot.volume) %>%
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


## Linear mixed effects model
# This model tries to explain the change by time and supplement, accounting for potential differences in baseline values and that the same participants
# are measured at multiple time points. 
# It produces results on both the time effect and the difference between the groups at any timepoint. We are interested in the difference between groups.

m1 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat)
plot(m1)

summary(m1)

## Fold-change estimated means
# Gets estimated means from the model, these are average increase at pre = 0 (the average pre value).
# These are log-fold change values (changeble with the mutate function)

confint.m1 <- confint(emmeans(m1, specs = ~"supplement|time")) %>%
        data.frame()

# Figure
pos <- position_dodge(width = 0.2) # creates a position dodge easier to implement

totvolfold.plot <- confint.m1 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, exp(emmean), group = supplement, fill = supplement)) +
        annotate("text", x = "change.2", y = 1.32, label = "p = NS", size = 2.5) +
        geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "Session 2", "change.3" = "Session 3",
                                  "change.4" = "Session 4", "change.5" = "Session 5", 
                                  "change.6" = "Session 6")) +
        labs(x = "", y = "Training volume \n(Fold change)\n", fill = "Supplement") +
        theme_classic() +
        theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1))
        #theme(plot.background = element_rect(fill = "gray80")) 

# Code for annotating each time point
# annotate("text", x = c("change.1", "change.2", "change.3", "change.4", "change.5", "change.6"),
  #       y = c(1.05, 1.2, 1.18, 1.3, 1.27, 1.32), label = "p > 0.05", size = 2.5) +

saveRDS(totvolfold.plot, "./data/derivedData/totvolfold.plot.RDS")

