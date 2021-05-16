## Training intensity as  absolute and fold change analysis

# This script provides log-fold change score calculation of training intensity as %1RM, and illustration of absolute and fold change from baseline to 
# session 6. First, data is handled, then summarized to mean ± SD to illustrate changes in absolute values. Subsequently, log-fold change scores 
# are calculated, and a linear mixed effects models is used together with emmeans to analyze the interaction between time- and/or supplement on 
# changes in training intensity. 

# Packages
library(readxl);library(tidyverse);library(knitr);library(lme4);library(broom);library(emmeans)

# Data
tot.load <- read_excel("./data/training/ribose.training.xlsx", na = "NA") %>%
        select(subject, timepoint, supplement, exercise, load)

tot.rm <- read_excel("./data/tests/ribose.1rm.xlsx", na = "NA") %>%
        select(subject, supplement, exercise, rm) 

# ## Handling the data by creating a new factor called time from timepoint. This factor combines any observation at T1 and T2 to baseline, etc. 
# The code also sorts the order of the factor time, from baseline to session 6, using time = factor(time, levels c()), and sets placebo to be compared to 
# glucose via supplement = factor(supplement, levels = c()).

# Total 
tot.loadh <- tot.load %>%
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
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        print()


## Baseline analysis - comparison of the two legs
# A baseline analysis comparing training intensity at baseline sessions between the two legs via a paired t.test, and providing a summary of mean training
# intensity and sd

base.joined <- tot.loadh %>%
        full_join(tot.rm) %>%
        mutate(p.rm = (load/rm)*100) %>%
        filter(time == "baseline") %>%
        select(subject, time, supplement, p.rm, exercise) %>%
        group_by(subject, supplement, exercise) %>%
        summarise(m = mean(p.rm)) %>%
        pivot_wider(names_from = supplement,
                    values_from = m) %>%
        print()

rm.ttest <- t.test(base.joined$glucose, base.joined$placebo, paired = TRUE)

rm.summary <- tot.loadh %>%
        full_join(tot.rm) %>%
        mutate(p.rm = (load/rm)*100) %>%
        filter(time == "baseline") %>%
        select(subject, time, supplement, p.rm, exercise) %>%
        group_by(supplement, exercise) %>%
        summarise(m = mean(p.rm),
               s = sd(p.rm)) %>%
        print()

## Join data by full_join() and calulates %1RM as a measure of training intensity

tot.joined <- tot.loadh %>%
        full_join(tot.rm) %>%
        mutate(p.rm = (load/rm)*100) %>%
        print()
        
mean.joined <- tot.joined %>%
        select(subject, supplement, time, p.rm) %>%
        group_by(supplement, time) %>%
        summarise(mean.prm = mean(p.rm),
                  sd.prm = sd(p.rm)) %>%
        print()


## Absolute data - summarising mean and SD for barplot
# Mean and SD are calculated from "absolute" numbers of training intensity, creating a basis for a barplot. 

tot.barplot <- ggplot(mean.joined, aes(fill = supplement, y = mean.prm, x = time)) +
        annotate("text", x = c("session2", "session3", "session4", "session5", "session6"),
                 y = c(82, 86, 90, 92, 94), label = "†") +
        geom_bar(position = "dodge", stat = "identity") +
        geom_errorbar(aes( ymin=mean.prm-sd.prm, ymax=mean.prm+sd.prm),
                     width = 0.2,
                     position = position_dodge(width = 1)) +
        labs(x = "", y = "Training intensity \n(%1RM)\n", fill = "") +
        scale_x_discrete(labels=c("baseline" = "Baseline", "session2" = "Session 2", "session3" = "Session 3",
                                  "session4" = "Session 4", "session5" = "Session 5", 
                                  "session6" = "Session 6")) +
        scale_y_continuous(limits = c(0, 100), breaks = c(0, 10, 20, 30, 40, 50, 60,
                                                          70, 80, 90, 100),
                           expand = expansion(0)) +
        theme_classic() +
        theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1))

saveRDS(tot.barplot, "./data/derivedData/tot.barplot.RDS")


## Change data
# The code beneath summarizes the mean values at each time, grouped by subject, time and supplement, creating a wider data set with observations of 
# participants glucose measurements per time point.
# Then, mutate() is used to calculate change scores, where each timepoint is log-transformed and compared to baseline. baseline = baseline - mean(baseline,
# na.rm = TRUE) mean centers the baseline values. Subject, supplement, baseline and change scores are then selected and pivoted for modeling.

tot.change <- tot.joined %>%
        dplyr::select(subject, time, p.rm, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(mean.prm = mean(p.rm)) %>%
        pivot_wider(names_from = time,
                    values_from = mean.prm) %>%
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

# Total 
m1 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = tot.change)
plot(m1)

summary(m1)

## Fold-change estimated means
# Gets estimated means from the model, these are average increase at pre = 0 (the average pre value).
# These are log-fold change values (changeble with the mutate function), reverse transformed with exp() for illustration

confint.m1 <- confint(emmeans(m1, specs = ~"supplement|time")) %>%
        data.frame() %>%
        print()

# Emmeans figure 
pos <- position_dodge(width = 0.2) # creates a position dodge

totalint.plot <- confint.m1 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, exp(emmean), group = supplement, fill = supplement)) +
        annotate("text", x = "change.2", y = 1.29, label = "p = NS", size = 2.5) +
        geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "Session 2", "change.3" = "Session 3",
                                  "change.4" = "Session 4", "change.5" = "Session 5", 
                                  "change.6" = "Session 6")) +
        scale_y_continuous(limits = c(1, 1.3), breaks = c(1.05, 1.10, 1.15, 1.20, 1.25, 1.30),
                           expand = expansion(0)) +
        labs(x = "", y = "Training intensity \n(Fold change)\n", fill = "") +
        theme_classic() +
        #theme(plot.background = element_rect(fill = "gray80")) +
        theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1))

saveRDS(totalint.plot, "./data/derivedData/totalint.plot.RDS")

## Annotate code
# Over one time point: annotate("text", x = "change.2", y = 1.3, label = "Between legs: p > 0.05", size = 2.5) +
# Over all time points:
# annotate("text", x = c("change.2", "change.3", "change.4", "change.5", "change.6"),
#y = c(1.13, 1.17, 1.21, 1.24, 1.28), label = "p > 0.05", size = 2.5) +
