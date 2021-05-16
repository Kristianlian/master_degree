#### RPE analysis

# This script analyzes log-fold change in RPE through the intervention per supplement. a linear mixed effects models is used together with emmeans to 
# analyze the interaction between time- and/or supplement on changes in training intensity. 

## Packages
library(readxl);library(tidyverse);library(knitr);library(lme4);library(broom);library(emmeans)

## Handling the data by creating a new factor called time from timepoint. This factor combines any observation at T1 and T2 to baseline, etc. 
# The code also sorts the order of the factor time, from baseline to session 6, using time = factor(time, levels c()), and sets placebo to be compared to 
# glucose via supplement = factor(supplement, levels = c()).

ss.dat <- read_excel("./data/training/ribose.rpe.xlsx", na = "NA") %>%
        select(subject, timepoint, session.score, supplement) %>%
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


## Absolute data - summarising mean and SD for barplot
# Mean and SD are calculated from absolute (kg) total volume, creating a basis for a barplot. 

ss.exp <- ss.dat %>%
        select(subject, supplement, time, session.score) %>%
        group_by(supplement, time) %>%
        summarise(mean.ss = mean(session.score),
                  sd.ss = sd(session.score)) %>%
        print()

ss.barplot <- ggplot(ss.exp, aes(fill = supplement, y = mean.ss, x = time)) +
        geom_bar(position = "dodge", stat = "identity") +
        #geom_errorbar(aes(x = time, ymin = mean.prm - sd.prm, ymax = mean.prm + sd.prm),
        #             width = 0.2)
        labs(x = "", y = "RPE \n(1-10)\n", fill = "Supplement") +
        scale_x_discrete(labels=c("baseline" = "", "session2" = "", "session3" = "",
                                  "session4" = "", "session5" = "", 
                                  "post" = "")) +
        scale_y_continuous(limits = c(0, 10), breaks = c(0, 1, 2, 3, 4, 5, 6,
                                                          7, 8, 9, 10),
                           expand = expansion(0)) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

saveRDS(ss.barplot, "./data/derivedData/ss.barplot.RDS")

## Change-data
# The code beneath summarizes the mean values at each time, grouped by subject, time and supplement, creating a wider data set with observations of 
# participants glucose measurements per time point.
# Then, mutate() is used to calculate change scores, where each timepoint is log-transformed and compared to baseline. baseline = baseline - mean(baseline,
# na.rm = TRUE) mean centers the baseline values. Subject, supplement, baseline and change scores are then selected and pivoted for modeling.

change_dat.ss <- ss.dat %>%
        dplyr::select(subject, time, session.score, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(session.score = mean(session.score, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = session.score) %>%
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



## Linear mixed effects model
# This model tries to explain the change by time and supplement, accounting for potential differences in baseline values and that the same participants
# are measured at multiple time points. 
# It produces results on both the time effect and the difference between the groups at any timepoint. We are interested in the difference between groups.

m1 <- lmerTest::lmer(change ~ baseline + time + supplement:time + (1|subject),
                     data = change_dat.ss)
plot(m1)

summary(m1)

## Fold-change estimated means
# Gets estimated means from the model, these are average increase at pre = 0 (the average pre value).
# These are log-fold change values (changeble with the mutate function)

confint.m1 <- confint(emmeans(m1, specs = ~"supplement|time")) %>%
        data.frame() %>%
        print()

# Figure
# Plots fold change of estimated marginal means per supplement and saves the figures to be used in the thesis

pos <- position_dodge(width = 0.2) # creates a position dodge

# Session score figure
ss.plot <- confint.m1 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, exp(emmean), group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "", "change.2" = "", "change.3" = "",
                                  "change.4" = "", "change.5" = "", 
                                  "change.post" = "")) +
        labs(x = "", y = "RPE (1-10) \n(Fold-change)\n", fill = "Supplement") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

saveRDS(ss.plot, "./data/derivedData/ss.plot.RDS")


