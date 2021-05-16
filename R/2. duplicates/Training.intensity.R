## Training intensity as %RM analysis

library(readxl);library(tidyverse);library(knitr);library(lme4);library(broom);library(emmeans)

lp.load <- read_excel("./data/training/ribose.training.xlsx", na = "NA") %>%
        select(subject, timepoint, supplement, exercise, load) %>%
        filter(exercise != "KE") 

lp.rm <- read_excel("./data/tests/ribose.1rm.xlsx", na = "NA") %>%
        select(subject, supplement, exercise, rm) %>%
        filter(exercise != "KE") 

ke.load <- read_excel("./data/training/ribose.training.xlsx", na = "NA") %>%
        select(subject, timepoint, supplement, exercise, load) %>%
        filter(exercise != "LP")

ke.rm <- read_excel("./data/tests/ribose.1rm.xlsx", na = "NA") %>%
        select(subject, supplement, exercise, rm) %>%
        filter(exercise != "LP") 

# Data handling
# Load
lp.loadh <- lp.load %>%
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

ke.loadh <- ke.load %>%
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

## Join data and calculate load as %1RM
lp.joined <- lp.loadh %>%
        full_join(lp.rm) %>%
        mutate(p.rm = (load/rm)*100)

ke.joined <- ke.loadh %>%
        full_join(ke.rm) %>%
        mutate(p.rm = (load/rm)*100) 
        
# Alternative code, just for mean and SD:
#select(subject, supplement, time, p.rm) %>%
#       group_by(supplement, time) %>%
#      summarise(Mean.prm = mean(p.rm),
#               sd.prm = sd(p.rm)) %>%

## Change data

# Leg press
lp.change <- lp.joined %>%
        dplyr::select(subject, time, p.rm, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(mean.prm = mean(p.rm)) %>%
        pivot_wider(names_from = time,
                    values_from = mean.prm) %>%
        ungroup() %>%
        mutate(change.2 = session2-baseline,
               change.3 = session3-baseline,
               change.4 = session4-baseline,
               change.5 = session5-baseline,
               change.6 = session6-baseline,
               baseline = baseline - mean(baseline, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4, change.5, change.6) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.6)) %>%
        print()
        
# Knee extension
ke.change <- ke.joined %>%
        dplyr::select(subject, time, p.rm, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(mean.prm = mean(p.rm)) %>%
        pivot_wider(names_from = time,
                    values_from = mean.prm) %>%
        ungroup() %>%
        mutate(change.2 = session2-baseline,
               change.3 = session3-baseline,
               change.4 = session4-baseline,
               change.5 = session5-baseline,
               change.6 = session6-baseline,
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
                     data = lp.change)
plot(m1)

summary(m1)

# Knee extension
m2 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = ke.change)
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


# Leg press
pos <- position_dodge(width = 0.2)

lp.plot <- confint.m1 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "Session 2", "change.3" = "Session 3",
                                  "change.4" = "Session 4", "change.5" = "Session 5", 
                                  "change.6" = "Session 6")) +
        labs(x = "", y = "Leg press training intensity \n(%1RM change)\n", fill = "Supplement") +
        theme_classic() +
        #theme(plot.background = element_rect(fill = "gray80")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))


ke.plot <- confint.m2 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "Session 2", "change.3" = "Session 3",
                                  "change.4" = "Session 4", "change.5" = "Session 5", 
                                  "change.6" = "Session 6")) +
        labs(x = "", y = "Knee extension training intensity \n(%1RM change)\n", fill = "Supplement") +
        theme_classic() +
        #theme(plot.background = element_rect(fill = "gray80")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Table
tab_model(m1, m2,
          string.pred = "Coefficient",
          string.ci = " Conf.Int (95%)",
          string.p = "P-value",
          dv.labels = c("LP %1RM change", "KE %1RM change"),
          show.re.var = FALSE,
          show.icc = FALSE)

# What does it mean that there are significant effects in the baseline coefficient?
