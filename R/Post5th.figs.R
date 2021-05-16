### Humac figure post fifth RT session

# This script is equal to Post5th.change, except from not log-transforming the data to illustrate absolute change scores for figure illustration of 
# change from pre-post per humac test per supplement. For discription of codes, see Post5th.change.

# Packages
library(readxl);library(tidyverse);library(nlme);library(lme4);library(broom);library(knitr);library(emmeans)

## Data handling

humac <- read_excel("./data/tests/ribose.humac.xlsx", na = "NA") %>%
        
        
        mutate(time = if_else(timepoint == "D-1", 
                              "baseline", 
                              if_else(timepoint %in% c("D4", "D5"), 
                                      "test1", 
                                      if_else(timepoint %in% c("D8", "D9"), 
                                              "test2", 
                                              if_else(timepoint %in% c("T3", "T4") & acute %in% c("rest", "post30min", "post2h"),
                                                      "test3",
                                                      if_else(acute == "post23h", "test4", timepoint)))))) %>%
        mutate(time = factor(time, levels = c("baseline", "test1", "test2", "test3", "test4")), 
               acute = factor(acute, levels = c("rest", "post30min", "post2h", "post23h")), 
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        print()

rest.dat <- humac %>%
        filter(acute == "rest" ) %>%
        print()

## Change data

# Isometric
isom.dat <- rest.dat %>%
        filter(test == "isom") %>%
        print()

change_dat <- isom.dat %>%
        dplyr::select(subject, time, supplement, peak.torque) %>%
        group_by(subject, time, supplement) %>%
        summarise(peak.torque = mean(peak.torque, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = peak.torque) %>%
        
        ungroup() %>%
        mutate(change.2 = test1-baseline,
               change.3 = test2-baseline,
               change.4 = test3-baseline,
               baseline = baseline - mean(baseline, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.4)) %>%
        print()

# Isok.60

isok60.dat <- rest.dat %>%
        filter(test == "isok.60") %>%
        print()

change_dat2 <- isok60.dat %>%
        dplyr::select(subject, time, supplement, peak.torque) %>%
        group_by(subject, time, supplement) %>%
        summarise(peak.torque = mean(peak.torque, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = peak.torque) %>%
        
        ungroup() %>%
        mutate(change.2 = test1-baseline,
               change.3 = test2-baseline,
               change.4 = test3-baseline,
               baseline = baseline - mean(baseline, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.4)) %>%
        print()

## Isok.240

isok240.dat <- rest.dat %>%
        filter(test == "isok.240") %>%
        print()

change_dat3 <- isok240.dat %>%
        dplyr::select(subject, time, supplement, peak.torque) %>%
        group_by(subject, time, supplement) %>%
        summarise(peak.torque = mean(peak.torque, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = peak.torque) %>%
        
        ungroup() %>%
        mutate(change.2 = test1-baseline,
               change.3 = test2-baseline,
               change.4 = test3-baseline,
               baseline = baseline - mean(baseline, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.2, change.3, change.4) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.2:change.4)) %>%
        print()


## Linear mixed effects model

# Isometric
m1 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat)
plot(m1)

summary(m1)

# Isok.60
m2 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat2)
plot(m2)

summary(m2)

# Isok.240
m3 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat3)
plot(m3)

summary(m3)

## Emmeans

# Isometric

confint.m1 <- confint(emmeans(m1, specs = ~"supplement|time")) %>%
        data.frame()

# Isok.60

confint.m2 <- confint(emmeans(m2, specs = ~"supplement|time")) %>%
        data.frame()

# Isok.240

confint.m3 <- confint(emmeans(m3, specs = ~"supplement|time")) %>%
        data.frame()

# Figure 
# Plots fold change of estimated marginal means per supplement and saves the figures to be used in the thesis

pos <- position_dodge(width = 0.2)

# Isom
isom5.plot <- confint.m1 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
        annotate("text", x = "change.1", y = 49, label = "p = NS", size = 3) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "Test 1", "change.3" = "Test 2",
                                  "change.4" = "Test 3")) +
        labs(x = "", y = "Isometric \n(Nm change)\n", fill = "Supplement") +
        theme_classic() +
        theme(axis.line.x = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())

saveRDS(isom5.plot, "./data/derivedData/isom5.plot.RDS")

# Isok 60
isok605.plot <- confint.m2 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "Test 1", "change.3" = "Test 2",
                                  "change.4" = "Test 3")) +
        labs(x = "", y = "Isokinetic 60 \n(Nm change)\n", fill = "Supplement") +
        theme_classic() +
        theme(axis.line.x = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())

saveRDS(isok605.plot, "./data/derivedData/isok605.plot.RDS")

# Isok 240
isok2405.plot <- confint.m3 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
       # annotate("text", x = "change.2", y = 49, label = "Between legs: p > 0.05", size = 2) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "Session 2", "change.3" = "Session 4",
                                  "change.4" = "Session 6")) +
        labs(x = "", y = "Isokinetic 240 \n(Nm change)\n", fill = "Supplement") +
        theme_classic() 
       # theme(axis.text.x = element_text(angle = 45, hjust = 1))

saveRDS(isok2405.plot, "./data/derivedData/isok2405.plot.RDS")


