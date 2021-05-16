#### RPE analysis

## Author: Kristian Lian
## Project: Ribose

# Purpose: This script plots and analyses mean rpe pre vs. post and from session 1-5 per supplement.

## Time-points
# T1 and T2: First session of either leg #1 (T1) and leg #2 (T2), aka baseline
# D3&4: Session 2
# D5&6: Session 3
# D7&8: Session 4
# D9&10: Session 5
# T3&T4: Post testing + session 6

## Data
# Subject
# Leg: left or right leg
# Session score: 1-10 point scale of perceived exertion of the entire session, taken 15min post session
# rpe: Feeling in the leg that exercised preceding day, 1-10 point scale
# vas: 
# Supplement: glucose or placebo

## Packages
library(readxl);library(tidyverse);library(knitr);library(lme4);library(broom);library(emmeans);library(lm.beta)
library(sjPlot)

## Data handling

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

rpe.dat <- read_excel("./data/training/ribose.rpe.xlsx", na = "NA") %>%
        select(subject, timepoint, rpe, supplement) %>%
        filter(timepoint != "T1") %>%
        mutate(time = if_else(timepoint %in% c("T2", "D3"),
                              "baseline",
                              if_else(timepoint %in% c("D4", "D5"),
                                      "session2",
                                      if_else(timepoint %in% c("D6", "D7"),
                                              "session3",
                                              if_else(timepoint %in% c("D8", "D9"),
                                                      "session4",
                                                      if_else(timepoint %in% c("D10", "T3"),
                                                              "session5",
                                                              if_else(timepoint == "T4",
                                                                      "post", timepoint))))))) %>%
        mutate(time = factor(time, levels = c("baseline", "session1", "session2", "session3", "session4", "session5", "post")),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        mutate(rpe = as.numeric(rpe)) %>%
        print()


## Change data
# Session score

ss.tab <- ss.dat %>%
        dplyr::select(subject, time, session.score, supplement) %>%
        group_by(time, supplement) %>%
        summarise(session.score = mean(session.score, na.rm = TRUE),
                  SD = sd(session.score, na.rm = TRUE),
                  n = n()) %>%
        #pivot_wider(names_from = time,
         #            values_from = session.score) %>%
        print()
        

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

# RPE

rpe.tab <- rpe.dat %>%
        dplyr::select(subject, time, rpe, supplement) %>%
        group_by(time, supplement) %>%
        summarise(rpe = mean(rpe, na.rm = TRUE),
                  SD = sd(rpe, na.rm = TRUE),
                  n = n()) %>%
        #pivot_wider(names_from = time,
        #            values_from = session.score) %>%
        print()

change_dat.rpe <- rpe.dat %>%
        dplyr::select(subject, time, rpe, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(rpe = mean(rpe, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = rpe) %>%
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

# Session score model
m1 <- lmerTest::lmer(change ~ baseline + time + supplement:time + (1|subject),
                     data = change_dat.ss)
plot(m1)

summary(m1)


# RPE model
m2 <- lmerTest::lmer(change ~ baseline + time + supplement:time + (1|subject),
                     data = change_dat.rpe)
plot(m2)

summary(m2)

### Get estimated means from the model, these are average increase at 
# pre = 0 (the average pre value)
# remember that these are absolute (kg) values (changeble in the mutate function)

# Session score
confint.m1 <- confint(emmeans(m1, specs = ~"supplement|time")) %>%
        data.frame() %>%
        print()

## Emmeans figure
pos <- position_dodge(width = 0.2)

# Session score figure
confint.m1 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0,
                upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0,
                upper.CL = 0, .before =2) %>%
       # mutate(change = (exp(emmean)-1)*100,
        #       lower.ci = (exp(lower.CL)-1)*100,
         #      upper.ci = (exp(upper.CL)-1)*100) %>%
       # print()
        ggplot(aes(time, exp(emmean), group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "Session 2", "change.3" = "Session 3",
                                  "change.4" = "Session 4", "change.5" = "Session 5", 
                                  "change.post" = "Session 6")) +
        labs(x = "Time-point", y = "Mean change session score \n(change scores)\n", fill = "Supplement") +
        theme_classic() +
        theme(plot.background = element_rect(fill = "gray80"))

saveRDS(confint.m1, "./data/derivedData/ss.fig.RDS")
        
        

# Session score table

# Table 1

cbind(coef(summary(m1)), data.frame(confint(m1))[1:11, ])

tab_model(m1,
          string.pred = "Coefficient",
          string.ci = " Conf.Int (95%)",
          string.p = "P-value",
          dv.labels = "Change",
          show.re.var = FALSE,
          show.icc = FALSE)



# RPE figure

confint.m2 <- confint(emmeans(m2, specs = ~"supplement|time")) %>%
        data.frame()

confint.m2 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0,
                upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0,
                upper.CL = 0, .before =2) %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
     #   geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
      #                position = pos,
       #               width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.1" = "Baseline", "change.2" = "Session 2", "change.3" = "Session 3",
                                  "change.4" = "Session 4", "change.5" = "Session 5", 
                                  "change.post" = "Session 6")) +
        labs(x = "Time-point", y = "Mean change RPE \n(change scores)\n", fill = "Supplement") +
        theme_classic() +
        theme(plot.background = element_rect(fill = "gray80"))

# RPE table

cbind(coef(summary(m2)), data.frame(confint(m2))[1:11, ])

# Table
tab_model(m2,
          string.pred = "Coefficient",
          string.ci = " Conf.Int (95%)",
          string.p = "P-value",
          dv.labels = "Session score (1-10",
          show.re.var = FALSE,
          show.icc = FALSE)

# Comments: 
# Ask Daniel about the plots, there seem to be patterns forming regardless of log-transforming or not
# Is there a better way to show these data? Its to be used as a background variable, along with training volume

