## Session score table

# This script creates a table of the RPE analysis

## Packages
library(readxl);library(tidyverse);library(knitr);library(lme4);library(broom);library(emmeans)

## Data handling

# Session score
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


## Mean and SD

# Session score
ss.exp <- ss.dat %>%
        select(subject, supplement, time, session.score) %>%
        group_by(supplement, time) %>%
        summarise(mean.ss = mean(session.score),
                  sd.ss = sd(session.score)) %>%
        mutate(stat = paste0(round(mean.ss, 1), " (", round(sd.ss, 1), ")")) %>%
        select(supplement, time, stat) %>%
        pivot_wider(names_from = time,
                    values_from = stat) %>%
        print()

ss.exp %>%
        kable(col.names = c("Supplement", "Baseline", "Session 2", "Session 3", "Session 4", "Session 5",
                            "Session 6"))

        




