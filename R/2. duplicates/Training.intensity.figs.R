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

tot.load <- read_excel("./data/training/ribose.training.xlsx", na = "NA") %>%
        select(subject, timepoint, supplement, exercise, load)

tot.rm <- read_excel("./data/tests/ribose.1rm.xlsx", na = "NA") %>%
        select(subject, supplement, exercise, rm) 


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


## Join data and calculate load as %1RM
lp.joined <- lp.loadh %>%
        full_join(lp.rm) %>%
        mutate(p.rm = (load/rm)*100) %>%
        select(subject, supplement, time, p.rm) %>%
               group_by(supplement, time) %>%
              summarise(mean.prm = mean(p.rm),
                       sd.prm = sd(p.rm)) %>%
        print()

ke.joined <- ke.loadh %>%
        full_join(ke.rm) %>%
        mutate(p.rm = (load/rm)*100) %>%
        select(subject, supplement, time, p.rm) %>%
        group_by(supplement, time) %>%
        summarise(mean.prm = mean(p.rm),
                  sd.prm = sd(p.rm)) %>%
        print()

# Total

tot.joined <- tot.loadh %>%
        full_join(tot.rm) %>%
        mutate(p.rm = (load/rm)*100) %>%
        select(subject, supplement, time, p.rm) %>%
        group_by(supplement, time) %>%
        summarise(mean.prm = mean(p.rm),
                  sd.prm = sd(p.rm)) %>%
        print()

## Figures

# Leg press
lp.plot <- ggplot(lp.joined, aes(fill = supplement, y = mean.prm, x = time)) +
        geom_bar(position = "dodge", stat = "identity") +
        #geom_errorbar(aes(x = time, ymin = mean.prm - sd.prm, ymax = mean.prm + sd.prm),
         #             width = 0.2)
        labs(x = "Time-point", y = "LP Exercise intensity \n(%1RM)\n", fill = "Supplement") +
        scale_y_continuous(limits = c(0, 100), breaks = c(0, 10, 20, 30, 40, 50, 60,
                                                          70, 80, 90, 100),
                           expand = expand_scale(0)) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(legend.position = "none")
        
# Knee extension
ke.plot <- ggplot(ke.joined, aes(fill = supplement, y = mean.prm, x = time)) +
        geom_bar(position = "dodge", stat = "identity") +
        #geom_errorbar(aes(x = time, ymin = mean.prm - sd.prm, ymax = mean.prm + sd.prm),
        #             width = 0.2)
        labs(x = "Time-point", y = "KE exercise intensity \n(%1RM)\n", fill = "Supplement") +
        scale_y_continuous(limits = c(0, 100), breaks = c(0, 10, 20, 30, 40, 50, 60,
                                                        70, 80, 90, 100),
                           expand = expand_scale(0)) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Joined figure

lpvol.plot <- readRDS("./data/derivedData/lpvol.plot.RDS")

kevol.plot <- readRDS("./data/derivedData/kevol.plot.RDS")

library(cowplot)

plot_grid(lpvol.plot, kevol.plot, lp.plot, ke.plot, labels = c("A", "B", "C", "D"), ncol = 2,
          nrow = 2, rel_widths = c(1, 1.3, 1, 1.3))




# Alternative code, just for mean and SD:
#select(subject, supplement, time, p.rm) %>%
#       group_by(supplement, time) %>%
#      summarise(Mean.prm = mean(p.rm),
#               sd.prm = sd(p.rm)) %>%

