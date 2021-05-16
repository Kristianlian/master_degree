#### Blood glucose change fig

## Author: Kristian Lian
## Project: Ribose

## This script analyses is equal to mean-bloodglucose-change, except from not calculating log-fold change. Instead, absolute change scores are calculated
# to illustrate the changes in absolute plasma glucose values. For description of the codes, see mean-bloodglucose-change.

# Packages
library(dplyr)
library(tidyverse)
library(readxl)
library(nlme)
library(lme4)
library(emmeans)
library(tidyr)


## Data handling

gluc.dat <- read_excel("./data/glucose/fingerstick.xlsx", na = "NA")
gluc.dat$sample_time <- as.character(gluc.dat$sample_time)

glu.dat <- gluc.dat %>%
        mutate(time = if_else(sample_time == "0",
                              "baseline",
                              if_else(sample_time == "45",
                                      "min45",
                                      if_else(sample_time == "90",
                                              "min90",
                                              if_else(sample_time == "120",
                                                      "min120",
                                                      if_else(sample_time == "135",
                                                              "min135",
                                                              if_else(sample_time == "150",
                                                                      "min150",
                                                                      if_else(sample_time == "270",
                                                                              "min270", sample_time))))))),
               time = factor(time, levels = c("baseline", "min45", "min90", "min120", "min135", "min150", "min270"))) %>%
        mutate(glu = as.numeric(glu),
               lak = as.numeric(lak)) %>%
        print()

## Change data

change_dat.glu <- glu.dat %>%
        # filter(subject != "107") %>%
        dplyr::select(subject, time, glu, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(glu = mean(glu, na.rm = TRUE)) %>%
        pivot_wider(names_from = time, 
                    values_from = glu) %>%
        #print()
        
        ungroup() %>%
        mutate(change.45 = min45-baseline,
               change.90 = min90-baseline,
               change.120 = min120-baseline,
               change.135 = min135-baseline,
               change.150 = min150-baseline,
               change.270 = min270-baseline,
               baseline = baseline - mean(baseline, na.rm = TRUE),
               supplement = factor(supplement, levels = c("placebo", "glucose"))) %>%
        select(subject, supplement, baseline, change.45, change.90, change.120, change.135, 
               change.150, change.270) %>%
        pivot_longer(names_to = "time",
                     values_to = "change",
                     cols = (change.45:change.270)) %>%
        print()

## Linear mixed effects model

m1 <- lmerTest::lmer(change ~ 0 + baseline + time + supplement:time + (1|subject),
                     data = change_dat.glu)
plot(m1)

summary(m1)


## Emmeans

confint.m1 <- confint(emmeans(m1, specs = ~"supplement|time")) %>%
        data.frame()


## Figure
# Plots fold change of estimated marginal means per supplement and saves the figures to be used in the thesis

pos <- position_dodge(width = 0.2)

gluc.fig <- confint.m1 %>%
        data.frame() %>%
        add_row(supplement = "placebo", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before =1) %>%
        add_row(supplement = "glucose", time = "change.1", emmean = 0, SE = 0, df = 0, lower.CL = 0, upper.CL = 0, .before = 2) %>%
        mutate(time.c = gsub("change.", "", time), 
               time.c = as.numeric(time.c)) %>%
        #mutate(time = factor(time, levels = c("change.45", "change.90", "change.120", "change.135", "change.150", "change.270"))) %>%
        #mutate(time = as.numeric(time)) %>%
        ggplot(aes(time.c, emmean, group = supplement, fill = supplement)) +
        annotate("text", x = c(120, 135, 150), y = c(2.5, 2.1, 2.1), label = "*") +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      size = 0.5,
                      position = pos) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 2) +
        #scale_x_discrete(labels=c("change.45" = "45", "change.90" = "90",
        #                         "change.120" = "120", "change.135" = "135", 
        #                        "change.150" = "150", "change.270" = "270")) +
        labs(x = "Time-point", y = "Plasma glucose levels \n(mmol/L change)\n", fill = "") +
        theme_classic() +
        # theme(plot.background = element_rect(fill = "gray80")) +
        scale_x_continuous(limits = c(0, 300), breaks = c(0, 45, 90, 120, 135, 150, 270),
                           expand = expansion(0), labels = c("0" = "-120", "45" = "-90", "90" = "-30", "120" = "0", "135" = "15",
                                                             "150" = "30", "270" = "120")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
        

saveRDS(gluc.fig, "./data/derivedData/gluc.fig.RDS")



