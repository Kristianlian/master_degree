#### Change pre vs. post volume and load


# Packages
library(readxl);library(tidyverse);library(knitr);library(lme4);library(broom);library(emmeans)

# Load data

tot.load <- readRDS("./data/derivedData/tot.load.RDS")
tot.vol <- readRDS("./data/derivedData/tot.vol.RDS")

chang.load <- readRDS("./data/derivedData/change_dat.load.RDS")
change.vol <- readRDS("./data/derivedData/change_dat.vol.RDS")

confint.load <- readRDS("./data/derivedData/confint.load.RDS")
confint.vol <- readRDS("./data/derivedData/confint.vol.RDS")

m.load <- readRDS("./data/derivedData/m1.load.RDS")
m.vol <- readRDS("./data/derivedData/m1.vol.RDS")

## Join data

joined.dat <- confint.load %>%
        full_join(confint.vol) %>%
        print()


## Figure
pos <- position_dodge(width = 2)

joined.dat %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = pos,
                      width = 0.2) +
        geom_line(position = pos) +
        geom_point(shape = 21, position = pos, size = 3) +
        scale_x_discrete(labels=c("change.2" = "Session 2", "change.3" = "Session 3",
                                  "change.4" = "Session 4", "change.5" = "Session 5", 
                                  "change.post" = "Session 6")) +
        labs(x = "Time-point", y = "Absolute mean change in total session load \n(log-kg)\n", fill = "Supplement") +
        theme_classic() +
        theme(plot.background = element_rect(fill = "gray80"))


# Comments:
# How to I join these two data sets in one ggplot correctly?