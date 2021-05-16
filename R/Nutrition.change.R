### Analysis of nutrition 

# This script analyses the macro-nutrient intake of the subjects per day per supplement

# Packages

library(knitr)
library(rmarkdown)
library(tidyverse)
library(readxl)
library(readr)
library(tinytex)
library(tidyr)
library(knitr)
library(broom)
library(arsenal) 
library(lme4)
library(dplyr)
library(lmerTest)
library(emmeans)
library(dbplyr)
library(magrittr)
library(dabestr)

# Subject 106 and 116 did not deliver their nutrient data, as such they are excluded from the analysis

# Summarise macro nutrients
# Supplemental protein (50g whey protein isolate = 42g protein) and glucose (90g) was added by investigators 
# prior to analysis

fp.weight <-read_excel("data/ribose_dxa.xlsx") %>%
        select(weight, subject) %>%
        print()

nutha <- read_excel("data/Ribose_nutrition_result.xlsx")%>%
        select(timepoint, subject, group, calories, fat, carbohydrates, protein, sup_pro, sup_gluc, kcal_glu) %>%
        mutate(timepoint = if_else(timepoint %in% c("T1","T2"), 
                                   "Day 1",
                                   if_else(timepoint %in% c("3", "4"),
                                           "Day 2",
                                           if_else(timepoint %in% c("5", "6"),
                                                   "Day 3",
                                                   if_else(timepoint %in% c("7", "8"),
                                                           "Day 4",
                                                           if_else(timepoint %in% c("9", "10"),
                                                                   "Day 5",
                                                                   if_else(timepoint %in% c("T3", "T4"),
                                                                           "Day 6",
                                                                           "na"))))))) %>%
        group_by(timepoint, subject, group) %>%
        inner_join(fp.weight) %>%
        summarise(protein = sum(protein + sup_pro),
                  fat = sum(fat),
                  calories = sum(calories + kcal_glu),
                  carbohydrates = sum(carbohydrates + sup_gluc),
                  weight = sum(weight)) %>%
        mutate(proprkg = (protein/weight)) %>%
        print()




#change data
#Protein

pro_change <- nutha %>%
        select(timepoint, subject, group, protein,) %>%
        print()
#ungroup() %>%
#pivot_wider(names_from = group,
#           values_from = c(protein)) %>%
#mutate(change = glucose - placebo) %>%
#print()


prolm <- lmer(protein ~ timepoint + timepoint:group + (1|subject), data = pro_change)

plot(prolm)

summary(prolm)

#Fat
fat_change <- nutha %>%
        select(timepoint, subject, group, fat,) %>%
        print()
#ungroup() %>%
#pivot_wider(names_from = group,
#           values_from = c(fat)) %>%
#mutate(change = glucose - placebo) %>%
#print()

fatlm <- lmer(fat  ~ timepoint + timepoint:group +(1|subject), data = fat_change)

plot(fatlm)

summary(fatlm)

#Carbohydrates
carb_change <- nutha %>%
        select(timepoint, subject, group, carbohydrates,) %>%
        print()
#ungroup() %>%
#pivot_wider(names_from = group,
#           values_from = c(carbohydrates)) %>%
#mutate(change = glucose - placebo) %>%
#print()

carblm <- lmer(carbohydrates  ~ timepoint + timepoint:group + (1|subject), data = carb_change)

plot(carblm)

summary(carblm)

#calories
cal_change <- nutha %>%
        select(timepoint, subject, group, calories,) %>%
        print()
#ungroup() %>%
#pivot_wider(names_from = group,
#           values_from = c(calories)) %>%
#mutate(change = glucose - placebo) %>%
#print()


callm <- lmer(calories  ~ 0 + timepoint + timepoint:group + (1|subject), data = cal_change)

plot(callm)

summary(callm)

#Protein pr bodyweigh



prowe <- nutha %>%
        select(timepoint, subject, group, proprkg) %>%
        print()

prokglm <- lmer(proprkg  ~ 0 + timepoint + timepoint:group + (1|subject), data = prowe)

plot(prokglm)

summary(prokglm)

