## Analysis of baseline DXA results.

# This script analyses the baseline DXA values and creates a table

# Packages

library(readxl)
library(tidyverse)
library(nlme)
library(lme4)
library(broom)
library(knitr)
library(emmeans)
library(knitr)

#1 descritpive table of subject
#
#count gender to be put in table
gendercount <- read_excel("data/ribose_dxa.xlsx") %>%
        count(sex) %>%
        pivot_wider(names_from = sex,
                    values_from = n) %>%
        print()

#Calculates grams to kg
table.dxa <- read_excel("data/ribose_dxa.xlsx") %>%
        select(sex, subject, age, height, weight, leanmasskg, fatmasskg, ffmg) %>%
        #make grams into kg
        mutate(ffkg = ffmg/1000) %>%
        select(sex, age, height, weight, leanmasskg, fatmasskg, ffkg) %>%
        print()


#dxa table pr. gender               
dxa.table <- table.dxa %>%
        pivot_longer(names_to = "variable",
                     values_to = "values", col= age:ffkg) %>%
        select(sex, variable, values) %>%
        group_by(sex, variable) %>%
        summarise(m = mean(values, na.rm = TRUE),
                  s= sd(values, na.rm = TRUE)) %>%
        mutate(stat = paste0(round(m, 1), " (", round(s, 1), ")")) %>%
        select(sex, variable, stat) %>%
        pivot_wider(names_from = sex,
                    values_from = stat) %>%
        mutate(variable = if_else(variable %in% c("age"),
                                  "Age",
                                  if_else(variable %in% c("fatmasskg"),
                                          "Fat mass",
                                          if_else(variable %in% c("ffkg"),
                                                  "Fat free mass",
                                                  if_else(variable %in% c("height"),
                                                          "Height",
                                                          if_else(variable %in% c("leanmasskg"),
                                                                  "Lean mass",
                                                                  if_else(variable %in% c("weight"),
                                                                          "Weight",
                                                                          ""))))))) %>%
        print()


dxa.table %>%
        kable()

# dxa table on total group.
dxa.table2 <- table.dxa %>%
        pivot_longer(names_to = "variable",
                     values_to = "values", col= age:ffkg) %>%
        select(variable, values) %>%
        group_by(variable) %>%
        summarise(m = mean(values, na.rm = TRUE),
                  s= sd(values, na.rm = TRUE)) %>%
        mutate(stat = paste0(round(m, 1), " (", round(s, 1), ")")) %>%
        select(variable, stat) %>%
        mutate(variable = if_else(variable %in% c("age"),
                                  "Age",
                                  if_else(variable %in% c("weight"),
                                          "Weight",
                                          if_else(variable %in% c("ffkg"),
                                                  "Fat free mass",
                                                  if_else(variable %in% c("height"),
                                                          "Height",
                                                          if_else(variable %in% c("leanmasskg"),
                                                                  "Lean mass",
                                                                  if_else(variable %in% c("fatmasskg"),
                                                                          "Fatt mass",
                                                                          ""))))))) %>%
        print()

dxa.table2 %>%
        kable()

saveRDS(dxa.table, file = "./data/derivedData/dxatable.RDS")


#2. take out kg of subject and merge it with nutrition to protein pr kg pr group


saveRDS(fp.weight, file = "./data/derivedData/fp_weight.RDS")

#DXA pr. leg.

#Load dxa data


dxaleg <-  read_excel("data/ribose_dxaprleg.xlsx") %>%
        select(subject, sex, leg, leanmassg, fatmassg, totalmasskg) %>%
        #make grams into kg
        mutate(leanmass = leanmassg/1000,
               fatmass = fatmassg/1000)%>%
        select(subject, sex, leg, leanmass, fatmass, totalmasskg) %>%
        print()
#Make analysis bewteen groups. 
dxa.stat <- dxaleg %>%
        select(subject, leg, leanmass, fatmass, totalmasskg) %>%
        print()

#lean mass between legs, within subject.
#no difference bewtween legs

leanmass.lm <- lmer(leanmass ~ leg + leg:subject + (1|subject), data = dxa.stat)

plot(leanmass.lm)

summary(leanmass.lm)

emmeans(leanmass.lm, specs = ~ "leg") %>%
        confint()

# fat mass between legs, within subject
#no difference bewtween legs

fat.lm <- lmer(fatmass ~ leg + leg:subject +(1|subject), data = dxa.stat)

plot(fat.lm)

summary(fat.lm)

emmeans(fat.lm, specs = ~ "leg") %>%
        confint()

#Totalmass between legs, within subject.
#no difference bewtween legs

total.lm <- lmer(totalmasskg ~leg + leg:subject + (1|subject), data = dxa.stat)

plot(total.lm)

summary(total.lm)

supsup <- emmeans(total.lm, specs = ~ "leg") %>%
        confint() %>%
        print()


#Table showing leanmass, fatmass and total mass in kg pr. leg
dxalegtable <- dxaleg %>%
        select(leg, leanmass, fatmass, totalmasskg) %>%
        pivot_longer(names_to = "variable",
                     values_to = "values", col= leanmass:totalmasskg) %>%
        group_by(leg, variable) %>%
        summarise(m = mean(values, na.rm = TRUE),
                  s= sd(values, na.rm = TRUE)) %>%
        mutate(stat = paste0(round(m, 1), " (", round(s, 1), ")")) %>%
        select(leg, variable, stat) %>%
        pivot_wider(names_from = leg,
                    values_from = stat) %>%
        print()

#table in markdown
dxalegtable %>%
        kable()


saveRDS(dxalegtable, file = "./data/derivedData/dxalegtable.RDS")
