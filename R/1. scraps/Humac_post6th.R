#### Humac analysis pre vs. post sixth RT session


## Author: Kristian Lian
##: Project: Ribose

# Purpose: This script plots mean torque per supplement (both through intervention and pre vs. post) results from the ribose project, 
# and analyses the data per test (isometric, isokinetic 60, isokinetic 240) in a linear model.

## Time-points
# D-1: Baseline, before any supplementation or training
# D4, D5, D8 and D9: Day 4, 5, 8 and 9 of the intervention, humac testing of the leg that performed
# RT the preceding day
# T3: Post testing leg #1 (leg that started the intervention). Leg #1 is tested four times at T3/T4:
# Test 1 leg 1: 1.5hrs after protein ingestion, 45min before RT (T3)
# Test 2 leg 1: 30min after RT (T3)
# Test 3 leg 1: 2hrs after RT (T3)
# Test 4 leg 1: ~23hrs after RT (T4)
# Test 1 serve as a post test for the 5 RT sessions and pre test before the sixth session, test 2,
# 3, and 4 serve as post test following sixth session
# T4 and 13 follow the same design for leg #2

## Data
# Date of testing
# Subject
# Test type: isok.60 (isokinetic 60), isok.240 (isokinetic 240), isom (isometric)
# Peak.torque: Highest peak torque from each test
# Leg: left or right leg
# Supplement: glucose or placebo

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




## Figure

pos <- position_dodge(width = 0.2)


humac %>%
        filter(acute %in% c("rest", "post23h") )%>%
      
        ggplot(aes(time, peak.torque, color = supplement)) +
        geom_boxplot() +
        labs(x = "Time-point", y = "Peak torque (nm)", color = "Supplement",
             title = "Mean peak torque developement per supplement",
             subtitle = "Different patterns per supplement",
             caption = "Values are mean") +
        facet_wrap(~test, scales = "free")
        

### Models

# select each test

  
# Isometric test

m1 <- lmer(peak.torque ~ acute * supplement + (1|subject), 
           data = filter(humac, time %in% c("test3","test4"),
                         test == "isom"))



plot(m1)

summary(m1)

confint(m1)

# Get mean estimates from model isometric

emmeans(m1, specs = ~ "acute|supplement") %>%
        data.frame() %>%
        ggplot(aes(acute, emmean, color = supplement, group = supplement)) +
        geom_errorbar(aes(ymin= lower.CL, ymax = upper.CL), width = 0.2) +
  geom_point() +      
  geom_line() 



pairs(emmeans(m1, specs = ~ "acute|supplement"), reverse = TRUE) %>%
        confint()



cbind(coef(summary(m1)), data.frame(confint(m1))[3:10, ])

  

# Isokinetic 60 test

m2 <- lmer(peak.torque ~ acute * supplement + (1|subject), 
           data = filter(humac, time %in% c("test3","test4"),
                         test == "isok.60"))



plot(m2)

summary(m2)

confint(m2)

# Get mean estimates from model isokinetic 60



emmeans(m2, specs = ~ "acute|supplement") %>%
  data.frame() %>%
  ggplot(aes(acute, emmean, color = supplement, group = supplement)) +
  
  geom_errorbar(aes(ymin= lower.CL, ymax = upper.CL), width = 0.2) +
  geom_line() +
  geom_point() 


pairs(emmeans(m2, specs = ~ "acute|supplement"), reverse = TRUE) %>%
  confint()



cbind(coef(summary(m2)), data.frame(confint(m2))[3:10, ])


# Isokinetic 240 test

m3 <- lmer(peak.torque ~ acute * supplement + (1|subject), 
           data = filter(humac, time %in% c("test3","test4"),
                         test == "isok.240"))



plot(m3)

summary(m3)

confint(m3)

# Get mean estimates from model isokinetic 240



emmeans(m3, specs = ~ "acute|supplement") %>%
  data.frame() %>%
  ggplot(aes(acute, emmean, color = supplement, group = supplement)) +
  
  geom_errorbar(aes(ymin= lower.CL, ymax = upper.CL), width = 0.2) +
  geom_line() +
  geom_point() 


pairs(emmeans(m3, specs = ~ "acute|supplement"), reverse = TRUE) %>%
  confint()



cbind(coef(summary(m3)), data.frame(confint(m3))[3:10, ])

# Fortsettelse: Legg til layer med iso testene per supplement i emmeans figuren
