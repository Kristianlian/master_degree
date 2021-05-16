## Total RNA analysis

# This script imports cleaned raw RNA data and analyzes baseline values via t-testing, calculates log-fold change score from pre to post, analyzes 
# change scores via a linear mixed effects model and emmeans, and produces the figure illustrating fold change in total RNA from pre to post

## Load packages 
library(tidyverse); library(readxl);library(nlme);library(lme4);library(emmeans); library(knitr);library(broom);
library(kableExtra);library(sjPlot);library(sjmisc)


## Load data 

rna_dat <- readRDS("./data/derivedData/rna_complete.RDS")

rna.ratio <- read_excel("./data/rna/RNA.raw.xlsx", na = "NA")
  

# 260/280 Ratio for estimation of sample purity
mean.rat <- rna.ratio %>%
  select(sample, RATIO_R260280_1, RATIO_R260280_2, RATIO_R260280_3, RATIO_R260280_4) %>%
  pivot_longer(names_to = "rep",
               values_to = "ratio",
               cols = RATIO_R260280_1:RATIO_R260280_4) %>%
  group_by(sample) %>%
  #print()
  mutate(mean.ratio = mean(ratio),
         sd.ratio = sd(ratio)) %>%
  print()

# Loads code_key containing ID list of which leg received which supplement

code <- read_excel("./data/code.key.xlsx")

# Sorts baseline data (T1 and T2) from the two legs into pre, and post biopsy data (T3 and T4) into post, calculates total RNA per wet muscle weight

rna_dat2 <- rna_dat %>%
        inner_join(code) %>%
        mutate(time = if_else(time %in% c("T1", "T2"), "pre", "post"), 
               time = factor(time, levels = c("pre", "post"))) %>%
        # filter(outlier == "in") %>%
        group_by(subject, time, rep, outlier, technical, biopsy, supplement) %>%
        summarise(weight = mean(weight), 
                  RNA = mean(RNA_tot)) %>%
        ungroup() %>%
        mutate(weight.mc = weight/mean(weight)) %>%
        mutate(RNA.weight = RNA / weight) %>%
        print()

# Illustration of individual subjects change in total RNA

rna_dat2 %>%
          filter(technical == "in") %>%
        group_by(subject, time,  supplement) %>%
        summarise(weight = mean(weight), 
                  RNA = mean(RNA)) %>%
        mutate(RNA = RNA/weight) %>%
        ggplot(aes(time, RNA, group=paste(subject, supplement), 
                   color = supplement)) + 
        scale_y_continuous(limits = c(400, 1600)) +
        geom_line() +
        geom_point()

## Baseline analysis: calculation of mean and SD of baseline values, t-test on baseline values

base.rna <- rna_dat2 %>%
  filter(time == "pre") %>%
  select(subject, time, supplement, RNA) %>%
  group_by(subject, supplement) %>%
  summarise(m = mean(RNA)) %>%
  pivot_wider(names_from = supplement,
              values_from = m) %>%
  print()

rna.ttest <- t.test(base.rna$GLUCOSE, base.rna$PLACEBO, paired = TRUE)

rna.sum <- rna_dat2 %>%
  filter(time == "pre") %>%
  select(subject, time, supplement, RNA) %>%
  group_by(supplement) %>%
  summarise(m = mean(RNA),
            s = sd(RNA)) %>%
  print()



# Log-fold change score calculation: We are interested in if change in GLUCOSE is different from change in 
# PLACEBO. 

change_dat <- rna_dat2 %>%
   #  filter(technical == "in", 
   #        outlier == "in") %>%
    dplyr::select(subject:rep, supplement, RNA.weight) %>%
    
    group_by(subject, time, supplement) %>%
    summarise(RNA.weight = mean(RNA.weight, na.rm = TRUE)) %>%
    
    pivot_wider(names_from = time, 
                values_from = RNA.weight) %>%
    ungroup() %>%
    mutate(change = log(post) - log(pre), 
           pre = pre - mean(pre, na.rm = TRUE), 
           supplement = factor(supplement, levels = c("PLACEBO", "GLUCOSE")))

## Linear mixed effects model 
# The data set is fitted in an ANCOVA model, explaining change based on the condition (supplement), accounting for differential baseline values and
# multiple sampling of the same subjects. The effect of condition will answer the question.

m1 <- lmerTest::lmer(change ~ pre + supplement + (1|subject), 
           data = change_dat)
    
plot(m1)

summary(m1)

## Fold-change estimated means
# Gets estimated means from the model, these are average increase at pre = 0 (the average pre value).
# These are log-fold change values (changeble with the mutate function), reverse transformed with exp() for illustration

totrna.fig <- confint(emmeans(m1, specs = ~"supplement")) %>%
  data.frame() %>%
  mutate(time = "post") %>%
  add_row(supplement = c("PLACEBO", "GLUCOSE"), time = "pre", emmean = 0, SE = 0, lower.CL = 0, upper.CL = 0) %>%
  mutate(time = factor(time, levels = c("pre", "post"), labels = c("Baseline", "Post"))) %>%
  # print()
  ggplot(aes(time, exp(emmean), group = supplement, fill = supplement)) +
  annotate("text", x = "Post", y = 1.5, label = "p = 0.499", size = 3) +
  geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)),
                width = 0.1,
                position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2)) +
  geom_point(position = position_dodge(width = 0.2), shape = 21, size = 3) +
  labs(x = "Time-point", y = "Total RNA per mg muscle tissue \n(fold change)\n", fill = "") +
  theme_classic()  
  #theme(legend.position = "none")

saveRDS(totrna.fig, "./data/derivedData/totrna.fig.RDS")




                 







