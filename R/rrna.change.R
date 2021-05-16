## qPCR change analysis
#### qPCR-analysis



# This script creates and save the data frame qpcr_data2, calculates normalization factor (based on lambda), and
# plots mean expression per treatment

# Packages

library(tidyverse)
library(readxl)
library(nlme)

# Load data

qpcr_data  <- readRDS("./data/derivedData/qpcr-data.RDS")


samples <- read_excel("./data/rna/RNA.raw.xlsx") %>%
        dplyr::select(subject, time_rep, sample, weight) %>%
        print()

code_key <- read_excel("./data/code.key.xlsx") %>%
        dplyr::select(subject, supplement, time) %>%
        print()



# Calculate normalization factor based on lambda

nf <- qpcr_data %>%
        mutate(sample = as.numeric(sample)) %>% 
        inner_join(samples) %>%
        separate(time_rep, into = c("time","rep"), sep = "rna") %>%
        mutate(rep = paste0("cdna", rep)) %>%
        inner_join(code_key) %>% 
        mutate(time = if_else(time %in% c("T1", "T2"), "Pre", "Post"), 
               time = factor(time, levels = c("Pre", "Post")), 
               expr = eff ^ -cq) %>%
        
        filter(target %in% c("Lambda F2R2", "Lambda F3R3"), 
               cq < 35) %>%
        separate(target, into = c("target", "nf_primer")) %>%
        
        mutate(nf.w = expr * weight, 
               nf.w = nf.w / max(nf.w)) %>%
        
        
        dplyr::select(sample, subject, time, rep, supplement, nf.w, nf_primer) %>%
        print()


qpcr_data2 <- qpcr_data %>%
        mutate(sample = as.numeric(sample)) %>% 
        inner_join(samples) %>%
        inner_join(nf) %>%
        separate(time_rep, into = c("time","rep"), sep = "rna") %>%
        mutate(rep = paste0("cdna", rep)) %>%
        inner_join(code_key) %>% 
        mutate(time = if_else(time %in% c("T1", "T2"), "Pre", "Post"), 
               time = factor(time, levels = c("Pre", "Post")), 
               expr = eff ^ -cq) %>%
        
        filter(!(target %in% c("Lambda F2R2", "Lambda F3R3"))) %>%   
        inner_join(nf) %>%
        print()

saveRDS(qpcr_data2, "./data/derivedData/qpcr_data2.RDS")


# Plotting mean expression per treatment for rRNA

qpcr_data2 %>%
        mutate(norm.expr = log(expr / nf.w)) %>%
        group_by(time, supplement, target) %>%
        summarise(m = mean(norm.expr, na.rm = TRUE))  %>%
        #print()
        
        
        ggplot(aes(time, m, color = supplement, group = supplement)) + 
        geom_line() +
        geom_point() +
        facet_wrap( ~ target, scales = "free")

qdat <- qpcr_data2 %>%
        mutate(target = gsub("rRNA", "", target),
               target = gsub("  ", " ", target),
               target = paste0("trg_", target)) %>%
        separate(target, into = c("target", "primer"), sep = " ") %>%
        
        # Create a weight-normalized variable
        mutate(nf.expr = log(expr / nf.w), 
               nf.w = scale(nf.w),
               # Technical random effect
               technical = paste(subject, time, supplement, rep, sep = "_"),
               biological = paste0("S", sample)) %>%
        filter(!(target %in% c("trg_MHC1", "trg_MHC2A", "trg_MHC2X"))) %>%
        
        
        print()

# Calculate mean Ct and primer efficiency for methods chapter (primer sequences)

rrna.ct <- qdat %>%
        select(target, primer, cq) %>%
        filter(target %in% c("trg_18s", "trg_28s", "trg_5.8s", "trg_5s", "trg_47s")) %>%
        group_by(target, primer) %>%
        summarise(cq.m = mean(cq),
                  sd.m = sd(cq)) %>%
        print()

rrna.e <- qdat %>%
        select(target, primer, eff) %>%
        filter(target %in% c("trg_18s", "trg_28s", "trg_5.8s", "trg_5s", "trg_47s")) %>%
        group_by(target, primer) %>%
        summarise(E = mean(eff)) %>%
        print()

lambda.ct <- qpcr_data %>%
        select(target, cq) %>%
        filter(target %in% c("Lambda F2R2", "Lambda F3R3"), 
               cq < 35) %>% # Filter out values above 35, likely to be polluted/erroneous
        group_by(target) %>%
        summarise(cq.m = mean(cq),
                  cq.s = sd(cq)) %>%
        print()
        
lambda.e <- qpcr_data %>%
        select(target, cq, eff) %>%
        filter(target %in% c("Lambda F2R2", "Lambda F3R3"), 
               cq < 35) %>%
        select(target, eff) %>%
        group_by(target) %>%
        summarise(E = mean(eff)) %>%
        print()
        
## Change data - preliminary handling, filtering out each rRNA for separate analysis

rrna18 <- qdat %>%
        filter(target == "trg_18s") %>%
        print()

rrna28 <- qdat %>%
        filter(target == "trg_28s") %>%
        print()

rrna5.8 <- qdat %>%
        filter(target == "trg_5.8s") %>%
        print()

rrna5 <- qdat %>%
        filter(target == "trg_5s") %>%
        print()

rrna47 <- qdat %>%
        filter(target == "trg_47s") %>%
        print()

## Baseline analysis and change scores per rRNA
# A baseline analysis comparing rRNA expression at baseline between the two legs via a paired t.test, and providing a summary of mean rRNA expression
# and sd
## Change scores per rRNA
# The code beneath summarizes the mean values at each time, grouped by subject, time and supplement, creating a wider data set with observations of 
# participants glucose measurements per time point.
# Then, mutate() is used to calculate change scores, where each timepoint is log-transformed and compared to baseline. baseline = baseline - mean(baseline,
# na.rm = TRUE) mean centers the baseline values. Subject, supplement, baseline and change scores are then selected and pivoted for modeling.

# 18S
change.18 <- rrna18 %>%
        dplyr::select(subject, time, rep, nf.expr, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(nf.expr = mean(nf.expr, na.rm = TRUE)) %>%
        pivot_wider(names_from = time,
                    values_from = nf.expr) %>%
        ungroup() %>%
        # print()
        mutate(change = Post-Pre,
               pre = Pre - mean(Pre, na.rm = TRUE),
               supplement = factor(supplement, levels = c("PLACEBO", "GLUCOSE"))) %>%
        print()

base.18 <- change.18 %>%
        select(subject, supplement, Pre) %>%
        pivot_wider(names_from = supplement,
                    values_from = Pre) %>%
        print()

ttest.18 <- t.test(base.18$GLUCOSE, base.18$PLACEBO, paired = TRUE)
        
# 28S
change.28 <- rrna28 %>%
        dplyr::select(subject, time, rep, nf.expr, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(nf.expr = mean(nf.expr, na.rm = TRUE)) %>%
        pivot_wider(names_from = time,
                    values_from = nf.expr) %>%
        ungroup() %>%
        # print()
        mutate(change = Post-Pre,
               pre = Pre - mean(Pre, na.rm = TRUE),
               supplement = factor(supplement, levels = c("PLACEBO", "GLUCOSE"))) %>%
        print()

base.28 <- change.28 %>%
        select(subject, supplement, Pre) %>%
        pivot_wider(names_from = supplement,
                    values_from = Pre) %>%
        print()

ttest.28 <- t.test(base.28$GLUCOSE, base.28$PLACEBO, paired = TRUE)

# 5.8S
change.5.8 <- rrna5.8 %>%
        dplyr::select(subject, time, rep, nf.expr, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(nf.expr = mean(nf.expr, na.rm = TRUE)) %>%
        pivot_wider(names_from = time,
                    values_from = nf.expr) %>%
        ungroup() %>%
        # print()
        mutate(change = Post-Pre,
               pre = Pre - mean(Pre, na.rm = TRUE),
               supplement = factor(supplement, levels = c("PLACEBO", "GLUCOSE"))) %>%
        print()

base.5.8 <- change.5.8 %>%
        select(subject, supplement, Pre) %>%
        pivot_wider(names_from = supplement,
                    values_from = Pre) %>%
        print()

ttest.5.8 <- t.test(base.5.8$GLUCOSE, base.5.8$PLACEBO, paired = TRUE)

# 5S
change.5 <- rrna5 %>%
        dplyr::select(subject, time, rep, nf.expr, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(nf.expr = mean(nf.expr, na.rm = TRUE)) %>%
        pivot_wider(names_from = time,
                    values_from = nf.expr) %>%
        ungroup() %>%
        # print()
        mutate(change = Post-Pre,
               pre = Pre - mean(Pre, na.rm = TRUE),
               supplement = factor(supplement, levels = c("PLACEBO", "GLUCOSE"))) %>%
        print()

base.5 <- change.5 %>%
        select(subject, supplement, Pre) %>%
        pivot_wider(names_from = supplement,
                    values_from = Pre) %>%
        print()

ttest.5 <- t.test(base.5$GLUCOSE, base.5$PLACEBO, paired = TRUE)

# 47S

change.47 <- rrna47 %>%
        dplyr::select(subject, time, rep, nf.expr, supplement) %>%
        group_by(subject, time, supplement) %>%
        summarise(nf.expr = mean(nf.expr, na.rm = TRUE)) %>%
        pivot_wider(names_from = time,
                    values_from = nf.expr) %>%
        ungroup() %>%
        # print()
        mutate(change = Post-Pre,
               pre = Pre - mean(Pre, na.rm = TRUE),
               supplement = factor(supplement, levels = c("PLACEBO", "GLUCOSE"))) %>%
        print()

base.47 <- change.47 %>%
        select(subject, supplement, Pre) %>%
        pivot_wider(names_from = supplement,
                    values_from = Pre) %>%
        print()

ttest.47 <- t.test(base.47$GLUCOSE, base.47$PLACEBO, paired = TRUE)

# Linear mixed effects model
# This model sets an intercept per participant (mixed model), and controls for differences in pre/baseline values. It then tries to explain the changes in
# groups by the supplement provided, for each rRNA.

# 18S
m1 <- lmerTest::lmer(change ~ pre + supplement + (1|subject), 
                     data = change.18)

plot(m1)

summary(m1)

# 28S
m2 <- lmerTest::lmer(change ~ pre + supplement + (1|subject), 
                     data = change.28)

plot(m2)

summary(m2)

# 5.8S
m3 <- lmerTest::lmer(change ~ pre + supplement + (1|subject), 
                     data = change.5.8)

plot(m3)

summary(m3)

# 5S
m4 <- lmerTest::lmer(change ~ pre + supplement + (1|subject), 
                     data = change.5)

plot(m4)

summary(m4)

# 47S
m5 <- lmerTest::lmer(change ~ pre + supplement + (1|subject), 
                     data = change.47)

plot(m5)

summary(m5)


## Fold-change estimated means
# Gets estimated means from the model, these are average increase at pre = 0 (the average pre value).
# These are log-fold change values (changeble with the mutate function)

# 18S
confint(emmeans(m1, specs = ~"supplement")) %>%
        data.frame()

# 28S
confint(emmeans(m2, specs = ~"supplement")) %>%
        data.frame()

# 5.8S
confint(emmeans(m3, specs = ~"supplement")) %>%
        data.frame()

# 5S
confint(emmeans(m4, specs = ~"supplement")) %>%
        data.frame()

# 47S
confint(emmeans(m5, specs = ~"supplement")) %>%
        data.frame()


