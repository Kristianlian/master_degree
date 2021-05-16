## qPCR change analysis
#### qPCR-analysis



# This script creates and save the data frame qpcr_data2, calculates normalization factor (based on lambda), and
# plots mean expression per treatment


# The data frame contains:
# 1. Sample name (sample)
# 2. Target genes (target)
# 3. Cq-values (Cq)
# 4. Amplification efficiencies (eff)
# 5. Subject
# 6. Time point (time)
# 7. Replicate number (rep)
# 8. Biopsy weight (weight)
# 9. Supplement
# 10. Expression (expr)
# 11. Normalization factor weigh (nf.w)

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

## Change data
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

# Change scores

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


# Create model: 
# Needs to have an intercept per participant (mixed model)
# Control for pre values.

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


### Get estimated means from the model, these are average increase at 
# pre = 0 (the average pre value)
# remember that these are originally log-fold change values, transformed back using exp(emmean) etc

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

# Tabel

tab_model(m1, m2, m3, m4, m5,
          pred.labels = c("Glucose", "Placebo", "Placebo vs. Glucose"),
          string.pred = "Coefficient",
          string.ci = " Conf.Int (95%)",
          string.p = "P-value",
          dv.labels = c("18S", "28S", "5.8S", "5S", "47S"),
          show.re.var = FALSE,
          show.icc = FALSE)
