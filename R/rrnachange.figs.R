## qPCR change analysis
#### qPCR-analysis



# This script is equal to rrna.change, except it produces figures illustrating reverse transformed fold change from pre-post per rRNA per supplement. 
# For discription of codes, see rrna.change

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


# Linear mixed effects model

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


### Get estimated means from the model

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
        data.frame() %>%
        print()

## Figures
# Plots fold change of estimated marginal means per supplement and saves the figures to be used in the thesis

# 18S
plot.18 <- confint(emmeans(m1, specs = ~"supplement")) %>%
        data.frame() %>%
        mutate(time = "post") %>%
        add_row(supplement = c("PLACEBO", "GLUCOSE"), time = "pre", emmean = 0, lower.CL = 0, upper.CL = 0) %>%
        mutate(time = factor(time, levels = c("pre", "post"), labels = c("Baseline", "Post"))) %>%
        ggplot(aes(time, exp(emmean), group = supplement, fill = supplement)) +
        annotate("text", x = "Post", y = 2.1, label = "p = 0.585", size = 3) +
        geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)),
                      width = 0.1,
                      position = position_dodge(width = 0.2)) +
        geom_line(position = position_dodge(width = 0.2)) +
        geom_point(shape = 21, size = 3, position = position_dodge(width = 0.2)) +
        labs(x = "", y = "18S rRNA \n(fold change)\n", fill = "") +
        theme_classic()

saveRDS(plot.18, "./data/derivedData/plot.18.RDS")

# 28S
plot.28 <- confint(emmeans(m2, specs = ~"supplement")) %>%
        data.frame() %>%
        mutate(time = "post") %>%
        add_row(supplement = c("PLACEBO", "GLUCOSE"), time = "pre", emmean = 0, lower.CL = 0, upper.CL = 0) %>%
        mutate(time = factor(time, levels = c("pre", "post"), labels = c("Baseline", "Post"))) %>%
        ggplot(aes(time, exp(emmean), group = supplement, fill = supplement)) +
        annotate("text", x = "Post", y = 2.1, label = "p = 0.740", size = 3) +
        geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)),
                      width = 0.1,
                      position = position_dodge(width = 0.2)) +
        geom_line(position = position_dodge(width = 0.2)) +
        geom_point(shape = 21, size = 3, position = position_dodge(width = 0.2)) +
        labs(x = "", y = "28S rRNA \n(fold change)\n", fill = "") +
        theme_classic()

saveRDS(plot.28, "./data/derivedData/plot.28.RDS")

# 5.8S
plot.5.8 <- confint(emmeans(m3, specs = ~"supplement")) %>%
        data.frame() %>%
        mutate(time = "post") %>%
        add_row(supplement = c("PLACEBO", "GLUCOSE"), time = "pre", emmean = 0, lower.CL = 0, upper.CL = 0) %>%
        mutate(time = factor(time, levels = c("pre", "post"), labels = c("Baseline", "Post"))) %>%
        ggplot(aes(time, exp(emmean), group = supplement, fill = supplement)) +
        annotate("text", x = "Post", y = 2.1, label = "p = 0.935", size = 3) +
        geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)),
                      width = 0.1,
                      position = position_dodge(width = 0.2)) +
        geom_line(position = position_dodge(width = 0.2)) +
        geom_point(shape = 21, size = 3, position = position_dodge(width = 0.2)) +
        labs(x = "", y = "5.8S rRNA \n(fold change)\n", fill = "") +
        theme_classic()

saveRDS(plot.5.8, "./data/derivedData/plot.5.8.RDS")

# 5S
plot.5 <- confint(emmeans(m4, specs = ~"supplement")) %>%
        data.frame() %>%
        mutate(time = "post") %>%
        add_row(supplement = c("PLACEBO", "GLUCOSE"), time = "pre", emmean = 0, lower.CL = 0, upper.CL = 0) %>%
        mutate(time = factor(time, levels = c("pre", "post"), labels = c("Baseline", "Post"))) %>%
        ggplot(aes(time, exp(emmean), group = supplement, fill = supplement)) +
        annotate("text", x = "Post", y = 1.9, label = "p = 0.790", size = 3) +
        geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)),
                      width = 0.1,
                      position = position_dodge(width = 0.2)) +
        geom_line(position = position_dodge(width = 0.2)) +
        geom_point(shape = 21, size = 3, position = position_dodge(width = 0.2)) +
        labs(x = "", y = "5S rRNA \n(fold change)\n", fill = "") +
        theme_classic()

saveRDS(plot.5, "./data/derivedData/plot.5.RDS")

# 47S
plot.47 <- confint(emmeans(m5, specs = ~"supplement")) %>%
        data.frame() %>%
        mutate(time = "post") %>%
        add_row(supplement = c("PLACEBO", "GLUCOSE"), time = "pre", emmean = 0, lower.CL = 0, upper.CL = 0) %>%
        mutate(time = factor(time, levels = c("pre", "post"), labels = c("Baseline", "Post"))) %>%
        ggplot(aes(time, exp(emmean), group = supplement, fill = supplement)) +
        annotate("text", x = "Post", y = 3.5, label = "p = 0.502", size = 3) +
        geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)),
                      width = 0.1,
                      position = position_dodge(width = 0.2)) +
        geom_line(position = position_dodge(width = 0.2)) +
        geom_point(shape = 21, size = 3, position = position_dodge(width = 0.2)) +
        labs(x = "Time-point", y = "47S pre-rRNA \n(fold change)\n", fill = "") +
        theme_classic()

saveRDS(plot.47, "./data/derivedData/plot.47.RDS")


