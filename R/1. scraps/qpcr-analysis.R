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
        
        
        ggplot(aes(time, m, color = supplement, group = supplement)) + 
        geom_line() +
        geom_point() +
                       facet_wrap( ~ target, scales = "free")


### Fitting models #################

qpcr_data2 <- readRDS("./data/derivedData/qpcr_data2.RDS")




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







# Fixed and random effects formulas for the first step of model building
# Compared to Method paper (Hammarström et al 2018), the fixed effects are reduced
# to only contain gene-specific time + time:sets
fixed <- nf.expr ~ 0  + target + target:time + target:supplement + target:time:supplement
random <- list(subject = ~1, technical = ~1, nf_primer = ~ 1)





## m1 is the model assuming homoscedastic errors.
m1 <- lme(fixed, random = random, data = qdat,
          control=list(msMaxIter=120,
                       opt = "nloptwrap",msVerbose=TRUE), 
          method = "REML", 
          na.action = na.exclude) # This allows to check progress


#### Models allowing for heteroscedasticity in the residuals. ####
# The variance is believed to be different from gene to gene due to different expression, primer design etc.
# One can also expect that residual variance to be related to cq-values as higher cq-values 
# implies higher meassurement error due to stochastic variation in the beginning of amplification.


# Varioance functions with cq values
varfun1 <- varIdent(form = ~1 | target) # Heterogeneity per gene
varfun2 <- varExp(form = ~ cq | target) # Variance changes as a exponential of cq-values



m2 <- update(m1, weights = varfun1)
m3 <- update(m1, weights = varfun2)



anova(m1, m2, m3)

# Preliminary results, including varExp per target is better than no variace function and 
# different variance per target.


## Check if var-cov matrix is ok...

intervals(m3)
# All random effects are estimated, however, 
# the random effect for technical has large CI


anova(m1, m5)


### Model diagnostics ##### 

plot(m3, resid(., type = "n")~fitted(.)|target)
qqnorm(m3, ~resid(., type = "n")|target, abline = c(0,1))

# Looks good



#### Get estimates from each gene #####################


summary(m1)




library(emmeans)


estimated_means <- emmeans(m3, specs = ~ time * supplement | target)



## Figures for presentation

# Compute foldchange from log-ratio values
logratio2foldchange <- function(logratio, base=2)
{
  retval <- base^(logratio)
  retval <- ifelse(retval < 1, -1/retval, retval)
  retval
}

rrna.pre <- estimated_means %>%
  data.frame() %>%
  filter(time == "Pre") %>%
  pivot_wider(names_from = "time",
              values_from = "emmean") %>%
  select(Pre)

rrna.post <- estimated_means %>%
  data.frame() %>%
  filter(time == "Post") %>%
  pivot_wider(names_from = "time",
              values_from = "emmean") %>%
  select(Post)

prepost <- c(rrna.pre, rrna.post) %>%
  print()

matrix.pp <- matrix(unlist(prepost), ncol = 2, byrow = TRUE)

matrix.pp


  logratio2foldchange(matrix.pp, base = 2)

        



 estimated_means %>%
   data.frame() %>%
  filter(target %in% c("trg_18s", "trg_28s", "trg_5.8s", "trg_5s")) %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.05,
                      position = position_dodge(width = 0.2)) +
        geom_line(position = position_dodge(width = 0.2)) +
        geom_point(shape = 21, position = position_dodge(width = 0.2),
                   size = 3) +
        labs(x = "Time-point", y = "Ribosomal RNA \n(log-abundance per tissue weight)\n", fill = "Supplement") +
        theme_classic() +
        theme(plot.background = element_rect(fill = "gray80")) +
        facet_wrap(~ target, scales = "free")

estimated_means %>%
        data.frame() %>%
        filter(target %in% c("trg_47s")) %>%
        ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      position = position_dodge(width = 0.2), 
                      width = 0.05) +
        geom_line(position = position_dodge(width = 0.2)) +
        geom_point(shape = 21, position = position_dodge(width = 0.2),
                   size = 3) +
        labs(x = "Time-point", y = "Ribosomal RNA \n(log-abundance per tissue weight)\n", fill = "Supplement") +
        theme_classic() +
        theme(plot.background = element_rect(fill = "gray80"))

summary(m3)


# Contemporary presentation, keep this or change to fold-change?
# Make a table containing rRNA outcomes and p-values/CI


