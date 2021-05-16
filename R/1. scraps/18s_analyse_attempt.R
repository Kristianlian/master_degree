#### 18s rRNA analysis

# Packages
library(tidyverse);library(knitr);library(nlme);library(emmeans);library(lme4)

# Outlier tech
outlier_tech <- data.frame(subject =
                                   c(101, 102, 102, 109, 109, 112, 112, 113, 113),
                           time = c("T3", "T3", "T4", "T1", "T2", "T1", "T2", "T2", "T2"),
                           rep = c("2", "1", "1", "1", "2", "1", "2", "1", "2"),
                           technical = "out")

# Data
rrna18.dat <- readRDS("./data/derivedData/qpcr_data2.RDS") %>%
        mutate(norm.expr = log(expr / nf.w)) %>%
        filter(target == "18s rRNA F2R2") %>%
        mutate(time = factor(time, levels = c("pre", "post"))) #%>%
        #left_join(outlier_tech) %>%
       # mutate(technical = if_else(is.na(technical), "in", technical)) %>%
       # print()


# Model

lm1 <- lm(norm.expr ~ time * supplement, data = rrna18.dat)
summary(lm1)
plot(lm1)

## Notes on model:
# Intercept: Mean in the glucose grp at baseline
# Timepost: Mean in the glucose grp post intervention
# SupplementPLACEBO: The difference at baseline between groups 
# Timepost:supplementPLACEBO: Difference in placebo compared to glucose at time point post
## Assumption checks
# Residuals vs. Fitted: 97, 9, 193
# Normal Q-Q: 97, 9, 193
# Scale-location: 97, 9, 193
# 

lmer1 <- lmer(norm.expr ~ time * supplement + (1|subject), data = rrna18.dat)
summary(lmer1)
confint(lmer1)
plot(lmer1)
