#### rRNA analysis

## This script creates and save the data frame rrna_data and plots mean rRNA expression
# of each rRNA target gene per supplement

# Packages
library(tidyverse)
library(knitr)

# Load data
rrna_data <- readRDS("./data/derivedData/qpcr_data2.RDS") %>%
        filter(target %in% c("47s rRNA F1R1", "18s rRNA F2R2", "28s rRNA F2R2", "5.8s rRNA F2R2", "5s rRNA F3R3"))


# Exploratory statistics 
rrna_summary <- rrna_data %>%
        mutate(norm.expr = log(expr / nf.w)) %>%
        group_by(target, time, supplement) %>%
        summarise(Mean = mean(norm.expr, na.rm = TRUE),
                  SD = sd(norm.expr, na.rm = TRUE),
                  n = n())
        
rrna_summary %>%
        kable(digits = c(NA, NA, NA, 4, 4, 0),
              col.names = c("Target", "Time-point", "Supplement", "Mean", "SD", "n"))

# Plotting mean expression per rRNA gene per supplement

rrna_data %>%
        mutate(norm.expr = log(expr / nf.w)) %>%
        group_by(time, supplement, target) %>%
        summarise(m = mean(norm.expr, na.rm = TRUE))  %>%
        
        
        ggplot(aes(time, m, color = supplement, group = supplement)) + 
        geom_line() +
        geom_point() +
        facet_wrap( ~ target, scales = "free")

# Fitting model

rrna.model <- rrna_data %>%
        mutate(norm.expr = log(expr / nf.w)) %>%
        group_by(time, supplement, target)

lm1 <- lm(norm.expr ~ supplement * time, data = rrna.model)
summary(lm1)





