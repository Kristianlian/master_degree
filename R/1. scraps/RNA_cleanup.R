##-------------------------------------
## RNA-cleanup.R
##
## Title: RNA cleanup
## Purpose: Filter away bad RNA estimates based on RNA/weight relationship
## Author: DH/KL
##
##
##
##-------------------------------------
## Notes:
# Filter technical outliers by model log(rna) ~ muscle weight.
# Set threshold for filtering by adjusting threshold for standardized
# residuals. A sensible threshold is 3 sd's from the modelled relationship
# between RNA and muscle weight, taking training status into account.
#
#
#
#
#
## ------------------------------------

# Packages
library(tidyverse);library(readxl);library(nlme)




# A data frame with technical outliers

outlier_tech <- data.frame(subject =
                                   c(101, 102, 102, 109, 109, 112, 112, 113, 113),
                           time = c("T3", "T3", "T4", "T1", "T2", "T1", "T2", "T2", "T2"),
                           rep = c("2", "1", "1", "1", "2", "1", "2", "1", "2"),
                           technical = "out")

# Load data

dat <- read_excel("./data/RNA_raw.xlsx") %>%
        pivot_longer(names_to = "variable", 
                     values_to = "value", 
                     cols = RATIO_R260280_1:RNA_conc_4) %>%
        separate(variable, into = c("Type", "subtype", "measurement")) %>%
        mutate(Type = if_else(Type == "RNA", Type, subtype)) %>%
        select(-subtype) %>%
        pivot_wider(names_from = Type, 
                    values_from = value) %>%
        mutate(RNA_tot = RNA * 2 * eluate) %>%
        separate(time_rep, into = c("time", "rep"), sep = "rna") %>%
        mutate(trained = if_else(time %in% c("T3", "T4"), "trained", "untrained"), 
               biopsy = paste0(subject, time)) %>%
        left_join(outlier_tech) %>%
        mutate(technical = if_else(is.na(technical), "in", technical)) %>%
        print()










