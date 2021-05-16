##-------------------------------------
## analyze-rna.R
##
## Title:
## Purpose: 
## Author:
##
##
##
##-------------------------------------
## Notes:
#
#
#
#
#
#
#
#
## ------------------------------------

## Load packages 
library(tidyverse); library(readxl); library(lme4);library(emmeans); library(knitr);library(broom);
library(kableExtra);library(sjPlot);library(sjmisc)


## Load data 

rna_dat <- readRDS("./data/derivedData/rna_complete.RDS")

# Load code_key
code <- read_excel("./data/code.key.xlsx")

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



#### Basic model 

# We are interested in if change in GLUCOSE is different from change in 
# PLACEBO. We will use an ANCOVA model with:

# change ~ pre + condition

# The effect of condition will answer the question.



change_dat <- rna_dat2 %>%
     filter(technical == "in", 
           outlier == "in") %>%
    dplyr::select(subject:rep, supplement, RNA.weight) %>%
    
    group_by(subject, time, supplement) %>%
    summarise(RNA.weight = mean(RNA.weight, na.rm = TRUE)) %>%
    
    pivot_wider(names_from = time, 
                values_from = RNA.weight) %>%
    ungroup() %>%
    mutate(change = log(post) - log(pre), 
           pre = pre - mean(pre, na.rm = TRUE), 
           supplement = factor(supplement, levels = c("PLACEBO", "GLUCOSE"))) %>%
        print()
   # print()
  #  pivot_longer(names_to = "time",
   #              values_to = "change2",
    #             cols = (pre:post)) %>%
    #select(subject, supplement, time, change2) %>%
    #print()




# Create model: 
# Needs to have an intercept per participant (mixed model)
# Control for pre values.

m1 <- lmerTest::lmer(change ~ 0 + pre + supplement + (1|subject), 
           data = change_dat)

    
#m2 <- lmerTest::lmer(change2 ~ time + supplement:time + (1|subject),
 #                    data = change_dat)
plot(m1)

summary(m1)




### Get estimated means from the model, these are average increase at 
# pre = 0 (the average pre value)
# remember that these are originally log-fold change values, transformed back using exp(emmean) etc

confint(emmeans(m1, specs = ~"supplement")) %>%
    data.frame() %>%
        mutate(time = "post") %>%
        add_row(supplement = c("PLACEBO", "GLUCOSE"), time = "pre", emmean = 0, SE = 0, lower.CL = 0, upper.CL = 0) %>%
        mutate(time = factor(time, levels = c("pre", "post"), labels = c("Baseline", "Post"))) %>%
       # print()
    ggplot(aes(time, exp(emmean), group = supplement, fill = supplement)) +
    geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)),
                  width = 0.1,
                  position = position_dodge(width = 0.2)) +
        geom_line(position = position_dodge(width = 0.2)) +
    geom_point(position = position_dodge(width = 0.2), shape = 21, size = 3) +
    labs(x = "Supplement", y = "Total RNA by biopsy weight \n(fold change)\n", fill = "Groups",
         caption = "Values are fold change and confidence limits") +
    theme_classic()




# Legge inn baseline for å vise endringen
# Bachelor: absolutte tall pre-post, angi tidseffekter
# volum, intensitet (% 1RM), glukose, total RNA, humac - absolutte verdier 
# rRNA: 
# Table of summary
tab_model(m1)

tab_model(m1,
          pred.labels = c("Glucose", "Placebo", "Placebo vs. Glucose"),
          string.pred = "Coefficient",
          string.ci = " Conf.Int (95%)",
          string.p = "P-value",
          dv.labels = "Change",
          show.re.var = FALSE,
          show.icc = FALSE)

##### OLD CODE BELOW ####################


# Calculate mean biopsy weight for plotting models
biopsy_size <- mean(rna_dat2$weight)

# Model total RNA per weight


m <- lme(log(RNA/weight) ~ time * supplement, 
         random = list(subject = ~ 1, 
                       biopsy = ~ 1), 
         data = rna_dat2)

plot(m)


m2 <- lme(log(RNA) ~ weight.mc + time * supplement, 
          random = list(subject = ~ 1, 
                        biopsy = ~ 1), 
          data = rna_dat2)


plot(m2)


m3 <- lme(log(RNA) ~ weight.mc + technical + time * supplement, 
          random = list(subject = ~ 1), 
          data = rna_dat2)

plot(m3)



summary(m3)

emmeans(m3, specs = ~ time | supplement) %>%
        data.frame() %>%
        ggplot(aes(time, exp(emmean) / biopsy_size, group = supplement, fill = supplement)) + 
        geom_errorbar(aes(ymin = exp(lower.CL) / biopsy_size, ymax = exp(upper.CL) / biopsy_size), 
                      position = position_dodge(width = 0.2), 
                      width = 0.1) +
       # geom_line(position = position_dodge(width = 0.2)) +
        geom_point(shape = 21, position = position_dodge(width = 0.2),
                   size = 3) +
        labs(x = "Time-point", y = "Mean total RNA per biopsy size", fill = "Supplement",
             title = "Mean total RNA pre vs. post",
             subtitle = "Different patterns per supplement",
             caption = "Values are mean and SD") +
        theme_bw() +
        theme(plot.background = element_rect(fill = "gray80"))

# Contemporary presentation, keep this or change to fold-change?
# Make a table containing rRNA outcomes and p-values/CI

### Get fold-changes from mixed-model

placebo.pre <- 
        
        
        
        
        
        
        
        
        
        summary(m)


rna_dat2 %>%
        #    filter(technical == "in") %>%
        
        group_by(subject, time, supplement) %>%
        mutate(RNA = RNA/weight) %>%
        summarise(RNA = mean(RNA)) %>%
        
        pivot_wider(names_from = time, 
                    values_from =  RNA) %>%
        
        mutate(change = post / pre, 
               pre = pre - mean(pre)) %>%
        
        dplyr::select(subject, supplement, change) %>%
        pivot_wider(names_from = supplement, values_from = change) %>%
        
        mutate(diff = GLUCOSE / PLACEBO, 
               incr = if_else(diff > 1, "favor_glucose", "favor_placebo")) %>%
        group_by(incr)  %>%
        summarise(n = n()) %>%
        
        print()



group_by(supplement, incr) %>%
        summarise(n = n(), 
                  fc = mean(change, na.rm = TRUE)) %>%
        print()






lme(log(post) ~ supplement, 
    random = list(subject = ~ 1), 
    data = ., 
    na.action = na.omit) %>%
        summary()
print()

                 







