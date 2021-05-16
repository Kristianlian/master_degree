## Preliminary total RNA analysis

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