#### Glucose percentage change

## Packages
library(tidyverse);library(readxl);library(nlme);library(lme4);library(knitr);library(broom);library(emmeans);library(dplyr)

## Data handling

# Intervention data
gluc.dat <- read_excel("./data/glucose/fingerstick.xlsx", na = "NA")
gluc.dat$sample_time <- as.character(gluc.dat$sample_time)

glu.dat <- gluc.dat %>%
        mutate(time = if_else(sample_time == "0",
                              "baseline",
                              if_else(sample_time == "45",
                                      "45min",
                                      if_else(sample_time == "90",
                                              "90min",
                                              if_else(sample_time == "120",
                                                      "120min",
                                                      if_else(sample_time == "135",
                                                              "135min",
                                                              if_else(sample_time == "150",
                                                                      "150min",
                                                                      if_else(sample_time == "270",
                                                                              "270min", sample_time))))))),
               time = factor(time, levels = c("baseline", "45min", "90min", "120min", "135min", "150min", "270min"))) %>%
        mutate(glu = as.numeric(glu),
               lak = as.numeric(lak))

# Control data
gluc.resp <- read_excel("./data/glucose/ribose.glucresponsetest.xlsx", na = "NA")
gluc.resp$sample_time <- as.character(gluc.resp$sample_time)

glu.con <- gluc.resp %>%
        mutate(time = if_else(sample_time == "0",
                              "baseline",
                              if_else(sample_time == "45",
                                      "45min",
                                      if_else(sample_time == "90",
                                              "90min",
                                              if_else(sample_time == "120",
                                                      "120min",
                                                      if_else(sample_time == "135",
                                                              "135min",
                                                              if_else(sample_time == "150",
                                                                      "150min",
                                                                      if_else(sample_time == "270",
                                                                              "270min", sample_time))))))),
               time = factor(time, levels = c("baseline", "45min", "90min", "120min", "135min", "150min", "270min"))) %>%
        mutate(glu = as.numeric(glu))



# All subjects percentage change
glu.joined <- glu.dat %>% 
        full_join(glu.con) %>%
        select(subject, timepoint, time, glu, supplement) %>%
        na.omit() %>%
        group_by(supplement, time) 
  print()
        #summarise(Mean = mean(glu, na.rm = TRUE)) %>% #,
        # SD = sd(glu, na.rm = TRUE)) 
       # pivot_wider(names_from = time,
              #      values_from = Mean) %>%
        #print()

## Fit model

m1 <- lmer(glu ~ time * supplement + (1|subject), data = glu.joined)
plot(m1)
summary(m1)
confint(m1)


# Get estimated means

emmeans(m1, specs = ~ "time|supplement") %>%
  data.frame() %>%
  ggplot(aes(time, emmean, group = supplement, fill = supplement)) +
  geom_errorbar(aes(ymin= lower.CL, ymax = upper.CL), position = position_dodge(width = 0.2), 
                width = 0.05) +
  geom_line(position = position_dodge(width = 0.2)) +
  geom_point(shape = 21, position = position_dodge(width = 0.2),
             size = 3) +
  labs(x = "Time-point", y = "Blood glucose \n(mmol/L)\n", color = "Supplement",
       caption = "Values are mean and CL") +
  theme_bw() +
  theme(plot.background = element_rect(fill = "gray80"))



# Figure: mean percentage change in blood glucose levels


## To do:
# Scale timeline to actual distance between measuring 
# Add spread to figure (CL)
# Add arrows for protein ingestion, RT (as e.g. a shaded area), and statistics to figure (power.point)

# Contemporary figure
glu.change <- read_excel("./data/glucose/ribose.gluchange.xlsx", na = "NA") %>%
  select(supplement, cbaseline, c45min, c90min, c120min, c135min, c150min, c270min) %>%
  mutate(c270min = as.numeric(c270min)) %>%
  pivot_longer(names_to = "time",
               values_to = "change",
               cols = cbaseline:c270min) %>%
  print()

glu.change %>%
  mutate(time = factor(time, levels = c("cbaseline", "c45min", "c90min", "c120min", "c135min", "c150min", "c270min"))) %>%
  ggplot(aes(time, change, group = supplement, fill = supplement)) +
  geom_line() +
  geom_point(shape = 21, size = 3) +
  scale_x_discrete(breaks = c("cbaseline", "c45min", "c90min", "c120min", "c135min", "c150min", "c270min"),
                   labels = c("Baseline", "45min post \nprotein\n", "90min post \nprotein\n",
                              "Pre RT", "Middle of RT", "Post RT", "2hrs \npost RT\n")) +
  labs(x = "Time-point", y = "Blood glucose \n(percentage change from baseline)\n", color = "Supplement",
       caption = "Values are mean percentage change per supplement") +
  theme_bw() +
  theme(plot.background = element_rect(fill = "gray80"))