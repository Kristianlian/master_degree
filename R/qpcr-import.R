### QPCR IMPORT #################



### Import and cleaning of qPRC raw data

#This script imports and saves a data frame with cq-values per sample and average efficiencies.
# Data are saved in data/derivedData/qpcr-data.RDS


library(qpcR)
library(qpcrpal)
library(parallel)
library(tidyverse)

## Load batch 
batch <- prepare_batch("./data/qpcr_raw/", skip = 19, equipment = "quant")

# Model all runs
models <- model_qpcr(batch)


model.test <- test_models(models)

model.test$figure

best_fit_models <- model.test$results %>%
        group_by(target) %>%
        filter(n == max(n)) %>%
        data.frame()

effi <- list()
cqs <- list()


for(i in 1:nrow(best_fit_models)) {
        
        
        batch_subset <- batch %>%
                filter(target == best_fit_models[i,1])
        
        models <- model_qpcr(batch_subset, model = best_fit_models[i,2])
        
        
        effi[[i]] <- analyze_efficiency(models, method = "cpD2", model = "linexp")
        cqs[[i]] <- analyze_models(models)
        
}

 
results_cq <- bind_rows(cqs)
results_eff <- bind_rows(effi) 



results_eff %>%
        separate(ID, sep = "_", into = c("sample", "n1", "n2", "target")) %>%
        filter(eff < 100) %>%
        ggplot(aes(target, eff)) + geom_point()

efficiencies <- results_eff %>%
        separate(ID, sep = "_", into = c("sample", "n1", "n2", "target")) %>%
        filter(eff < 100) %>%
        group_by(target) %>%
        summarise(eff = mean(eff, na.rm = TRUE)) %>%
        print()

import <- results_cq %>%
        separate(ID, sep = "_", into = c("sample", "n1", "n2", "target")) %>%
        dplyr::select(sample, target, cq = cpD2) %>%
        inner_join(efficiencies) %>%
        print()



saveRDS(import, "./data/derivedData/qpcr-data.RDS")
