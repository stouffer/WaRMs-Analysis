####
# sample workflow for data analysis from manuscript with the BRMS library
# this code fits the baseline, ambient(no treatment effect), and all the treatment effect models to a sample dataset
####

# libraries not loaded with other script
require('tidyverse')

# load the sample data frame
load('example.Rdata')

# read in the model-fitting function
source('brmsModelFunctions.R')

######## fit using poisson distribution #####################
# data for poisson models is called CA_low_stan

# fit null or baseline model 
null.model <- model.fit(CA_low_stan, "Null")

# Save model outputs if needed
# saveRDS(null.model, "nullmodel_CA_Low.rds")

# fit ambient or no treatment effect model
ambient.model <- model.fit(CA_low_stan, "Ambient")

# fit single treatment effects model
warming.model <- model.fit(CA_low_stan, "Warm")
removal.model <- model.fit(CA_low_stan, "Removal")

# fit both treatments
removalpluswarming.model <- model.fit(CA_low_stan, "Removal_plus_warming")

# fit full model with both treatments and interaction
removaltimeswarming.model <- model.fit(CA_low_stan, "Removal_times_warming")

# model comparison with WAIC
null.model <- add_criterion(null.model, "waic")
ambient.model <- add_criterion(ambient.model, "waic")
removal.model <- add_criterion(removal.model, "waic") 
warming.model <- add_criterion(warming.model, "waic")
removalpluswarming.model <- add_criterion(removalpluswarming.model, "waic")
removaltimeswarming.model <- add_criterion(removaltimeswarming.model, "waic")


CA_low_waic <- loo_compare(null.model,  ambient.model, removal.model, warming.model, 
                           removalpluswarming.model, removaltimeswarming.model, criterion = "waic")

# model comparison with WAIC weights
model_weights(null.model,  ambient.model, removal.model, warming.model, 
              removalpluswarming.model, removaltimeswarming.model,  
              weights = "waic") %>%
  as_tibble() %>% 
  rename(weight = value) %>% 
  mutate(model  = c("Null", "Ambient", "Removal", "Warm", "Removal_plus_warming", "Removal_times_warming"),
         weight = weight %>% round(digits = 2)) %>% 
  select(model, weight) %>% 
  arrange(desc(weight)) %>% 
  knitr::kable()

