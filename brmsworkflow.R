####
# sample workflow for data analysis from manuscript translated into the BRMS library
# names of models remain the same for both poisson and beta distribution
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
both.warmandremoval.model <- model.fit(CA_low_stan, "Removal_plus_warming")

# fit full model with both treatments and interaction
withinteraction.model <- model.fit(CA_low_stan, "Removal_times_warming")

# model comparison with WAIC
null.model <- add_criterion(null.model, "waic")
ambient.model <- add_criterion(ambient.model, "waic")
removal.model <- add_criterion(removal.model, "waic") 
warming.model <- add_criterion(warming.model, "waic")
both.warmandremoval.model <- add_criterion(both.warmandremoval.model, "waic")
withinteraction.model <- add_criterion(withinteraction.model, "waic")


CA_low_waic <- loo_compare(null.model,  ambient.model, removal.model, warming.model, 
                       both.warmandremoval.model, withinteraction.model, criterion = "waic")

# model comparison with WAIC weights
model_weights(null.model,  ambient.model, removal.model, warming.model, 
              both.warmandremoval.model, withinteraction.model,  
              weights = "waic") %>%
  as_tibble() %>% 
  rename(weight = value) %>% 
  mutate(model  = c("null", "amb", "removal", "warm", "both", "full"),
         weight = weight %>% round(digits = 2)) %>% 
  select(model, weight) %>% 
  arrange(desc(weight)) %>% 
  knitr::kable()

######## fit using beta distribution #####################
# need beta data CA_low_stan_beta 
# this is the same as the above data but now abundance measures are between 1 and 0
# results not presented in manuscript

# fit null or baseline model 
null.model.beta <- model.fit.beta(CA_low_stan_beta, "Null")

# fit ambient or no treatment effect model
ambient.model.beta  <- model.fit.beta(CA_low_stan_beta, "Ambient")

# fit single treatment effects model
warming.model.beta  <- model.fit.beta(CA_low_stan_beta, "Warm")
removal.model.beta  <- model.fit.beta(CA_low_stan_beta, "Removal")

# fit both treatments
both.warmandremoval.model.beta  <- model.fit.beta(CA_low_stan_beta, "NoInteraction")

# fit full model with both treatments and interaction
withinteraction.model.beta <- model.fit.beta(CA_low_stan_beta, "Full")

# the same comparison workflow as above can be used on the beta models. 

