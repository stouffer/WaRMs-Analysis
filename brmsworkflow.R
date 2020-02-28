####
# sample workflow for data analysis from manuscript translated into the BRMS library
# names of models remain the same for both poisson and beta distribution
# this code fits the baseline, ambient(no treatment effect), and all the treatment effect models to a sample dataset
####

# load the sample data frame
load('example.Rdata')

# read in the model-fitting function
source('brmsModelFunctions.R')
# functions take dataset, location, elevation, wont work with family bc need a different formula

# data for poisson models CA_low_stan

# fit null or baseline model 
null.model <- null.mod(CA_low_stan)

# Save model outputs if needed
# saveRDS(null.model, "nullmodel_CA_Low.rds")


# fit ambient or no treatment effect model
ambient.model <- Amb.mod(CA_low_stan)

# fit single treatment effects model
warming.model <- Warm.mod(CA_low_stan)
removal.model <- Rem.mod(CA_low_stan)

# fit both treatments
both.warmandremoval.model <- NoInteraction.mod(CA_low_stan)

# fit full model with both treatments and interaction
withinteraction.model <- Full.mod(CA_low_stan)

# model comparison
null.model <- add_criterion(null.model, "waic")
ambient.model <- add_criterion(ambient.model, "waic")
removal.model <- add_criterion(removal.model, "waic") 
warming.model <- add_criterion(warming.model, "waic")
both.warmandremoval.model <- add_criterion(both.warmandremoval.model, "waic")
withinteraction.model <- add_criterion(withinteraction.model, "waic")


CA_low_waic <- loo_compare(null.model,  ambient.model, removal.model, warming.model, 
                       both.warmandremoval.model, withinteraction.model, criterion = "waic")

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

### fit using beta distribution
# model functions have same name
source('brmsBetaFunctions.R')

# need beta data CA_low_stan_beta 
# this is the same as the above data but now abundance measures are between 1 and 0

# fit null or baseline model 
null.model.beta <- null.mod(CA_low_stan)

# fit ambient or no treatment effect model
ambient.model.beta  <- Amb.mod(CA_low_stan)

# fit single treatment effects model
warming.model.beta  <- Warm.mod(CA_low_stan)
removal.model.beta  <- Rem.mod(CA_low_stan)

# fit both treatments
both.warmandremoval.model.beta  <- NoInteraction.mod(CA_low_stan)

# fit full model with both treatments and interaction
withinteraction.model.beta <- Full.mod(CA_low_stan)



