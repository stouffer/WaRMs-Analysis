####
# functions to post process a list of model fits
####

# model fit statistics are computed with brms
require(brms)

# contain functions used in preparation of the results table
require(tidyverse)
require(knitr)

# given a list of model fits, compile a table with the model comparison information
results.table <- function(model.fits){
  # compute WAICs for all model fits
  model.waic.summaries <- lapply(
    model.fits,
    waic
  )

  # extract the estimated WAIC for each model
  model.waics <- sapply(
    model.waic.summaries,
    function(x){ x$estimates["waic","Estimate"] }
  )

  # extract the effective number of parameters
  model.params <- sapply(
    model.waic.summaries,
    function(x){ x$estimates["p_waic","Estimate"] }
  )

  # compute model weights based on WAIC values
  model.weights <- exp(-(model.waics - min(model.waics)) / 2)
  model.weights <- model.weights / sum(model.weights)

  # put everything together into a data frame
  results <- data.frame(
    row.names = names(model.fits),
    WAIC = model.waics,
    pWAIC = model.params,
    Weight = model.weights
  )

  # tidy up the results table
  results <- results %>%
    as_tibble() %>%
    arrange(desc(Weight)) %>%
    mutate(Model = rownames(results)) %>%
    mutate(WAIC = WAIC %>% round(digits = 2)) %>%
    mutate(pWAIC = pWAIC %>% round(digits = 2)) %>%
    mutate(Weight = Weight %>% round(digits = 2)) %>%
    select(Model, WAIC, pWAIC, Weight) %>%
    kable()

  return(results)
}
