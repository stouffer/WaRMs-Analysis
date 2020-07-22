####
# code that fits all candidate models to an example dataset and generates a model comparison table
####

# model fits are performed with brms
library(brms)

# allow for parallelization of MCMC chains
options(mc.cores = parallel::detectCores())

# read in the model-fitting functions
source('brmsBetaFunctions.R')

# load the sample data frame
example_data <- readRDS('../example.Rdata')

# define the different models which are to be fit to the data
models <- c(
  "Null",
  "Ambient",
  "Warming",
  "Removal",
  "Removal_+_Warming",
  "Removal_x_Warming"
)

# fit each of the above specified models to the data
model.fits <- sapply(
  models,
  function(model.name,df){
    fit <- model.fit.beta(df, model.name)
    return(fit)
  },
  df = example_data,
  simplify = FALSE
)

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

# put everything together into a data frame
results <- data.frame(
  row.names = models,
  WAIC = model.waics,
  pWAIC = model.params,
  Weight = model.weights
)

# # model comparison with WAIC weights
# model_weights(null.model,  ambient.model, removal.model, warming.model, 
#               removalpluswarming.model, removaltimeswarming.model,  
#               weights = "waic") %>%
#   as_tibble() %>% 
#   rename(weight = value) %>% 
#   mutate(model  = c("Null", "Ambient", "Removal", "Warm", "Removal_plus_warming", "Removal_times_warming"),
#          weight = weight %>% round(digits = 2)) %>% 
#   select(model, weight) %>% 
#   arrange(desc(weight)) %>% 
#   knitr::kable()

