####
# code that samples from the posterior to demonstrate plausibility of models when not fit to data
####

# read in the model-fitting functions
source('brmsBetaFunctions.R')

# load the sample data frame
# data frame names are labelled as location_elevations_stan.rds
example_data <- readRDS('../data/CA_low_stan.rds')

# define the different models which are to be fit to the data
models <- c(
  "Null",
  "Recruitment",
   "Ambient",
   "Warming",
   "Removal",
   "Removal_+_Warming",
  "Removal_x_Warming"
)

# generate samples purely from the prior to visualize prior predictions
model.fits.priors <- sapply(
  models,
  function(model.name,df){
    fit <- model.fit.beta(df, model.name, sample_prior="only")
    return(fit)
  },
  df = example_data,
  simplify = FALSE
)

# use conditional effects to plot prior predictive distributions for the null model
conditional_effects(model.fits[["Null"]])

# use conditional effects to plot prior predictive distributions for the full model
conditional_effects(model.fits[["Removal_x_Warming"]])
