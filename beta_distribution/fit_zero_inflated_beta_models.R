#' @title Population dynamics model fitting workflow
#' @description Code that fits all candidate models to an example dataset and generates a model comparison table
#' @author Michelle 
#' @param models can be one or a list from the following options "Null", "Recruitment", "Ambient", "Warming", "Removal", "Removal_+_Warming", "Removal_x_Warming" required for model.formula function 
#' @param example_data data frame being analyzed, this is a single site as an example but other sites are provided in the same format
#' @param formula out put of model.formula function
#' @return model.fits applies the model.fit.beta function to the list of model names and returns one or a list of brms formulas
#' @return results.table() returns a table of waic weights of the models run



# read in the model-fitting functions
source('brmsBetaFunctions.R')

# read in the post-processing functions
source('resultsFunctions.R')

# load the sample data frame
# data frame names are codes as location_elevations_stan.rds
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

# generate a pretty table of model comparison statistics
# note that despite the same seed settings, there may be slight differences in model WAIC and weights to that reported in the manuscript
results <- results.table(model.fits)
