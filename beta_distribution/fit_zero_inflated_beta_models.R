####
# code that fits all candidate models to an example dataset and generates a model comparison table
####

# read in the model-fitting functions
source('brmsBetaFunctions.R')

# read in the post-processing functions
source('resultsFunctions.R')

# load the sample data frame
example_data <- readRDS('../example.Rdata')

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
results <- results.table(model.fits)
