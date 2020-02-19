# Stan Models as function
library(brms)
options(mc.cores = parallel::detectCores())

# functions require dataset, location, and elevation

null.mod <- function(df){
  null  <- (brm(
    data = df, 
    family = zero_inflated_beta(link="identity"),
    # use a function to tell brms it is a multi-level model
    bf( 
      abundance ~  (1/(1+exp(-p))) + prev_abund, 
      p ~ 1 + (1 | ID| focal),
      nl=T #allows intercept (immigration) to vary with focal
    ),
    #iter = 6000, warmup = 1000, chains = 2, cores = 2, # default is 4 chains
    seed = 12 #set seed
  ))
  return(null)
}


Amb.mod <- function(df){
  formula <- #multilevel formulat
    bf(abundance ~(1/(1+exp(-p))) + prev_abund*exp(growth),
                p ~ 1 + (1 | ID1 | focal), #allows intercept (immigration) to vary with focal
                growth ~ 1 + (1 | ID1 | focal), nl = T) #allows growth intercept (ambient plots) to vary with focal
# here ID is used to link the two levels of models
  Amb <- brm(
    formula = formula,
    family = zero_inflated_beta(link="identity"),
    data = df,
    #iter = 6000, warmup = 1000, chains = 2, cores = 2, # default is 4 chains
    control = list(adapt_delta=0.95, max_treedepth = 20),
    seed = 12 #set seed
  )
  return(Amb)
}

## single treatment models
# here, under growth, each treatment and the intercept are allowed to vary by focal
Rem.mod <- function(df){
  formula <- bf(abundance ~ (1/(1+exp(-p))) + prev_abund*exp(growth),
                p ~ 1 + (1 | ID1 | focal),
                growth ~ 1 +removal + (1 + removal| ID1 | focal), nl = T)
  Rem <- brm(
    formula = formula,
    family = zero_inflated_beta(link="identity"),
    data = df,
    #iter = 4000, warmup = 1000, chains = 1, cores = 1, # default is 4 chains
    control = list(adapt_delta=0.99, max_treedepth = 20),
    seed = 12 #set seed
  )
  return(Rem)
}

Warm.mod <- function(df){
  formula <- bf(abundance ~ (1/(1+exp(-p))) + prev_abund*exp(growth),
                p ~ 1 + (1 | ID | focal),
                growth ~ 1 + warming + (1 +warming| ID | focal), nl = T)
  
  Warm <- brm(
    formula = formula,
    family = zero_inflated_beta(link="identity"),
    data = df,
  # iter = 6000, warmup = 1000, chains = 2, cores = 2, # default is 4 chains
   control = list(adapt_delta=0.99, max_treedepth = 20),
    seed = 12 #set seed
  )
  return(Warm)
}

## multitreatment models
# here, under growth, both treatments (and their interaction) along with ambient vary by focal
NoInteraction.mod <- function(df){
  formula <- bf(abundance ~ (1/(1+exp(-p))) + prev_abund*exp(growth),
                p ~ 1 + (1 | ID | focal),
                growth ~ 1 + warming+ removal + (1 + warming+ removal | ID | focal), nl = T)
  NoInt <- brm(
    formula = formula,
    family = zero_inflated_beta(link="identity"),
    data = df,
    #   iter = 6000, warmup = 1000, chains = 2, cores = 2, #will up later
    control = list(adapt_delta=0.99, max_treedepth = 20),
    seed = 12 #set seed
  )
  return(NoInt)
}


Full.mod <- function(df){
  formula <- bf(abundance ~ (1/(1+exp(-p))) + prev_abund*exp(growth),
                p ~ 1 + (1 | ID | focal),
                growth ~ 1 + warming+ removal + removal:warming + (1 + warming+ removal + removal:warming| ID | focal), nl = T)
  Full <- brm(
    formula = formula,
    family = zero_inflated_beta(link="identity"),
    data = df,
   # iter = 6000, warmup = 1000, chains = 2, cores = 2, #default is 4 chains
    control = list(adapt_delta=0.99, max_treedepth = 20),
    seed = 12 #set seed
  )
  return(Full)
}


