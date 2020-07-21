####
# functions to fit population dynamics models within a bayesian hierarchical framework
####

##################### beta distribution ##################### 
# slightly different model formula as things are now on 0 to 1 scale

library(brms)
options(mc.cores = parallel::detectCores())

model.formula <- function(model.name){
  if(model.name == "Null"){
    formula <- bf(
      abundance ~ (1/(1+exp(-q))) + prev_abund, 
      q ~ 1 + (1 | ID | focal), #allows influx intercept (ambient plots) to vary with focal
      nl=T
    )
  }
  
  if(model.name == "Ambient"){
    # define the model structure
    formula <- bf(
      abundance ~ (1/(1+exp(-q))) + prev_abund*exp(g),
      q ~ 1 + (1 | ID | focal),
      g ~ 1 + (1 | ID | focal), #allows growth intercept (ambient plots) to vary with focal
      nl = T # here ID is used to link the two levels of models
    )
  }
  
  
  if(model.name == "Removal"){
    # define the model structure
    formula <- bf(
      abundance ~ (1/(1+exp(-q))) + prev_abund*exp(g),
      q ~ 1 + (1 | ID | focal),
      g ~ 1 + removal + (1 + removal | ID | focal), # density-dependent growth varies  with removal treatment
      nl = T
    )
  }
  
  if(model.name == "Warm"){
    # define the model structure
    formula <- bf(
      abundance ~ (1/(1+exp(-q))) + prev_abund*exp(g),
      q ~ 1 + (1 | ID | focal),
      g ~ 1 + warming + (1 + warming | ID | focal),# density-dependent growth varies with warming treatment
      nl = T
    )
  }
  
  
  if(model.name == "Removal_plus_warming"){ # R plus W in Table 1
    # define the model structure
    formula <- bf(
      abundance ~ (1/(1+exp(-q))) + prev_abund*exp(g),
      q ~ 1 + (1 | ID | focal),
      g ~ 1 + warming + removal + (1 + warming + removal | ID | focal), # contains both removal and warming treatment but not their interaction
      nl = T
    )
  }
  
  if(model.name == "Removal_times_warming"){ # R times W in Table 1
    # define the model structure
    formula <- bf(
      abundance ~ (1/(1+exp(-q))) + prev_abund*exp(g),
      q ~ 1 + (1 | ID | focal),
      g ~ 1 + warming + removal + removal:warming + (1 + warming + removal + removal:warming | ID | focal), # both treatments and the interaction
      nl = T
    )
  }
  
  return(formula)
}


model.fit.beta <- function(df, model.name){
  formula <- model.formula(model.name)
  
  # fits the appropriate model with brms
  mf <- brm(
    formula = formula,
    family = zero_inflated_beta(link="identity"),
    #prior = prior, # used default priors
    data = df,
    # iterations and chains used in manuscript can adjust as necessary
    iter = 6000,
    warmup = 1000,
    chains = 2, # default is 4 chains
    cores = 2, # specify cores used
    control = list(adapt_delta=0.99, max_treedepth = 20),
    seed = 12 # set seed
  )
  
  return(mf)
}


