####
# functions to fit population dynamics models within a bayesian hierarchical framework
####

library(brms)
options(mc.cores = parallel::detectCores())

model.formula <- function(model.name){
  if(model.name == "Null"){
    formula <- bf(
      abundance ~ (100/(1+exp(-q))) + prev_abund,
      q ~ 1 + (1 | ID | focal),
      nl=T
    )
  }

  if(model.name == "Ambient"){
    # define the model structure
    formula <- bf(
      abundance ~ (100/(1+exp(-q))) + prev_abund*exp(g),
      q ~ 1 + (1 | ID | focal),
      g ~ 1 + (1 | ID | focal),
      nl = T
    )
  }

  if(model.name == "Removal"){
    # blah blah blah
  }

  return(formula)
}

model.prior <- function(df, formula, model.name){
  # used to specify the prior for the random effects  
  prior <- get_prior(formula, data = df, family = poisson())
  if(model.name != "Null"){
    # define cor prior as lkj 2
    prior$prior[2] <- "lkj(2)"

    # sd from focal is cauchy
    prior$prior[8] <- "cauchy(0, 2)"
  }

  return(prior)
}

model.fit <- function(df, model.name){
  formula <- model.formula(model.name)
  prior <- model.prior(df, formula, model.name)

  # fits the appropriate model with brms
  mf <- brm(
    formula = formula,
    family = poisson(link ="identity"),
    prior = prior,
    data = df,
    iter = 6000,
    warmup = 1000,
    chains = 4,
    cores = 4, 
    control = list(adapt_delta=0.99, max_treedepth = 20),
    seed = 12 #set seed
  )

  return(mf)
}

# influx varies by focal species
# density-dependent growth varies by focal species and with removal treatment
Rem.mod <- function(df){
  formula <- bf(
    abundance ~ (100/(1+exp(-q))) + prev_abund*exp(g),
    q ~ 1 + (1 | ID | focal),
    g ~ 1+ removal +(1 + removal | ID | focal),
    nl = T
  )
  
  prior <- get_prior(formula, data = df, family = poisson())
  prior$prior[2] <- "lkj(2)"  # define cor prior as lkj 2
  prior$prior[8] <- "cauchy(0, 2)" # sd from focal is cauchy
  
  Rem <- brm(
    formula = formula,
    family = poisson(link ="identity"),
    prior = prior,
    data = df,
    iter = 6000, warmup = 1000, chains = 4, cores = 4, 
    control = list(adapt_delta=0.99, max_treedepth = 20),
    seed = 12 #set seed
  )
  return(Rem)
}

# influx varies by focal species
# density-dependent growth varies by focal species and with warming treatment
Warm.mod <- function(df){
  formula <- bf(
    abundance ~  (100/(1+exp(-q)))+ prev_abund*exp(g),
    q ~ 1 + (1 | ID | focal),
    g ~ 1 + warming + (1 + warming | ID | focal),
    nl = T
  )
  
  prior <- get_prior(formula, data = df, family = poisson())
  prior$prior[2] <- "lkj(2)"   # define cor prior as lkj 2
  prior$prior[8] <- "cauchy(0, 2)" # sd from focal is cauchy
  
  Warm <- brm(
    formula = formula,
    family = poisson(link ="identity"),
    prior = prior,
    data = df,
    iter = 6000, warmup = 1000, chains = 4, cores = 4, 
    control = list(adapt_delta=0.99, max_treedepth = 20),
    seed = 12 #set seed
  )
  return(Warm)
}

# influx varies by focal species
# density-dependent growth varies by focal species and with removal and warming treatments but no interaction
NoInteraction.mod <- function(df){
  formula <- bf(
    abundance ~ (100/(1+exp(-q))) + prev_abund*exp(g),
    q ~ 1 + (1 | ID | focal),
    g ~ 1 + warming + removal + (1 + warming + removal | ID | focal),
    nl = T
  )
  
  prior <- get_prior(formula, data = df, family = poisson())
  prior$prior[2] <- "lkj(2)"  # define cor prior as lkj 2
  prior$prior[8] <- "cauchy(0, 2)" # sd from focal is cauchy
  
  NoInt <- brm(
    formula = formula,
    family = poisson(link ="identity"),
    prior = prior,
    data = df,
    iter = 6000, warmup = 1000, chains = 4, cores = 4, 
    control = list(adapt_delta=0.99, max_treedepth = 20),
    seed = 12 #set seed
  )
  return(NoInt)
}

Full.mod <- function(df){
  formula <- bf(
    abundance ~ (100/(1+exp(-q))) + prev_abund*exp(g),
    q ~ 1 + (1 | ID | focal),
    g ~ 1 + warming + removal + removal:warming + (1 + warming + removal + removal:warming | ID | focal),
    nl = T
  )
  
  prior <- get_prior(formula, data = df, family = poisson())
  prior$prior[2] <- "lkj(2)"  # define cor prior as lkj 2
  prior$prior[8] <- "cauchy(0, 2)" # sd from focal is cauchy
  Full <- brm(
    formula = formula,
    family = poisson(link ="identity"),
    prior = prior,
    data = df,
    iter = 6000, warmup = 1000, chains = 4, cores = 4, 
    control = list(adapt_delta=0.99, max_treedepth = 20),
    seed = 12 #set seed
  )
  return(Full)
}


