####
# functions to fit population dynamics models within a bayesian hierarchical framework
####

library(brms)
options(mc.cores = parallel::detectCores())

# funtions require dataset
null.mod <- function(df){
  null  <- brm(
    data = df, 
    family = poisson(link ="identity"),
    # use a function to tell brms it is a multi-level model
    bf( 
      abundance ~  (100/(1+exp(-p))) + prev_abund, 
      p ~ 1 + (1 | ID| focal),
      nl=T #allows intercept (immigration) to vary with focal
    ),
    iter = 6000, warmup = 1000, chains = 4, cores = 4, 
    control = list(adapt_delta=0.99, max_treedepth = 20),
    seed = 12 #set seed
  )
  return(null)
}


Amb.mod <- function(df){
  formula <- bf(abundance ~ (100/(1+exp(-p))) + prev_abund*exp(growth),
                p ~ 1 + (1 | ID | focal),
                growth ~ 1 + (1 | ID | focal), nl = T) #allows growth intercept (ambient plots) to vary with focal
  # here ID is used to link the two levels of models
  
  prior <- get_prior(formula, data = df, family = poisson())
  prior$prior[2] <- "lkj(2)"  # define cor prior as lkj 2
  prior$prior[8] <- "cauchy(0, 2)" # sd from focal is cauchy
  
  Amb <- brm(
    formula = formula,
    family = poisson(link ="identity"),
    prior = prior,
    data = df,
    iter = 6000, warmup = 1000, chains = 4, cores = 4, 
    control = list(adapt_delta=0.99, max_treedepth = 20),
    seed = 12 #set seed
  )
  return(Amb)
}

## single treatment models
# here, under growth, each treatment and the intercept are allowed to vary by focal
Rem.mod <- function(df){
  formula <- bf(abundance ~ (100/(1+exp(-p))) + prev_abund*exp(growth),
                p ~ 1 + (1 | ID | focal),
                growth ~ 1+removal +(1 +removal| ID | focal), nl = T)
  
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

Warm.mod <- function(df){
  formula <- bf(abundance ~  (100/(1+exp(-p)))+ prev_abund*exp(growth),
                p ~ 1 + (1 | ID | focal),
                growth ~ 1 +warming+(1 +warming| ID | focal), nl = T)
  
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

## multitreatment models
# here, under growth, both treatments (and their interaction) along with ambient vary by focal
NoInteraction.mod <- function(df){
  formula <- bf(abundance ~  (100/(1+exp(-p))) + prev_abund*exp(growth),
                p ~ 1 + (1 | ID | focal),
                growth ~ 1 + warming+ removal +(1 + warming+ removal + removal:warming| ID | focal), nl = T)
  
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
  formula <- bf(abundance ~ (100/(1+exp(-p))) + prev_abund*exp(growth),
                p ~ 1 + (1 | ID | focal),
                growth ~ 1 + warming+ removal + removal:warming+(1 + warming+ removal + removal:warming| ID | focal), nl = T)
  
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


