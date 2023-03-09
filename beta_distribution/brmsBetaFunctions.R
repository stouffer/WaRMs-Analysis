#' @title Population dynamics model functions
#' @description These functions fit population dynamics models within a Bayesian hierarchical framework. There is a function for the model 
#' @author Michelle 
#' @param model.name can be one or a list from the following options "Null", "Recruitment", "Ambient", "Warming", "Removal", "Removal_+_Warming", "Removal_x_Warming" required for model.formula function 
#' @param df data frame being analyzed
#' @param formula out put of model.formula function
#' @return model.formula() returns one or a list of brms formulas
#' @return model.prior() returns priors for specified models
#' @return model.fit.beta() returns brms object of model fits


# model fits are performed with brms
require(brms)

# allow for parallelization of MCMC chains
options(mc.cores = parallel::detectCores())

# the non-linear population-dynamics model depends on what treatment effects we include
model.formula <- function(model.name){
  family <- zero_inflated_beta(link="identity",link_phi="log",link_zi="logit")

  # future cover mirrors past cover
  if(model.name == "Null"){
    formula <- bf(
      cover ~ prev_cover,
      family = family
    )
  }

  # plots experience recruitment and
  # no density-dependent growth
  if(model.name == "Recruitment"){
    formula <- bf(
      cover ~ (0.5/(1+exp(-q))) + prev_cover,
      q ~ 1 + (1 | ID | focal),
      nl=T,
      family = family
    )
  }

  # plots experience recruitment and
  # density-dependent growth does not vary by treatment
  if(model.name == "Ambient"){
    formula <- bf(
      cover ~ (0.5/(1+exp(-q))) + prev_cover*exp(g),
      q ~ 1 + (1 | ID | focal),
      g ~ 1 + (1 | ID | focal),
      nl = T,
      family = family
    )
  }

  # plots experience recruitment and
  # density-dependent growth varies with removal treatment
  if(model.name == "Removal"){
    formula <- bf(
      cover ~ (0.5/(1+exp(-q))) + prev_cover*exp(g),
      q ~ 1 + (1 | ID | focal),
      g ~ 1 + removal + (1 + removal | ID | focal),
      nl = T,
      family = family
    )
  }

  # plots experience recruitment and
  # density-dependent growth varies with warming treatment
  if(model.name == "Warming"){
    formula <- bf(
      cover ~ (0.5/(1+exp(-q))) + prev_cover*exp(g),
      q ~ 1 + (1 | ID | focal),
      g ~ 1 + warming + (1 + warming | ID | focal),
      nl = T,
      family = family
    )
  }

  # plots experience recruitment and
  # density-dependent growth varies with removal and warming treatments
  # (called Removal + Warming in Table 1)
  if(model.name == "Removal_+_Warming"){
    formula <- bf(
      cover ~ (0.5/(1+exp(-q))) + prev_cover*exp(g),
      q ~ 1 + (1 | ID | focal),
      g ~ 1 + warming + removal + (1 + warming + removal | ID | focal),
      nl = T,
      family = family
    )
  }

  # plots experience recruitment and 
  # density-dependent growth varies with removal and warming treatments and their interaction
  if(model.name == "Removal_x_Warming"){ 
    formula <- bf(
      cover ~ (0.5/(1+exp(-q))) + prev_cover*exp(g),
      q ~ 1 + (1 | ID | focal),
      g ~ 1 + warming + removal + removal:warming + (1 + warming + removal + removal:warming | ID | focal),
      nl = T,
      family = family
    )
  }

  return(formula)
}

# function that specifies the prior for various population effects
model.prior <- function(df, formula){
  prior <- get_prior(formula, data = df)
  # set b prior just once
  # this prior is broad though the model specification restricts recruitment to be [0,0.5]
  # as this becomes exp(G) we restrict the range to avoid unbiological predictions (e.g., growth from low percent cover to percent cover > 1)
  prior$prior[which(prior$class=="b")] <- "normal(0,5)"


  # set all sd priors
  prior$prior[which(prior$class=="sd")] <- "student_t(3, 0, 10)"
  # set all corr priors
  prior$prior[which(prior$class=="cor")] <- "lkj(1)"
  return(prior)
}

# given a dataset (df) and a model name (model.name), fit the appropriate model with brms
# the option sample_prior allows us to generate predictions from the prior only to check their biological realism
model.fit.beta <- function(df, model.name, sample_prior = "no"){
  formula <- model.formula(model.name)
  prior <- model.prior(df, formula)
  
  # fits the appropriate model with brms and options as described in
  mf <- brm(
    formula = formula,
    prior = prior,
    data = df,
    iter = 4000,
    warmup = 1000,
    chains = 2,
    cores = 2,
    control = list(adapt_delta=0.90, max_treedepth = 20),
    seed = 12,
    sample_prior = sample_prior
    )
}


