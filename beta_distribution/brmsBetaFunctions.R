####
# functions to fit population dynamics models within a bayesian hierarchical framework
####

# the non-linear population-dynamics model depends on what treatment effects we include
model.formula <- function(model.name){
  family <- zero_inflated_beta(link="identity",link_phi="log",link_zi="logit")

  # plots experience recruitment
  # no density-dependent growth
  if(model.name == "Null"){
    formula <- bf(
      cover ~ (0.5/(1+exp(-q))) + prev_cover,
      q ~ 1 + (1 | ID | focal),
      nl=T,
      family = family
    )
  }

  # plots experience recruitment
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

  # plots experience recruitment
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

  # plots experience recruitment
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

  # plots experience recruitment
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

  # plots experience recruitment
  # density-dependent growth varies with removal and warming treatments and their interaction
  # (called Removal x Warming in Table 1)
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

  # set the prior for recruitment so that we can sample from the prior
  # this prior is broad though the model specification restricts recruitment to be [0,0.5]
  prior$prior[which(prior$class=="b" & prior$coef=="" & prior$group=="" & prior$nlpar=="q" )] <- "normal(0,5)"

  # set the prior for the population effects effecting growth G
  # as this becomes exp(G) we restrict the range to avoid unbiological predictions (e.g., growth from low percent cover to percent cover > 1)
  prior$prior[which(prior$class=="b" & prior$coef=="" & prior$group=="" & prior$nlpar=="g" )] <- "normal(0,0.5)"

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
    control = list(adapt_delta=0.90),
    seed = 12,
    sample_prior = sample_prior
  )
  
  return(mf)
}
