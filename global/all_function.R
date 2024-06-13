
case_compare <- function(state, observed, pars = NULL) {
  exp_noise <- 1e4
  
  # incidence based on model's "n_SI_daily" from gen_sir$new(pars = list(), time = 0, n_particles = 1L)$info()
  incidence_modelled <- state[5, , drop = TRUE] # n_SI_daily is located at state[5, , ]
  
  # incidence based on data
  incidence_observed <- observed$cases # daily new cases
  
  if (is.na(observed$cases)) {
    loglik_cases <- numeric(n)
    
  } else {
    n <- ncol(state)
    lamb <- incidence_modelled + rexp(n, exp_noise)
    loglik_cases <- dpois(x = incidence_observed, lambda = lamb)
    
  }
  return(loglik_cases)
}

# That transform function
# https://github.com/mrc-ide/mcstate/blob/da9f79e4b5dd421fd2e26b8b3d55c78735a29c27/tests/testthat/test-if2.R#L40
# https://github.com/mrc-ide/mcstate/issues/184
parameter_transform <- function(pars) {
  I_ini <- pars[["I_ini"]]
  just_beta <- pars[["just_beta"]]
  just_sigma <- pars[["just_sigma"]]
  
  list(I_ini = I_ini,
       just_beta = just_beta,
       just_sigma = just_sigma)
  
}

transform <- function(pars) {
  parameter_transform(pars)
}

prepare_parameters <- function(initial_pars, priors, proposal, transform) {
  
  mcmc_pars <- mcstate::pmcmc_parameters$new(
    list(mcstate::pmcmc_parameter("I_ini", 0.1, min = 0, max = 0.5,
                                  prior = function(s) log(1e-10)),
         mcstate::pmcmc_parameter("just_beta", 0.5, min = 0, max = 0.8,
                                  prior = priors$just_beta),
         mcstate::pmcmc_parameter("just_sigma", 0.01, min = 0, max = 1,
                                  prior = priors$just_sigma)
    ),
    proposal = proposal,
    transform = transform)
  
}

prepare_priors <- function(pars) {
  priors <- list()
  
  # priors$I_ini <- function(s) {
  #   dunif(s, min = 0, max = 0.5, log = TRUE)
  # } # assume I_ini draws from uniform distribution
  priors$just_beta <- function(s) {
    dgamma(s, shape = 1, scale = 0.1, log = TRUE)
  }
  priors$just_sigma <- function(s) {
    dgamma(s, shape = 1, scale = 1, log = TRUE)
  }
  priors
}

pmcmc_further_process <- function(n_steps, pmcmc_result) {
  processed_chains <- mcstate::pmcmc_thin(pmcmc_result, burnin = n_steps/2, thin = 2)
  parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
  parameter_mean_hpd
  
  mcmc1 <- coda::as.mcmc(cbind(pmcmc_result$probabilities, pmcmc_result$pars))
  mcmc1
}

ess_calculation <- function(mcmc1){
  calc <- list(ess = coda::effectiveSize(mcmc1),
               acceptance_rate = 1 - coda::rejectionRate(mcmc1))
  calc
}

pmcmc_trace <- function(mcmc1) {
  plot(mcmc1) # to save the figures into pdf
  # png("pictures/mcmc1.png", res = 1200)
  # plot(mcmc1)
  # dev.off()
}
