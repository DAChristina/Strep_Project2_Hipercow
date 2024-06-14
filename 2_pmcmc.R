# 2. Data Fitting ##############################################################
library(mcstate)
library(coda)
library(odin.dust)
library(dust)

source("global/all_function.R") # Collected functions stored here!

# The anatomy of an mcstate particle filter, as noted above, consists of three main components: \n 
# 1. A set of observations to fit the model to, generated using mcstate::particle_filter_data(). \n 
# 2. A model to fit, which must be a dust generator, either dust::dust() or odin.dust::odin_dust(). \n 
# 3. A comparison function, which is an R function which calculates the likelihood of the state given the data at one time point.

# There is a calibration function in mcstate to fit our model to data.
# https://mrc-ide.github.io/mcstate/articles/sir_models.html

# To make my life easier I compile the Serotype 1 cases into a new object called sir_data
# data is fed as an input to mcstate::particle_filter_data
incidence <- read.csv("inputs/incidence.csv")

dt <- 1 # rate must be an integer; 0.25 to make it 4 days, I make it 1
sir_data <- mcstate::particle_filter_data(data = incidence,
                                          time = "day",
                                          rate = 1 / dt,
                                          initial_time = 0) # Initial time makes t0 start from 0 (not 1)

# Annotate the data so that it is suitable for the particle filter to use
rmarkdown::paged_table(sir_data)


## 2a. Model Load ##############################################################
# The model below is stochastic, closed system SADR model that I have created before
# I updated the code, filled the parameters with numbers;
# e.g.dt <- user(0) because if dt <- user() generates error during MCMC run
gen_sir <- odin.dust::odin_dust("inputs/sir_stochastic.R")

# This is part of sir odin model:
pars <- list(I_ini = (-3), # in toy data the real value = 0.0015*S_ini (100 people)
             just_beta = 0.5, # in toy data the real value = 0.3
             just_sigma = 0.05 # in toy data the real value = 0.01
)

# https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/the-model-over-time
n_particles <- 50 # Trial n_particles = 50
filter <- mcstate::particle_filter$new(data = sir_data,
                                       model = gen_sir, # Use odin.dust input
                                       n_particles = n_particles,
                                       compare = case_compare,
                                       seed = 1L)

filter$run(pars)

# Variance and particles estimation (as suggested by Rich)
# parallel::detectCores() # I have 4 cores
# x <- replicate(30, filter$run(pars))
# x
# 
# var(x)
# # [1] 0.03251836
# 500 / 0.03251836 # Trial 500 particles for 0.03251836 var; how many particles are needed?
# # [1] 15375.93


priors <- prepare_priors(pars)
proposal_matrix <- diag(1, 3)
# rownames(proposal_matrix) <- c("just_beta", "just_sigma")
# colnames(proposal_matrix) <- c("just_beta", "just_sigma")

mcmc_pars <- prepare_parameters(initial_pars = pars, priors = priors, proposal = proposal_matrix, transform = transform)

# n_steps <- 100 #1e6

# Use deterministic model by add filter_deterministic
# https://mrc-ide.github.io/mcstate/articles/deterministic.html
# Index function is optional when only a small number of states are used in comparison function.
filter_deterministic <- mcstate::particle_deterministic$new(data = sir_data,
                                                            model = gen_sir,
                                                            compare = case_compare
                                                            # index = index
)

# I change pmcmc_run into a function that involve control inside:
# pmcmc_run <- mcstate::pmcmc(mcmc_pars, filter_deterministic, control = control)

# Directory for saving the outputs
dir.create("outputs", FALSE, TRUE)

pmcmc_run <- function(n_particles, n_steps){
  filter <- mcstate::particle_filter$new(data = sir_data,
                                         model = gen_sir, # Use odin.dust input
                                         n_particles = n_particles,
                                         compare = case_compare,
                                         seed = 1L)
  
  control <- mcstate::pmcmc_control(n_steps = n_steps,
                                    rerun_every = 50,
                                    rerun_random = TRUE,
                                    progress = TRUE)
  
  # The pmcmc
  pmcmc_result <- mcstate::pmcmc(mcmc_pars, filter_deterministic, control = control)
  pmcmc_result
  saveRDS(pmcmc_result, "outputs/pmcmc_result.rds")
  
  new_proposal_mtx <- cov(pmcmc_result$pars)
  write.csv(new_proposal_mtx, "outputs/new_proposal_mtx.csv", row.names = TRUE)
  
  lpost_max <- which.max(pmcmc_result$probabilities[, "log_posterior"])
  write.csv(as.list(pmcmc_result$pars[lpost_max, ]),
            "outputs/initial.csv", row.names = FALSE)
  
  # Further processing for thinning chains
  mcmc1 <- pmcmc_further_process(n_steps, pmcmc_result)
  write.csv(mcmc1, "outputs/mcmc1.csv", row.names = TRUE)
  
  # Calculating ESS & Acceptance Rate
  calc_ess <- ess_calculation(mcmc1)
  write.csv(calc_ess, "outputs/calc_ess.csv", row.names = TRUE)
  
  # Figures! (still failed, margin error)
  fig <- pmcmc_trace(mcmc1)
  
}

# Test
# pmcmc_run_result <- pmcmc_run(10, 100)

# 3. Tuning the pMCMC ##########################################################
pmcmc_tuning <- function(n_particles, n_steps){
  # New proposal matrix
  new_proposal_matrix <- as.matrix(read.csv("outputs/new_proposal_mtx.csv"))
  new_proposal_matrix <- new_proposal_matrix[, -1]
  new_proposal_matrix <- apply(new_proposal_matrix, 2, as.numeric)
  new_proposal_matrix <- (new_proposal_matrix + t(new_proposal_matrix)) / 2
  colnames(new_proposal_matrix) <- c("I_ini", "just_beta", "just_sigma")
  rownames(new_proposal_matrix) <- c("I_ini", "just_beta", "just_sigma")
  # isSymmetric(new_proposal_matrix)
  
  tune_mcmc_pars <- prepare_parameters(initial_pars = pars, priors = priors, proposal = new_proposal_matrix, transform = transform)
  
  tune_control <- mcstate::pmcmc_control(n_steps = n_steps,
                                    n_chains = 4,
                                    rerun_every = 50,
                                    rerun_random = TRUE,
                                    progress = TRUE)
  
  # The pmcmc
  tune_pmcmc_result <- mcstate::pmcmc(tune_mcmc_pars, filter_deterministic, control = tune_control)
  tune_pmcmc_result
  saveRDS(tune_pmcmc_result, "outputs/tune_pmcmc_result.rds")
  
  # new_proposal_mtx <- cov(pmcmc_result$pars)
  # write.csv(new_proposal_mtx, "outputs/new_proposal_mtx.csv", row.names = TRUE)
  
  tune_lpost_max <- which.max(tune_pmcmc_result$probabilities[, "log_posterior"])
  write.csv(as.list(tune_pmcmc_result$pars[, ]),#lpost_max, ]),
            "outputs/tune_initial.csv", row.names = FALSE)
  
  # Further processing for thinning chains
  mcmc2 <- tuning_pmcmc_further_process(n_steps, tune_pmcmc_result)
  write.csv(mcmc2, "outputs/mcmc2.csv", row.names = TRUE)
  
  # Calculating ESS & Acceptance Rate
  tune_calc_ess <- ess_calculation(mcmc2)
  write.csv(tune_calc_ess, "outputs/tune_calc_ess.csv", row.names = TRUE)
  
  # Figures! (still failed, margin error)
  fig <- pmcmc_trace(mcmc2)
  
}

# Test
# pmcmc_run_result <- pmcmc_run(10, 100)
# tuning_run_result <- pmcmc_tuning(10, 100)
