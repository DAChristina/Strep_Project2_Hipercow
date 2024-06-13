# 3. Tuning the pMCMC ##########################################################
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
pars <- list(I_ini = 0.001, # in toy data the real value = 0.0015*S_ini (100 people)
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

priors <- prepare_priors(pars)

proposal_matrix <- as.matrix(read.csv("outputs/new_proposal_mtx.csv")) # change proposal matrix for tuning
mcmc_pars <- prepare_parameters(initial_pars = pars, priors = priors, proposal = proposal_matrix, transform = transform)



