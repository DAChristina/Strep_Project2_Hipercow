# A simple SIR model
freq <- user(1)
dt <- 1/freq
initial(time) <- 0

# 1. PARAMETERS ################################################################
S_ini <- user(1e5) # required in mcState
I_ini <- user(0.0015) # required in mcState
just_beta <- user(0.5)
just_sigma <- user(0.01) # required in mcState

# 2. INITIAL VALUES ############################################################
initial(S) <- S_ini
initial(I) <- I_ini*S_ini
initial(R) <- 0
initial(n_SI_daily) <- 0
initial(n_SI_cumul) <- 0

# 3. UPDATES ###################################################################
N <- S + I + R
just_lambda <- just_beta*I/N

# Individual probabilities of transition
p_SI <- 1- exp(-just_lambda * dt)
p_IR <- 1- exp(-just_sigma * dt)

# Draws for numbers changing between compartments
n_SI <- rbinom(S, p_SI)
n_IR <- rbinom(I, p_IR)

# The transitions
update(time) <- (step + 1) * dt
update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR
# that "little trick" previously explained in https://github.com/mrc-ide/dust/blob/master/src/sir.cpp for cumulative incidence:
# based on tutorial: https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/the-model
update(n_SI_daily) <- if (step %% freq == 0) n_SI else n_SI_daily + n_SI
update(n_SI_cumul) <- n_SI_cumul + n_SI