# To change the model into stochastic fuction, use the probabilistic function:
## Individual probabilities of transition:
## Definition of the time-step and output as "time"
freq <- user(1)
dt <- 1/freq
initial(time) <- 0

# 1. PARAMETERS ################################################################
N <- user(6.7e7) # FIXED England's pop size is roughly 67,000,000
log_A_ini <- user(0) # S_ini*10^(log10(-5.69897)) = 120 people; change A_ini into log10(A_ini)
D_ini <- user(0) 
time_shift_1 <- user(0)
time_shift_2 <- user(0)
beta_0 <- user(0)
beta_1 <- user(0)
beta_2 <- user(0)

max_wane <- user(-5) # FIXED, scaled waning immunity
min_wane <- user(-10) # FIXED, scaled waning immunity
scaled_wane <- user(0)

# Vaccination:
# https://webarchive.nationalarchives.gov.uk/ukgwa/20211105111851mp_/https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/540290/hpr2416_ppv.pdf
# https://fingertips.phe.org.uk/search/PPV#page/4/gid/1/pat/159/par/K02000001/ati/15/are/E92000001/iid/30313/age/27/sex/4/cat/-1/ctp/-1/yrr/1/cid/4/tbm/1
# vacc_elderly <- 0.7*0.57 # FIXED PPV23 vaccination coverage * efficacy
# ratio of vaccinated elderly for >64 y.o. people, averaged 69.7243% ~ 70%
vacc <- 0.9*0.862*0.02 # FIXED PCV13 vaccination coverage * efficacy * proportion of kids below 2 y.o.
# ratio of vaccinated kids, averaged 90%
# vacc <- (vacc_elderly + vacc_kids)/2 # FIXED, average

# Country calibration:
# Children: 1.07638532472038 (it is called delta in the spreadsheet)
# Adults: 0.536936186788821 (basically gamma)
# Average: 0.8066608
UK_calibration <- user(0.8066608) # FIXED (Lochen et al., 2022)

log_delta <- user(0) # required in mcState
sigma_1 <- user(1/15.75) # FIXED per-day, carriage duration (95% CI 7.88-31.49) (Serotype 1) (Chaguza et al., 2021)
sigma_2 <- user(0)
mu_0 <- 1/(80.70*365) # background mortality FIXED based on the inverse of life expectancy
mu_1 <- user(192/(4064*4745)) # FIXED disease-associated mortality; ratio 192/4064 in 4745 days
pi <- user(3.141593) # FIXED

# 2. INITIAL VALUES ############################################################
A_ini <- 10^(log_A_ini)*N
initial(A) <- A_ini
initial(D) <- D_ini
initial(S) <- N - A_ini - D_ini
initial(R) <- 0
initial(n_AD_daily) <- 0
initial(n_AD_cumul) <- 0

# 3. UPDATES ###################################################################
beta_temporary <- beta_0*((1+beta_1*cos(2*pi*((time_shift_1*365)+time)/365)) + (1+beta_2*sin(2*pi*((time_shift_2*365)+time)/365)))
# Infant vaccination coverage occurs when PCV13 introduced in April 2010 (day 2648 from 01.01.2003)
# https://fingertips.phe.org.uk/search/vaccination#page/4/gid/1/pat/159/par/K02000001/ati/15/are/E92000001/iid/30306/age/30/sex/4/cat/-1/ctp/-1/yrr/1/cid/4/tbm/1/page-options/tre-do-0
# https://cran.r-project.org/web/packages/finalsize/vignettes/varying_contacts.html
beta <- if (time >= 2648) beta_temporary*(1-vacc) else beta_temporary

# lambda <- beta*(A+D)/N # infectious state from Asymtomatic & Diseased individuals
lambda <- if ((A+D) > 0) beta*(A+D)/N else 0
delta <- (10^(log_delta))*UK_calibration

log_wane <- scaled_wane*(max_wane-min_wane)+min_wane # scaled_wane*(max_wane−min_wane)+min_wane; rescaled using (wane-wane_min)/(wane_max-wane_min)
wane <- 10^(log_wane)
# wane <- 0

# Individual probabilities of transition
p_SA <- 1- exp(-(lambda+mu_0) * dt)
p_Asym <- 1- exp(-(delta+mu_0+sigma_1) * dt)
p_Dis <- 1- exp(-(sigma_2+mu_0+mu_1) * dt)
p_RS <- 1- exp(-(wane+mu_0) * dt)

# Draws for numbers changing between compartments
# Leaving S
n_Suscep <- rbinom(S, p_SA)
n_SA <- rbinom(n_Suscep, lambda/(lambda+mu_0))
n_S_dead <- n_Suscep - n_SA

# Leaving A
n_Asym <- rbinom(A, p_Asym) # n_Asym <- n_AD + n_AR cause cyclic dependency error
n_AD <- rbinom(n_Asym, delta/(delta+mu_0+sigma_1))
n_AR <- rbinom((n_Asym - n_AD), sigma_1/(mu_0+sigma_1))
n_A_dead <- n_Asym - n_AD - n_AR

# Leaving D
n_Dis <- rbinom(D, p_Dis)
n_DR <- rbinom(n_Dis, sigma_2/(sigma_2+mu_0+mu_1))
n_Dd <- rbinom((n_Dis - n_DR), mu_1/(mu_0+mu_1))
n_D_dead <- n_Dis - n_DR - n_Dd

# Leaving R
n_Resist <- rbinom(R, p_RS)
n_RS <- rbinom(n_Resist, wane)
n_R_dead <- n_Resist - n_RS

# Closed system: births = deaths; all born susceptible
n_S_born <- n_S_dead + n_Dd + n_D_dead + n_A_dead + n_R_dead

# The transitions
update(time) <- (step + 1) * dt
update(S) <- S - n_SA + n_RS + n_S_born - n_S_dead
update(A) <- A + n_SA - (n_AD + n_AR) - n_A_dead
update(D) <- D + n_AD - (n_DR + n_Dd) - n_D_dead
update(R) <- R + n_AR + n_DR - n_RS - n_R_dead
# that "little trick" previously explained in https://github.com/mrc-ide/dust/blob/master/src/sir.cpp for cumulative incidence:
# based on tutorial: https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/the-model
update(n_AD_daily) <- if (step %% freq == 0) n_AD else n_AD_daily + n_AD
update(n_AD_cumul) <- n_AD_cumul + n_AD # no interest in asymptomatic cases that've recovered