
library(odin.dust)
gen_sir <- odin.dust::odin_dust("inputs/sir_stochastic.R")

# Running the SIR model with dust
pars <- list(I_ini = 0.0015,
             just_beta = 0.3,
             just_sigma = 0.01
)

sir_model <- gen_sir$new(pars = pars,
                         time = 1,
                         n_particles = 15L,
                         n_threads = 4L,
                         seed = 1L)

# sir_model$state() # test array OR matrix state

# update_state is required "every single time" to run & produce matrix output (don't know why)
sir_model$update_state(pars = pars,
                       time = 0) # make sure time is 0

# all_date <- incidence$day
# all_date <- data.frame(col = integer(4745))
n_times <- 200 # 500 for trial
n_particles <- 15
x <- array(NA, dim = c(sir_model$info()$len, n_particles, n_times))

# R0 estimation
R0 <- pars$just_beta/pars$just_sigma
R0

for (t in seq_len(n_times)) {
  x[ , , t] <- sir_model$run(t)
}
time <- x[1, 1, ] # because in the position of [1, 1, ] is time
x <- x[-1, , ] # compile all matrix into 1 huge df, delete time (position [-1, , ])

# Some viz
par(mfrow = c(2, 1), mar = c(5.1, 5.1, 0.5, 0.5), mgp = c(3.5, 1, 0), las = 1)
cols <- c(S = "#8c8cd9", I = "darkred", R = "#999966", n_SI_daily = "#cc0099", n_SI_cumul = "green")
matplot(time, t(x[1, , ]), type = "l",
        xlab = "Time", ylab = "Number of individuals",
        col = cols[["S"]], lty = 1, ylim = c(0,max(x[1,,])),
        main = "All model")
matlines(time, t(x[2, , ]), col = cols[["I"]], lty = 1)
matlines(time, t(x[3, , ]), col = cols[["R"]], lty = 1)
matlines(time, t(x[4, , ]), col = cols[["n_SI_daily"]], lty = 1)
matlines(time, t(x[5, , ]), col = cols[["n_SI_cumul"]], lty = 1)
legend("right", lwd = 1, col = cols, legend = names(cols), bty = "n")

matplot(time, t(x[4, , ]), type = "l",
        xlab = "Time", ylab = "Number of individuals",
        col = cols[["n_SI_daily"]], lty = 1, ylim = c(0,max(x[4,,])),
        main = "Focused on Daily Case")

legend("right", lwd = 1, col = cols, legend = names(cols), bty = "n")
max(x[5,,]) # Check max n_SI_daily

# Toy data creation ############################################################
# glimpse(x)
new_toyData <- as.data.frame(x[4, 1, 1:n_times])
colnames(new_toyData) <- "cases"
# glimpse(new_toyData)

library(tidyverse)
incidence <- tibble(day = 1:n_times) %>% 
  bind_cols(new_toyData)

write.csv(incidence, file="inputs/incidence.csv", row.names =F)
