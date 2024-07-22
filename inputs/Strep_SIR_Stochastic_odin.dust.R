
library(tidyverse)
library(odin.dust)
gen_sir <- odin.dust::odin_dust("inputs/sir_stochastic.R")

# Running the SIR model with dust (parameters consisting of value, lo_CI, hi_CI)
results <- read.csv("outputs/main/seasonality_waning[-10, 0]_nice_final_with_vacc/tune_initial_with_CI.csv")
pars <- list(log_A_ini = results[1,2],
             time_shift_1 = results[2,2],
             time_shift_2 = results[3,2],
             beta_0 = results[4,2],
             beta_1 = results[5,2],
             beta_2 = results[6,2],
             scaled_wane = results[7,2],
             log_delta = results[8,2],
             # Other fixed values:
             max_wane = (0),
             min_wane = (-10),
             vacc = (0.9*0.862*0.02),
             sigma_1 = (1/15.75),
             sigma_2 = (1)
) # Serotype 1 is categorised to have the lowest carriage duration

pars_lo_CI <- list(log_A_ini = results[1,3],
                   time_shift_1 = results[2,3],
                   time_shift_2 = results[3,3],
                   beta_0 = results[4,3],
                   beta_1 = results[5,3],
                   beta_2 = results[6,3],
                   scaled_wane = results[7,3],
                   log_delta = results[8,3]
)

pars_hi_CI <- list(log_A_ini = results[1,4],
                   time_shift_1 = results[2,4],
                   time_shift_2 = results[3,4],
                   beta_0 = results[4,4],
                   beta_1 = results[5,4],
                   beta_2 = results[6,4],
                   scaled_wane = results[7,4],
                   log_delta = results[8,4]
)

n_times <- 4745 # 4745 or similar to the number of date range (of the provided data), or try 500 for trial
n_particles <- 15L

sir_model <- gen_sir$new(pars = pars,
                         time = 1,
                         n_particles = n_particles,
                         n_threads = 4L,
                         seed = 1L)

# sir_model$state() # test array OR matrix state
sir_lo_CI <- gen_sir$new(pars = pars_lo_CI, time = 1,
                         n_particles = n_particles, n_threads = 4L, seed = 1L)
sir_hi_CI <- gen_sir$new(pars = pars_hi_CI, time = 1,
                         n_particles = n_particles, n_threads = 4L, seed = 1L)


model <- array(NA, dim = c(sir_model$info()$len, n_particles, n_times))
# Save arrays for low CI and high CI separately:
lo_CI <- array(NA, dim = c(sir_lo_CI$info()$len, n_particles, n_times))
hi_CI <- array(NA, dim = c(sir_hi_CI$info()$len, n_particles, n_times))

# Beta check
time <- seq(1, n_times, 1)
# beta_temporary <- beta_0*((1+beta_1*cos(2*pi*((time_shift_1*365)+time)/365)) + (1+beta_2*sin(2*pi*((time_shift_2*365)+time)/365)))
beta_temporary <- pars$beta_0*((1+pars$beta_1*cos(2*pi*((pars$time_shift_1*365)+time)/365)) + (1+pars$beta_2*sin(2*pi*((pars$time_shift_2*365)+time)/365)))
beta <- ifelse(time >= 2648, beta_temporary*(1-(0.9*0.862*0.02)), beta_temporary)
print(c(max(beta), max(beta_temporary))) # beta_temporary simulated with no infant vaccination
print(c(min(beta), min(beta_temporary)))

# R0 estimation (R0 changes due to seasonality)
# R0_no_vacc <- beta_temporary/(10^(pars$log_delta)+pars$sigma_1+pars$))
# R0_vacc <- (beta/(pars$log_delta+pars$sigma_1)) +  ((pars$log_delta)*(beta)) / ((pars$log_delta + 192/(4064*4745))*(pars$sigma_2 + 192/(4064*4745)))
# max(R0)
# min(R0)

for (t in seq_len(n_times)) {
  model[ , , t] <- sir_model$run(t)
  lo_CI[ , , t] <- sir_model$run(t)
  hi_CI[ , , t] <- sir_model$run(t)
}

time <- model[1, 1, ] # because in the position of [1, 1, ] is time
model <- model[-1, , ] # compile all matrix into 1 huge df, delete time (position [-1, , ])
time <- lo_CI[1, 1, ] # because in the position of [1, 1, ] is time
lo_CI <- lo_CI[-1, , ] # compile all matrix into 1 huge df, delete time (position [-1, , ])
time <- hi_CI[1, 1, ] # because in the position of [1, 1, ] is time
hi_CI <- hi_CI[-1, , ] # compile all matrix into 1 huge df, delete time (position [-1, , ])


## 1. Data Load ################################################################
incidence <- read.csv("inputs/incidence.csv")
daily_incidence_data <- incidence %>% 
  dplyr::rename(time = day,
                value_data = cases) %>% 
  dplyr::mutate(compartment = "Diseased",
                replicate = 1)

daily_incidence_modelled <- 
  reshape2::melt(model) %>% 
  dplyr::rename(index = Var1,     # Var1 = dimension that stored SADR values
                replicate = Var2, # Var2 = particles
                time = Var3       # Var3 = time
  ) %>% 
  dplyr::filter(index < 5) %>% 
  dplyr::mutate(compartment = 
                  dplyr::case_when(index == 1 ~ "Asymptomatic",
                                   index == 2 ~ "Diseased",
                                   index == 3 ~ "Susceptible",
                                   index == 4 ~ "Recovered"
                  )) %>% 
  dplyr::rename(value_model = value) %>% 
  dplyr::select(-index)

daily_low_CI <- 
  reshape2::melt(lo_CI) %>% 
  dplyr::rename(index = Var1,     # Var1 = dimension that stored SADR values
                replicate = Var2, # Var2 = particles
                time = Var3       # Var3 = time
  ) %>% 
  dplyr::filter(index < 5) %>% 
  dplyr::mutate(compartment = 
                  dplyr::case_when(index == 1 ~ "Asymptomatic",
                                   index == 2 ~ "Diseased",
                                   index == 3 ~ "Susceptible",
                                   index == 4 ~ "Recovered"
                  )) %>% 
  dplyr::rename(value_low_CI = value) %>% 
  dplyr::select(-index)

daily_high_CI <- 
  reshape2::melt(hi_CI) %>% 
  dplyr::rename(index = Var1,     # Var1 = dimension that stored SADR values
                replicate = Var2, # Var2 = particles
                time = Var3       # Var3 = time
  ) %>% 
  dplyr::filter(index < 5) %>% 
  dplyr::mutate(compartment = 
                  dplyr::case_when(index == 1 ~ "Asymptomatic",
                                   index == 2 ~ "Diseased",
                                   index == 3 ~ "Susceptible",
                                   index == 4 ~ "Recovered"
                  )) %>% 
  dplyr::rename(value_high_CI = value) %>% 
  dplyr::select(-index)

daily_joined <- daily_incidence_modelled %>% 
  dplyr::left_join(daily_low_CI, by =  c("time", "compartment", "replicate")) %>% 
  dplyr::left_join(daily_high_CI, by =  c("time", "compartment", "replicate")) %>% 
  dplyr::left_join(daily_incidence_data, by =  c("time", "compartment", "replicate"))

write.csv(daily_joined, "inputs/daily_joined.csv", row.names = F)

Sys.sleep(10) # wait 10 secs before generate pictures

## 2. Generate pics! ###########################################################
daily_joined <- read.csv("inputs/daily_joined.csv")

# Viz daily simulation
vaccine_UK <- data.frame(
  day = c(1343, 2648),
  week = c(1343/7, 2648/7),
  y_crd = c(log1p(4e7),log1p(4e7)),
  vaccine = c("PCV7", "PCV13")) # PCV7 start 4 September 2006, PCV13 from April 2010
# SOURCE: https://www.gov.uk/government/publications/pneumococcal-disease-caused-by-strains-in-prevenar-13-and-not-in-prevenar-7-vaccine/pneumococcal-disease-infections-caused-by-serotypes-in-prevenar-13-and-not-in-prevenar-7

col_compartment <- c(Susceptible = "#8c8cd9",
                     Asymptomatic = "darkred",
                     Diseased = "#cc0099",
                     Recovered = "#999966",
                     Data_dis = "deepskyblue3"
)

png("pictures/data_plus_model_daily.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(daily_joined, aes(x = time, y = value_model,
                         group = interaction(compartment,replicate),
                         colour = compartment)) +
  geom_line() +
  geom_line(data = daily_joined %>% 
              dplyr::filter(compartment == "Diseased",
                            !is.na(value_data)) %>% 
              dplyr::mutate(compartment = 
                              case_when(compartment == "Diseased" ~ "Data_dis")),
            aes(y = value_data, colour = compartment),
            show.legend = T) +
  geom_vline(data = vaccine_UK, aes(xintercept = day,
                                    colour = "black"),
             linetype = "dashed") +
  # geom_label(data = vaccine_UK, aes(x = day, y = y_crd, label = vaccine),
  #            colour = "black") +
  scale_color_manual(values = c(col_compartment),
                     name = "States",
                     breaks = c("Susceptible", "Asymptomatic", "Diseased", "Recovered", "Data_dis"),
                     labels = c("Susceptible", "Asymptomatic", "Diseased", "Recovered", "Diseased Data")
  ) +
  scale_y_continuous(trans = "log1p") +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in x-axis
  ggtitle("Serotype 1 Cases (Aggregated by Days)") +
  xlab("Time (in Day)") +
  ylab("Number of People") +
  theme_bw()
dev.off()

## 1.2. Weekly incidence
daily_joined <- daily_joined %>% 
  dplyr::mutate(weekly = ceiling(time/7))

weekly_joined <- daily_joined %>% 
  dplyr::group_by(replicate, weekly, compartment) %>% 
  dplyr::summarise(value_model = sum(value_model, na.rm = T),
                   value_low_CI = sum(value_low_CI, na.rm = T),
                   value_high_CI = sum(value_high_CI, na.rm = T),
                   value_data = sum(value_data, na.rm = T),
                   .groups = "drop") %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(value_data = case_when(replicate > 1 ~ NA,
                                       .default = as.numeric(value_data)))

png("pictures/data_plus_model_weekly.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(weekly_joined, aes(x = weekly, y = value_model,
                          group = interaction(compartment,replicate),
                          colour = compartment)) +
  geom_line() +
  geom_line(data = weekly_joined %>% 
              dplyr::filter(compartment == "Diseased",
                            !is.na(value_data)) %>% 
              dplyr::mutate(compartment = 
                              case_when(compartment == "Diseased" ~ "Data_dis")),
            aes(y = value_data, colour = compartment),
            show.legend = T) +
  geom_vline(data = vaccine_UK, aes(xintercept = week,
                                    colour = "black"),
             linetype = "dashed") +
  # geom_label(data = vaccine_UK, aes(x = week, y = y_crd, label = vaccine),
  #            colour = "black") +
  scale_color_manual(values = c(col_compartment),
                     name = "States",
                     breaks = c("Susceptible", "Asymptomatic", "Diseased", "Recovered", "Data_dis"),
                     labels = c("Susceptible", "Asymptomatic", "Diseased", "Recovered", "Diseased Data")
  ) +
  scale_y_continuous(trans = "log1p") +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in x-axis
  ggtitle("Serotype 1 Cases (Aggregated by Weeks)") +
  xlab("Time (in Week)") +
  ylab("Number of People") +
  theme_bw()
dev.off()

# Weekly incidence focused only on diseased people
weekly_D_only <- weekly_joined %>% 
  dplyr::filter(compartment == "Diseased")

png("pictures/data_plus_model_weekly_diseased_only.png", width = 17, height = 12, unit = "cm", res = 1200)
col_imD_weekly <- c(cases_weekly = "deepskyblue3",
                    model_weekly = "#cc0099",
                    CIs = "yellow")
ggplot(weekly_D_only, aes(x = weekly,
                          group = interaction(compartment,replicate),
                          colour = compartment)) +
 
  geom_line(aes(y = value_low_CI, colour = "CIs")) + # --> low CI
  geom_line(aes(y = value_high_CI, colour = "CIs")) + # --> high CI
  
  geom_line(aes(y = value_model, colour = "model_weekly")) + # --> model
  geom_line(aes(y = value_data, colour = "cases_weekly")) + # --> data
  
  geom_vline(data = vaccine_UK, aes(xintercept = week,
                                    colour = "black"),
             linetype = "dashed") +
  geom_label(aes(x = 1343/7, y = 45, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7
  geom_label(aes(x = 2648/7, y = 45, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13
  scale_x_continuous() +
  # scale_y_continuous(trans = "log1p") +
  scale_color_manual(values = col_imD_weekly,
                     name = "Cases",
                     breaks = c("cases_weekly", "model_weekly"),
                     labels = c("Data", "Model")
  ) +
  ggtitle("The Comparison of Model Output and Counts of Serotype 1 in England") +
  xlab("Week") +
  ylab("Serotype 1 Cases (Aggregated by Week)") +
  theme_bw()
dev.off()

