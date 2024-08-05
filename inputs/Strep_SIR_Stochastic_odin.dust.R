
library(tidyverse)
library(odin.dust)
gen_sir <- odin.dust::odin_dust("inputs/sir_stochastic.R")

# Running the SIR model with dust (parameters consisting of value, lo_CI, hi_CI)
results <- read.csv("outputs/main/seasonality_waning[-10, -5]_nice_final_with_vacc_modiv_beta_trial2/tune_initial_with_CI.csv")
pars <- list(log_A_ini = results[1,2],
             time_shift_1 = results[2,2],
             time_shift_2 = results[3,2],
             beta_0 = results[4,2],
             beta_1 = results[5,2],
             beta_2 = results[6,2],
             scaled_wane = results[7,2],
             log_delta = results[8,2],
             # Other fixed values:
             max_wane = (-5),
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

n_times <- 4745*4 # 4745 or similar to the number of date range (of the provided data), or try 500 for trial
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
mu_0 <- 1/(80.70*365) # FIXED background mortality
mu_1 <- 192/(4064*4745) # FIXED disease-associated mortality; ratio 192/4064 in 4745 days
# beta_temporary <- beta_0*((1+beta_1*cos(2*pi*((time_shift_1*365)+time)/365)) + (1+beta_2*sin(2*pi*((time_shift_2*365)+time)/365)))
beta_temporary <- pars$beta_0*((1+pars$beta_1*cos(2*pi*((pars$time_shift_1*365)+time)/365)) + (1+pars$beta_2*sin(2*pi*((pars$time_shift_2*365)+time)/365)))
beta <- ifelse(time >= 2648, beta_temporary*(1-(0.9*0.862*0.02)), beta_temporary)
print(c(max(beta), max(beta_temporary))) # beta_temporary simulated with no infant vaccination
print(c(min(beta), min(beta_temporary)))

# R0 estimation (R0 changes due to seasonality)
R0_no_vacc <- beta_temporary/((mu_0+(10^(pars$log_delta))+pars$sigma_1)*((pars$sigma_2)+(mu_0+mu_1)))
R0_vacc <- beta/((mu_0+(10^(pars$log_delta))+pars$sigma_1)*((pars$sigma_2)+(mu_0+mu_1)))

plot(time, R0_no_vacc)
plot(time, R0_vacc)




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

write.csv(daily_joined, "inputs/daily_joined_simulated4745times4.csv", row.names = F)
write.csv(daily_incidence_modelled, "raw_data/daily_joined_longer_simulated4745times4.csv", row.names = F)

Sys.sleep(10) # wait 10 secs before generate pictures

## 2. Generate pics! ###########################################################
daily_joined <- read.csv("inputs/daily_joined.csv") %>% # save data for 4745 days only
  dplyr::mutate(date = as.Date("2003-01-01") + days(time - 1))
daily_incidence_modelled <- read.csv("raw_data/daily_joined_longer_simulated4745times4.csv") %>% # save data for longer days
  dplyr::mutate(date = as.Date("2003-01-01") + days(time - 1))


# Viz daily simulation
vaccine_UK <- data.frame(
  date = c("2006-09-04", "2010-04-01"),
  day = c(1343, 2648),
  week = c(1343/7, 2648/7),
  y_crd = c(log1p(4e7),log1p(4e7)),
  vaccine = c("PCV7", "PCV13")) # PCV7 start 4 September 2006, PCV13 from April 2010
# SOURCE: https://www.gov.uk/government/publications/pneumococcal-disease-caused-by-strains-in-prevenar-13-and-not-in-prevenar-7-vaccine/pneumococcal-disease-infections-caused-by-serotypes-in-prevenar-13-and-not-in-prevenar-7

col_compartment <- c("Susceptible" = "#8c8cd9",
                     "Asymptomatic" = "orange",
                     "Diseased" = "#cc0099",
                     "Recovered" = "#999966",
                     "Data_dis" = "lightblue3",
                     # Vaccination
                     "PCV7" = "gray70",
                     "PCV13" = "gray20"
)

png("pictures/data_plus_model_daily_simulated4745times4.png", width = 24, height = 12, unit = "cm", res = 1200)
ggplot(daily_incidence_modelled, aes(x = date, y = value_model,
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
  geom_vline(data = vaccine_UK, aes(xintercept = as.Date(vaccine_UK$date),
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(col_compartment),
                     name = "States",
                     breaks = c("Susceptible", "Asymptomatic", "Diseased", "Recovered", "Data_dis", "PCV7", "PCV13"),
                     labels = c("Susceptible", "Asymptomatic", "Diseased", "Recovered", "Diseased Data", "PCV7", "PCV13")
  ) +
  scale_y_continuous(trans = "log1p") +
  scale_x_date(date_breaks = "5 years",
               # date_minor_breaks = "1 month",
               date_labels = "%Y") + # try "%b %Y"
  xlim(as.Date(c('1/1/2003', '12/12/2030'), format="%d/%m/%Y")) + # edit this for data fitting: 20/10/2015
  ggtitle("Simulation Model for Serotype 1 Cases (Aggregated by Days)") +
  xlab("Time") +
  ylab("Number of People") +
  theme_bw()
dev.off()

## 1.2. Weekly incidence
daily_joined <- daily_joined %>% 
  dplyr::mutate(weekly = ceiling(time/7))

weekly_joined <- daily_joined %>% 
  dplyr::group_by(replicate, weekly, compartment) %>% 
  dplyr::summarise(value_model = sum(value_model, na.rm = T),
                   date = max(date),
                   value_low_CI = sum(value_low_CI, na.rm = T),
                   value_high_CI = sum(value_high_CI, na.rm = T),
                   value_data = sum(value_data, na.rm = T),
                   .groups = "drop") %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(value_data = case_when(replicate > 1 ~ NA,
                                       .default = as.numeric(value_data)))

weekly_incidence_modelled <- daily_incidence_modelled %>% 
  dplyr::mutate(weekly = ceiling(time/7)) %>% 
  dplyr::group_by(replicate, weekly, compartment) %>% 
  dplyr::summarise(value_model = sum(value_model, na.rm = T),
                   date = max(date),
                   # value_low_CI = sum(value_low_CI, na.rm = T),
                   # value_high_CI = sum(value_high_CI, na.rm = T),
                   # value_data = sum(value_data, na.rm = T),
                   .groups = "drop") %>% 
  dplyr::ungroup() #%>% 
  # dplyr::mutate(value_data = case_when(replicate > 1 ~ NA,
  #                                      .default = as.numeric(value_data)))

png("pictures/data_plus_model_weekly_simulated4745times4.png", width = 24, height = 12, unit = "cm", res = 1200)
ggplot(weekly_incidence_modelled, aes(x = date, y = value_model,
                          group = interaction(compartment,replicate),
                          colour = compartment)) +
  geom_line() +
  geom_line(data = weekly_joined %>% 
              dplyr::filter(compartment == "Diseased",
                            !is.na(value_data)) %>% 
              dplyr::mutate(compartment = 
                              case_when(compartment == "Diseased" ~ "Data_dis")),
            aes(y = value_data, colour = compartment),
            size = 1.1,
            show.legend = T) +
  geom_vline(data = vaccine_UK, aes(xintercept = as.Date(vaccine_UK$date),
                                    colour = vaccine),
             linetype = "dashed") +
  # geom_label(data = vaccine_UK, aes(x = week, y = y_crd, label = vaccine),
  #            colour = "black") +
  scale_color_manual(values = c(col_compartment),
                     name = "States",
                     breaks = c("Susceptible", "Asymptomatic", "Diseased", "Recovered", "Data_dis", "PCV7", "PCV13"),
                     labels = c("Susceptible", "Asymptomatic", "Diseased", "Recovered", "Diseased Data", "PCV7", "PCV13")
  ) +
  scale_y_continuous(trans = "log1p") +
  scale_x_date(date_breaks = "5 years",
               # date_minor_breaks = "1 month",
               date_labels = "%Y") + # try "%b %Y"
  xlim(as.Date(c('1/1/2003', '12/12/2030'), format="%d/%m/%Y")) + # edit this for data fitting: 20/10/2015
  ggtitle("Simulation Model for Serotype 1 Cases (Aggregated by Days)") +
  ggtitle("Serotype 1 Cases (Aggregated by Weeks)") +
  xlab("Time (in Week)") +
  ylab("Number of People") +
  theme_bw()
dev.off()

# Weekly incidence focused only on diseased people
weekly_D_only <- weekly_joined %>% 
  dplyr::filter(compartment == "Diseased")

weekly_modelled_longer <- weekly_incidence_modelled %>% 
  dplyr::filter(compartment == "Diseased")

png("pictures/data_plus_model_weekly_diseased_only_simulated4745times4.png", width = 24, height = 12, unit = "cm", res = 1200)
col_imD_weekly <- c(cases_weekly = "lightblue3",
                    model_weekly = "#cc0099",
                    CIs = "yellow",
			  # Vaccination
                     "PCV7" = "gray70",
                     "PCV13" = "gray20"
)
ggplot(weekly_modelled_longer, aes(x = date,
                          group = interaction(compartment,replicate),
                          colour = compartment)) +
 
  # geom_line(aes(y = value_low_CI, colour = "CIs")) + # --> low CI
  # geom_line(aes(y = value_high_CI, colour = "CIs")) + # --> high CI
  
  geom_line(aes(y = value_model, colour = "model_weekly")) + # --> model
  # geom_line(aes(y = value_data, colour = "cases_weekly")) + # --> data
  geom_line(data = weekly_joined %>% 
              dplyr::filter(compartment == "Diseased",
                            !is.na(value_data)) %>% 
              dplyr::mutate(compartment = 
                              case_when(compartment == "Diseased" ~ "cases_weekly")),
            aes(y = value_data, colour = compartment),
            show.legend = T,
            size = 1.1) +
  geom_vline(data = vaccine_UK, aes(xintercept = as.Date(vaccine_UK$date),
                                    colour = vaccine),
             linetype = "dashed") +
  scale_x_date(date_breaks = "5 years",
               # date_minor_breaks = "1 month",
               date_labels = "%Y") + # try "%b %Y"
  xlim(as.Date(c('1/1/2003', '12/12/2015'), format="%d/%m/%Y")) + # edit this for data fitting: 20/10/2015
  scale_color_manual(values = col_imD_weekly,
                     name = "Cases",
                     breaks = c("cases_weekly", "model_weekly", "PCV7", "PCV13"),
                     labels = c("Data", "Model", "PCV7", "PCV13")
  ) +
  ggtitle("The Comparison of Model Output and Daily Epidemiological Data of Serotype 1 in England") +
  xlab("Year") +
  ylab("Serotype 1 Cases (Aggregated by Week)") +
  theme_bw()
dev.off()

