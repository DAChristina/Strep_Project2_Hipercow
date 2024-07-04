# Data preparation
library(tidyverse)
library(readxl)
library(cowplot)

# New updated data with meningitis (25.04.2024)
# All df are stored in raw_data
dat <- read_excel("raw_data/serotype1_UKHSA_imperial_date_age_region_MOLIS_withdeath_meningitis_clean.xlsx") #%>% 
# glimpse()

dat <- dat %>% 
  dplyr::rename(Earliest.specimen.date = Earliestspecimendate,
         current.region.name = currentregionname)

dat_G <- dat %>% 
  dplyr::mutate(AGEYR = ifelse(AGEYR >= 90, 90, as.numeric(AGEYR)), # For incidence calculation, data grouped for people aged 90+
         year = year(Earliest.specimen.date),
         month = month(Earliest.specimen.date),
         vacc = case_when(
           year < 2006 ~ "Pre-PCV7",
           year >= 2006 & year < 2011 ~ "PCV7",
           year >= 2011 ~ "PCV13",
           TRUE ~ NA_character_
         ),
         ageGroup = case_when( # edit 5 age bands
           AGEYR < 5 ~ "<5",
           AGEYR >= 5 & AGEYR < 19 ~ "5-18",
           AGEYR >= 19 & AGEYR < 31 ~ "19-30",
           AGEYR >= 31 & AGEYR < 65 ~ "31-64",
           AGEYR >= 65 ~ "65+",
           is.na(AGEYR) ~ "Unknown" # 16 IDs have no AGEYR
           # TRUE ~ "Unknown" 
         ),
         current.region.name = ifelse(current.region.name == "EASTERN", "EAST", current.region.name), # Wrong perception of "EASTERN" that should be "EAST"
         current.region.name = case_when(
           current.region.name == "E MIDS" ~ "East Midlands",
           current.region.name == "EAST" ~ "East of England",
           current.region.name == "LONDON" ~ "London",
           current.region.name == "N EAST" ~ "North East",
           current.region.name == "N WEST" ~ "North West",
           current.region.name == "S EAST" ~ "South East",
           current.region.name == "S WEST" ~ "South West",
           current.region.name == "W MIDS" ~ "West Midlands",
           current.region.name == "YORK&HUM" ~ "Yorkshire and The Humber",
           TRUE ~ current.region.name
         ),
         ageLabel = ifelse(AGEYR >= 90, 90, as.numeric(AGEYR)), # For incidence calculation, data grouped for people aged 90+
  ) #%>% 
  # glimpse()

# Basic case count data without age structure or regions
# Create all hypothetical recorded disease date

# Separate data based on ageGroup
temp_list <- list()
for(i in  unique(dat_G$ageGroup)){
  temp_list[i] <- dat_G[dat_G$ageGroup == i, "Earliest.specimen.date"]
  
}

# Manually extract data per-ageGroup (coz' for loop failed for list)
df_1_toddler <- data.frame(date = temp_list$`<5`)
df_1_toddler <- df_1_toddler %>% 
  dplyr::group_by(date) %>% 
  dplyr::summarise(cases_1_toddler = n()) %>% 
  ungroup()

df_2_518 <- data.frame(date = temp_list$`5-18`)
df_2_518 <- df_2_518 %>% 
  dplyr::group_by(date) %>% 
  dplyr::summarise(cases_2_518 = n()) %>% 
  ungroup()

df_3_1930 <- data.frame(date = temp_list$`19-30`)
df_3_1930 <- df_3_1930 %>% 
  dplyr::group_by(date) %>% 
  dplyr::summarise(cases_3_1930 = n()) %>% 
  ungroup()

df_4_3164 <- data.frame(date = temp_list$`31-64`)
df_4_3164 <- df_4_3164 %>% 
  dplyr::group_by(date) %>% 
  dplyr::summarise(cases_4_3164 = n()) %>% 
  ungroup()

df_5_65plus <- data.frame(date = temp_list$`65+`)
df_5_65plus <- df_5_65plus %>% 
  dplyr::group_by(date) %>% 
  dplyr::summarise(cases_5_65plus = n()) %>% 
  ungroup()


# ALL dates!
dat_G$Earliest.specimen.date <- as.Date(dat_G$Earliest.specimen.date)
all_date <- data.frame(allDate = seq.Date(from = min(dat_G$Earliest.specimen.date),
                                          to = max(dat_G$Earliest.specimen.date), 
                                          by = 1))
all_date$day <- 1:nrow(all_date)
# Coz the incidence only requires 2 columns called "counts" and "Day" in NUMBERS
# The counts (but in 0 counts the date are not recorded)

compiled <- all_date %>% 
  dplyr::left_join(df_1_toddler, by =c("allDate" = "date")) %>% 
  dplyr::left_join(df_2_518, by =c("allDate" = "date")) %>% 
  dplyr::left_join(df_3_1930, by =c("allDate" = "date")) %>% 
  dplyr::left_join(df_4_3164, by =c("allDate" = "date")) %>% 
  dplyr::left_join(df_5_65plus, by =c("allDate" = "date")) %>% 
  replace(is.na(.), 0) #%>% # NA means no data of meningitis or 30 days death, changed them to 0
# glimpse()


# Total population data by age, year for each region
# SOURCE: https://www.nomisweb.co.uk/
# pop <- read_excel("nomis_2024_04_15_124553_DCedit.xlsx") %>% 
# glimpse()
# I don't think I need total population for now,
# Examples on https://github.com/mrc-ide/mcstate/blob/master/inst/sir_incidence.csv
# Requires case count per aligned day only

# Viz per-day counts by base R plot
# png("pictures/daily_cases.png", width = 17, height = 30, unit = "cm", res = 1200)
par(mfrow= c(5, 1), bty = "n", mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0))
col_imD <- c(counts_Ser1 = "deepskyblue3")
for(i in c("cases_1_toddler", "cases_2_518", "cases_3_1930", "cases_4_3164", "cases_5_65plus")){
  
  plot(compiled$allDate, compiled[[i]], type = "b",
       xlab = "Date (year)", ylab = "Counts",
       ylim = c(0, max(compiled[[i]]+1)),
       col = col_imD[1], pch = 20,
       main = paste0(i))
  
  legend("topleft", names(col_imD), fill = col_imD, bty = "n")
}
# dev.off()

# Viz per-week counts by ggplot
Natm_n_imD <- compiled
Natm_n_imD$weeks <- cut(Natm_n_imD[,"allDate"], breaks="week")
Nat_weekly <- Natm_n_imD %>% 
  dplyr::group_by(weeks) %>% 
  dplyr::summarise(cases_1_toddler_weekly = sum(cases_1_toddler),
                   cases_2_518_weekly = sum(cases_2_518),
                   cases_3_1930_weekly = sum(cases_3_1930),
                   cases_4_3164_weekly = sum(cases_4_3164),
                   cases_5_65plus_weekly = sum(cases_5_65plus)) %>% 
  dplyr::ungroup() #%>% 

# png("pictures/weekly_cases.png", width = 17, height = 18, unit = "cm", res = 1200)
cases_1 <- ggplot(Nat_weekly, aes(as.Date(weeks))) +
  geom_line(aes(y = cases_1_toddler_weekly), colour = "deepskyblue3") +
  labs(x = NULL, y = NULL) +
  ggtitle("<5")

cases_2 <- ggplot(Nat_weekly, aes(as.Date(weeks))) +
  geom_line(aes(y = cases_2_518_weekly), colour = "deepskyblue3") +
  labs(x = NULL, y = NULL) +
  ggtitle("5-18")

cases_3 <- ggplot(Nat_weekly, aes(as.Date(weeks))) +
  geom_line(aes(y = cases_3_1930_weekly), colour = "deepskyblue3") +
  labs(x = NULL, y = NULL) +
  ggtitle("19-30")

cases_4 <- ggplot(Nat_weekly, aes(as.Date(weeks))) +
  geom_line(aes(y = cases_4_3164_weekly), colour = "deepskyblue3") +
  labs(x = NULL, y = NULL) +
  ggtitle("31-64")

cases_5 <- ggplot(Nat_weekly, aes(as.Date(weeks))) +
  geom_line(aes(y = cases_5_65plus_weekly), colour = "deepskyblue3") +
  labs(x = NULL, y = NULL) +
  ggtitle("65+")

cases <- cowplot::plot_grid(cases_1, cases_2, cases_3, cases_4, cases_5,
                   # labels = c("A", "B", "C", "D", "E"),
                   nrow = 5,
                   scale = 0.95)

x_label <- cowplot::ggdraw() +
  cowplot::draw_label("Year", x=0.53, y=  0, vjust=-0.5, angle= 0)
y_label <- cowplot::ggdraw() +
  cowplot::draw_label("Serotype 1 Cases (Aggregated by Week)", x=0, y=-50, vjust= 1.5, angle=90)

# Combine all together
final_plot <- cowplot::plot_grid(
  y_label, cases, x_label,
  ncol = 1,
  rel_heights = c(0.01, 1, 0.01),
  rel_widths = c(1)
)
final_plot
# dev.off()

## 2. Data Fitting #############################################################
# The anatomy of an mcstate particle filter, as noted above, consists of three main components: \n 
# 1. A set of observations to fit the model to, generated using mcstate::particle_filter_data(). \n 
# 2. A model to fit, which must be a dust generator, either dust::dust() or odin.dust::odin_dust(). \n 
# 3. A comparison function, which is an R function which calculates the likelihood of the state given the data at one time point.

# There is a calibration function in mcstate to fit our model to data.
# https://mrc-ide.github.io/mcstate/articles/sir_models.html
incidence <- compiled %>% 
  dplyr::select(-allDate) # ignore allDate

dir.create("inputs")
write.csv(incidence, "inputs/incidence.csv", row.names = FALSE)
write.csv(compiled, "inputs/incidence_weekly.csv", row.names = FALSE)


