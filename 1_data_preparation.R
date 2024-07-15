# Data preparation
library(tidyverse)
library(readxl)
library(epitools)
library(BactDating)

## 1. Data Viz and Analysis! ###################################################
# New updated data with meningitis (25.04.2024)
# All df are stored in raw_data
dat <- readxl::read_excel("raw_data/serotype1_UKHSA_imperial_date_age_region_MOLIS_withdeath_meningitis_clean.xlsx") #%>% 
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
         ),
         ageGroup7 = case_when(
           AGEYR < 2 ~ "<2",
           AGEYR >= 2 & AGEYR < 5 ~ "2-4",
           AGEYR >= 5 & AGEYR < 15 ~ "5-14",
           AGEYR >= 15 & AGEYR < 31 ~ "15-30", # Edit the Age-band into 15-30 & 31-44
           AGEYR >= 31 & AGEYR < 45 ~ "31-44", # Edit the Age-band into 15-30 & 31-44
           AGEYR >= 45 & AGEYR < 65 ~ "45-64",
           AGEYR >= 65 ~ "65+",
           is.na(AGEYR) ~ "Unknown" # 16 IDs have no AGEYR
         ),
         ageGroup2 = case_when(
           AGEYR < 15 ~ "children",
           AGEYR >= 15 ~ "adults",
           is.na(AGEYR) ~ "Unknown" # 16 IDs have no AGEYR
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

# EpiDescription based on incidences and CI
# Total population data by age, year for each region
# SOURCE: https://www.nomisweb.co.uk/
pop <- readxl::read_excel("raw_data/nomis_2024_04_15_124553_DCedit.xlsx") #%>% 
  # glimpse()

pop_l <- pop %>% 
  tidyr::pivot_longer(cols = `2001`:`2022`,
               names_to = "Year",
               values_to = "PopSize") %>% 
  dplyr::mutate(Age = gsub("Age ", "", Age),
         Age = ifelse(Age == "Aged 90+", 90, as.numeric(Age)), # For incidence calculation, data grouped for people aged 90+
         ageGroup = case_when( # edit 5 age bands
           Age < 5 ~ "<5",
           Age >= 5 & Age < 19 ~ "5-18",
           Age >= 19 & Age < 31 ~ "19-30",
           Age >= 31 & Age < 65 ~ "31-64",
           Age >= 65 ~ "65+",
           is.na(Age) ~ "Unknown" # 16 IDs have no Age
           # TRUE ~ "Unknown" 
         ),
         ageGroup7 = case_when(
           Age < 2 ~ "<2",
           Age >= 2 & Age < 5 ~ "2-4",
           Age >= 5 & Age < 15 ~ "5-14",
           Age >= 15 & Age < 31 ~ "15-30", # Edit the Age-band into 15-30 & 31-44
           Age >= 31 & Age < 45 ~ "31-44", # Edit the Age-band into 15-30 & 31-44
           Age >= 45 & Age < 65 ~ "45-64",
           Age >= 65 ~ "65+",
           is.na(Age) ~ "Unknown" # 16 IDs have no AGEYR
         ),
         ageGroup2 = case_when(
           Age < 15 ~ "children",
           Age >= 15 ~ "adults",
           is.na(Age) ~ "Unknown" # 16 IDs have no AGEYR
         ),
         Year = as.numeric(Year)) %>% 
  glimpse()

# Vaccination programme:
# SOURCE: https://www.gov.uk/government/publications/pneumococcal-the-green-book-chapter-25
vaccine_UK <- data.frame(
  year = c(2006, 2011),
  vaccine = c("PCV7", "PCV13")
)

# Simple counts & incidence per year
all_year <- dat_G %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise(counts = n()) %>% 
  dplyr::ungroup()

pop_year <- pop_l %>% 
  dplyr::group_by(Year) %>% 
  dplyr::summarise(PopSize = sum(PopSize)) %>% 
  dplyr::ungroup()

all_combined <- merge(all_year, pop_year,
                       by.x = c("year"),
                       by.y = c("Year")) %>%
  dplyr::mutate(Conf_Int = epitools::binom.exact(counts, PopSize),
                incid_Ser1 = Conf_Int$proportion) # per-100,000 population

# Colour names:
# https://www.datanovia.com/en/blog/awesome-list-of-657-r-color-names/
col_map <- c(# 5 age bands 
             "<5" = "indianred4",
             "5-18" = "orange",
             "19-30" = "seagreen4",
             "31-64" = "steelblue",
             "65+" = "purple3",
             "Unknown" = "black",
             # 7 age bands
             "<2" = "indianred4", 
             "2-4" = "indianred2", 
             "5-14" = "orange",
             "15-30" = "seagreen1", # Edit the Age-band into 15-30 & 31-44
             "31-44" = "seagreen4", # Edit the Age-band into 15-30 & 31-44
             "45-64" = "steelblue",
             "65+" = "purple3",
             # 2 age bands
             "children" = "darkred",
             "adults" = "darkblue"
)

vacc_map <- c("PCV7" = "gray80",
              "PCV13" = "gray20")

col_imD <- c(incid_Ser1 = "deepskyblue3",
             incid_m = "green",
             incid_D = "maroon")

# Viz counts
png("pictures/counts_allages.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(all_year, aes(x = year, y = counts)) +
  geom_line(size = 1.5) +
  geom_vline(data = vaccine_UK, aes(xintercept = year),
             linetype = "dashed") +
  # scale_color_manual(values = "black"
  # ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  geom_label(aes(x = 2006, y = 150, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = 2011, y = 150, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("The Counts of Serotype 1 in England") +
  xlab("Year") +
  ylab("Serotype 1 Cases")
dev.off()

# Viz incidence
png("pictures/incidence_allages.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(all_combined, aes(x = year, y = Conf_Int$proportion*100000)) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = Conf_Int$lower*100000, ymax = Conf_Int$upper*100000), # It doesn't matter whether I add the CI or not because the Pop data is quite huge, I suppose (?)
                width = .1) +
  geom_vline(data = vaccine_UK, aes(xintercept = year),
             linetype = "dashed") +
  scale_color_manual(values = "black"
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  scale_linetype_manual(values = c(vacc_map),
                        name = "Vaccine",
                        labels = c("PCV7", "PCV13")) +
  geom_label(aes(x = 2006, y = 0.15, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = 2011, y = 0.15, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("The Incidence of Serotype 1 in England \n(per 100,000)") +
  xlab("Year") +
  ylab("Serotype 1 Incidence")
dev.off()


# CI calculations for children-adults
ageGroup2 <- dat_G %>% 
  dplyr::group_by(year, ageGroup2) %>% 
  dplyr::summarise(counts = n()) %>% 
  dplyr::ungroup()

pop_ageGroup2 <- pop_l %>% 
  dplyr::group_by(Year, ageGroup2) %>% 
  dplyr::summarise(PopSize = sum(PopSize)) %>% 
  dplyr::ungroup()

all_ageGroup2 <- merge(ageGroup2, pop_ageGroup2,
               by.x = c("year","ageGroup2"),
               by.y = c("Year", "ageGroup2")) %>%
  dplyr::mutate(Conf_Int = epitools::binom.exact(counts, PopSize),
         incid_Ser1 = Conf_Int$proportion) # per-100,000 population

write.csv(all_ageGroup2, "raw_data/incidence_CI_per_year_2_ageGroup.csv", row.names = FALSE)

# Viz counts
png("pictures/counts_2ageGroups.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(all_ageGroup2, aes(x = year, y = counts, group = ageGroup2,
                color = ageGroup2)) +
  geom_line(size = 1.5) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(col_map),
                     name = "Demographic",
                     breaks = c("children", "adults"),
                     labels = c("Children (< 15)", "Adults")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  geom_label(aes(x = 2006, y = 150, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = 2011, y = 150, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("The Counts of Serotype 1 in England by Demographic Groups") +
  xlab("Year") +
  ylab("Serotype 1 Cases")
dev.off()

# Viz incidence
png("pictures/incidence_2ageGroups.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(all_ageGroup2, aes(x = year, y = Conf_Int$proportion*100000, group = ageGroup2,
                  color = ageGroup2)) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = Conf_Int$lower*100000, ymax = Conf_Int$upper*100000), # It doesn't matter whether I add the CI or not because the Pop data is quite huge, I suppose (?)
                width = .1) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(col_map),
                     name = "Demographic",
                     breaks = c("children", "adults"),
                     labels = c("Children (< 15)", "Adults")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  scale_linetype_manual(values = c(vacc_map),
                        name = "Vaccine",
                        labels = c("PCV7", "PCV13")) +
  geom_label(aes(x = 2006, y = 0.15, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = 2011, y = 0.15, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("The Incidence of Serotype 1 in England \nby Demographic Groups (per 100,000)") +
  xlab("Year") +
  ylab("Serotype 1 Incidence")
dev.off()

# CI calculations for 5 ageGroups
ageGroup5 <- dat_G %>% 
  dplyr::group_by(year, ageGroup) %>% 
  dplyr::summarise(counts = n()) %>% 
  dplyr::ungroup()

pop_ageGroup5 <- pop_l %>% 
  dplyr::group_by(Year, ageGroup) %>% 
  dplyr::summarise(PopSize = sum(PopSize)) %>% 
  dplyr::ungroup()

all_ageGroup5 <- merge(ageGroup5, pop_ageGroup5,
                       by.x = c("year","ageGroup"),
                       by.y = c("Year", "ageGroup")) %>%
  dplyr::mutate(Conf_Int = epitools::binom.exact(counts, PopSize),
                incid_Ser1 = Conf_Int$proportion) # per-100,000 population

write.csv(all_ageGroup5, "raw_data/incidence_CI_per_year_5_ageGroup.csv", row.names = FALSE)

# Viz counts
png("pictures/counts_5ageGroups.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(all_ageGroup5, aes(x = year, y = counts, group = ageGroup,
                          color = ageGroup)) +
  geom_line(size = 1.5) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(col_map),
                     name = "Demographic",
                     breaks = c("<5", "5-18", "19-30", "31-64", "65+", "Unknown"),
                     labels = c("<5", "5-18", "19-30", "31-64", "65+", "Unknown")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  geom_label(aes(x = 2006, y = 175, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = 2011, y = 250, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("The Counts of Serotype 1 in England by Demographic Groups") +
  xlab("Year") +
  ylab("Serotype 1 Cases")
dev.off()

# Viz incidence
png("pictures/incidence_5ageGroups.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(all_ageGroup5, aes(x = year, y = Conf_Int$proportion*100000, group = ageGroup,
                  color = ageGroup)) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = Conf_Int$lower*100000, ymax = Conf_Int$upper*100000), # It doesn't matter whether I add the CI or not because the Pop data is quite huge, I suppose (?)
                width = .1) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(col_map),
                     name = "Demographic",
                     breaks = c("<5", "5-18", "19-30", "31-64", "65+", "Unknown"),
                     labels = c("<5", "5-18", "19-30", "31-64", "65+", "Unknown")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  scale_linetype_manual(values = c(vacc_map),
                        name = "Vaccine",
                        labels = c("PCV7", "PCV13")) +
  geom_label(aes(x = 2006, y = 0.15, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = 2011, y = 0.15, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("The Incidence of Serotype 1 in England \nby Demographic Groups (per 100,000)") +
  xlab("Year") +
  ylab("Serotype 1 Incidence")
dev.off()

# CI calculations for 7 ageGroup7s
ageGroup7 <- dat_G %>% 
  dplyr::group_by(year, ageGroup7) %>% 
  dplyr::summarise(counts = n()) %>% 
  dplyr::ungroup()

pop_ageGroup7 <- pop_l %>% 
  dplyr::group_by(Year, ageGroup7) %>% 
  dplyr::summarise(PopSize = sum(PopSize)) %>% 
  dplyr::ungroup()

all_ageGroup7 <- merge(ageGroup7, pop_ageGroup7,
                       by.x = c("year","ageGroup7"),
                       by.y = c("Year", "ageGroup7")) %>%
  dplyr::mutate(Conf_Int = epitools::binom.exact(counts, PopSize),
                incid_Ser1 = Conf_Int$proportion) # per-100,000 population

write.csv(all_ageGroup7, "raw_data/incidence_CI_per_year_7_ageGroup.csv", row.names = FALSE)

# Viz counts
png("pictures/counts_7ageGroups.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(all_ageGroup7, aes(x = year, y = counts, group = ageGroup7,
                          color = ageGroup7)) +
  geom_line(size = 1.5) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(col_map),
                     name = "Demographic",
                     breaks = c("<2", "2-4", "5-14", "15-30", "31-44", "45-64", "65+", "Unknown"),
                     labels = c("<2", "2-4", "5-14", "15-30", "31-44", "45-64", "65+", "Unknown")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  geom_label(aes(x = 2006, y = 150, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = 2011, y = 150, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("The Counts of Serotype 1 in England by Demographic Groups") +
  xlab("Year") +
  ylab("Serotype 1 Cases")
dev.off()

# Viz incidence
png("pictures/incidence_7ageGroups.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(all_ageGroup7, aes(x = year, y = Conf_Int$proportion*100000, group = ageGroup7,
                          color = ageGroup7)) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = Conf_Int$lower*100000, ymax = Conf_Int$upper*100000), # It doesn't matter whether I add the CI or not because the Pop data is quite huge, I suppose (?)
                width = .1) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(col_map),
                     name = "Demographic",
                     breaks = c("<2", "2-4", "5-14", "15-30", "31-44", "45-64", "65+", "Unknown"),
                     labels = c("<2", "2-4", "5-14", "15-30", "31-44", "45-64", "65+", "Unknown")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  scale_linetype_manual(values = c(vacc_map),
                        name = "Vaccine",
                        labels = c("PCV7", "PCV13")) +
  geom_label(aes(x = 2006, y = 2.5, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = 2011, y = 2.5, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("The Incidence of Serotype 1 in England \nby Demographic Groups (per 100,000)") +
  xlab("Year") +
  ylab("Serotype 1 Incidence")
dev.off()

# Basic case count data without age structure or regions
# Create all hypothetical recorded disease date
dat_G$Earliest.specimen.date <- as.Date(dat_G$Earliest.specimen.date)
all_date <- data.frame(allDate = seq.Date(from = min(dat_G$Earliest.specimen.date),
                                          to = max(dat_G$Earliest.specimen.date), 
                                          by = 1))
all_date$day <- 1:nrow(all_date)

# Coz the incidence only requires 2 columns called "counts" and "Day" in NUMBERS
# The counts (but in 0 counts the date are not recorded)
Natm_ni <- dat_G %>% 
  dplyr::group_by(Earliest.specimen.date) %>% 
  dplyr::summarise(counts_Ser1 = n()) %>% 
  dplyr::ungroup() #%>% 
# glimpse()

Natm_nmeningitis <- dat_G %>% 
  dplyr::filter(MeningitisFlag == "Y") %>% 
  dplyr::group_by(Earliest.specimen.date) %>% 
  dplyr::summarise(counts_meningitis = n()) %>% 
  dplyr::ungroup() #%>% 
# glimpse()

Natm_n30DDeath <- dat_G %>% 
  dplyr::filter(`30daydeath` == "D") %>% 
  dplyr::group_by(Earliest.specimen.date) %>% 
  dplyr::summarise(counts_30DDeath = n()) %>% 
  dplyr::ungroup() #%>% 
# glimpse()


# Create a new df based on counts per day for Serotype 1, meningitis, and 30 days death
Natm_n_i <- dplyr::full_join(all_date, Natm_ni,
                      by = c("allDate" = "Earliest.specimen.date"))

Natm_n_im <- dplyr::full_join(Natm_n_i, Natm_nmeningitis,
                       by = c("allDate" = "Earliest.specimen.date"))

Natm_n_imD <- dplyr::full_join(Natm_n_im, Natm_n30DDeath,
                        by = c("allDate" = "Earliest.specimen.date")) %>% 
  replace(is.na(.), 0) #%>% # NA means no data of meningitis or 30 days death, changed them to 0
  # glimpse()


# Examples on https://github.com/mrc-ide/mcstate/blob/master/inst/sir_incidence.csv
# Requires case count per aligned day only

# Viz per-day counts by base R plot
png("pictures/daily_cases.png", width = 17, height = 12, unit = "cm", res = 1200)
par(bty = "n", mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0))
col_imD <- c(counts_Ser1 = "deepskyblue3",
             counts_meningitis = "green",
             counts_30DDeath = "maroon")
plot(Natm_n_imD$allDate, Natm_n_imD$counts_Ser1, type = "b",
     xlab = "Date (year)", ylab = "Counts",
     ylim = c(0, max(Natm_n_imD$counts_Ser1)+2),
     col = col_imD[1], pch = 20)

lines(Natm_n_imD$allDate, Natm_n_imD$counts_meningitis,
      type = "b", col = col_imD[2], pch = 20)
lines(Natm_n_imD$allDate, Natm_n_imD$counts_30DDeath,
      type = "b", col = col_imD[3], pch = 20)
legend("topleft", names(col_imD), fill = col_imD, bty = "n")
dev.off()

# Viz per-week counts by base R plot
# https://stackoverflow.com/questions/30431444/plotting-by-week-with-ggplot-in-r
# https://stackoverflow.com/questions/3777174/plotting-two-variables-as-lines-using-ggplot2-on-the-same-graph
Natm_n_imD$weeks <- cut(Natm_n_imD[,"allDate"], breaks="week")
Nat_weekly <- Natm_n_imD %>% 
  dplyr::group_by(weeks) %>% 
  dplyr::summarise(counts_Ser1_weekly = sum(counts_Ser1),
                   counts_meningitis_weekly = sum(counts_meningitis),
                   counts_30DDeath_weekly = sum(counts_30DDeath)) %>% 
  dplyr::ungroup() #%>% 

dir.create("inputs")

incidence_weekly <- Nat_weekly %>% 
  dplyr::select(weeks, counts_Ser1_weekly) %>% 
  dplyr::rename(cases = counts_Ser1_weekly) # That annoying name

write.csv(incidence_weekly, "inputs/incidence_weekly.csv", row.names = FALSE)

png("pictures/weekly_cases.png", width = 17, height = 12, unit = "cm", res = 1200)
col_imD_weekly <- c(counts_Ser1_weekly = "deepskyblue3",
             counts_meningitis_weekly = "green",
             counts_30DDeath_weekly = "maroon")
ggplot(Nat_weekly, aes(as.Date(weeks))) +
  geom_line(aes(y = counts_Ser1_weekly, colour = "counts_Ser1_weekly")) +
  geom_line(aes(y = counts_meningitis_weekly, colour = "counts_meningitis_weekly")) +
  geom_line(aes(y = counts_30DDeath_weekly, colour = "counts_30DDeath_weekly")) +
  scale_x_date() +
  scale_color_manual(values = col_imD_weekly,
                     name = "Cases",
                     breaks = c("counts_Ser1_weekly", "counts_meningitis_weekly", "counts_30DDeath_weekly"),
                     labels = c("Serotype 1", "Meningitis", "30 Days Death")
  ) +
  ggtitle("The Counts of Serotype 1 in England") +
  xlab("Year") +
  ylab("Serotype 1 Cases (Aggregated by Week)")
dev.off()

## 2. Data Fitting #############################################################
# The anatomy of an mcstate particle filter, as noted above, consists of three main components: \n 
# 1. A set of observations to fit the model to, generated using mcstate::particle_filter_data(). \n 
# 2. A model to fit, which must be a dust generator, either dust::dust() or odin.dust::odin_dust(). \n 
# 3. A comparison function, which is an R function which calculates the likelihood of the state given the data at one time point.

# There is a calibration function in mcstate to fit our model to data.
# https://mrc-ide.github.io/mcstate/articles/sir_models.html
incidence <- Natm_n_imD %>% 
  dplyr::select(day, counts_Ser1) %>% 
  dplyr::rename(cases = counts_Ser1) # That annoying name

dir.create("inputs")
write.csv(incidence, "inputs/incidence.csv", row.names = FALSE)

png("pictures/hist_daily_cases.png", width = 17, height = 12, unit = "cm", res = 1200)
hist(incidence$cases,
     main = "Histogram of Daily Cases",
     xlab = "Daily Incidence") # huge zero daily cases occur
dev.off()

## 3. Serotypes and GPSCs Determination ########################################
# Load dat_G first.
link_ID <- readxl::read_excel("raw_data/gubbins/ukhsa_assemblies_02_07_24.xlsx")

# Load All Serotype 1: n739 samples of GPSC31 and n7 samples of GPSC2
all_serotype1 <- dplyr::bind_rows(
  read.table("raw_data/gubbins/remove_GPSC31_list.txt") %>%
    dplyr::mutate(GPSC = 31),
  read.table("raw_data/gubbins/remove_GPSC2_list.txt") %>%
    dplyr::mutate(GPSC = 2)
) %>% 
  dplyr::mutate(serotype = 1)

# Load n703 samples of GPSC31; filtered based on N50 >= 30kb
link_GPSC31_n703 <- read.table("raw_data/gubbins/remove_GPSC31_choosen_n703_list.txt") %>% 
  dplyr::mutate(n703_choosen_GPSC31 = 1)

joined_serotype_GPSC <- link_ID %>% 
  dplyr::full_join(all_serotype1, by = c("assembly_name" = "V1")) %>% 
  dplyr::full_join(link_GPSC31_n703, by = c("assembly_name" = "V1"))

dat_G <- dat_G %>% 
  dplyr::full_join(joined_serotype_GPSC, by = "ID")

## 4. AMR Analysis #############################################################
# Load dat_G first.
resistance_smx <- 
  dplyr::bind_rows(
    read.csv("raw_data/gubbins/n739/resistance_folp_smx.csv") %>% 
      dplyr::rename(resistance_smx = Resistance),
    read.csv("raw_data/gubbins/GPSC2_n7/resistance_folp_smx.csv") %>% 
      dplyr::rename(resistance_smx = Resistance)
    )

resistance_tmp <- 
  dplyr::bind_rows(
    read.csv("raw_data/gubbins/n739/resistance_dhfr_tmp.csv") %>% 
      dplyr::rename(resistance_tmp = Resistance),
    read.csv("raw_data/gubbins/GPSC2_n7/resistance_dhfr_tmp.csv") %>% 
      dplyr::rename(resistance_tmp = Resistance)
  )

resistance_smx_tmp <- 
  dplyr::full_join(resistance_smx, resistance_tmp, by = "isolate_id") %>% 
  dplyr::mutate(isolate_id = paste0(isolate_id, ".fasta"))

dat_G <- dat_G %>% 
  dplyr::full_join(resistance_smx_tmp, by = c("assembly_name" = "isolate_id")) %>% 
  dplyr::mutate(ngsid = substr(assembly_name, 1, 8)) # correction for ngsid including those that sequenced but have no EpiData

## 5. Clade Analysis from Microreact ###########################################
# Load dat_G first.

# Context: Run microreact first and determine some clades.
clade_assignment_df <-
  dplyr::bind_rows(
    read.table("raw_data/gubbins/n703/n703_clade1.txt") %>% dplyr::mutate(clade = "clade1"),
    read.table("raw_data/gubbins/n703/n703_clade2.txt") %>% dplyr::mutate(clade = "clade2"),
    read.table("raw_data/gubbins/n703/n703_clade3.txt") %>% dplyr::mutate(clade = "clade3"),
  ) %>% 
  dplyr::mutate(V1 = paste0(V1, ".fasta"))

joined_clades <- 
  dplyr::full_join(dat_G, clade_assignment_df, by = c("assembly_name" = "V1"))

# Clade viz!
freq_clades <- joined_clades %>%
  dplyr::group_by(year) %>%
  dplyr::mutate(Sample_size = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(clade,year) %>%
  dplyr::mutate(Frequency = n()/Sample_size) %>%
  dplyr::mutate(Count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(clade),
                !is.na(year)) %>% # Possibly NA in year is 2015
  dplyr::mutate(Conf_Int = epitools::binom.exact(Count, Sample_size),
                Prop = Conf_Int$proportion) %>% # per-100,000 population
  dplyr::select(clade,year,Frequency, Count, Prop, Conf_Int, current.region.name) %>%
  dplyr::distinct()

ggplot(freq_clades,
       aes(x = year, y = Frequency,
           colour = clade,
           fill = clade)) +
  geom_line() +
  geom_errorbar(aes(ymin = Conf_Int$lower, ymax = Conf_Int$upper),
                width = .05) +
  geom_vline(data = vaccine_UK, aes(xintercept = year),
             linetype = "dashed") +
  geom_label(aes(label = Count),
             # fill = "white",
             color = "black",
             nudge_y = 0.001,
             nudge_x = 0.05) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  geom_label(aes(x = 2006, y = 0.15, label = "PCV7"),
             fill = "white", color = "black") +
  geom_label(aes(x = 2011, y = 0.15, label = "PCV13"),
             fill = "white", color = "black") +
  theme_bw()

max_freq_region <- freq_clades %>% 
  dplyr::group_by(year, clade) %>%
  dplyr::arrange(desc(Frequency)) %>%
  dplyr::slice(1)

ggplot(max_freq_region,
       aes(x = year, y = Frequency,
           colour = clade,
           fill = clade)) +
  geom_line() +
  geom_errorbar(aes(ymin = Conf_Int$lower, ymax = Conf_Int$upper),
                width = .05) +
  geom_vline(data = vaccine_UK, aes(xintercept = year),
             linetype = "dashed") +
  geom_label(aes(label = paste0(Count, current.region.name)),
             # fill = "white",
             color = "black",
             nudge_y = 0.001,
             nudge_x = 0.05) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  geom_label(aes(x = 2006, y = 0.15, label = "PCV7"),
             fill = "white", color = "black") +
  geom_label(aes(x = 2011, y = 0.15, label = "PCV13"),
             fill = "white", color = "black") +
  theme_bw()

## 6. Data Preparation for Microreact ##########################################
# Load dat_G first.
tre <- BactDating::loadGubbins("raw_data/gubbins/n739/n739_")

tre_names <- as.data.frame(tre$tip.label)
tre_names$ID_contigs <- substr(tre$tip.label, 1, 8)
tre_names <- 
  dplyr::full_join(tre_names, joined_clades, by = c("ID_contigs" = "ngsid"))

write.csv(tre_names, "raw_data/gubbins/n739/phandango_microreact_check/microreact_tre_names.csv", row.names = FALSE)
write.csv(tre_names, "raw_data/gubbins/n703/phandango_microreact_check/microreact_tre_names.csv", row.names = FALSE)
write.csv(tre_names, "raw_data/gubbins/GPSC2_n7/phandango_microreact_check/microreact_tre_names.csv", row.names = FALSE)


## 7. Some Stats Analysis Between Clades #######################################
joined_clades <- read.csv("raw_data/gubbins/n739/phandango_microreact_check/microreact_tre_names.csv")

# GPSCs vs. ages
anova_GPSCs <- aov(ageLabel ~ GPSC, data = joined_clades)
summary(anova_GPSCs)
# > summary(anova_GPSCs)
# Df Sum Sq Mean Sq F value Pr(>F)
# GPSC          1     57    57.0   0.097  0.755
# Residuals   731 429044   586.9               
# 3342 observations deleted due to missingness
# No statistically significant differences in the mean ages across GPSCs (Pr = 0.324).


# GPSC31 Clade vs. ages
anova_clade_GPSC31 <- aov(ageLabel ~ clade, data = joined_clades)
summary(anova_clade_GPSC31)
# > summary(anova_clade_GPSC31)
# Df Sum Sq Mean Sq F value Pr(>F)
# clade         2   1321   660.4   1.129  0.324
# Residuals   688 402409   584.9               
# 3384 observations deleted due to missingness
# No statistically significant differences in the mean ages across clades (Pr = 0.324).


# GPSCs vs. regions 
chi_GPSC <- chisq.test(table(joined_clades$current.region.name, joined_clades$GPSC))
print(chi_GPSC)
# data:  table(joined_clades$current.region.name, joined_clades$GPSC)
# X-squared = 13.979, df = 8, p-value = 0.08231

# GPSC31 Clade vs. regions 
chi_region <- chisq.test(table(joined_clades$current.region.name, joined_clades$clade))
print(chi_region)
# data:  table(joined_clades$current.region.name, joined_clades$clade)
# X-squared = 48.53, df = 16, p-value = 3.92e-05

# Post-hoc
std_resid_region <- chi_region$stdres
print(std_resid_region)
mosaicplot(table(joined_clades$current.region.name, joined_clades$clade),
           main = "Mosaic Plot of Region and GPSC31 Clades", shade = TRUE, legend = TRUE)

# Create barplot of region and clades:
filtered_NA <- joined_clades %>% 
  dplyr::filter(!is.na(current.region.name),
                !is.na(clade),
                !is.na(year))
filtered_NA_counts <- filtered_NA %>%
  group_by(year, current.region.name, clade) %>%
  summarise(count = n()) %>%
  ungroup()

ggplot(filtered_NA_counts, aes(x = year, y = count, colour = clade)) +
  geom_line() +
  labs(title = "Count of Clades by Region",
       x = "Region",
       y = "Counts",
       fill = "Clades of GPSC31") +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  facet_wrap(~ current.region.name, scales = "free") +
  # annotate("text", x = 2006, y = 5, label = "PCV7") +
  # annotate("text", x = 2011, y = 5, label = "PCV13") +
  annotate("segment", x = 2006, xend = 2006, y = -Inf, yend = Inf, linetype = "dashed") +
  annotate("segment", x = 2011, xend = 2011, y = -Inf, yend = Inf, linetype = "dashed") +
  theme_minimal()

# GPSCs vs. vaccination era
chi_vacc_GPSCs <- chisq.test(table(joined_clades$vacc, joined_clades$GPSC))

# Post-hoc
std_resid_vacc_GPSCs <- chi_vacc_GPSCs$stdres
print(std_resid_vacc_GPSCs)
mosaicplot(table(joined_clades$vacc, joined_clades$GPSC),
           main = "Mosaic Plot of Region and GPSC31 Clades", shade = TRUE, legend = TRUE)

# GPSC31 Clade vs. vaccination era
chi_vacc_clades <- chisq.test(table(joined_clades$vacc, joined_clades$clade))

# Post-hoc
std_resid_vacc_clades <- chi_vacc_clades$stdres
print(std_resid_vacc_clades)
mosaicplot(table(joined_clades$vacc, joined_clades$clade),
           main = "Mosaic Plot of Region and GPSC31 Clades", shade = TRUE, legend = TRUE)

# Other not-significant results:
# GPSCs vs. ages
chisq.test(table(joined_clades$ageGroup2, joined_clades$GPSC))
chisq.test(table(joined_clades$ageGroup, joined_clades$GPSC))
chisq.test(table(joined_clades$ageGroup7, joined_clades$GPSC))

# GPSCs vs. meningitis
chisq.test(table(joined_clades$MeningitisFlag, joined_clades$GPSC))
# GPSCs vs. death
# Highly correlated!(significance but nonsense coz' all death driven by GPSC31)
chi_death_GPSCs <- chisq.test(table(joined_clades$X30daydeath, joined_clades$GPSC))

# Post-hoc
std_resid_death_GPSCs <- chi_death_GPSCs$stdres
print(std_resid_death_GPSCs)
mosaicplot(table(joined_clades$X30daydeath, joined_clades$GPSC),
           main = "Mosaic Plot of Region and GPSC31 Clades", shade = TRUE, legend = TRUE)


# GPSC31 Clade vs. ages
chisq.test(table(joined_clades$ageGroup2, joined_clades$clade))
chisq.test(table(joined_clades$ageGroup, joined_clades$clade))
chisq.test(table(joined_clades$ageGroup7, joined_clades$clade))

# GPSC31 vs. meningitis
chisq.test(table(joined_clades$MeningitisFlag, joined_clades$clade))
# GPSC31 vs. death
chisq.test(table(joined_clades$X30daydeath, joined_clades$clade))

