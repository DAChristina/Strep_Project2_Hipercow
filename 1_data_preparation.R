# Data preparation
library(tidyverse)
library(readxl)
library(epitools)
library(BactDating)

## 1. Data Viz and Analysis! ###################################################
# New updated data with meningitis (25.04.2024)
# All df are stored in raw_data
dat <- dplyr::full_join(
  readxl::read_excel("raw_data/serotype1_UKHSA_imperial_date_age_region_MOLIS_withdeath_meningitis_clean.xlsx") %>% 
    dplyr::rename(Earliest.specimen.date = Earliestspecimendate,
                  current.region.name = currentregionname) %>% 
    dplyr::select(-sequenced),
  # The link to sequenced IDs
  readxl::read_excel("raw_data/gubbins/ukhsa_assemblies_02_07_24.xlsx"),
  by = "ID"
)


# Additional data from 2015 and other sequenced samples (n = 9 + n = 2):
add_n11_n20 <- dplyr::bind_rows(
  readxl::read_excel("raw_data/nine_missing_specimen_date_serotype1_extended.xlsx"),
  readxl::read_excel("raw_data/two_missing_specimen_date_serotype1_extended.xlsx"),
  # Additional 20 sequenced data from 2017 to 2023 (UPDATE 17 Oct. 2024)
  readxl::read_excel("raw_data/serotype1_UKHSA_imperial_routine.xlsx") %>% 
    dplyr::mutate(Earliestspecimendate = as.Date(Earliestspecimendate, format = "%Y-%m-%d"),
                  AGEYR = as.numeric(AGEYR)),
) %>% 
  dplyr::mutate(RevisedOPIEID = NA) %>% 
  dplyr::rename(Earliest.specimen.date = Earliestspecimendate,
                current.region.name = currentregionname) %>% 
  dplyr::select(ID, RevisedOPIEID, Earliest.specimen.date, AGEYR,
                current.region.name, MeningitisFlag, `30daydeath`,
                ngsid, assembly_name)

dat_G <- dplyr::bind_rows(dat, add_n11_n20) %>% 
  # dat %>% 
  dplyr::mutate(ngsid = as.numeric(ngsid),
                AGEYR = ifelse(AGEYR >= 90, 90, as.numeric(AGEYR)), # For incidence calculation, data grouped for people aged 90+
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
  )

# Mannually delete 3 rows with duplicated IDs 
# because distinct() & arrange()-filter() failed to produce what I really want
write.csv(dat_G, "raw_data/serotype1_UKHSA_imperial_date_age_region_MOLIS_sequenced_postThesis.csv", row.names = FALSE)

dat_G <- read.csv("raw_data/serotype1_UKHSA_imperial_date_age_region_MOLIS_sequenced_postThesis_cleaned.csv")

# Temporary data prep for microreact
temporary_microreact_check <- dat_G %>% 
  dplyr::mutate(microreact_ID = stringr::str_remove(assembly_name, ".fasta"),
                Earliest.specimen.date = as.Date(Earliest.specimen.date))
write.csv(temporary_microreact_check, "raw_data/temporary_microreact_check.csv", row.names = FALSE)


# EpiDescription based on incidences and CI
# Total population data by age (rows), year (columns) for each region
# I update population data from 2001 to 2023 (NOMISWEB Update 2024-07-24)
# DOWNLOADED 17 Oct. 2024
# SOURCE: https://www.nomisweb.co.uk/
# pop <- readxl::read_excel("raw_data/nomis_2024_04_15_124553_DCedit.xlsx") # ver.1 2001-2022

pop <- readxl::read_excel("raw_data/nomis_2024_10_17_DCedit.xlsx") %>%  # ver.2 2001-2023
  dplyr::mutate(`2024` = `2023`) # Temporary for 2024 population; ONS hasn't released the data yet!


pop_l <- pop %>% 
  tidyr::pivot_longer(cols = `2001`:`2024`,
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
col_map <- c(
  # Vaccination
  "Pre-PCV7" = "lightblue",
  "PCV7" = "gray70",
  "PCV13" = "gray20",
  # 5 age bands 
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
  "adults" = "darkblue",
  # Regions (From north to south)
  "North West" = "indianred4",
  "North East" = "steelblue",
  "Yorkshire and The Humber" = "seagreen4",
  "East Midlands" = "purple3",
  "West Midlands" = "orange",
  "East of England" = "indianred2",
  "London" = "seagreen1",
  "South East" = "deepskyblue",
  "South West" = "gold1",
  # Cases vs sequenced
  "Serotype 1 Case" = "gray75",
  "Sequenced" = "deepskyblue3",
  "Meningitis" = "green",
  "30 Day Death" = "maroon"
)

vacc_map <- c("PCV7" = "gray70",
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
  ylab("Serotype 1 Cases") +
  theme_bw()
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
  ylab("Serotype 1 Incidence") +
  theme_bw()
dev.off()

png("pictures/combined_allages.png", width = 28, height = 10, unit = "cm", res = 1200)
England_allYear <-
  cowplot::plot_grid(plotlist = list(ggplot(all_year, aes(x = year, y = counts)) +
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
                                       ggtitle("Case Counts of IPD Caused by Serotype 1\nin England") +
                                       xlab("Year") +
                                       ylab("Serotype 1 Cases") +
                                       theme_bw(),
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
                                       ggtitle("Incidence of IPD Caused by Serotype 1\nin England (per 100,000 population)") +
                                       xlab("Year") +
                                       ylab("Serotype 1 Incidence") +
                                       theme_bw()))

England_allYear
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
  ylab("Serotype 1 Cases") +
  theme_bw()
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
  ylab("Serotype 1 Incidence") +
  theme_bw()
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
  ylab("Serotype 1 Cases") +
  theme_bw()
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
  ylab("Serotype 1 Incidence") +
  theme_bw()
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

# Kruskal test
kruskal.test(Conf_Int$proportion ~ ageGroup7, data = all_ageGroup7)

# Post-hoc
FSA::dunnTest(Conf_Int$proportion ~ ageGroup7, data = all_ageGroup7, method = "bonferroni")


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
  ylab("Serotype 1 Cases") +
  theme_bw()
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
  ylab("Serotype 1 Incidence") +
  theme_bw()
dev.off()

png("pictures/combined_7ageGroups.png", width = 28, height = 10, unit = "cm", res = 1200)
England_7ageGroups <-
  cowplot::plot_grid(plotlist = list(ggplot(all_ageGroup7, aes(x = year, y = counts, group = ageGroup7,
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
                                       ggtitle("Case Counts of IPD Caused by Serotype 1\nin England by Demographic Groups") +
                                       xlab("Year") +
                                       ylab("Serotype 1 Cases") +
                                       theme_bw() +
                                       theme(legend.position="none"),
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
                                       ggtitle("Incidence of IPD Caused by Serotype 1\nin England by Demographic Groups (per 100,000 population)") +
                                       xlab("Year") +
                                       ylab("Serotype 1 Incidence") +
                                       theme_bw() +
                                       theme(legend.position="none")
                                     ))

England_7ageGroups

legend <- cowplot::get_legend(ggplot(all_ageGroup7, aes(x = year, y = counts, group = ageGroup7,
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
                                ggtitle("Case Counts of IPD Caused by Serotype 1\nin England by Demographic Groups") +
                                xlab("Year") +
                                ylab("Serotype 1 Cases") +
                                theme_bw() + 
                                theme(legend.box.margin = margin(0, 0, 0, 12))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
cowplot::plot_grid(England_7ageGroups, legend, rel_widths = c(2, 0.3))
dev.off()

# CI calculations for 9 Regions
region <- dat_G %>% 
  dplyr::group_by(year, current.region.name) %>% 
  dplyr::summarise(counts = n()) %>% 
  dplyr::ungroup()

pop_region <- pop_l %>% 
  dplyr::group_by(Year, Region) %>% 
  dplyr::summarise(PopSize = sum(PopSize)) %>% 
  dplyr::ungroup()

all_reg <- merge(region, pop_region,
                 by.x = c("year","current.region.name"),
                 by.y = c("Year", "Region")) %>%
  dplyr::mutate(Conf_Int = epitools::binom.exact(counts, PopSize),
                incid_Ser1 = Conf_Int$proportion) # per-100,000 population

write.csv(all_reg, "raw_data/incidence_CI_per_year_region.csv", row.names = FALSE)

# Kruskal test
kruskal.test(Conf_Int$proportion ~ current.region.name, data = all_reg)

# Post-hoc
FSA::dunnTest(Conf_Int$proportion ~ current.region.name, data = all_reg, method = "bonferroni")

# Viz counts
png("pictures/counts_yearRegion.png", width = 20, height = 12, unit = "cm", res = 1200)
ggplot(all_reg, aes(x = year, y = counts, group = current.region.name,
                    color = current.region.name)) +
  geom_line(size = 1.5) +
  geom_vline(data = vaccine_UK, aes(xintercept = year),
             colour = "black",
             linetype = "dashed") +
  scale_color_manual(values = c(col_map),
                     name = "Region",
                     breaks = c("North West", "North East", "Yorkshire and The Humber",
                                "East Midlands", "West Midlands", "East of England",
                                "London", "South East", "South West"),
                     labels =  c("North West", "North East", "Yorkshire and The Humber",
                                 "East Midlands", "West Midlands", "East of England",
                                 "London", "South East", "South West")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  geom_label(aes(x = 2006, y = 100, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = 2011, y = 100, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("The Counts of Serotype 1 in England by Region") +
  xlab("Year") +
  ylab("Serotype 1 Cases") +
  theme_bw()
dev.off()

# Viz incidence
png("pictures/incidence_yearRegion.png", width = 20, height = 12, unit = "cm", res = 1200)
ggplot(all_reg, aes(x = year, y = Conf_Int$proportion*100000, group = current.region.name,
                    color = current.region.name)) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = Conf_Int$lower*100000, ymax = Conf_Int$upper*100000), # It doesn't matter whether I add the CI or not because the Pop data is quite huge, I suppose (?)
                width = .1) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(col_map),
                     name = "Region",
                     breaks = c("North West", "North East", "Yorkshire and The Humber",
                                "East Midlands", "West Midlands", "East of England",
                                "London", "South East", "South West"),
                     labels =  c("North West", "North East", "Yorkshire and The Humber",
                                 "East Midlands", "West Midlands", "East of England",
                                 "London", "South East", "South West")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  scale_linetype_manual(values = c(vacc_map),
                        name = "Vaccine",
                        labels = c("PCV7", "PCV13")) +
  geom_label(aes(x = 2006, y = 2.5, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = 2011, y = 2.5, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("The Incidence of Serotype 1 in England \nby Region (per 100,000)") +
  xlab("Year") +
  ylab("Serotype 1 Incidence") +
  theme_bw()
dev.off()

png("pictures/combined_region.png", width = 30, height = 12, unit = "cm", res = 1200)
England_region <-
  cowplot::plot_grid(plotlist = list(ggplot(all_reg, aes(x = year, y = counts, group = current.region.name,
                                                         color = current.region.name)) +
                                       geom_line(size = 1.5) +
                                       geom_vline(data = vaccine_UK, aes(xintercept = year),
                                                  colour = "black",
                                                  linetype = "dashed") +
                                       scale_color_manual(values = c(col_map),
                                                          name = "Region",
                                                          breaks = c("North West", "North East", "Yorkshire and The Humber",
                                                                     "East Midlands", "West Midlands", "East of England",
                                                                     "London", "South East", "South West"),
                                                          labels =  c("North West", "North East", "Yorkshire and The Humber",
                                                                      "East Midlands", "West Midlands", "East of England",
                                                                      "London", "South East", "South West")
                                       ) +
                                       scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
                                       geom_label(aes(x = 2006, y = 100, label = "PCV7"),
                                                  fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
                                       geom_label(aes(x = 2011, y = 100, label = "PCV13"),
                                                  fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
                                       ggtitle("Case Counts of IPD Caused by Serotype 1\nin England by Region") +
                                       xlab("Year") +
                                       ylab("Serotype 1 Cases") +
                                       theme_bw() +
                                       theme(legend.position="none"),
                                     ggplot(all_reg, aes(x = year, y = Conf_Int$proportion*100000, group = current.region.name,
                                                         color = current.region.name)) +
                                       geom_line(size = 1.5) +
                                       geom_errorbar(aes(ymin = Conf_Int$lower*100000, ymax = Conf_Int$upper*100000), # It doesn't matter whether I add the CI or not because the Pop data is quite huge, I suppose (?)
                                                     width = .1) +
                                       geom_vline(data = vaccine_UK, aes(xintercept = year,
                                                                         colour = vaccine),
                                                  linetype = "dashed") +
                                       scale_color_manual(values = c(col_map),
                                                          name = "Region",
                                                          breaks = c("North West", "North East", "Yorkshire and The Humber",
                                                                     "East Midlands", "West Midlands", "East of England",
                                                                     "London", "South East", "South West"),
                                                          labels =  c("North West", "North East", "Yorkshire and The Humber",
                                                                      "East Midlands", "West Midlands", "East of England",
                                                                      "London", "South East", "South West")
                                       ) +
                                       scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
                                       scale_linetype_manual(values = c(vacc_map),
                                                             name = "Vaccine",
                                                             labels = c("PCV7", "PCV13")) +
                                       geom_label(aes(x = 2006, y = 2.5, label = "PCV7"),
                                                  fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
                                       geom_label(aes(x = 2011, y = 2.5, label = "PCV13"),
                                                  fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
                                       ggtitle("Incidence of IPD Caused by Serotype 1\nin England by Region (per 100,000 population)") +
                                       xlab("Year") +
                                       ylab("Serotype 1 Incidence") +
                                       theme_bw() +
                                       theme(legend.position="none")
  ))

England_region

legend <- cowplot::get_legend(ggplot(all_reg, aes(x = year, y = counts, group = current.region.name,
                                                  color = current.region.name)) +
                                geom_line(size = 1.5) +
                                geom_vline(data = vaccine_UK, aes(xintercept = year),
                                           colour = "black",
                                           linetype = "dashed") +
                                scale_color_manual(values = c(col_map),
                                                   name = "Region",
                                                   breaks = c("North West", "North East", "Yorkshire and The Humber",
                                                              "East Midlands", "West Midlands", "East of England",
                                                              "London", "South East", "South West"),
                                                   labels =  c("North West", "North East", "Yorkshire and The Humber",
                                                               "East Midlands", "West Midlands", "East of England",
                                                               "London", "South East", "South West")
                                ) +
                                scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
                                geom_label(aes(x = 2006, y = 100, label = "PCV7"),
                                           fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
                                geom_label(aes(x = 2011, y = 100, label = "PCV13"),
                                           fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
                                ggtitle("Case Counts of IPD Caused by Serotype 1\nin England by Region") +
                                xlab("Year") +
                                ylab("Serotype 1 Cases") +
                                theme_bw() +
                                theme(legend.box.margin = margin(0, 0, 0, 12))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
cowplot::plot_grid(England_region, legend, rel_widths = c(2, 0.4))
dev.off()

# 1.4. Choropleth! #############################################################
# 1.4.1. Three maps based on the vaccination era ###############################

# Workflow:
# 1. Create/load the df
# 2. Filter to the specific:
# 2.1. Year range
# 2.2. Regions of interest?
# 3. Data preparation for GADM map
# 4. Combine G_BF_ESPEN_Admin1_forGADM to *.shp data (CANNOT BE RUN VICE-VERSA)!!!

# load all_reg first
all_reg <- read.csv("raw_data/incidence_CI_per_year_region.csv")

selected_map <- all_reg %>% 
  dplyr::select(year, current.region.name, counts, PopSize) %>% 
  dplyr::mutate(
    vacc = case_when(
      year < 2006 ~ "Pre-PCV7",
      year >= 2006 & year < 2011 ~ "PCV7",
      year >= 2011 ~ "PCV13",
      TRUE ~ NA_character_
    )) %>% 
  dplyr::group_by(vacc, current.region.name) %>% 
  dplyr::mutate(ave_N = mean(PopSize)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(vacc, current.region.name) %>% 
  dplyr::mutate(n = sum(counts),
                Incid = n/ave_N) %>% 
  distinct()


# 1.4.1.2. The geometry data ###################################################
# 1. Download the *.shp data at:
# https://geoportal.statistics.gov.uk/datasets/ons::regions-december-2023-boundaries-en-bfc-2/about
GBR_shp_zip = "raw_data/Regions_December_2023_Boundaries_EN_BFC_-7302455802062307841.zip"
GBR_shp_out = "raw_data/GBR_shp"
unzip(GBR_shp_zip, exdir=GBR_shp_out)

# 2. Load the file
library(sf) # tidyverse cannot edit geospatial data.

GBR_spdf <- sf::read_sf(dsn = "raw_data/GBR_shp/RGN_DEC_2023_EN_BFC.shp")

# summary(GBR_spdf) # tells you the max and min coordinates, the kind of projection in use
# length(GBR_spdf) # how many regions you have

# Coz' this is an sf object
head(sf::st_drop_geometry(GBR_spdf))
glimpse(GBR_spdf)
sort(unique(GBR_spdf$RGN23NM))
sort(unique(dat_G$current.region.name))

# NEXT:
# U have to combine *.csv OR *.xlsx data to *.shp (CANNOT BE RUN VICE-VERSA)!!!
# Convert the 'SpatialPolygonsDataFrame' to 'sf' object first!

# Combine the df to *.shp data (CANNOT BE RUN VICE-VERSA)!!!
# Convert the 'SpatialPolygonsDataFrame' to 'sf' object first:
# Recall GBR_spdf <- st_read(dsn = GBR_shp_path_LINUX)

# Notes: DO NOT Filter sf Data!!!!
GBR_spdf_sf <- sf::st_as_sf(GBR_spdf, coords = c("longitude", "latitude"), crs = '4326')
glimpse(GBR_spdf_sf)

dat_1PrePCV7 <- selected_map %>% 
  dplyr::filter(vacc == "Pre-PCV7")
dat_2PCV7 <- selected_map %>% 
  dplyr::filter(vacc == "PCV7")
dat_3PCV13 <- selected_map %>% 
  dplyr::filter(vacc == "PCV13")

Comm_merged_1PrePCV7 <- merge(GBR_spdf_sf, dat_1PrePCV7, by.x = 'RGN23NM', by.y = 'current.region.name', all.x = TRUE) %>% 
  dplyr::mutate(Incid = as.numeric(Incid)) %>% 
  glimpse()
Comm_merged_2PCV7 <- merge(GBR_spdf_sf, dat_2PCV7, by.x = 'RGN23NM', by.y = 'current.region.name', all.x = TRUE) %>% 
  dplyr::mutate(Incid = as.numeric(Incid)) %>% 
  glimpse()
Comm_merged_3PCV13 <- merge(GBR_spdf_sf, dat_3PCV13, by.x = 'RGN23NM', by.y = 'current.region.name', all.x = TRUE) %>% 
  dplyr::mutate(Incid = as.numeric(Incid)) %>% 
  glimpse()


# Max-min counts for plot colours:
max_c <- max(selected_map$Incid)
min_c <- min(selected_map$Incid)

# Trial plot of the combined data:
# Compiled_geom is not required.
# Base plot trial source: https://r-charts.com/spatial/choropleth-map/?utm_content=cmp-true
plot(GBR_spdf_sf$geometry,
     main = "The 9 Regions of England")

# Centroids (failed):
# https://r-graph-gallery.com/169-170-basic-manipulation-of-shapefiles.html#centroid
# centroids <- st_centroid(GBR_spdf_sf)#, of_largest_polygon = F)
# centers <- cbind(centroids, st_coordinates(centroids)) # Small manipulation to add coordinates as columns

# Combined plot:
breaks_seq <- seq(min_c, max_c, length.out = 200)

plot(Comm_merged_1PrePCV7[, "Incid"],
     breaks = breaks_seq,
     # nbreaks = 4,
     pal = colorRampPalette(c("white", "maroon"))(199), # fill with (length breaks - 1)
     main = "The 9 Regions of England During the Pre-PCV7 Era")
# text(centers$X, centers$Y, Comm_merged_1PrePCV7$counts, cex = .9, col = "black")
# Weird success

# The real plot
DataMerged <- c("Comm_merged_1PrePCV7", "Comm_merged_2PCV7", "Comm_merged_3PCV13")
DataTitle <- list(Comm_merged_1PrePCV7 = "Pre-PCV7 Era (2003-2005)",
                  Comm_merged_2PCV7 = "PCV7 Era (2006-2010)",
                  Comm_merged_3PCV13 = "PCV13 Era (2011-2015)")

# par(mfrow = c(1,3)) # mfrow failed to load for loop
for (f in DataMerged) {
  df <- get(f)
  file_path <- file.path("pictures/", paste0(f, ".png"))
  
  mainn <- paste("\n", DataTitle[[f]])
  
  png(file = file_path, width = 10, height = 10, units = "cm", res = 800)
  plot(df[, "Incid"],
       breaks = breaks_seq,
       pal = colorRampPalette(c("white", "maroon"))(199), # fill with (length breaks - 1)
       main = mainn)
  dev.off()
}


# Basic table for region-ageGroup counts
regAgeGroup_details <- dat_G %>% 
  dplyr::mutate(
    Meningitis = case_when(
      MeningitisFlag == "Y" ~ 1,
      MeningitisFlag == "N" ~ 0,
      TRUE ~ NA_real_ # basically numeric
    ),
    Death = case_when(
      X30daydeath == "D" ~ 1,
      is.na(X30daydeath) ~ 0
      # TRUE ~ NA_real_
    )
  ) %>% 
  dplyr::group_by(current.region.name, ageGroup7) %>% 
  dplyr::summarise("Serotype 1 Case" = n(),
                   "Meningitis" = sum(Meningitis),
                   "30 Day Death" = sum(Death)) %>% 
  dplyr::ungroup()

reg_byRegionOnly <- regAgeGroup_details %>% 
  dplyr::group_by(current.region.name) %>% 
  dplyr::summarise("Serotype 1 Case" = sum(`Serotype 1 Case`),
                   "Meningitis" = sum(Meningitis),
                   "30 Day Death" = sum(`30 Day Death`)) %>% 
  dplyr::ungroup()


# Facet wrap region, year and ageGroup7 ########################################
ageGroup7_region <- dat_G %>% 
  dplyr::group_by(current.region.name, year, ageGroup7) %>% 
  dplyr::summarise(counts = n()) %>% 
  dplyr::ungroup()

pop_ageGroup7_region <- pop_l %>% 
  dplyr::group_by(Region, Year, ageGroup7) %>% 
  dplyr::summarise(PopSize = sum(PopSize)) %>% 
  dplyr::ungroup()

all_ageGroup7_region <- merge(ageGroup7_region, pop_ageGroup7_region,
                              by.x = c("current.region.name", "year","ageGroup7"),
                              by.y = c("Region", "Year", "ageGroup7")) %>%
  dplyr::mutate(Conf_Int = epitools::binom.exact(counts, PopSize),
                incid_Ser1 = Conf_Int$proportion) # per-100,000 population


png("pictures/incidence_7ageGroups_region.png", width = 24, height = 12, unit = "cm", res = 1200)
ggplot(all_ageGroup7_region, aes(x = year, y = Conf_Int$proportion*100000, group = ageGroup7,
                                 color = ageGroup7)) +
  geom_line(size = 0.7) +
  geom_errorbar(aes(ymin = Conf_Int$lower*100000, ymax = Conf_Int$upper*100000), # It doesn't matter whether I add the CI or not because the Pop data is quite huge, I suppose (?)
                width = .1) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(col_map),
                     name = "Demographic",
                     breaks = c("<2", "2-4", "5-14", "15-30", "31-44", "45-64", "65+", "Unknown", "PCV7", "PCV13"),
                     labels = c("<2", "2-4", "5-14", "15-30", "31-44", "45-64", "65+", "Unknown", "PCV7", "PCV13")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  # scale_linetype_manual(values = c(vacc_map),
  #                       name = "Vaccine",
  #                       labels = c("PCV7", "PCV13")) +
  # geom_label(aes(x = 2006, y = 2.5, label = "PCV7"),
  #            fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  # geom_label(aes(x = 2011, y = 2.5, label = "PCV13"),
  #            fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("The Incidence of Serotype 1 in England \nby Demographic Groups (per 100,000)") +
  xlab("Year") +
  ylab("Serotype 1 Incidence") +
  facet_wrap(~ current.region.name, scales = "free_y") +
  theme_bw()
dev.off()


# Basic case count data without age structure or regions #######################
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
  dplyr::filter(X30daydeath == "D") %>% 
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
  ylab("Serotype 1 Cases (Aggregated by Week)") +
  theme_bw()
dev.off()

# Additional pic for only Serotype 1 case
ggplot(incidence_weekly, aes(as.Date(weeks))) +
  geom_line(aes(y = cases, colour = "counts_Ser1_weekly")) +
  scale_x_date() +
  scale_color_manual(values = col_imD_weekly,
                     name = "Cases",
                     breaks = c("counts_Ser1_weekly"),
                     labels = c("Serotype 1")
  ) +
  ggtitle("The Counts of Serotype 1 in England") +
  xlab("Year") +
  ylab("Serotype 1 Cases (Aggregated by Week)") +
  theme_bw()

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

# Load MLST result (not included new surveillance 2017-2024)
mlst_results <- dplyr::bind_rows(
  read.table("raw_data/gubbins/n739/mlst_results_GPSC31_n739.csv"),
  read.table("raw_data/gubbins/GPSC2_n7/mlst_results_GPSC2_n7.csv") %>% 
    dplyr::mutate(V3 = as.character(V3))
) %>% 
  dplyr::select(V1, V3) %>% 
  dplyr::rename(MLST_ST = V3)


# Load All Serotype 1: better load the data from PopPUNK result!
# New data is updated from popPUNK output
all_serotype1 <- dplyr::bind_rows(
  read.table("raw_data/gubbins/remove_GPSC31_list.txt") %>% # pre-thesis data
    dplyr::mutate(GPSC = 31),
  read.table("raw_data/gubbins/remove_GPSC2_list.txt") %>% # pre-thesis data
    dplyr::mutate(GPSC = 2),
  read.csv("raw_data/gubbins/outputs_PopPUNK_external_clusters_n20.csv") %>% # post-thesis data
    dplyr::mutate(V1 = paste0(sample, ".fasta")) %>% 
    dplyr::select(V1, GPSC)
) %>% 
  dplyr::mutate(serotype = 1)


# Load filtered samples based on N50 >= 30kb (data should have been updated before!)
stats <- dplyr::bind_rows(
  read.csv("raw_data/gubbins/sanger_stats_compiled_GPSC31.csv"),
  read.csv("raw_data/gubbins/sanger_stats_compiled_GPSC2.csv")
) %>% 
  dplyr::mutate(my_ID = stringr::str_remove(my_ID, "stats_"),
                V1 = paste0(my_ID, ".fasta")) %>% 
  dplyr::rename(stats_sumBase = sum,
                stats_nsumBase = n_sum,
                stats_N50 = N50) %>% 
  dplyr::select(V1, stats_sumBase, stats_nsumBase, stats_N50)

# New choosen data (data should have been updated before!)
link_GPSC31_n712 <- read.table("outputs/genomics/remove_GPSC31_choosen_list.txt") %>% 
  dplyr::mutate(n712_choosen_GPSC31 = 1)

joined_serotype_GPSC <- dat_G %>% 
  dplyr::full_join(stats, by = c("assembly_name" = "V1")) %>% 
  dplyr::full_join(all_serotype1, by = c("assembly_name" = "V1")) %>% 
  dplyr::full_join(mlst_results, by = c("assembly_name" = "V1")) %>% 
  dplyr::full_join(link_GPSC31_n712, by = c("assembly_name" = "V1"))


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

joined_AMR <- joined_serotype_GPSC %>% 
  dplyr::full_join(resistance_smx_tmp, by = c("assembly_name" = "isolate_id")) %>% 
  dplyr::mutate(ngsid = substr(assembly_name, 1, 8)) # correction for ngsid including those that sequenced but have no EpiData

# load joined_clades first
AMR_summary <- joined_clades %>% 
  dplyr::select(GPSC, MLST_ST, resistance_smx, resistance_tmp) %>% 
  dplyr::group_by(GPSC, MLST_ST, resistance_smx) %>% 
  dplyr::mutate(smx = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(GPSC, MLST_ST, resistance_tmp) %>% 
  dplyr::mutate(tmp = n()) %>% 
  dplyr::ungroup() %>% 
  distinct() %>% 
  view()

write.csv(AMR_summary, "raw_data/AMR.csv")

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
  dplyr::full_join(joined_AMR, clade_assignment_df, by = c("assembly_name" = "V1"))

# Clade viz!
freq_clades <- joined_clades %>%
  dplyr::filter(n703_choosen_GPSC31 == 1) %>% 
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

clades_map <- c("Serotype 1 Case" = "gray75",
                "clade1" = "steelblue",
                "clade2" = "darkgreen",
                "clade3" = "red",
                # Vaccination
                "PCV7" = "gray70",
                "PCV13" = "gray20")

png("pictures/GPSC31_clades.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(freq_clades, aes(x = year, y = Frequency, group = clade,
                        color = clade)) +
  geom_line(size = 1.5) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_colour_manual(values = c(clades_map),
                      name = "GPSC31 Clades",
                      breaks = c("clade1", "clade2", "clade3"),
                      labels = c("Clade 1", "Clade 2", "Clade 3")) +
  geom_errorbar(aes(ymin = Conf_Int$lower, ymax = Conf_Int$upper),
                width = .05) +
  geom_label(aes(label = Count),
             # fill = "white",
             color = "black",
             nudge_y = 0.001,
             nudge_x = 0.05) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  geom_label(aes(x = 2006, y = 0.95, label = "PCV7"),
             fill = "white", color = "black") +
  geom_label(aes(x = 2011, y = 0.95, label = "PCV13"),
             fill = "white", color = "black") +
  ggtitle("The Frequency of GPSC31 Clades in England") +
  xlab("Year") +
  ylab("Frequency") +
  theme_bw()
dev.off()

max_freq_region <- freq_clades %>% 
  dplyr::group_by(year, clade) %>%
  dplyr::arrange(desc(Frequency)) %>%
  dplyr::slice(1)

# Edit name & manually add the label later
max_freq_region <- max_freq_region %>% 
  dplyr::mutate(current.region.name = case_when(
    year == 2014 ~ "0",
    TRUE ~ current.region.name
  ))

png("pictures/GPSC31_clades_with_maxregion.png", width = 25, height = 12, unit = "cm", res = 1200)
ggplot(max_freq_region, aes(x = year, y = Frequency, group = clade,
                            color = clade)) +
  geom_line(size = 1.5) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(clades_map),
                     name = "GPSC31 Clades",
                     breaks = c("clade1", "clade2", "clade3"),
                     labels = c("Clade 1", "Clade 2", "Clade 3")) +
  geom_errorbar(aes(ymin = Conf_Int$lower, ymax = Conf_Int$upper),
                width = .05) +
  geom_label(aes(label = paste0(Count, " (", current.region.name, ")")),
             # fill = "white",
             size = 2.3,
             color = "black",
             # nudge_y = 0.001,
             # nudge_x = 0.25,
             position=position_jitter()) +
  # Mannually add label for 2014 (overlapped data):
  geom_label(label="4 (South West)\n4 (East Midlands)", 
             x=2014, y=0.0689,
             label.padding = unit(0.55, "lines"),
             size = 2.3,
             color = "black",
             nudge_y = 0.001,
             nudge_x = 0.55
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  geom_label(aes(x = 2006, y = 0.9, label = "PCV7"),
             fill = "white", color = "black") +
  geom_label(aes(x = 2011, y = 0.9, label = "PCV13"),
             fill = "white", color = "black") +
  ggtitle("The Frequency of GPSC31 Clades in England") +
  xlab("Year") +
  ylab("Frequency") +
  theme_bw()
dev.off()

# Clade viz for each region!
freq_clades_region <- joined_clades %>%
  dplyr::filter(n703_choosen_GPSC31 == 1) %>% 
  dplyr::group_by(year, current.region.name) %>%
  dplyr::mutate(Sample_size = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(clade, year, current.region.name) %>%
  dplyr::mutate(Frequency = n()/Sample_size) %>%
  dplyr::mutate(Count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(clade),
                !is.na(year),
                !is.na(current.region.name)) %>% # Possibly NA in year is 2015
  # dplyr::mutate(Conf_Int = epitools::binom.exact(Count, Sample_size),
  #               Prop = Conf_Int$proportion) %>% # per-100,000 population
  dplyr::select(clade,year,current.region.name,Frequency, Count, Sample_size) %>% #, Conf_Int, current.region.name) %>%
  dplyr::distinct()

png("pictures/GPSC31_clades_region.png", width = 24, height = 12, unit = "cm", res = 1200)
ggplot(freq_clades_region, aes(x = year, y = Frequency, group = clade,
                               color = clade)) +
  geom_line(size = 1.5) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_colour_manual(values = c(clades_map),
                      name = "GPSC31 Clades",
                      breaks = c("clade1", "clade2", "clade3", "PCV7", "PCV13"),
                      labels = c("Clade 1", "Clade 2", "Clade 3", "PCV7", "PCV13")) +
  # geom_errorbar(aes(ymin = Conf_Int$lower, ymax = Conf_Int$upper),
  #               width = .05) +
  geom_label(aes(label = Count),
             # fill = "white",
             color = "black",
             size = 2,
             nudge_y = 0.001,
             nudge_x = 0.05) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  facet_wrap(~ current.region.name, scales = "free") +
  ggtitle("The Frequency of GPSC31 Clades in England") +
  xlab("Year") +
  ylab("Frequency") +
  theme_bw()
dev.off()

# I'm just kinda curious (2)...
aggregated_clades <- joined_clades %>% 
  dplyr::mutate(
    Earliest.specimen.date = as.Date(Earliest.specimen.date),
    weeks = cut(Earliest.specimen.date, breaks="week"),
    months = cut(Earliest.specimen.date, breaks="month")
  ) %>% 
  # dplyr::group_by(year, current.region.name) %>% 
  # dplyr::mutate("Serotype 1 Case" = n()) %>% 
  # dplyr::ungroup() %>% 
  dplyr::filter(!is.na(clade)) %>% 
  dplyr::group_by(year, current.region.name) %>% 
  dplyr::mutate(Sample_size = n()) %>% 
  dplyr::group_by(clade, year, current.region.name) %>%
  dplyr::mutate(n_per_clade = n(),
                Frequency_clades = n()/Sample_size) %>%
  dplyr::ungroup() %>% 
  dplyr::select(year, current.region.name,
                clade, n_per_clade, Sample_size, Frequency_clades) %>%
  dplyr::filter(!is.na(current.region.name)) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(Conf_Int = epitools::binom.exact(n_per_clade, Sample_size))

png("pictures/GPSC31_clades_FREQ_year_region.png", width = 22, height = 14, unit = "cm", res = 1200)
ggplot(aggregated_clades, aes(x = year, y = Frequency_clades, group = clade, colour = clade)) +
  # geom_area(stat = "bin", fill = "gray75")+
  geom_line(size = 1.5) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(clades_map),
                     name = "GPSC31 Clades",
                     breaks = c("Serotype 1 Case", "clade1", "clade2", "clade3", "PCV7", "PCV13"),
                     labels = c("Serotype 1 Case", "Clade 1", "Clade 2", "Clade 3", "PCV7\n(2006)\n", "PCV13\n(2011)")) +
  geom_errorbar(aes(ymin = Conf_Int$lower, ymax = Conf_Int$upper),
                width = .02) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  scale_y_continuous(trans = "log1p") +
  facet_wrap(~ current.region.name, scales = "free") +
  labs(title = "The Frequency of GPSC31 Clades per Year in England \nby Regions (per 100,000)",
       x = "Year",
       y = "Frequency per Year, Region",
       fill = "Clades of GPSC31") +
  theme_bw()
dev.off()


## 6. Data Preparation for Microreact ##########################################
# Load dat_G first.
tre <- BactDating::loadGubbins("raw_data/gubbins/n739/n739_")

tre_names <- as.data.frame(tre$tip.label)
tre_names$ID_contigs <- substr(tre$tip.label, 1, 8)
tre_names <- 
  dplyr::full_join(tre_names, joined_clades, by = c("ID_contigs" = "ngsid"))

# Write cleaned EpiData plus genomic analysis
write.csv(joined_clades, "raw_data/serotype1_UKHSA_imperial_date_age_region_MOLIS_sequenced_cleaned_joined.csv", row.names = FALSE)

write.csv(tre_names, "raw_data/gubbins/n739/phandango_microreact_check/microreact_tre_names.csv", row.names = FALSE)
write.csv(tre_names, "raw_data/gubbins/n703/phandango_microreact_check/microreact_tre_names.csv", row.names = FALSE)
write.csv(tre_names, "raw_data/gubbins/GPSC2_n7/phandango_microreact_check/microreact_tre_names.csv", row.names = FALSE)


## 7. Other Data Viz ###########################################################
# Basic table for region-ageGroup counts
joined_clades <- read.csv("raw_data/serotype1_UKHSA_imperial_date_age_region_MOLIS_sequenced_cleaned_joined.csv")

aggregated_data <- joined_clades %>% 
  dplyr::mutate(
    Meningitis = case_when(
      MeningitisFlag == "Y" ~ 1,
      MeningitisFlag == "N" ~ 0,
      TRUE ~ NA_real_ # basically numeric
    ),
    Death = case_when(
      X30daydeath == "D" ~ 1,
      is.na(X30daydeath) ~ 0
      # TRUE ~ NA_real_
    ),
  ) %>% 
  dplyr::group_by(Earliest.specimen.date, year, current.region.name, ageGroup7) %>% 
  dplyr::mutate("Serotype 1 Case" = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(Earliest.specimen.date, year, current.region.name, ageGroup7) %>% 
  dplyr::mutate("Sequenced" = sum(serotype, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(Earliest.specimen.date, year, current.region.name, ageGroup7) %>% 
  dplyr::mutate("Meningitis" = sum(Meningitis, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(Earliest.specimen.date, year, current.region.name, ageGroup7) %>% 
  dplyr::mutate("30 Day Death" = sum(Death, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Earliest.specimen.date, year, current.region.name, ageGroup7,
                `Serotype 1 Case`, Sequenced, Meningitis, `30 Day Death`) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(Earliest.specimen.date = as.Date(Earliest.specimen.date),
                weeks = cut(Earliest.specimen.date, breaks="week"))

# Some tables related to sequenced data
table_region <- aggregated_data %>% 
  dplyr::group_by(current.region.name) %>% 
  dplyr::summarise("Serotype 1 Case" = sum(`Serotype 1 Case`),
                   "Sequenced" = sum(Sequenced), # coz' serotype is 1
                   "Meningitis" = sum(Meningitis),
                   "30 Day Death" = sum(`30 Day Death`)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(Freq_sequenced = format((Sequenced/`Serotype 1 Case`*100), digits = 5),
                "Sequenced Samples" = paste0(Freq_sequenced, "%", " =(", Sequenced, "/", `Serotype 1 Case`, ")")) %>% 
  dplyr::select(current.region.name, `Sequenced Samples`) %>% 
  view()

table_ageGroup <- aggregated_data %>% 
  dplyr::group_by(ageGroup7) %>% 
  dplyr::summarise("Serotype 1 Case" = sum(`Serotype 1 Case`),
                   "Sequenced" = sum(Sequenced), # coz' serotype is 1
                   "Meningitis" = sum(Meningitis),
                   "30 Day Death" = sum(`30 Day Death`)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(Freq_sequenced = format((Sequenced/`Serotype 1 Case`*100), digits = 5),
                "Sequenced Samples" = paste0(Freq_sequenced, "%", " =(", Sequenced, "/", `Serotype 1 Case`, ")")) %>% 
  dplyr::select(ageGroup7, `Sequenced Samples`) %>% 
  view()


cases_year <- aggregated_data %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise("Serotype 1 Case" = sum(`Serotype 1 Case`),
                   "Sequenced" = sum(Sequenced), # coz' serotype is 1
                   "Meningitis" = sum(Meningitis),
                   "30 Day Death" = sum(`30 Day Death`)) %>% 
  dplyr::ungroup() %>% 
  reshape2::melt(id.vars = c("year"),
                 measure.vars = c("Serotype 1 Case", "Sequenced", "Meningitis", "30 Day Death"))

png("pictures/counts_allages_withSeq.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(cases_year, aes(x = year, y = value, group = variable,
                       color = variable)) +
  geom_line(size = 1.5) +
  geom_vline(data = vaccine_UK, aes(xintercept = year),
             linetype = "dashed") +
  scale_color_manual(values = c(col_map),
                     name = "Cases",
                     breaks = c("Serotype 1 Case", "Sequenced", "Meningitis", "30 Day Death"),
                     labels = c("Serotype 1 Cases", "Sequenced", "Meningitis", "30 Day Death")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  geom_label(aes(x = 2006, y = 150, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = 2011, y = 150, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("The Counts of Serotype 1 in England") +
  xlab("Year") +
  ylab("Serotype 1 Cases") +
  theme_bw()
dev.off()

cases_weeks <- aggregated_data %>% 
  dplyr::group_by(weeks) %>% 
  dplyr::summarise("Serotype 1 Case" = sum(`Serotype 1 Case`),
                   "Sequenced" = sum(Sequenced), # coz' serotype is 1
                   "Meningitis" = sum(Meningitis),
                   "30 Day Death" = sum(`30 Day Death`)) %>% 
  dplyr::ungroup() %>% 
  reshape2::melt(id.vars = c("weeks"),
                 measure.vars = c("Serotype 1 Case", "Sequenced", "Meningitis", "30 Day Death")) %>% 
  dplyr::filter(variable != "Meningitis" & variable != "30 Day Death")

png("pictures/counts_allages_withSeq_weeks.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(cases_weeks, aes(x = as.Date(weeks), y = value, group = variable,
                        color = variable)) +
  geom_line(size = 0.75) +
  geom_vline(data = vaccine_UK, aes(xintercept = as.Date(c("2006-09-04", "2010-10-01"))),
             linetype = "dashed") +
  scale_color_manual(values = c(col_map),
                     name = "Cases",
                     breaks = c("Serotype 1 Case", "Sequenced", "Meningitis", "30 Day Death"),
                     labels = c("Serotype 1 Cases", "Sequenced", "Meningitis", "30 Day Death")
  ) +
  # scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  geom_label(aes(x = as.Date("2006-09-04"), y = 41, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = as.Date("2010-10-01"), y = 41, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("Serotype 1 Cases in England (Aggregated by Week)") +
  xlab("Time") +
  ylab("Serotype 1 Cases") +
  theme_bw()
dev.off()

# I'm just being curious...
aggregated_clades <- joined_clades %>% 
  dplyr::mutate(
    Clade1 = case_when(
      clade == "clade1" ~ 1,
      clade == "clade2" ~ 0,
      clade == "clade3" ~ 0,
      TRUE ~ NA_real_ # basically numeric
    ),
    Clade2 = case_when(
      clade == "clade1" ~ 0,
      clade == "clade2" ~ 1,
      clade == "clade3" ~ 0,
      TRUE ~ NA_real_ # basically numeric
    ),
    Clade3 = case_when(
      clade == "clade1" ~ 0,
      clade == "clade2" ~ 0,
      clade == "clade3" ~ 1,
      TRUE ~ NA_real_ # basically numeric
    )) %>% 
  dplyr::group_by(Earliest.specimen.date, year, current.region.name, ageGroup7) %>% 
  dplyr::mutate("Serotype 1 Case" = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(Earliest.specimen.date, year, current.region.name, ageGroup7) %>% 
  dplyr::mutate(clade1 = sum(Clade1, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(Earliest.specimen.date, year, current.region.name, ageGroup7) %>% 
  dplyr::mutate(clade2 = sum(Clade2, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(Earliest.specimen.date, year, current.region.name, ageGroup7) %>% 
  dplyr::mutate(clade3 = sum(Clade3, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Earliest.specimen.date, year, current.region.name, ageGroup7,
                `Serotype 1 Case`, clade1, clade2, clade3) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(Earliest.specimen.date = as.Date(Earliest.specimen.date),
                weeks = cut(Earliest.specimen.date, breaks="week"),
                months = cut(Earliest.specimen.date, breaks="month"))


clades_year <- aggregated_clades %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise("Serotype 1 Case" = sum(`Serotype 1 Case`),
                   clade1 = sum(clade1), # coz' serotype is 1
                   clade2 = sum(clade2),
                   clade3 = sum(clade3)) %>% 
  dplyr::ungroup() %>% 
  reshape2::melt(id.vars = c("year"),
                 measure.vars = c("Serotype 1 Case", "clade1", "clade2", "clade3"))

png("pictures/counts_allages_withSeq_clades.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(clades_year, aes(x = year, y = value, group = variable,
                        color = variable)) +
  geom_line(size = 1.5) +
  geom_vline(data = vaccine_UK, aes(xintercept = year),
             linetype = "dashed") +
  scale_color_manual(values = c(clades_map),
                     name = "Cases (GPSC31)",
                     breaks = c("Serotype 1 Case", "clade1", "clade2", "clade3"),
                     labels = c("Serotype 1 Cases", "Clade 1", "Clade 2", "Clade 3")
  ) +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  geom_label(aes(x = 2006, y = 150, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = 2011, y = 150, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("Case Counts of IPD Caused by Serotype 1\nin England") +
  xlab("Year") +
  ylab("Serotype 1 Cases") +
  theme_bw()
dev.off()


clades_months <- aggregated_clades %>% 
  dplyr::group_by(months) %>% 
  dplyr::summarise("Serotype 1 Case" = sum(`Serotype 1 Case`),
                   clade1 = sum(clade1), # coz' serotype is 1
                   clade2 = sum(clade2),
                   clade3 = sum(clade3)) %>% 
  dplyr::ungroup() %>% 
  reshape2::melt(id.vars = c("months"),
                 measure.vars = c("Serotype 1 Case", "clade1", "clade2", "clade3"))

png("pictures/counts_allages_withSeq_months_clades.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(clades_months, aes(x = as.Date(months), y = value, group = variable,
                          color = variable)) +
  geom_line(size = 0.75) +
  geom_vline(data = vaccine_UK, aes(xintercept = as.Date(c("2006-09-04", "2010-10-01"))),
             linetype = "dashed") +
  scale_color_manual(values = c(clades_map),
                     name = "Cases (GPSC31)",
                     breaks = c("Serotype 1 Case", "clade1", "clade2", "clade3"),
                     labels = c("Serotype 1 Cases", "Clade 1", "Clade 2", "Clade 3")
  ) +
  # scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  geom_label(aes(x = as.Date("2006-09-04"), y = 90, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = as.Date("2010-10-01"), y = 90, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("Serotype 1 Cases in England (Aggregated by Month)") +
  xlab("Time") +
  ylab("Serotype 1 Cases") +
  theme_bw()
dev.off()


clades_weeks <- aggregated_clades %>% 
  dplyr::group_by(weeks) %>% 
  dplyr::summarise("Serotype 1 Case" = sum(`Serotype 1 Case`),
                   clade1 = sum(clade1), # coz' serotype is 1
                   clade2 = sum(clade2),
                   clade3 = sum(clade3)) %>% 
  dplyr::ungroup() %>% 
  reshape2::melt(id.vars = c("weeks"),
                 measure.vars = c("Serotype 1 Case", "clade1", "clade2", "clade3"))

png("pictures/counts_allages_withSeq_weeks_clades.png", width = 17, height = 12, unit = "cm", res = 1200)
ggplot(clades_weeks, aes(x = as.Date(weeks), y = value, group = variable,
                         color = variable)) +
  geom_line(size = 0.75) +
  geom_vline(data = vaccine_UK, aes(xintercept = as.Date(c("2006-09-04", "2010-10-01"))),
             linetype = "dashed") +
  scale_color_manual(values = c(clades_map),
                     name = "Cases (GPSC31)",
                     breaks = c("Serotype 1 Case", "clade1", "clade2", "clade3"),
                     labels = c("Serotype 1 Cases", "Clade 1", "Clade 2", "Clade 3")
  ) +
  # scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  geom_label(aes(x = as.Date("2006-09-04"), y = 41, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = as.Date("2010-10-01"), y = 41, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("Serotype 1 Cases in England (Aggregated by Week)") +
  xlab("Time") +
  ylab("Serotype 1 Cases") +
  theme_bw()
dev.off()


## 8. Some Stats Analysis ######################################################
joined_clades <- read.csv("raw_data/serotype1_UKHSA_imperial_date_age_region_MOLIS_sequenced_cleaned_joined.csv")
# Meningitis-death odd ratio for all data
# (Death+ & Meningitis+)
all_Dm <- joined_clades %>% 
  filter(MeningitisFlag == "Y" & !is.na(X30daydeath)) %>% # not NA equals to "D"
  nrow()

# (Death+ & Meningitis-)
all_Dmmin <- joined_clades %>% 
  filter(MeningitisFlag == "N" & !is.na(X30daydeath)) %>% # not NA equals to "D"
  nrow()

# (Death- & Meningitis+)
all_Dminm <- joined_clades %>% 
  filter(MeningitisFlag == "Y" & is.na(X30daydeath)) %>% 
  nrow()

# (Death- & Meningitis-)
all_Dminmmin <- joined_clades %>% 
  filter(MeningitisFlag == "N" & is.na(X30daydeath)) %>% 
  nrow()

# Matrix (D+M+, D-M+, D+M-, D-M-)
Outcome <- c("Death", "Alive")
Meningitis <- c("Meningitis", "Not Meningitis")


Input_mD <- matrix(c(all_Dm, all_Dminm,
                     all_Dmmin, all_Dminmmin),
                   nrow = 2, ncol = 2, byrow = T)
dimnames(Input_mD) <- list('Meningitis'=Meningitis, 'Outcome'=Outcome)

All_OD <- epitools::oddsratio(Input_mD)
All_OD

# GPSCs vs. regions 
chisq.test(table(joined_clades$current.region.name, joined_clades$GPSC))
# data:  table(joined_clades$current.region.name, joined_clades$GPSC)
# X-squared = 13.979, df = 8, p-value = 0.08231

# GPSC31 Clade vs. regions 
chi_region <- chisq.test(table(joined_clades$current.region.name, joined_clades$clade))
print(chi_region)
prop.table(chi_region)
# data:  table(joined_clades$current.region.name, joined_clades$clade)
# X-squared = 48.53, df = 16, p-value = 3.92e-05
# Post-hoc
std_resid_region <- chi_region$stdres
print(std_resid_region)

# Whether a mosaic plot is the best way to viz the distribution of clades (not related to year)
mosaicplot(table(joined_clades$current.region.name, joined_clades$clade),
           main = "Mosaic Plot of Region and GPSC31 Clades", shade = TRUE, legend = TRUE)

# Create barplot of region vs. clades:
filtered_NA <- joined_clades %>% 
  dplyr::filter(!is.na(current.region.name),
                !is.na(clade),
                !is.na(year))
filtered_NA_counts <- filtered_NA %>%
  dplyr::group_by(year, current.region.name, clade) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup()
clades_combined <- merge(filtered_NA_counts, pop_year,
                         by.x = c("year"),
                         by.y = c("Year")) %>%
  # Epitools take long time to compile!
  # dplyr::mutate(Conf_Int = epitools::binom.exact(count, PopSize),
  #               incid_clades_year = Conf_Int$proportion) # per-100,000 population
  dplyr::mutate(incid_clades_year = count/PopSize)

# Counts (filtered)
png("pictures/GPSC31_counts_by_region.png", width = 25, height = 17, unit = "cm", res = 1200)
ggplot(filtered_NA_counts, aes(x = year, y = count, colour = clade)) +
  geom_line(size = 1.5) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(clades_map),
                     name = "GPSC31 Clades",
                     breaks = c("clade1", "clade2", "clade3", "PCV7", "PCV13"),
                     labels = c("Clade 1", "Clade 2", "Clade 3", "PCV7", "PCV13")) +
  # annotate("segment", x = 2006, xend = 2006, y = -Inf, yend = Inf, linetype = "dashed") +
  # annotate("segment", x = 2011, xend = 2011, y = -Inf, yend = Inf, linetype = "dashed") +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  facet_wrap(~ current.region.name, scales = "free") +
  labs(title = "Count of GPSC31 Clades by Region",
       x = "Year",
       y = "Counts",
       fill = "Clades of GPSC31") +
  theme_bw()
dev.off()

# Clade frequency per year only (filtered)
png("pictures/GPSC31_freq_by_year_only.png", width = 25, height = 17, unit = "cm", res = 1200)
ggplot(clades_combined, aes(x = year, y = incid_clades_year*100000, colour = clade)) +
  geom_line(size = 1.5) +
  # Epitools takes long time to compile!
  # geom_errorbar(aes(ymin = Conf_Int$lower*100000, ymax = Conf_Int$upper*100000), # It doesn't matter whether I add the CI or not because the Pop data is quite huge, I suppose (?)
  #               width = .1) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(clades_map),
                     name = "GPSC31 Clades",
                     breaks = c("clade1", "clade2", "clade3", "PCV7", "PCV13"),
                     labels = c("Clade 1", "Clade 2", "Clade 3", "PCV7", "PCV13")) +
  # annotate("segment", x = 2006, xend = 2006, y = -Inf, yend = Inf, linetype = "dashed") +
  # annotate("segment", x = 2011, xend = 2011, y = -Inf, yend = Inf, linetype = "dashed") +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  facet_wrap(~ current.region.name, scales = "free") +
  labs(title = "The Frequency of GPSC31 Clades per Year in England \nby Regions (per 100,000)",
       x = "Year",
       y = "Frequency per Year",
       fill = "Clades of GPSC31") +
  theme_bw()
dev.off()

# Clade frequency per year, region (filtered)
pop_reg <- pop_l %>% 
  dplyr::group_by(Year, Region) %>% 
  dplyr::summarise(PopSize = sum(PopSize), .groups="keep") %>% 
  dplyr::ungroup()
clades_combined_region <- merge(filtered_NA_counts, pop_reg,
                                by.x = c("year", "current.region.name"),
                                by.y = c("Year", "Region")) %>%
  # Epitools takes long time to compile!
  # dplyr::mutate(Conf_Int = epitools::binom.exact(count, PopSize),
  #               incid_clades_year_reg = Conf_Int$proportion) # per-100,000 population
  dplyr::mutate(incid_clades_year_reg = count/PopSize)

# Population Growth CheckPoint!
# ggplot(pop_reg, aes(x = Year, y = PopSize, colour = Region)) +
#   geom_line() +
#   facet_wrap(~ Region, scales = "free")
# pop_test <- pop_reg %>% 
#   dplyr::group_by(Year) %>% 
#   dplyr::summarise(PopSize = sum(PopSize)) %>% 
#   dplyr::ungroup() %>% 
# pop_year # compare

png("pictures/GPSC31_freq_by_year_region.png", width = 25, height = 17, unit = "cm", res = 1200)
ggplot(clades_combined_region, aes(x = year, y = incid_clades_year_reg*100000, colour = clade)) +
  geom_line(size = 1.5) +
  # Epitools take long time to compile!
  # geom_errorbar(aes(ymin = Conf_Int$lower*100000, ymax = Conf_Int$upper*100000), # It doesn't matter whether I add the CI or not because the Pop data is quite huge, I suppose (?)
  #               width = .1) +
  geom_vline(data = vaccine_UK, aes(xintercept = year,
                                    colour = vaccine),
             linetype = "dashed") +
  scale_color_manual(values = c(clades_map),
                     name = "GPSC31 Clades",
                     breaks = c("clade1", "clade2", "clade3", "PCV7", "PCV13"),
                     labels = c("Clade 1", "Clade 2", "Clade 3", "PCV7", "PCV13")) +
  # annotate("segment", x = 2006, xend = 2006, y = -Inf, yend = Inf, linetype = "dashed") +
  # annotate("segment", x = 2011, xend = 2011, y = -Inf, yend = Inf, linetype = "dashed") +
  scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  facet_wrap(~ current.region.name, scales = "free") +
  labs(title = "The Frequency of GPSC31 in England \nby Regions (per 100,000)",
       x = "Year",
       y = "Frequency per Year & Region",
       fill = "Clades of GPSC31") +
  theme_bw()
dev.off()

# GPSCs vs. vaccination era
joined_clades$year <- as.factor(joined_clades$year)

anova_year_GPSCs <- aov(year ~ GPSC, data = joined_clades)
summary(anova_year_GPSCs)
# Post-Hoc
TukeyHSD(anova_year_GPSCs)


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

# Other non-statistically significant results: #################################
# GPSCs vs. ages
anova_ages_GPSCs <- aov(ageLabel ~ GPSC, data = joined_clades)
summary(anova_ages_GPSCs)
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
# No statistically significant differences in the mean ages across clades.

# GPSCs vs. age groups
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
           main = "Mosaic Plot", shade = TRUE, legend = TRUE)

# GPSC31 Clade vs. ages
chisq.test(table(joined_clades$ageGroup2, joined_clades$clade))
chisq.test(table(joined_clades$ageGroup, joined_clades$clade))
chisq.test(table(joined_clades$ageGroup7, joined_clades$clade))

# GPSC31 vs. meningitis
chisq.test(table(joined_clades$MeningitisFlag, joined_clades$clade))
# GPSC31 vs. death
chisq.test(table(joined_clades$X30daydeath, joined_clades$clade))

