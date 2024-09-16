library(tidyverse)

# Join mleSky result, model result, epiData & GPSC31 genomics data
mlesky_table <- dplyr::left_join(
  mlesky_res <- read.csv("outputs/genomics/GPSC31_model_2_adaptation_true_rep_10000_mlesky.csv") %>% 
    dplyr::mutate(new_date = format(lubridate::date_decimal(Year), "%Y-%m-%d")) %>% # instead of "%Y-%m-%d"
    dplyr::mutate(new_date = as.Date(new_date)) %>%
    dplyr::mutate(new_yearMonth = format(lubridate::date_decimal(Year), "%Y-%m")), # year-month as characters
  joined_clades <- read.csv("raw_data/serotype1_UKHSA_imperial_date_age_region_MOLIS_sequenced_cleaned_joined.csv") %>% 
    dplyr::filter(n703_choosen_GPSC31 == 1) %>% 
    dplyr::mutate(Earliest.specimen.date = as.Date(Earliest.specimen.date)) %>% 
    dplyr::mutate(new_yearMonth = format(Earliest.specimen.date, "%Y-%m")) %>% # year-month as characters
    dplyr::mutate(new_yearMonth = as.Date(paste0(new_yearMonth, "-01"))) %>% 
    dplyr::group_by(new_yearMonth) %>% # or year or month
    dplyr::summarise(data_counts = n()) %>% # because data here is daily data per-IDs
    dplyr::ungroup(),
  by = "new_yearMonth") %>% #c("new_year" = "year")) #%>%
  dplyr::mutate(new_yearMonth = as.Date(paste0(new_yearMonth, "-01")))

vaccine_UK <- data.frame(
  date = c("2006-09-04", "2010-04-01"),
  vaccine = c("PCV7", "PCV13")) # PCV7 start 4 September 2006, PCV13 from April 2010
# SOURCE: https://www.gov.uk/government/publications/pneumococcal-disease-caused-by-strains-in-prevenar-13-and-not-in-prevenar-7-vaccine/pneumococcal-disease-infections-caused-by-serotypes-in-prevenar-13-and-not-in-prevenar-7

col_mlesky <- c(
  "Ne" = "maroon",
  "Model" = "deepskyblue3",
  "Data" = "steelblue",
  # Vaccination
  "PCV7" = "gray70",
  "PCV13" = "gray20",
  # Upper & lower ranges
  "lo" = "maroon",
  "hi" = "maroon")


# Expected winter at the end of December
winter <- mlesky_table %>%
  filter(format(new_yearMonth, "%m") == "12") %>%
  mutate(start_of_december = as.Date(paste0(format(new_yearMonth, "%Y"), "-12-15")),
         end_of_december = as.Date(paste0(format(new_yearMonth, "%Y"), "-12-31")))


# Data viz!
# 1. Focused on month-year:
png("pictures/genomics/mleSky_month_year_withRanges.png", width = 30, height = 12, unit = "cm", res = 1200)
ggplot(data = mlesky_table) +
  geom_line(linewidth = 1, group = 1, #data = mlesky_res,
            aes(x = new_yearMonth, y = Ne, colour = "Ne")) +
  geom_line(linetype = "dashed", aes(x = new_yearMonth, y = Lower_range,
                                     colour = "lo")) + # either Lower_bound or Lower_range
  geom_line(linetype = "dashed", aes(x = new_yearMonth, y = Upper_range,
                                     colour = "hi")) + # either Upper_bound or Upper_range
  geom_bar(stat = "identity", #data = all_day,
           aes(x = new_yearMonth, y = data_counts, fill = "Data")) +
  geom_vline(data = vaccine_UK, aes(xintercept = as.Date(date, format = "%Y-%m-%d"), colour = vaccine),
             linetype = "dashed") +
  scale_fill_manual(values = c(col_mlesky),
                    name = "Bar",
                    breaks = c("Data"),
                    labels = c("Data")
  ) +
  scale_color_manual(values = c(col_mlesky),
                     name = "Line",
                     breaks = c("Ne", "PCV7", "PCV13", "lo", "hi"),
                     labels = c("Ne", "PCV7", "PCV13", "Lower Range", "Upper Range")
  ) +
  geom_rect(data = winter, 
            aes(xmin = start_of_december, xmax = end_of_december, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.8) +
  # scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  scale_x_date(limits = as.Date(c("2003-06-01", "2015-10-31")),
               date_breaks = "1 year", date_labels = "%Y") +
  ggtitle("GPSC31 Population Size Aggregated Monthly in England") +
  xlab("Year") +
  ylab("Effective Population Size (Ne)") +
  theme_bw()
dev.off()


# 2. Aggregated by year only:
mlesky_table_year <- mlesky_table %>% 
  dplyr::mutate(year = lubridate::year(new_date)) %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise(data_counts = sum(data_counts),
                   Ne = sum(Ne),
                   Lower_range = sum(Lower_range),
                   Upper_range = sum(Upper_range),
                   Lower_bound = sum(Lower_bound),
                   Upper_bound = sum(Upper_bound)) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(year = as.Date(paste0(year, "-01-01")))

png("pictures/genomics/mleSky_year_2start2003_withRanges.png", width = 30, height = 12, unit = "cm", res = 1200)
ggplot(data = mlesky_table_year) +
  geom_line(linewidth = 1, group = 1, #data = mlesky_res,
            aes(x = year, y = Ne, colour = "Ne")) +
  geom_line(linetype = "dashed", aes(x = year, y = Lower_range,
                                     colour = "lo")) + # either Lower_bound or Lower_range
  geom_line(linetype = "dashed", aes(x = year, y = Upper_range,
                                     colour = "hi")) + # either Upper_bound or Upper_range
  geom_bar(stat = "identity", #data = all_day,
           aes(x = year, y = data_counts, fill = "Data")) +
  geom_vline(data = vaccine_UK, aes(xintercept = as.Date(date, format = "%Y-%m-%d"), colour = vaccine),
             linetype = "dashed") +
  scale_fill_manual(values = c(col_mlesky),
                    name = "Bar",
                    breaks = c("Data"),
                    labels = c("Data")
  ) +
  scale_color_manual(values = c(col_mlesky),
                     name = "Line",
                     breaks = c("Ne", "PCV7", "PCV13", "lo", "hi"),
                     labels = c("Ne", "PCV7", "PCV13", "Lower Range", "Upper Range")
  ) +
  # geom_rect(data = winter, 
  #           aes(xmin = start_of_december, xmax = end_of_december, ymin = -Inf, ymax = Inf),
  #           fill = "grey", alpha = 0.8) +
  # scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  scale_x_date(limits = as.Date(c("2003-01-01", "2015-10-31")),
               date_breaks = "10 years", date_labels = "%Y") +
  ggtitle("GPSC31 Population Size Aggregated Annually in England") +
  xlab("Year") +
  ylab("Effective Population Size (Ne)") +
  theme_bw()
dev.off()

# Monthly and annual table for Ne ##############################################
# DO NOT COMBINE DATA INTO ONE HUGE DF!!!
# model result + epiData
model_epi_daily <- read.csv("inputs/daily_joined.csv") %>% # save data for 4745 days only
  dplyr::mutate(date = as.Date("2003-01-01") + days(time - 1),
                value_epiData = value_data) %>% 
  dplyr::select(date, compartment, value_model, value_epiData) %>% # EpiData is available in value_data
  dplyr::filter(compartment %in% c("Asymptomatic", "Diseased")) %>% 
  # Generate weekly, monthly and annual data
  dplyr::mutate(time_start = as.numeric(date - min(date)),
                weekly = ceiling(time_start/7),
                month = lubridate::month(date),
                yearMonth = format(date, "%Y-%m"), # year-month as characters
                yearMonth = as.Date(paste0(yearMonth, "-01")),
                year = lubridate::year(date),
                year = as.Date(paste0(year, "-01-01")))

# mleSky result + genomics (GPSC31) data
mleSky_genomics_daily <- read.csv("outputs/genomics/GPSC31_model_2_adaptation_true_rep_10000_mlesky.csv") %>% 
  dplyr::mutate(new_date = format(lubridate::date_decimal(Year), "%Y-%m-%d")) %>% # instead of "%Y-%m-%d"
  dplyr::mutate(date = as.Date(new_date),
                compartment = "Asymptomatic") %>% # will be fitted to model result
  dplyr::select(date, compartment, Ne) %>% 
  # Generate weekly, monthly and annual data
  dplyr::mutate(time_start = as.numeric(date - min(date)),
                weekly = ceiling(time_start/7),
                month = lubridate::month(date),
                yearMonth = format(date, "%Y-%m"), # year-month as characters
                yearMonth = as.Date(paste0(yearMonth, "-01")),
                year = lubridate::year(date),
                year = as.Date(paste0(year, "-01-01")))

# weekly
combined_weekly <- dplyr::full_join(
  model_epi_daily %>% 
    dplyr::filter(compartment == "Asymptomatic") %>% 
    dplyr::group_by(weekly, date) %>% 
    dplyr::summarise(value_model = sum(value_model))%>% 
    dplyr::select(date, value_model) %>% 
    dplyr::distinct(),
  mleSky_genomics_daily %>% 
    dplyr::filter(compartment == "Asymptomatic") %>% 
    dplyr::group_by(weekly, date) %>% 
    dplyr::summarise(Ne = sum(Ne))%>% 
    dplyr::select(date, Ne) %>% 
    dplyr::distinct(),
  by = "date"
)

# Data viz for asymptomatic cases per-week (mleSky vs. model)
png("pictures/genomics/mleSky_vs_model_asymptomatic_weekly.png", width = 30, height = 12, unit = "cm", res = 1200)
ggplot(data = combined_weekly) +
  geom_line(linewidth = 1, group = 1, #data = mlesky_res,
            aes(x = date, y = Ne, colour = "Ne")) +
  # geom_line(linewidth = 1, group = 1, #data = mlesky_res,
  #           aes(x = date, y = value_model, colour = "Model")) +
  # geom_bar(stat = "identity", #data = all_day,
  #          aes(x = date, y = value_epiData, fill = "Data")) +
  geom_vline(data = vaccine_UK, aes(xintercept = as.Date(date, format = "%Y-%m-%d"), colour = vaccine),
             linetype = "dashed") +
  scale_fill_manual(values = c(col_mlesky),
                    name = "Bar",
                    breaks = c("Data"),
                    labels = c("Data")
  ) +
  scale_color_manual(values = c(col_mlesky),
                     name = "Line",
                     breaks = c("Ne", "Model", "PCV7", "PCV13", "lo", "hi"),
                     labels = c("Ne", "Model", "PCV7", "PCV13", "Lower Range", "Upper Range")
  ) +
  scale_x_date(limits = as.Date(c("2003-01-01", "2015-10-31")),
               date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(trans = "log1p") +
  ggtitle("GPSC31 Population Size Aggregated Annually in England") +
  xlab("Year") +
  ylab("Effective Population Size (Ne)") +
  theme_bw()
dev.off()

# monthly
combined_monthly <- dplyr::full_join(
  model_epi_daily %>% 
    dplyr::filter(compartment == "Asymptomatic") %>% 
    dplyr::group_by(yearMonth) %>% 
    dplyr::summarise(value_model = sum(value_model))%>% 
    dplyr::select(yearMonth, value_model) %>% 
    dplyr::distinct(),
  mleSky_genomics_daily %>% 
    dplyr::filter(compartment == "Asymptomatic") %>% 
    dplyr::group_by(yearMonth) %>% 
    dplyr::summarise(Ne = sum(Ne))%>% 
    dplyr::select(yearMonth, Ne) %>% 
    dplyr::distinct(),
  by = "yearMonth"
)

# Data viz for asymptomatic cases per-month (mleSky vs. model)
png("pictures/genomics/mleSky_vs_model_asymptomatic_monthly.png", width = 30, height = 12, unit = "cm", res = 1200)
ggplot(data = combined_monthly) +
  geom_line(linewidth = 1, group = 1, #data = mlesky_res,
            aes(x = yearMonth, y = Ne, colour = "Ne")) +
  geom_line(linewidth = 1, group = 1, #data = mlesky_res,
            aes(x = yearMonth, y = value_model, colour = "Model")) +
  # geom_bar(stat = "identity", #data = all_day,
  #          aes(x = date, y = value_epiData, fill = "Data")) +
  geom_vline(data = vaccine_UK, aes(xintercept = as.Date(date, format = "%Y-%m-%d"), colour = vaccine),
             linetype = "dashed") +
  scale_fill_manual(values = c(col_mlesky),
                    name = "Bar",
                    breaks = c("Data"),
                    labels = c("Data")
  ) +
  scale_color_manual(values = c(col_mlesky),
                     name = "Line",
                     breaks = c("Ne", "Model", "PCV7", "PCV13", "lo", "hi"),
                     labels = c("Ne", "Model", "PCV7", "PCV13", "Lower Range", "Upper Range")
  ) +
  scale_x_date(limits = as.Date(c("2003-01-01", "2015-10-31")),
               date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(trans = "log1p") +
  ggtitle("GPSC31 Population Size Aggregated Annually in England") +
  xlab("Year") +
  ylab("Effective Population Size (Ne)") +
  theme_bw()
dev.off()

# annually
combined_annually <- dplyr::full_join(
  model_epi_daily %>% 
    dplyr::filter(compartment == "Asymptomatic") %>% 
    dplyr::group_by(year) %>% 
    dplyr::summarise(value_model = sum(value_model))%>% 
    dplyr::select(year, value_model) %>% 
    dplyr::distinct(),
  mleSky_genomics_daily %>% 
    dplyr::filter(compartment == "Asymptomatic") %>% 
    dplyr::group_by(year) %>% 
    dplyr::summarise(Ne = sum(Ne))%>% 
    dplyr::select(year, Ne) %>% 
    dplyr::distinct(),
  by = "year"
)

# annual plot with epiData
epiData_year <- read.csv("raw_data/serotype1_UKHSA_imperial_date_age_region_MOLIS_sequenced_cleaned_joined.csv") %>% 
  dplyr::mutate(Earliest.specimen.date = as.Date(Earliest.specimen.date),
                year = year(Earliest.specimen.date)) %>% 
  dplyr::group_by(year) %>% # or year or month
  dplyr::summarise(epiData_counts = n()) %>% # because data here is daily data per-IDs
  dplyr::mutate(year = as.Date(paste0(year, "-01-01")))
  dplyr::ungroup()

# Data viz for asymptomatic cases per-year (mleSky vs. model)
png("pictures/genomics/mleSky_vs_model_asymptomatic_annually.png", width = 18, height = 12, unit = "cm", res = 1200)
ggplot(data = combined_annually) +
  geom_line(linewidth = 1, group = 1, #data = mlesky_res,
            aes(x = year, y = Ne, colour = "Ne")) +
  geom_line(linewidth = 1, group = 1, #data = mlesky_res,
            aes(x = year, y = value_model, colour = "Model")) +
  geom_bar(stat = "identity", data = epiData_year,
           aes(x = year, y = epiData_counts, fill = "Data")) +
  geom_vline(data = vaccine_UK, aes(xintercept = as.Date(date, format = "%Y-%m-%d"), colour = vaccine),
             linetype = "dashed") +
  scale_fill_manual(values = c(col_mlesky),
                    name = "Bar",
                    breaks = c("Data"),
                    labels = c("Epidemiological \nData")
  ) +
  scale_color_manual(values = c(col_mlesky),
                     name = "Line",
                     breaks = c("Ne", "Model", "PCV7", "PCV13", "lo", "hi"),
                     labels = c("Ne (GPSC31)", "Model", "PCV7", "PCV13", "Lower Range", "Upper Range")
  ) +
  scale_x_date(limits = as.Date(c("2003-01-01", "2015-10-31")),
               date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(trans = "log1p") +
  ggtitle("GPSC31 Population Size Aggregated Annually in England") +
  xlab("Year") +
  ylab("Effective Population Size (Ne)") +
  theme_bw()
dev.off()

























# Additional data
mlesky_table_epiData <- dplyr::left_join(
  mlesky_table,
  joined_clades <- read.csv("raw_data/serotype1_UKHSA_imperial_date_age_region_MOLIS_sequenced_cleaned_joined.csv") %>% 
    # dplyr::filter(n703_choosen_GPSC31 == 1) %>% # NO FILTER
    dplyr::mutate(Earliest.specimen.date = as.Date(Earliest.specimen.date)) %>% 
    dplyr::mutate(new_yearMonth = format(Earliest.specimen.date, "%Y-%m")) %>% # year-month as characters
    dplyr::mutate(new_yearMonth = as.Date(paste0(new_yearMonth, "-01"))) %>% 
    dplyr::group_by(new_yearMonth) %>% # or year or month
    dplyr::summarise(epiData_counts = n()) %>% # because data here is daily data per-IDs
    dplyr::ungroup(),
  by = "new_yearMonth") %>% 
  dplyr::mutate(new_yearMonth = as.Date(paste0(new_yearMonth, "-01"))) %>% 
  dplyr::select(new_yearMonth, Ne, data_counts, epiData_counts)

# Save mleSky monthly data
write.csv(mlesky_table_epiData, "raw_data/mlesky_table_epiData_monthly.csv", row.names = FALSE)

# Annual data
mlesky_table_epiData_annually <- mlesky_table_epiData %>% 
  dplyr::mutate(year = lubridate::year(new_yearMonth)) %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise(data_counts = sum(data_counts),
                   Ne = sum(Ne),
                   # Lower_range = sum(Lower_range),
                   # Upper_range = sum(Upper_range),
                   # Lower_bound = sum(Lower_bound),
                   # Upper_bound = sum(Upper_bound),
                   epiData_counts = sum(epiData_counts)) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(year = as.Date(paste0(year, "-01-01"))) %>% 
  dplyr::select(year, Ne, data_counts, epiData_counts)

# Save mleSky annual data
write.csv(mlesky_table_epiData_annually, "raw_data/mlesky_table_epiData_annually.csv", row.names = FALSE)
