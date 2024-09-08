library(tidyverse)

mlesky_res <- read.csv("outputs/genomics/GPSC31_model_2_adaptation_true_rep_10000_mlesky.csv") %>% 
  dplyr::mutate(new_date = format(lubridate::date_decimal(Year), "%Y-%m-%d")) %>% 
  dplyr::mutate(new_year = lubridate::year(new_date)) #%>% 
  # tidyr::replace_na(,list())

# load joined_clades first (or all_day???)
all_day <- joined_clades %>% 
  dplyr::filter(n703_choosen_GPSC31 == 1) %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise(counts = n()) %>% 
  dplyr::ungroup() #%>% 
  # dplyr::mutate(Earliest.specimen.date = as.Date(Earliest.specimen.date))

mlesky_combined <- merge(mlesky_res, all_day,
                      by.x = c("new_year"),
                      by.y = c("year"))

col_mlesky <- c(
  "Ne" = "black",
  "Data" = "steelblue")

png("pictures/genomics/mleSky.png", width = 30, height = 12, unit = "cm", res = 1200)
ggplot() +
  geom_line(linewidth = 1, data = mlesky_res, aes(x = new_year, y = Ne, colour = "Ne")) +
  # geom_line(linetype = "dashed", data = mlesky_res, aes(x = new_year, y = Lower_range)) +
  # geom_line(linetype = "dashed", data = mlesky_res, aes(x = new_year, y = Upper_range)) +
  geom_bar(stat = "identity", data = all_day, aes(x = year, y = counts, fill = "Data")) +
  geom_vline(data = vaccine_UK, aes(xintercept = year),
             linetype = "dashed") +
  scale_color_manual(values = c(col_mlesky),
                     name = "Line",
                     breaks = c("Ne"),
                     labels = c("Ne")
  ) +
  scale_fill_manual(values = c(col_mlesky),
                     name = "Bar",
                     breaks = c("Data"),
                     labels = c("Data")
  ) +
  # scale_x_continuous(breaks = ~ axisTicks(., log = FALSE)) + # delete weird decimals in Year
  xlim(1950, 2016) +
  geom_label(aes(x = 2006, y = 150, label = "PCV7"),
             fill = "white", color = "black") + # 2006 = PCV7 = "gray80"
  geom_label(aes(x = 2011, y = 150, label = "PCV13"),
             fill = "white", color = "black") + # 2011 = PCV13 = "gray20"
  ggtitle("GPSC31 Population Size in England") +
  xlab("Year") +
  ylab("Effective Population Size (Ne)") +
  theme_bw()
dev.off()
