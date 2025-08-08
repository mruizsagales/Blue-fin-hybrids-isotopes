#--------------------------------------------------------------------------------
# Baleen SIA in fin and blue whale hybrids 
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# 02 - Time estimation
#--------------------------------------------------------------------------------
setwd("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/02code")
source("01_Data_preparation.R")

combined_df <- NULL;

for (i in 1:length(unique(merge_good_ind$Whale))) { 
  df_i <- dplyr::filter(merge_good_ind, Whale == unique(merge_good_ind$Whale)[i]) # We select one fin whale individual
  growth_rate <- 16.07  # Baleen plate growth rate (centimeters per year; Ruiz-Sagalés et al., 2024)
  last_sample <- unique(ymd(df_i$Data_capt)) # Date of last keratin sample (catch date)
  df_i <- 
    df_i %>%
    dplyr::mutate(days = Cm * 365.25 / growth_rate) %>% # Number of days between samples 
    dplyr::mutate(days = days - days[1]) %>%
    dplyr::mutate(rev_days = days - days[length(days)]) %>% # Number of days starting at maximum Cm
    dplyr::mutate(sample_date = last_sample + rev_days) %>% # Estimated date
    dplyr::mutate(year = year(sample_date)) %>% # Estimated date year
    dplyr::mutate(year_rev = rev(sample_date)) # Estimated date retrospectively (real date)
  
  df_i$Year_from_sample_date <- format(as.Date(df_i$year_rev), "%Y") # year from the date estimated retrospectively
  df_i$Year_month <- format(as.Date(df_i$year_rev), "%Y/%m") # year and month from the date estimated retrospectively
  df_i$Month_from_sample_date <- format(as.Date(df_i$year_rev), "%m") # month from the date estimated retrospectively
  
  combined_df <- rbind(combined_df, df_i)
}

print(combined_df)
