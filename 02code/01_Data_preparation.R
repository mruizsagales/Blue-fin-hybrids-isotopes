
#--------------------------------------------------------------------------------
# Baleen SIA in fin and blue whale hybrids (Marc Ruiz-Sagalés, 6 de juny del 2025)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# 01 - Data preparation
#--------------------------------------------------------------------------------

# 1. Load libraries

{
  library(car)
  library(corrplot)
  library(devtools)
  library(dplyr)
  library(ecotraj)
  library(egg)
  library(ellipse)
  library(ggcorrplot)
  library(ggeffects)
  library(ggh4x)
  library(ggplot2)
  library(ggrepel)
  library(ggstatsplot)
  library(ggtext)
  library(here)
  library(lubridate)
  library(lme4)
  library(lmerTest)
  library(mgcv)
  library(nicheROVER)
  library(nichetools)
  library(patchwork)
  #library(pairwiseAdonis)
  library(purrr)
  library(RColorBrewer)
  library(readr)
  library(readxl)
  library(rstatix)
  library(stringr)
  library(SuessR)
  library(tidyr)
  library(tidyverse)
  library(TSA)
  library(vegan)
  library(zoo)
  
}

# 2. Import data

merge <- read_excel("/Users/marcruizisagales/Dropbox/00_papers/2023_hybrids_SIA_baleen/Blue_whale_results/PER_CAROL_BO_0_point_alignments_29_juny_total_-4_amb_sofre_blue_25_april_2024.xlsx")

# 3. Same number of measured points along the baleen (closest to the base for each baleen)
values_by_category <- list(
  F13065 = c(-4:25),
  F13066 = c(-3:26),
  F13068 = c(-3:26),
  F13073 = c(-1:28),
  F13076 = c(-4:25),
  F13129 = c(-4:25),
  F18022 = c(-4:25),
  F18098 = c(0:29),
  BMUS_2010 = c(-4:25),
  BMUS_1990 = c(-4:25)
)

filtered_data <- merge %>%
  filter(
    (Whale == "F13065" & Cm %in% values_by_category$F13065) |
      (Whale == "F13066" & Cm %in% values_by_category$F13066) |
      (Whale == "F13068" & Cm %in% values_by_category$F13068) |
      (Whale == "F13073" & Cm %in% values_by_category$F13073) |
      (Whale == "F13076" & Cm %in% values_by_category$F13076) |
      (Whale == "F13129" & Cm %in% values_by_category$F13129) |
      (Whale == "F18022" & Cm %in% values_by_category$F18022) |
      (Whale == "F18098" & Cm %in% values_by_category$F18098) |
      (Whale == "BMUS_2010" & Cm %in% values_by_category$BMUS_2010) |
      (Whale == "BMUS_1990" & Cm %in% values_by_category$BMUS_1990) 
  )

# 4. Linear interpolation for one cm (out of 300) without analysed isotope ratio

sum(is.na(filtered_data$dN))
sum(is.na(filtered_data$dC))
sum(is.na(filtered_data$dS))

filtered_data$dN <- na.approx(filtered_data$dN)
filtered_data$dC <- na.approx(filtered_data$dC)
filtered_data$dS <- na.approx(filtered_data$dS)

merge_good_ind<-filtered_data
merge_good_ind$Data_capt[241:300] <- as.Date(c(rep("2010/08/23", 30),rep("1987/08/23", 30))) # add year of stranding for the blue whale
