#--------------------------------------------------------------------------------
# Baleen SIA in fin and blue whale hybrids 
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# 01 - Data preparation
#--------------------------------------------------------------------------------

# Load libraries

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

# Import data

merge <- read_excel("/Users/marcruizisagales/Dropbox/00_papers/2023_hybrids_SIA_baleen/Blue_whale_results/PER_CAROL_BO_0_point_alignments_29_juny_total_-4_amb_sofre_blue_25_april_2024_rec.xlsx")

# 3. Same number of measured points along the baleen (closest to the base for each baleen)
values_by_category <- list(
  Bphy1 = c(-4:25),
  Bphy2 = c(-3:26),
  Bphy3 = c(-3:26),
  Bphy4 = c(-1:28),
  Bphy5 = c(-4:25),
  Hyb1 = c(-4:25),
  Hyb2 = c(-4:25),
  Hyb3 = c(0:29),
  Bmus1 = c(-4:25),
  Bmus2 = c(-4:25)
)

filtered_data <- merge %>%
  filter(
    (Whale == "Bphy1" & Cm %in% values_by_category$Bphy1) |
      (Whale == "Bphy2" & Cm %in% values_by_category$Bphy2) |
      (Whale == "Bphy3" & Cm %in% values_by_category$Bphy3) |
      (Whale == "Bphy4" & Cm %in% values_by_category$Bphy4) |
      (Whale == "Bphy5" & Cm %in% values_by_category$Bphy5) |
      (Whale == "Hyb1" & Cm %in% values_by_category$Hyb1) |
      (Whale == "Hyb2" & Cm %in% values_by_category$Hyb2) |
      (Whale == "Hyb3" & Cm %in% values_by_category$Hyb3) |
      (Whale == "Bmus1" & Cm %in% values_by_category$Bmus1) |
      (Whale == "Bmus2" & Cm %in% values_by_category$Bmus2) 
  )

# 4. Linear interpolation for one cm (out of 300) without analysed isotope ratio

sum(is.na(filtered_data$dN))
sum(is.na(filtered_data$dC))
sum(is.na(filtered_data$dS))

filtered_data$dN <- na.approx(filtered_data$dN)
filtered_data$dC <- na.approx(filtered_data$dC)
filtered_data$dS <- na.approx(filtered_data$dS)

merge_good_ind<-filtered_data
merge_good_ind$Data_capt[241:300] <- as.Date(c(rep("2010/08/23", 30),rep("1987/09/11", 30))) # add year of stranding for the blue whale


