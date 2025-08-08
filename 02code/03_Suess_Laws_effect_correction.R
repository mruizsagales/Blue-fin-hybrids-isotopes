#--------------------------------------------------------------------------------
# Baleen SIA in fin and blue whale hybrids 
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# 03 - Suess and Laws effect correction
#--------------------------------------------------------------------------------
setwd("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/02code")
source("02_Time_estimation.R")

combined_df$id <- seq.int(nrow(combined_df)) # add a sequence number to each row
combined_df <- combined_df %>% add_column(region = "Subpolar North Atlantic") # add the Subpolar North Atlantic region to each row
df1 <- combined_df # rename combined_df to df1
subset <- df1[c("id", "dC","Year_from_sample_date","region")] # select the merge3 column names (for RSuess) and name it subset
names(subset) <- c("id", "d13c","year","region") # rename subset columnames
subset$year <- as.numeric(subset$year) # define year as a numeric variable
subset <- as.data.frame(subset) # define subset as a dataframe
df2 <- SuessR(data=subset, correct.to = 2022) # correct the Suess and Laws effect to the year 2022
data_d13cor <- merge(df1,df2,by="id") # merge df1 and df2
df <- data_d13cor
names(df) <- c("id","Cm","d15n","dC","d34s","Whale_1","Whale","Sex","Status","Talla_fetus_(cm)","Sexe_fetus","Length","Talla_(peus)","Data_capt","Lat","Long","Edat","Year","species","days","rev_days","sample_date","year.x","year_rev","Year_from_sample_date","Year_month","Month_from_sample_date","region","year.y","d13c.uncor","Laws.cor","Suess.cor","net.cor","d13c")

#library(openxlsx) # save
#write.xlsx(df, file = "All_merged_time_calibrated_2013_to_2022_Suess_cor.xlsx")


