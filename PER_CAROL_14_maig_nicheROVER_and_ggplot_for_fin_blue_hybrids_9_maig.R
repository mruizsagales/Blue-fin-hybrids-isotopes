# ---- Bring in Packages ---- 
{
  library(dplyr)
  library(ellipse)
  library(ggplot2)
  library(ggtext)
  library(here)
  library(nicheROVER) 
  library(purrr)
  library(patchwork)
  library(readr)
  library(stringr)
  library(tidyr)
}

#--------------------------------------------------------------------------------
#Suess correction of time-calibrated stable isotope data
#--------------------------------------------------------------------------------

# 1. Load libraries
library(tidyverse)
library(lubridate)
library(ggrepel)
library(readxl)
library(ggplot2)
library(egg)
library(zoo)
library(TSA)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(SuessR)

# Trobar script fet per la ECS!

# Blue whale data



#SAME NUMBER OF SAMPLE POINTS
# 2. Import data
#merge <- read_excel("~/Desktop/Doctorat/Analisis_isotops_barbes/Projecte_barbes_clima/Environmental_vars/All_merged_time_calibrated_2013_to_2022_non_Suess_cor.xlsx")
merge <- read_excel("/Users/marcruizisagales/Desktop/Papers/Paper_hybrids/Blue_whale_results/PER_CAROL_BO_0_point_alignments_29_juny_total_-4_amb_sofre_blue_25_april_2024.xlsx")

values_by_category <- list(
  F13065 = c(-5:24),
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



filtered_data$dN <- na.approx(filtered_data$dN)
filtered_data$dC <- na.approx(filtered_data$dC)
filtered_data$dS <- na.approx(filtered_data$dS)

#seleccionar individus i fer algo d'interpolacio amb el na.approx?
#merge_good_ind <- merge %>% filter(Whale %in% c("F13065", "F13066", "F13068","F13073", "F13076", "F13129", "F18022", "F18098", "BMUS_2010","BMUS_1990"))

merge_good_ind<-filtered_data

merge_good_ind <-merge_good_ind[complete.cases(merge_good_ind$dC), ] # remove the rows with no d13C values in order to be able to correct for the Suess effect

merge_good_ind$Whale
merge_good_ind$Data_capt
merge_good_ind$Data_capt[240:299] <- as.Date(c(rep("2010/08/23", 30),rep("1987/08/23", 30))) #segona data inventada


#date
point_size <- 2
combined_df_blue <- NULL;

for (i in 1:length(unique(merge_good_ind$Whale))) { 
  df_i <- dplyr::filter(merge_good_ind, Whale == unique(merge_good_ind$Whale)[i]) # We select one fin whale individual from 2013
  growth_rate <- 15.5  # Baleen plate growth rate (centimeters per year)
  last_sample <- unique(ymd(df_i$Data_capt)) # Date of last keratin sample (catch date)
  df_i <- 
    df_i %>%
    dplyr::mutate(days = Cm * 365.25 / growth_rate) %>% # Number of days between samples starting at Cm 0
    dplyr::mutate(days = days - days[1]) %>% # Number of days between samples starting at Cm -4
    dplyr::mutate(rev_days = days - days[length(days)]) %>% # Number of days starting at maximum Cm
    dplyr::mutate(sample_date = last_sample + rev_days) %>% # Estimated date
    dplyr::mutate(year = year(sample_date)) %>% # Estimated date year
    dplyr::mutate(year_rev = rev(sample_date)) # Estimated date retrospectively (real date)
  
  df_i$Year_from_sample_date <- format(as.Date(df_i$year_rev), "%Y") # year from the date estimated retrospectively
  df_i$Year_month <- format(as.Date(df_i$year_rev), "%Y/%m") # year and month from the date estimated retrospectively
  df_i$Month_from_sample_date <- format(as.Date(df_i$year_rev), "%m") # month from the date estimated retrospectively
  
  # plot top and bottom axis
  isop <- ggplot(df_i, aes(x = rev(Cm))) + 
    theme_classic(base_size = 16) + 
    xlab("Baleen sample (cm from youngest sample)") + 
    scale_x_continuous(breaks = df_i$Cm, 
                       labels = rev(df_i$Cm),
                       sec.axis = sec_axis(~., 
                                           breaks = rev(df_i$Cm),
                                           labels = rev(df_i$sample_date),
                                           name = "Date"))
  
  #isop <- ggplot(df_i, aes(x = rev(sample_date))) + xlab("Date") + scale_x_continuous(breaks = seq(min(df_i$sample_date), max(df_i$sample_date), by = "1 month"),labels = format(seq(min(rev(df_i$sample_date)), max(rev(df_i$sample_date)), by = "1 month"), "%b-%Y"))
  
  # plot carbon stable isotope values
  isop <- isop + geom_line(aes(y = dC), col = "black") + 
    geom_point(aes(y = dC), size = point_size, shape = 21, fill = "black", col = "black") +
    ylab(expression(paste(delta^{13}, "C (\u2030)")))+ theme_article()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  a <- mean(na.omit(df_i$dN))# add the d15N data but need to scale it to d13C
  a
  b <- mean(na.omit(df_i$dC))# add the d15N data but need to scale it to d13C
  b
  v_shift <- a - b
  
  # plot nitrogen stable isotope values
  isop <- isop + geom_line(aes(y = dN - v_shift), col = "darkgrey") + 
    geom_point(aes(y = (dN - v_shift)), size = point_size, 
               shape = 21, fill = "grey", col = "darkgrey")
  
  # plot x and y axis
  isop <- isop + 
    scale_y_continuous(sec.axis = sec_axis(~.+v_shift, 
                                           name = expression(paste(delta^{15}, "N (\u2030)"))),limits = c(-30, -10)) + labs(title= unique(merge_good_ind$Whale)[i])
  new_df <- df_i
  combined_df_blue <- rbind(combined_df_blue, new_df)
  print(isop)
}


#merge1 <- merge %>% 
  dplyr::mutate(Species = ifelse(Whale %in% c("F13129", "F18022", "F18098"), "hybrid", "fin")) 

#merge2 <- dplyr::filter(merge1, Species == "fin") # remove blue-fin hybrids

length(unique(merge1$Whale)) # we have 30 fin whale individuals (being non-hybrids)

#merge3 <- dplyr::filter(merge2, Whale != "F18025") # select all the whales that are not F18025 (non well-sampled)



# 3. Suess and Laws effect correction

merge3 <-combined_df_blue

merge3$id <-seq.int(nrow(merge3)) # add a sequence number to each row

merge3 <- merge3 %>% add_column(region = "Subpolar North Atlantic") # add the Subpolar North Atlantic region to each row

df1 <- merge3 # rename merge3 to df1

subset <- df1[c("id", "dC","Year_from_sample_date","region")] # select the merge3 column names (for RSuess) and name it subset

names(subset) <- c("id", "d13c","year","region") # rename subset columnames

subset$year <- as.numeric(subset$year) # define year as a numeric variable

subset <- as.data.frame(subset) # define subset as a dataframe

df2 <- SuessR(data=subset, correct.to = 2022) # correct the Suess and Laws effect to the year 2022

fw.data_d13cor <- merge(df1,df2,by="id") # merge df1 and df2

library(openxlsx) # Save
#write.xlsx(fw.data_d13cor, file = "All_merged_time_calibrated_2013_to_2022_Suess_cor.xlsx")

# ---- Bring in example data; will need to replace with your own ----

df <- fw.data_d13cor
names(df) <- c("id","Cm","d15n","dC","d34s","Whale_1","Whale","Sex","Status","Talla_fetus_(cm)","Sexe_fetus","Length","Talla_(peus)",          
"Data_capt","Lat","Long","Edat","Year","species","days","rev_days","sample_date","year.x","year_rev","Year_from_sample_date","Year_month","Month_from_sample_date","region",                
"year.y","d13c.uncor","Laws.cor","Suess.cor","net.cor","d13c")

unique(df$Whale)


ggplot(data=df, aes(df$sample_date, df$d15n, color=Whale)) + geom_point() + geom_line() + ylab(expression(paste(delta^{15}, "N (\u2030)"))) + theme(aspect.ratio = 3/4) + xlab("Year") + theme_bw(base_size = 15) + facet_wrap(vars(spec))
ggplot(data=df, aes(df$Cm, df$d15n, color=Whale)) + geom_point() + geom_line() + ylab(expression(paste(delta^{15}, "N (\u2030)"))) + theme(aspect.ratio = 3/4) + xlab("Cm") + theme_bw(base_size = 15) + facet_wrap(vars(species))

ggplot(data=df, aes(df$sample_date, df$d13c, color=Whale)) + geom_point() + geom_line() + ylab(expression(paste(delta^{13}, "C (\u2030)"))) + theme(aspect.ratio = 3/4) + xlab("Year") + theme_bw(base_size = 15)
ggplot(data=df, aes(df$Cm, df$d13c, color=Whale)) + geom_point() + geom_line() + ylab(expression(paste(delta^{13}, "C (\u2030)"))) + theme(aspect.ratio = 3/4) + xlab("Cm") + theme_bw(base_size = 15) + facet_wrap(vars(species))

ggplot(data=df, aes(df$sample_date, df$d34s, color=Whale)) + geom_point() + geom_line() + ylab(expression(paste(delta^{34}, "S (\u2030)"))) + theme(aspect.ratio = 3/4) + xlab("Year") + theme_bw(base_size = 15)
ggplot(data=df, aes(df$Cm, df$d34s, color=Whale)) + geom_point() + geom_line() + ylab(expression(paste(delta^{34}, "S (\u2030)"))) + theme(aspect.ratio = 3/4) + xlab("Cm") + theme_bw(base_size = 15) + facet_wrap(vars(species))

ggplot(data=df, aes(df$sample_date, df$d13c, color=Whale)) + geom_line()

library(ggstatsplot)
set.seed(123)

ggbetweenstats(
  data  = df,
  x     = Whale,
  y     = d15n
)

ggbetweenstats(
  data  = df,
  x     = Whale,
  y     = d13c,
  pairwise.display = "s"
)

ggbetweenstats(
  data  = df,
  x     = Whale,
  y     = d34s
)


###################################################################################
#fer-ho amb nicheROVER
###################################################################################

library(nicheROVER)

fish <- df[,c(3,34,5,7,19)]
fish1 <-fish[complete.cases(fish), ]

#test for normality
library(rstatix)
fish1 %>% dplyr::group_by(species) %>% shapiro_test(d15n)
fish1 %>% dplyr::group_by(species) %>% shapiro_test(d13c)
fish1 %>% dplyr::group_by(species) %>% shapiro_test(d34s)

#test for homoscedasticity
bartlett.test(d15n ~ species, data = fish1)
bartlett.test(d13c ~ species, data = fish1)
bartlett.test(d34s ~ species, data = fish1)

data(fish) # 4 fish, 3 isotopes
aggregate(fish1[1:3], fish1[4], mean) # isotope means calculated for each speci
aggregate(fish1[1:3], fish1[4], sd) 
# fish data
data(fish)

# generate parameter draws from the "default" posteriors of each fish
nsamples <- 1e3
system.time({
  fish.par <- tapply(1:nrow(fish1), fish1$Whale,
                     function(ii) niw.post(nsamples = nsamples, X = fish1[ii,1:3]))
})

# various parameter plots
#clrs <- c("red", "blue", "orange") # colors for each species
#clrs <- brewer.pal(n = 10, name = "RdYlBu")

clrs <- c( "#313695","#313695", "#A50026" , "#A50026" , "#A50026", "#A50026" , "#A50026","#F46D43" , "#F46D43" , "#F46D43" )

# mu1 (del15N), mu2 (del13C), and Sigma12
par(mar = c(4, 4, .5, .1)+.1, mfrow = c(1,3))
niche.par.plot(fish.par, col = clrs, plot.index = 1)
niche.par.plot(fish.par, col = clrs, plot.index = 2)
niche.par.plot(fish.par, col = clrs, plot.index = 1:2)
legend("topleft", legend = names(fish.par), fill = clrs)

# all mu (del15N, del13C, del34S)
niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = FALSE)
legend("topleft", legend = names(fish.par), fill = clrs)

# all mu and Sigma
par(mar = c(4.2, 4.2, 2, 1)+.1)
niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = TRUE)
legend("topright", legend = names(fish.par), fill = clrs)

# 2-d projections of 10 niche regions
#clrs <- c("blue","red", "orange") # colors for each species
nsamples <- 1
fish.par <- tapply(1:nrow(fish1), fish1$Whale,
                   function(ii) niw.post(nsamples = nsamples, X = fish1[ii,1:3]))

# format data for plotting function
fish.data <- tapply(1:nrow(fish1), fish1$Whale, function(ii) X = fish1[ii,1:3])

niche.plot(niche.par = fish.par, niche.data = fish.data, pfrac = .05, alpha = 0.4,
           iso.names = expression(delta^{15}*N, delta^{13}*C, delta^{34}*S),
           col = clrs, xlab = expression("Isotope Ratio (per mil)"))

# niche overlap plots for 95% niche region sizes
nsamples <- 1000
fish.par <- tapply(1:nrow(fish1), fish1$Whale,
                   function(ii) niw.post(nsamples = nsamples, X = fish1[ii,1:3]))

# Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher accuracy.
# the variable over.stat can be supplied directly to the overlap.plot function

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1e3, alpha = c(.95, 0.99))

#The mean overlap metrics calculated across iteratations for both niche 
#region sizes (alpha = .95 and alpha = .99) can be calculated and displayed in an array.
over.mean <- apply(over.stat, c(1:2,4), mean)*100
round(over.mean, 2)

$ <- apply(over.stat*100, c(1:2, 4), quantile, prob = c(.025, .975), na.rm = TRUE)
round(over.cred[,,,1]) # display alpha = .95 niche region

# Overlap plot.Before you run this, make sure that you have chosen your 
#alpha level.
#clrs <- c("blue","red", "orange") # colors for each species
over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1e3, alpha = .95)
overlap.plot(over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")

# posterior distribution of (mu, Sigma) for each species
nsamples <- 1000
fish.par <- tapply(1:nrow(fish1), fish1$Whale,
                   function(ii) niw.post(nsamples = nsamples, X = fish1[ii,1:3]))

# posterior distribution of niche size by species
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})

# point estimate and standard error
rbind(est = colMeans(fish.size),
      se = apply(fish.size, 2, sd))

# boxplots
#clrs <- c("blue","red", "orange") # colors for each species
dev.new()
boxplot(fish.size, col = clrs, pch = 16, cex = .5,
        ylab = "Niche Size", xlab = "Species")




#----------------------------------------------------------------------------------------
### per sp

library(nicheROVER)

fish <- df[,c(3,34,5,19)]
fish1 <-fish[complete.cases(fish), ]

#diferences of SI ratios between sp
fish2 <- df[,c(3,34,5,19, 7)]
fish3 <-fish2[complete.cases(fish2), ]
fish3$species<- factor(fish3$species, 
                             levels = c('fin','hybrid','blue'))



#Nitrogen
library(lme4)
lmer_dN <- lmer(d15n ~ 1 + species + (1|Whale), data=fish3)
summary(lmer_dN)
plot(allEffects(lmer_dN))

# Model validation
plot_model(lmer_dN) # diagnostics
cv_m1 <- loo_cv(lmer_dN, df, "Whale", keep = "used") # cross-validation
accuracy(lmer_dN) # accuracy
compare_accuracy(lmer_dN, cv_m1) # compare accuracy
confint(lmer_dN) # confint
summary(lmer_dN)

# Model validation with boostraping
qqPlot(resid(lmer_dN)) # normality of the residuals
qqPlot(c(ranef(lmer_dN)$Whale$`(Intercept)`)) # normality of the random effects
plot(lmer_dN) # variance homogeneity

f = function(m) {res= fixef(m)[2]; names(res) = names(fixef(m)[2]); res}
boost <- lme4::bootMer(lmer_dN, f, nsim=200, type="semiparametric", use.u=TRUE)
plot(boost, index= 1)

#Carbon
fish3$species<- factor(fish3$species, 
                       levels = c('fin','hybrid','blue'))

lmer_dC <- lmer(d13c ~ 1 + species + (1|Whale), data=fish3)
summary(lmer_dC)
plot(allEffects(lmer_dC))

# Model validation
plot_model(lmer_dC) # diagnostics
cv_m1 <- loo_cv(lmer_dC, df, "Whale", keep = "used") # cross-validation
accuracy(lmer_dC) # accuracy
compare_accuracy(lmer_dC, cv_m1) # compare accuracy
confint(lmer_dC) # confint
summary(lmer_dC)

# Model validation with boostraping
qqPlot(resid(lmer_dC)) # normality of the residuals
qqPlot(c(ranef(lmer_dC)$Whale$`(Intercept)`)) # normality of the random effects
plot(lmer_dC) # variance homogeneity

f = function(m) {res= fixef(m)[2]; names(res) = names(fixef(m)[2]); res}
boost <- lme4::bootMer(lmer_dC, f, nsim=200, type="semiparametric", use.u=TRUE)
plot(boost, index= 1)

#Sulfur
fish3$species<- factor(fish3$species, 
                       levels = c('fin','hybrid','blue'))
lmer_dS <- lmer(d34s ~ 1 + species + (1|Whale), data=fish3)
summary(lmer_dS)
plot(allEffects(lmer_dS))

# Model validation
plot_model(lmer_dS) # diagnostics
cv_m1 <- loo_cv(lmer_dS, df, "Whale", keep = "used") # cross-validation
accuracy(lmer_dS) # accuracy
compare_accuracy(lmer_dS, cv_m1) # compare accuracy
confint(lmer_dS) # confint
summary(lmer_dS)

# Model validation with boostraping
qqPlot(resid(lmer_dS)) # normality of the residuals
qqPlot(c(ranef(lmer_dS)$Whale$`(Intercept)`)) # normality of the random effects
plot(lmer_dS) # variance homogeneity

f = function(m) {res= fixef(m)[2]; names(res) = names(fixef(m)[2]); res}
boost <- lme4::bootMer(lmer_dS, f, nsim=200, type="semiparametric", use.u=TRUE)
plot(boost, index= 1)

#

data(fish) # 4 fish, 3 isotopes
aggregate(fish1[1:3], fish1[4], mean) # isotope means calculated for each speci

# fish data
data(fish)

# generate parameter draws from the "default" posteriors of each fish
nsamples <- 1e3
system.time({
  fish.par <- tapply(1:nrow(fish1), fish1$species,
                     function(ii) niw.post(nsamples = nsamples, X = fish1[ii,1:3]))
})

# various parameter plots
clrs <- c( "#313695",  "#A50026", "#F46D43" )
#clrs <- c( "blue", "red", "orange") # colors for each species
#clrs <- brewer.pal(n = 10, name = "RdYlBu")

# mu1 (del15N), mu2 (del13C), and Sigma12
par(mar = c(4, 4, .5, .1)+.1, mfrow = c(1,3))
niche.par.plot(fish.par, col = clrs, plot.index = 1)
niche.par.plot(fish.par, col = clrs, plot.index = 2)
niche.par.plot(fish.par, col = clrs, plot.index = 1:2)
legend("topleft", legend = names(fish.par), fill = clrs)

# all mu (del15N, del13C, del34S)
niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = FALSE)
legend("topleft", legend = names(fish.par), fill = clrs)

# all mu and Sigma
par(mar = c(4.2, 4.2, 2, 1)+.1)
niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = TRUE)
legend("topright", legend = names(fish.par), fill = clrs)

# 2-d projections of 10 niche regions
#clrs <- c("blue","red", "orange") # colors for each species
nsamples <- 10
fish.par <- tapply(1:nrow(fish1), fish1$species,
                   function(ii) niw.post(nsamples = nsamples, X = fish1[ii,1:3]))

# format data for plotting function
fish.data <- tapply(1:nrow(fish1), fish1$species, function(ii) X = fish1[ii,1:3])

niche.plot(niche.par = fish.par, niche.data = fish.data, pfrac = .05, alpha= 0.4,
           iso.names = expression(delta^{15}*N, delta^{13}*C, delta^{34}*S),
           col = clrs, xlab = expression("Isotope Ratio (per mil)"))

# niche overlap plots for 95% niche region sizes
nsamples <- 1000
fish.par <- tapply(1:nrow(fish1), fish1$Whale,
                   function(ii) niw.post(nsamples = nsamples, X = fish1[ii,1:3]))

# Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher accuracy.
# the variable over.stat can be supplied directly to the overlap.plot function

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1e3, alpha = c(.95, 0.99))

#The mean overlap metrics calculated across iteratations for both niche 
#region sizes (alpha = .95 and alpha = .99) can be calculated and displayed in an array.
over.mean <- apply(over.stat, c(1:2,4), mean)*100
round(over.mean, 2)

over.cred <- apply(over.stat*100, c(1:2, 4), quantile, prob = c(.025, .975), na.rm = TRUE)
round(over.cred[,,,1]) # display alpha = .95 niche region

# Overlap plot.Before you run this, make sure that you have chosen your 
#alpha level.
#clrs <- c("blue","red", "orange") # colors for each species
over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1e3, alpha = .95)
overlap.plot(over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")

# posterior distribution of (mu, Sigma) for each species
nsamples <- 1000
fish.par <- tapply(1:nrow(fish1), fish1$species,
                   function(ii) niw.post(nsamples = nsamples, X = fish1[ii,1:3]))

# posterior distribution of niche size by species
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})

# point estimate and standard error
rbind(est = colMeans(fish.size),
      se = apply(fish.size, 2, sd))

# boxplots
#clrs <- c("blue","red", "orange") # colors for each species

boxplot(fish.size, col = clrs, pch = 16, cex = .5,
        ylab = "Niche Size", xlab = "Species")


## amb ggplot

# set the number of posterior samples that will be taken 
nsample <- 10000


# ---- Use map from purrr and niw.post from nicheROVER to estimate niches ----- 
fish_par <- fish1 %>% 
  split(.$species) %>% 
  map(~ dplyr::select(., d15n, d13c, d34s)) %>% 
  map(~niw.post(nsample = nsample, X = .))


# ---- Extract mu from nicheROVER object ----- 
df_mu <- map(fish_par, pluck, 1) %>% 
  imap(~ as_tibble(.x) %>% 
         dplyr::mutate( 
           metric = "mu", 
           species = .y
         )
  ) %>%
  dplyr::bind_rows() %>% 
  dplyr::mutate(
    species = factor(species, 
                     levels = c("fin", "hybrid", "blue"))
  ) %>% 
  dplyr::group_by(species) %>% 
  dplyr::mutate(
    sample_number = 1:10000
  ) %>% 
  dplyr::ungroup()




df_mu_long <- df_mu %>% 
  pivot_longer(cols = -c(metric, species, sample_number), 
               names_to = "isotope", 
               values_to = "mu_est") %>% 
  dplyr::mutate(
    element =  dplyr::case_when(
      isotope == "d15n" ~ "N",
      isotope == "d13c" ~ "C",
      isotope == "d34s" ~ "S",
    ), 
    neutron =  dplyr::case_when(
      isotope == "d15n" ~ 15,
      isotope == "d13c" ~ 13,
      isotope == "d34s" ~ 34,
    ) 
  )

# ---- Extract sigma from nicheROVER object ----- 

df_sigma <- map(fish_par, pluck, 2) %>%
  imap(~ as_tibble(.x) %>%
         dplyr::mutate(
           metric = "sigma",
           id = c("d15n", "d13c","d34s"),
           species = .y
         )
  ) %>%
  dplyr::bind_rows() %>%
  pivot_longer(cols = -c("id", "species", "metric"),
               names_to = "isotope",
               values_to = "post_sample"
  )  %>%
  separate(isotope, into = c("isotopes", "sample_number"), sep = "\\.")



df_sigma_cn <- df_sigma %>%
  dplyr::filter(id != isotopes)

# ---- Plot posterior samples for mu estimates from nicheROver ----- 

#change this color for a palettes 
posterior_plots <- df_mu_long %>%
  split(.$isotope) %>%
  imap(
    ~ ggplot(data = ., aes(x = mu_est)) +
      geom_density(aes(fill = species), alpha = 0.5) +
      scale_fill_brewer(palette = "Set1", name = "Species") +
      #scale_fill_viridis_d(begin = 0.25, end = 0.75, option = "D", name = "Species") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.title.x =  element_markdown(),
            axis.title.y =  element_markdown(),
            legend.position = "none"
      ) +
      labs(
        x = paste("\u00b5<sub>\U03B4</sub>", "<sub><sup>",
                  unique(.$neutron), "</sup></sub>",
                  "<sub>",unique(.$element), "</sub>", sep = ""),
        y = paste0("p(\u00b5 <sub>\U03B4</sub>","<sub><sup>",
                   unique(.$neutron), "</sub></sup>",
                   "<sub>",unique(.$element),"</sub>",
                   " | X)"), sep = "")
  )

posterior_plots_pl <- posterior_plots$d15n  +
  theme(legend.position = c(0.18, 0.84)) + 
  posterior_plots$d13c + posterior_plots$d34s 

ggsave("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/png/Posterior_plots.png", posterior_plots_pl, 
       device = png(width = 800, height = 300))

# ---- Prepare sigma dataframes for plotting ----



df_sigma_cn <- df_sigma_cn %>%
  mutate(
    element_id = case_when(
      id == "d15n" ~ "N",
      id == "d13c" ~ "C",
      id == "d34s" ~ "S",
    ),
    neutron_id = case_when(
      id == "d15n" ~ 15,
      id == "d13c" ~ 13,
      id == "d34s" ~ 34,
    ),
    element_iso = case_when(
      isotopes == "d15n" ~ "N",
      isotopes == "d13c" ~ "C",
      isotopes == "d34s" ~ "S",
    ),
    neutron_iso = case_when(
      isotopes == "d15n" ~ 15,
      isotopes == "d13c" ~ 13,
      isotopes == "d34s" ~ 34,
    )
  )

# ---- Plot posterior samples for mu estimates from nicheROver ----- 
sigma_plots <- df_sigma_cn %>%
  group_split(id, isotopes) %>%
  imap(
    ~ ggplot(data = ., aes(x = post_sample)) +
      geom_density(aes(fill = species), alpha = 0.5) +
      scale_fill_brewer(palette = "Set1", name = "Species") +
      #scale_fill_viridis_d(begin = 0.25, end = 0.75, option = "D", name = "Species") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.title.x =  element_markdown(),
            axis.title.y =  element_markdown(),
            legend.position = "none"
      ) +
      labs(
        x = paste("\U03A3","<sub>\U03B4</sub>",
                  "<sub><sup>", unique(.$neutron_id), "</sub></sup>",
                  "<sub>",unique(.$element_id),"</sub>"," ",
                  "<sub>\U03B4</sub>",
                  "<sub><sup>", unique(.$neutron_iso), "</sub></sup>",
                  "<sub>",unique(.$element_iso),"</sub>", sep = ""),
        y = paste("p(", "\U03A3","<sub>\U03B4</sub>",
                  "<sub><sup>", unique(.$neutron_id), "</sub></sup>",
                  "<sub>",unique(.$element_id),"</sub>"," ",
                  "<sub>\U03B4</sub>",
                  "<sub><sup>", unique(.$neutron_iso), "</sub></sup>",
                  "<sub>",unique(.$element_iso),"</sub>", " | X)", sep = ""),
      )
  )

sigma_plots_pl <- sigma_plots[[1]] +sigma_plots[[2]] +sigma_plots[[3]] +sigma_plots[[4]] +sigma_plots[[5]] +sigma_plots[[6]] +
  theme(legend.position = c(0.1, 0.82))

ggsave("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/png/Sigma_plots.png", sigma_plots_pl, 
       device = png(width = 800, height = 600))

# ---- Manipulate sigma dataframes for ellipse loops ----- 
df_sigma_wide <- df_sigma %>%
  pivot_wider(names_from = id,
              values_from = post_sample)


# ---- Set the confidence interval for ellipses ----
p.ell <- 0.4

# ---- create blank list to dump results and category object to loop for 
species_name <- unique(df_sigma_wide$species)



all_ellipses_NC <- list()


# ---- Loop to create ellipse ------ 
for (i in 1:length(species_name)) {
  
  sigma_species <- df_sigma_wide %>% 
    filter(species %in% species_name[i])
  
  sigma_species <- filter(sigma_species, isotopes != "d34s")
  
  mu_species <- df_mu %>% 
    filter(species %in% species_name[i])
  
  ell <- NULL
  post.id <- NULL
  
  for(j in 1:length(unique(sigma_species$sample_number))) {
    
    sigma_ind <- sigma_species %>%
      filter(sample_number %in% sample_number[j]) %>% 
      dplyr::select(d15n, d13c) 
    
    Sigma <- as.matrix(sigma_ind, 2, 2)
    row.names(Sigma) <- c("d15n", "d13c")
    
    mu <- mu_species %>%
      filter(sample_number %in% sample_number[j]) %>% 
      dplyr::select(sample_number, d15n, d13c) %>% 
      pivot_longer(cols = -sample_number, 
                   names_to = "isotope", 
                   values_to = "mu") %>% 
      .$mu
    
    out <- ellipse::ellipse(Sigma, centre = mu, which = c(1,2), level = p.ell)
    
    ell <- rbind(ell, out)
    post.id <- c(post.id, rep(j, nrow(out)))
  }
  ell <- as.data.frame(ell)
  ell$rep <- post.id
  all_ellipses_NC[[i]] <- ell
}

# ---- Combine ellipse list into dataframe and add species names back in -----
ellipse_df_NC <-  dplyr::bind_rows(all_ellipses_NC, .id = "id") %>% 
  dplyr::mutate(
    species = factor(
      dplyr::case_when(
        id == "1" ~ "fin",
        id == "2" ~ "hybrid",
        id == "3" ~ "blue"
      ), level = c("fin", "hybrid","blue")
    )
  ) %>% 
  as_tibble()

# ---- Randomly sample 10 ellipse for plotting ----- 

ellipse_df_NC %>% 
  dplyr::group_by(species, rep) %>% 
  nest() %>%
  dplyr::group_by(species) %>% 
  dplyr::slice_sample(n = 10, replace = TRUE) %>% 
  dplyr::ungroup() %>% 
  unnest(cols = c(data)) -> random_ellipse_NC 


# ---- Plot ellipse, biplot and each isotopes density ----- 

ellipse_plots_1 <- ggplot() + 
  geom_polygon(data = random_ellipse_NC,
               mapping = aes(x = d13c, y = d15n,
                             group = interaction(rep, species),
                             color = species),
               fill = NA,
               linewidth = 0.5) + 
  
  scale_color_brewer(palette = "Set1", name = "Species") +
  #scale_color_viridis_d(begin = 0.25, end = 0.75, option = "D", name = "Species") +
  #scale_x_continuous(breaks = rev(seq(-20, -40, -2))) +
  #scale_y_continuous(breaks = seq(6, 16, 2)) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta^{13}, "C (\u2030)")), 
       y = expression(paste(delta^{15}, "N (\u2030)")))


all_ellipses_NS <- list()


# ---- Loop to create ellipse ------ 
for (i in 1:length(species_name)) {
  
  sigma_species <- df_sigma_wide %>% 
    filter(species %in% species_name[i])
  
  sigma_species <- filter(sigma_species, isotopes != "d13c")
  
  mu_species <- df_mu %>% 
    filter(species %in% species_name[i])
  
  ell <- NULL
  post.id <- NULL
  
  for(j in 1:length(unique(sigma_species$sample_number))) {
    
    sigma_ind <- sigma_species %>%
      filter(sample_number %in% sample_number[j]) %>% 
      dplyr::select(d15n, d34s) 
    
    Sigma <- as.matrix(sigma_ind, 2, 2)
    row.names(Sigma) <- c("d15n", "d34s")
    
    mu <- mu_species %>%
      filter(sample_number %in% sample_number[j]) %>% 
      dplyr::select(sample_number, d15n, d34s) %>% 
      pivot_longer(cols = -sample_number, 
                   names_to = "isotope", 
                   values_to = "mu") %>% 
      .$mu
    
    out <- ellipse::ellipse(Sigma, centre = mu, which = c(1,2), level = p.ell)
    
    ell <- rbind(ell, out)
    post.id <- c(post.id, rep(j, nrow(out)))
  }
  ell <- as.data.frame(ell)
  ell$rep <- post.id
  all_ellipses_NS[[i]] <- ell
}

# ---- Combine ellipse list into dataframe and add species names back in -----
ellipse_df_NS <-  dplyr::bind_rows(all_ellipses_NS, .id = "id") %>% 
  dplyr::mutate(
    species = factor(
      dplyr::case_when(
        id == "1" ~ "fin",
        id == "2" ~ "hybrid",
        id == "3" ~ "blue"
      ), level = c("fin", "hybrid","blue")
    )
  ) %>% 
  as_tibble()

# ---- Randomly sample 10 ellipse for plotting ----- 

ellipse_df_NS %>% 
  group_by(species, rep) %>% 
  nest() %>%
  group_by(species) %>% 
  slice_sample(n = 10, replace = TRUE) %>% 
  ungroup() %>% 
  unnest(cols = c(data)) -> random_ellipse_NS 


# ---- Plot ellipse, biplot and each isotopes density ----- 

ellipse_plots_2 <- ggplot() + 
  geom_polygon(data = random_ellipse_NS,
               mapping = aes(x = d34s, y = d15n,
                             group = interaction(rep, species),
                             color = species),
               fill = NA,
               linewidth = 0.5) + 
  
  scale_color_brewer(palette = "Set1", name = "Species") +
  #scale_color_viridis_d(begin = 0.25, end = 0.75, option = "D", name = "Species") +
  #scale_x_continuous(breaks = rev(seq(-20, -40, -2))) +
  #scale_y_continuous(breaks = seq(6, 16, 2)) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta^{34}, "S (\u2030)")), 
       y = expression(paste(delta^{15}, "N (\u2030)")))

all_ellipses_CS <- list()


# ---- Loop to create ellipse ------ 
for (i in 1:length(species_name)) {
  
  sigma_species <- df_sigma_wide %>% 
    filter(species %in% species_name[i])
  
  sigma_species <- filter(sigma_species, isotopes != "d15n")
  
  mu_species <- df_mu %>% 
    filter(species %in% species_name[i])
  
  ell <- NULL
  post.id <- NULL
  
  for(j in 1:length(unique(sigma_species$sample_number))) {
    
    sigma_ind <- sigma_species %>%
      filter(sample_number %in% sample_number[j]) %>% 
      dplyr::select(d13c, d34s) 
    
    Sigma <- as.matrix(sigma_ind, 2, 2)
    row.names(Sigma) <- c("d13c", "d34s")
    
    mu <- mu_species %>%
      filter(sample_number %in% sample_number[j]) %>% 
      dplyr::select(sample_number, d13c, d34s) %>% 
      pivot_longer(cols = -sample_number, 
                   names_to = "isotope", 
                   values_to = "mu") %>% 
      .$mu
    
    out <- ellipse::ellipse(Sigma, centre = mu, which = c(1,2), level = p.ell)
    
    ell <- rbind(ell, out)
    post.id <- c(post.id, rep(j, nrow(out)))
  }
  ell <- as.data.frame(ell)
  ell$rep <- post.id
  all_ellipses_CS[[i]] <- ell
}

# ---- Combine ellipse list into dataframe and add species names back in -----
ellipse_df_CS <- bind_rows(all_ellipses_CS, .id = "id") %>% 
  mutate(
    species = factor(
      case_when(
        id == "1" ~ "fin",
        id == "2" ~ "hybrid",
        id == "3" ~ "blue"
      ), level = c("fin", "hybrid", "blue")
    )
  ) %>% 
  as_tibble()

# ---- Randomly sample 10 ellipse for plotting ----- 

ellipse_df_CS %>% 
  group_by(species, rep) %>% 
  nest() %>%
  group_by(species) %>% 
  slice_sample(n = 10, replace = TRUE) %>% 
  ungroup() %>% 
  unnest(cols = c(data)) -> random_ellipse_CS 


# ---- Plot ellipse, biplot and each isotopes density ----- 

ellipse_plots_3 <- ggplot() + 
  geom_polygon(data = random_ellipse_CS,
               mapping = aes(x = d34s, y = d13c,
                             group = interaction(rep, species),
                             color = species),
               fill = NA,
               linewidth = 0.5) + 
  
  scale_color_brewer(palette = "Set1", name = "Species") +
  #scale_color_viridis_d(begin = 0.25, end = 0.75, option = "D", name = "Species") +
  #scale_x_continuous(breaks = rev(seq(-20, -40, -2))) +
  #scale_y_continuous(breaks = seq(6, 16, 2)) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta^{34}, "S (\u2030)")), 
       y = expression(paste(delta^{13}, "C (\u2030)")))


iso_long <- fish1 %>% pivot_longer(cols = -species,
                                                         names_to = "isotope", 
                                                         values_to = "value") %>% 
  mutate(
    element = case_when(
      isotope == "d15n" ~ "N",
      isotope == "d13c" ~ "C",
      isotope == "d34s" ~ "S",
    ), 
    neutron = case_when(
      isotope == "d15n" ~ 15,
      isotope == "d13c" ~ 13,
      isotope == "d34s" ~ 34,
    )
  )

iso_density1 <- iso_long %>% 
  group_split(isotope) 
iso_density1 <- ggplot(data = iso_density1[[1]]) + 
  geom_density(aes(x = value, 
                   color = species), 
               alpha = 0.35, 
               linewidth = 0.8) +
  scale_color_brewer(palette = "Set1", name = "Species") +
  #scale_fill_viridis_d(begin = 0.25, end = 0.75, option = "D", name = "Species") +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = 'none', 
        legend.title.align = 0.5,
        legend.background = element_blank(), 
        axis.title.x = element_markdown(family = "sans")) + 
  labs(x =  paste("\U03B4",
                  "<sup>", unique(iso_density1[[1]]$neutron), "</sup>",unique(iso_density1[[1]]$element), " (\u2030)",
                  sep = ""), 
       y = "Density")

iso_density2 <- iso_long %>% 
  group_split(isotope) 
iso_density2 <- ggplot(data = iso_density2[[2]]) + 
  geom_density(aes(x = value, 
                   color = species), 
               alpha = 0.35, 
               linewidth = 0.8) +
  scale_color_brewer(palette = "Set1", name = "Species") +
  #scale_fill_viridis_d(begin = 0.25, end = 0.75, option = "D", name = "Species") +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.15, 0.75), 
        legend.title = element_blank(),
        #legend.title.align = 0.5,
        legend.background = element_blank(),
        legend.key.size = unit(0.5, "lines"),
        axis.title.x = element_markdown(family = "sans")) + 
  labs(x =  paste("\U03B4",
                  "<sup>", unique(iso_density2[[2]]$neutron), "</sup>",unique(iso_density2[[2]]$element), " (\u2030)",
                  sep = ""), 
       y = "Density")


iso_density3 <- iso_long %>% 
  group_split(isotope) 
iso_density3 <- ggplot(data = iso_density3[[3]]) + 
  geom_density(aes(x = value, 
                   color = species), 
               alpha = 0.35, 
               linewidth = 0.8) +
  scale_color_brewer(palette = "Set1", name = "Species") +
  #scale_fill_viridis_d(begin = 0.25, end = 0.75, option = "D", name = "Species") +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = 'none', 
        legend.title.align = 0.5,
        legend.background = element_blank(), 
        axis.title.x = element_markdown(family = "sans")) + 
  labs(x =  paste("\U03B4",
                  "<sup>", unique(iso_density3[[3]]$neutron), "</sup>",unique(iso_density3[[3]]$element), " (\u2030)",
                  sep = ""), 
       y = "Density")


d34s_density <- iso_density3

d15n_density <- iso_density2 

d13c_density <- iso_density1


iso_biplot_1 <- ggplot() + 
  geom_point(data = fish1, aes(x = d15n, y = d13c,
                                                     fill = species),
             shape = 21, colour = "black", 
             stroke = 0.2,
             size = 1.5, alpha = 0.70) +
  scale_color_brewer(palette = "Set1", name = "Species") +
  scale_fill_brewer(palette = "Set1", name = "Species") +
  #scale_color_viridis_d(begin = 0.25, end = 0.75, option = "D", name = "Species") +
  geom_smooth(data = fish1, aes(x = d15n, y = d13c, fill=species, col=species), method="lm", alpha=0.2) +
  #scale_x_continuous(breaks = rev(seq(-20, -39, -1))) +
  #scale_y_continuous(breaks = seq(5, 17, 1)) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(y = expression(paste(delta^{13}, "C (\u2030)")), 
       x = expression(paste(delta^{15}, "N (\u2030)")))

iso_biplot_2 <- ggplot() + 
  geom_point(data = fish1, aes(x = d15n, y = d34s,
                                                     fill = species),
             shape = 21, colour = "black", 
             stroke = 0.2,
             size = 1.5, alpha = 0.70) +
  scale_color_brewer(palette = "Set1", name = "Species") +
  scale_fill_brewer(palette = "Set1", name = "Species") +
  #scale_color_viridis_d(begin = 0.25, end = 0.75, option = "D", name = "Species") +
  geom_smooth(data = fish1, aes(x = d15n, y = d34s, fill=species, col=species), method="lm", alpha=0.2) +
  #scale_x_continuous(breaks = rev(seq(-20, -39, -1))) +
  #scale_y_continuous(breaks = seq(5, 17, 1)) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(y = expression(paste(delta^{34}, "S (\u2030)")), 
       x = expression(paste(delta^{15}, "N (\u2030)")))

iso_biplot_3 <- ggplot() + 
  geom_point(data = fish1, aes(x = d13c, y = d34s,
                                                     fill = species),
             shape = 21, colour = "black", 
             stroke = 0.2,
             size = 1.5, alpha = 0.70) +
  scale_color_brewer(palette = "Set1", name = "Species") +
  scale_fill_brewer(palette = "Set1", name = "Species") +
  #scale_color_viridis_d(begin = 0.25, end = 0.75, option = "D", name = "Species") +
  geom_smooth(data = fish1, aes(x = d13c, y = d34s, fill=species, col=species), method="lm", alpha=0.2) +
  #scale_x_continuous(breaks = rev(seq(-20, -39, -1))) +
  #scale_y_continuous(breaks = seq(5, 17, 1)) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(y = expression(paste(delta^{34}, "S (\u2030)")), 
       x = expression(paste(delta^{13}, "C (\u2030)")))


plot <- d15n_density + ellipse_plots_1 + ellipse_plots_2  + iso_biplot_1 + d13c_density + ellipse_plots_3 + iso_biplot_2 + iso_biplot_3 +  d34s_density +
  plot_annotation(tag_levels = "a", 
                  tag_suffix = ")") 


print(plot)
ggsave("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/png/Niche_rover_plots.png", plot, 
       device = png(width = 800, height = 600))

#####

# ---- Estimate niche similarities with nicheROVER ----- 
over_stat <- overlap(fish_par, nreps = nsample, nprob = 1000, 
                     alpha = 0.95)


# ---- Extract niche similarities and convert into a dataframe ---- 
over_stat_df <- over_stat %>% 
  as_tibble(rownames = "species_a") %>% 
  dplyr::mutate(
    id = 1:nrow(.), 
    species_a = factor(species_a, 
                       level = c("fin", "hybrid", "blue"))
  ) %>% 
  pivot_longer(cols = -c(id, species_a), 
               names_to = "species_b", 
               values_to = "mc_nr")  %>% 
  separate(species_b, into = c("species_c", "sample_number"), 
           sep = "\\.") %>% 
  dplyr::select(-id) %>% 
  dplyr::rename(species_b = species_c) %>% 
  dplyr::mutate(
    species_b =  factor(species_b, 
                        level = c("fin", "hybrid" ,"blue")
    ), 
    mc_nr_perc = mc_nr * 100
  )



over_sum <- over_stat_df %>% 
  dplyr::group_by(species_a, species_b) %>% 
  dplyr::summarise(
    mean_mc_nr = round(mean(mc_nr_perc), digits = 2),
    qual_2.5 = round(quantile(mc_nr_perc, probs = 0.025, na.rm = TRUE), digits = 2), 
    qual_97.5 = round(quantile(mc_nr_perc, probs = 0.975, na.rm = TRUE), digits = 2)
  ) %>% 
  dplyr::ungroup() %>% 
  pivot_longer(cols = -c(species_a, species_b, mean_mc_nr), 
               names_to = "percentage", 
               values_to = "mc_nr_qual") %>% 
  dplyr::mutate(
    percentage = as.numeric(str_remove(percentage, "qual_"))
  ) 

# ---- plot niche similarities ----- 
#we calculate the mean overlap metric between each species. Remember that the overlap metric is directional, 
#such that it represents the probability that an individual from Species A will be found in the niche of Species B
#Species A is along the rows and Species B is along columns. 
#The plots represent the posterior probability that an individual from the species indicated by the row will be found within the niche of the species indicated by the column header

Niche_overlap <- ggplot(data = over_stat_df, aes(x = mc_nr_perc)) + 
  geom_density(aes(fill = species_a)) + 
  geom_vline(data = over_sum, aes(xintercept = mean_mc_nr), 
             colour = "black", linewidth = 1) +
  geom_vline(data = over_sum, aes(xintercept = mc_nr_qual), 
             colour = "black", linewidth = 0.5, linetype = 2) +
  scale_fill_brewer(palette = "Set1", name = "Species") +
  ggh4x::facet_grid2(species_a ~ species_b, 
                     independent = "y",
                     scales = "free_y") + 
  theme_bw() + 
  theme(
    panel.grid = element_blank(), 
    axis.text = element_text(colour = "black"), 
    legend.background = element_blank(),
    strip.background = element_blank()
  ) +
  labs(x = paste("Overlap Probability (%)", "\u2013", 
                 "Niche Region Size: 95%"), 
       y = "p(Percent Overlap | X)")


ggsave("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/png/Niche_overlap.png", Niche_overlap, 
       device = png(width = 800, height = 600))

# ---- Estimate niche size ----- 
niche_size <- sapply(fish_par, function(spec) {
  apply(spec$Sigma, 3, niche.size)
})


# ---- Convert niche size into a dataframe ----- 
niche_size_df <- niche_size %>% 
  as_tibble() %>% 
  mutate(
    id = 1:nrow(.)
  ) %>% 
  pivot_longer(
    cols = -id, 
    names_to = "species", 
    values_to = "niche_size"
  ) %>% 
  mutate(
    id = 1:nrow(.), 
    species = factor(species,
                     level = c("fin", "hybrid", "blue"))
  )


# ---- Extract niche size mean, sd, and sem ----- 
niche_size_mean <- niche_size_df %>% 
  dplyr::group_by(species) %>% 
  dplyr::summarise(
    mean_niche = round(mean(niche_size), digits = 2), 
    sd_niche = round(sd(niche_size), digits = 2), 
    sem_niche = round(sd(niche_size) / sqrt(n()), digits = 2)
  ) %>% 
  dplyr::ungroup()


# ---- Plot niche sizes as violin plots ------ 

Niche_size <- ggplot(data = niche_size_df) + 
  geom_violin(
    aes(x = species, y = niche_size, fill=species),
    width = 0.2) + 
  scale_fill_brewer(palette = "Set1", name = "Species") +
  geom_point(data = niche_size_mean, aes(x = species, y = mean_niche)) +
  geom_errorbar(data = niche_size_mean, aes(x = species, 
                                            ymin = mean_niche  - sem_niche, 
                                            ymax = mean_niche  + sem_niche), 
                width = 0.05) + 
  theme_bw(base_size = 15) + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black"),
        legend.position = c(0.1, 0.9), aspect.ratio = 4/4, legend.title = element_blank()) + 
  labs(x = "", 
       y = paste("Niche Size (\u2030)"))

ggsave("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/png/Niche_size.png", Niche_size, 
       device = png(width = 800, height = 600))

#permanova Script Xon

#Instalación del paquete

install.packages("vegan")
library(vegan)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") #aquest ja esta instalat
library(pairwiseAdonis)


Y <- fish1[c(1,2,3)]
fish1

?permute::how

adonis(Y ~ fish1$species, permutations = 10000)   

#pairwise.adonis(mydata[,2:4],mydata$Species)

pairwise.adonis(Y, fish1$species, sim.method = "bray", p.adjust.m = "none",perm=10000)

#pairwise comparisons 

#hybrid vs fin 1 0.020948309 68.54684 0.2243415 9.999e-05  9.999e-05 ***
#hybrid vs blue  1 0.004720921 13.82104 0.0854094 9.999e-05  9.999e-05 ***
#fin vs blue  1 0.017188392 57.51735 0.2174426 9.999e-05  9.999e-05 ***

#Analysis of multivariate homogeneity of group dispersions (variances). Used as a measure of multivariate beta diversity. Bedadisper is a multivariate analogue of Levene’s test for homogeneity of variances. 
dev.off()
plot(dispersion, hull=FALSE, ellipse=TRUE, col= clrs) ##sd ellipse

#WITH ALL POINTS AVAILABLE

# 2. Import data
#merge <- read_excel("~/Desktop/Doctorat/Analisis_isotops_barbes/Projecte_barbes_clima/Environmental_vars/All_merged_time_calibrated_2013_to_2022_non_Suess_cor.xlsx")
merge <- read_excel("/Users/marcruizisagales/Desktop/Papers/Paper_hybrids/Blue_whale_results/PER_CAROL_BO_0_point_alignments_29_juny_total_-4_amb_sofre_blue_25_april_2024.xlsx")


#seleccionar individus i fer algo d'interpolacio amb el na.approx?
merge_good_ind <- merge %>% filter(Whale %in% c("F13065", "F13066", "F13068","F13073", "F13076", "F13129", "F18022", "F18098", "BMUS_2010","BMUS_1990"))


merge_good_ind <-merge_good_ind[complete.cases(merge_good_ind$dC), ] # remove the rows with no d13C values in order to be able to correct for the Suess effect

merge_good_ind$Whale
merge_good_ind$Data_capt
merge_good_ind$Data_capt[378:468] <- as.Date(c(rep("2010/08/23", 61),rep("1987/08/23", 30))) #segona 


#date
point_size <- 2
combined_df_blue <- NULL;

for (i in 1:length(unique(merge_good_ind$Whale))) { 
  df_i <- dplyr::filter(merge_good_ind, Whale == unique(merge_good_ind$Whale)[i]) # We select one fin whale individual from 2013
  growth_rate <- 15.5  # Baleen plate growth rate (centimeters per year)
  last_sample <- unique(ymd(df_i$Data_capt)) # Date of last keratin sample (catch date)
  df_i <- 
    df_i %>%
    dplyr::mutate(days = Cm * 365.25 / growth_rate) %>% # Number of days between samples starting at Cm 0
    dplyr::mutate(days = days - days[1]) %>% # Number of days between samples starting at Cm -4
    dplyr::mutate(rev_days = days - days[length(days)]) %>% # Number of days starting at maximum Cm
    dplyr::mutate(sample_date = last_sample + rev_days) %>% # Estimated date
    dplyr::mutate(year = year(sample_date)) %>% # Estimated date year
    dplyr::mutate(year_rev = rev(sample_date)) # Estimated date retrospectively (real date)
  
  df_i$Year_from_sample_date <- format(as.Date(df_i$year_rev), "%Y") # year from the date estimated retrospectively
  df_i$Year_month <- format(as.Date(df_i$year_rev), "%Y/%m") # year and month from the date estimated retrospectively
  df_i$Month_from_sample_date <- format(as.Date(df_i$year_rev), "%m") # month from the date estimated retrospectively
  
  # plot top and bottom axis
  isop <- ggplot(df_i, aes(x = rev(Cm))) + 
    theme_classic(base_size = 16) + 
    xlab("Baleen sample (cm from youngest sample)") + 
    scale_x_continuous(breaks = df_i$Cm, 
                       labels = rev(df_i$Cm),
                       sec.axis = sec_axis(~., 
                                           breaks = rev(df_i$Cm),
                                           labels = rev(df_i$sample_date),
                                           name = "Date"))
  
  #isop <- ggplot(df_i, aes(x = rev(sample_date))) + xlab("Date") + scale_x_continuous(breaks = seq(min(df_i$sample_date), max(df_i$sample_date), by = "1 month"),labels = format(seq(min(rev(df_i$sample_date)), max(rev(df_i$sample_date)), by = "1 month"), "%b-%Y"))
  
  # plot carbon stable isotope values
  isop <- isop + geom_line(aes(y = dC), col = "black") + 
    geom_point(aes(y = dC), size = point_size, shape = 21, fill = "black", col = "black") +
    ylab(expression(paste(delta^{13}, "C (\u2030)")))+ theme_article()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  a <- mean(na.omit(df_i$dN))# add the d15N data but need to scale it to d13C
  a
  b <- mean(na.omit(df_i$dC))# add the d15N data but need to scale it to d13C
  b
  v_shift <- a - b
  
  # plot nitrogen stable isotope values
  isop <- isop + geom_line(aes(y = dN - v_shift), col = "darkgrey") + 
    geom_point(aes(y = (dN - v_shift)), size = point_size, 
               shape = 21, fill = "grey", col = "darkgrey")
  
  # plot x and y axis
  isop <- isop + 
    scale_y_continuous(sec.axis = sec_axis(~.+v_shift, 
                                           name = expression(paste(delta^{15}, "N (\u2030)"))),limits = c(-30, -10)) + labs(title= unique(merge_good_ind$Whale)[i])
  new_df <- df_i
  combined_df_blue <- rbind(combined_df_blue, new_df)
  print(isop)
}


ggplot(data=combined_df_blue, aes(combined_df_blue$sample_date, combined_df_blue$dN, color=Whale)) + geom_line()
ggplot(data=combined_df_blue, aes(combined_df_blue$sample_date, combined_df_blue$dC, color=Whale)) + geom_line()
ggplot(data=combined_df_blue, aes(combined_df_blue$sample_date, combined_df_blue$dS, color=Whale)) + geom_line()

ggplot(data=df, aes(df$sample_date, df$d13c, color=Whale)) + geom_line()


#merge1 <- merge %>% 
#dplyr::mutate(Species = ifelse(Whale %in% c("F13129", "F18022", "F18098"), "hybrid", "fin")) 

#merge2 <- dplyr::filter(merge1, Species == "fin") # remove blue-fin hybrids

length(unique(merge1$Whale)) # we have 30 fin whale individuals (being non-hybrids)

#merge3 <- dplyr::filter(merge2, Whale != "F18025") # select all the whales that are not F18025 (non well-sampled)



# 3. Suess and Laws effect correction

merge3 <-combined_df_blue

merge3$id <-seq.int(nrow(merge3)) # add a sequence number to each row

merge3 <- merge3 %>% add_column(region = "Subpolar North Atlantic") # add the Subpolar North Atlantic region to each row

df1 <- merge3 # rename merge3 to df1

subset <- df1[c("id", "dC","Year_from_sample_date","region")] # select the merge3 column names (for RSuess) and name it subset

names(subset) <- c("id", "d13c","year","region") # rename subset columnames

subset$year <- as.numeric(subset$year) # define year as a numeric variable

subset <- as.data.frame(subset) # define subset as a dataframe

df2 <- SuessR(data=subset, correct.to = 2022) # correct the Suess and Laws effect to the year 2022

fw.data_d13cor <- merge(df1,df2,by="id") # merge df1 and df2

library(openxlsx) # Save
#write.xlsx(fw.data_d13cor, file = "All_merged_time_calibrated_2013_to_2022_Suess_cor.xlsx")

# ---- Bring in example data; will need to replace with your own ----

df <- fw.data_d13cor
names(df) <- c("id","Cm","d15n","dC","d34s","Whale_1","Whale","Sex","Status","Talla_fetus_(cm)","Sexe_fetus","Length","Talla_(peus)",          
               "Data_capt","Lat","Long","Edat","Year","species","days","rev_days","sample_date","year.x","year_rev","Year_from_sample_date","Year_month","Month_from_sample_date","region",                
               "year.y","d13c.uncor","Laws.cor","Suess.cor","net.cor","d13c")

unique(df$Whale)

###################################################################################
#fer-ho amb nicheROVER
###################################################################################

library(nicheROVER)

fish <- df[,c(3,34,5,7)]
fish1 <-fish[complete.cases(fish), ]

data(fish) # 4 fish, 3 isotopes
aggregate(fish1[1:3], fish1[4], mean) # isotope means calculated for each speci

# fish data
data(fish)

# generate parameter draws from the "default" posteriors of each fish
nsamples <- 1e3
system.time({
  fish.par <- tapply(1:nrow(fish1), fish1$Whale,
                     function(ii) niw.post(nsamples = nsamples, X = fish1[ii,1:3]))
})

# various parameter plots
#clrs <- c("red", "blue", "orange") # colors for each species
clrs <- brewer.pal(n = 10, name = "RdYlBu")

# mu1 (del15N), mu2 (del13C), and Sigma12
par(mar = c(4, 4, .5, .1)+.1, mfrow = c(1,3))
niche.par.plot(fish.par, col = clrs, plot.index = 1)
niche.par.plot(fish.par, col = clrs, plot.index = 2)
niche.par.plot(fish.par, col = clrs, plot.index = 1:2)
legend("topleft", legend = names(fish.par), fill = clrs)

# all mu (del15N, del13C, del34S)
niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = FALSE)
legend("topleft", legend = names(fish.par), fill = clrs)

# all mu and Sigma
par(mar = c(4.2, 4.2, 2, 1)+.1)
niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = TRUE)
legend("topright", legend = names(fish.par), fill = clrs)

# 2-d projections of 10 niche regions
#clrs <- c("blue","red", "orange") # colors for each species
nsamples <- 10
fish.par <- tapply(1:nrow(fish1), fish1$Whale,
                   function(ii) niw.post(nsamples = nsamples, X = fish1[ii,1:3]))

# format data for plotting function
fish.data <- tapply(1:nrow(fish1), fish1$Whale, function(ii) X = fish1[ii,1:3])

niche.plot(niche.par = fish.par, niche.data = fish.data, pfrac = .05, alpha = 0.4,
           iso.names = expression(delta^{15}*N, delta^{13}*C, delta^{34}*S),
           col = clrs, xlab = expression("Isotope Ratio (per mil)"))

# niche overlap plots for 95% niche region sizes
nsamples <- 1000
fish.par <- tapply(1:nrow(fish1), fish1$Whale,
                   function(ii) niw.post(nsamples = nsamples, X = fish1[ii,1:3]))

# Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher accuracy.
# the variable over.stat can be supplied directly to the overlap.plot function

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1e3, alpha = c(.95, 0.99))

#The mean overlap metrics calculated across iteratations for both niche 
#region sizes (alpha = .95 and alpha = .99) can be calculated and displayed in an array.
over.mean <- apply(over.stat, c(1:2,4), mean)*100
round(over.mean, 2)

over.cred <- apply(over.stat*100, c(1:2, 4), quantile, prob = c(.025, .975), na.rm = TRUE)
round(over.cred[,,,1]) # display alpha = .95 niche region

# Overlap plot.Before you run this, make sure that you have chosen your 
#alpha level.
#clrs <- c("blue","red", "orange") # colors for each species
over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1e3, alpha = .95)
overlap.plot(over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")

# posterior distribution of (mu, Sigma) for each species
nsamples <- 1000
fish.par <- tapply(1:nrow(fish1), fish1$Whale,
                   function(ii) niw.post(nsamples = nsamples, X = fish1[ii,1:3]))

# posterior distribution of niche size by species
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})

# point estimate and standard error
rbind(est = colMeans(fish.size),
      se = apply(fish.size, 2, sd))

# boxplots
#clrs <- c("blue","red", "orange") # colors for each species
boxplot(fish.size, col = clrs, pch = 16, cex = .5,
        ylab = "Niche Size", xlab = "Species")

#per sp

library(nicheROVER)

fish <- df[,c(3,34,5,19)]
fish1 <-fish[complete.cases(fish), ]

data(fish) # 4 fish, 3 isotopes
aggregate(fish1[1:3], fish1[4], mean) # isotope means calculated for each speci

# fish data
data(fish)

# generate parameter draws from the "default" posteriors of each fish
nsamples <- 1e3
system.time({
  fish.par <- tapply(1:nrow(fish1), fish1$species,
                     function(ii) niw.post(nsamples = nsamples, X = fish1[ii,1:3]))
})

# various parameter plots
clrs <- c( "blue", "red", "orange") # colors for each species
#clrs <- brewer.pal(n = 10, name = "RdYlBu")

# mu1 (del15N), mu2 (del13C), and Sigma12
par(mar = c(4, 4, .5, .1)+.1, mfrow = c(1,3))
niche.par.plot(fish.par, col = clrs, plot.index = 1)
niche.par.plot(fish.par, col = clrs, plot.index = 2)
niche.par.plot(fish.par, col = clrs, plot.index = 1:2)
legend("topleft", legend = names(fish.par), fill = clrs)

# all mu (del15N, del13C, del34S)
niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = FALSE)
legend("topleft", legend = names(fish.par), fill = clrs)

# all mu and Sigma
par(mar = c(4.2, 4.2, 2, 1)+.1)
niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = TRUE)
legend("topright", legend = names(fish.par), fill = clrs)

# 2-d projections of 10 niche regions
#clrs <- c("blue","red", "orange") # colors for each species
nsamples <- 10
fish.par <- tapply(1:nrow(fish1), fish1$species,
                   function(ii) niw.post(nsamples = nsamples, X = fish1[ii,1:3]))

# format data for plotting function
fish.data <- tapply(1:nrow(fish1), fish1$species, function(ii) X = fish1[ii,1:3])

niche.plot(niche.par = fish.par, niche.data = fish.data, pfrac = .05, alpha= 0.4,
           iso.names = expression(delta^{15}*N, delta^{13}*C, delta^{34}*S),
           col = clrs, xlab = expression("Isotope Ratio (per mil)"))

# niche overlap plots for 95% niche region sizes
nsamples <- 10
fish.par <- tapply(1:nrow(fish1), fish1$species,
                   function(ii) niw.post(nsamples = nsamples, X = fish1[ii,1:3]))

# Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher accuracy.
# the variable over.stat can be supplied directly to the overlap.plot function

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1e3, alpha = c(.95, 0.99))

#The mean overlap metrics calculated across iteratations for both niche 
#region sizes (alpha = .95 and alpha = .99) can be calculated and displayed in an array.
over.mean <- apply(over.stat, c(1:2,4), mean)*100
round(over.mean, 2)

over.cred <- apply(over.stat*100, c(1:2, 4), quantile, prob = c(.025, .975), na.rm = TRUE)
round(over.cred[,,,1]) # display alpha = .95 niche region

# Overlap plot.Before you run this, make sure that you have chosen your 
#alpha level.
#clrs <- c("blue","red", "orange") # colors for each species
over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1e3, alpha = .95)
overlap.plot(over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")

# posterior distribution of (mu, Sigma) for each species
nsamples <- 1000
fish.par <- tapply(1:nrow(fish1), fish1$species,
                   function(ii) niw.post(nsamples = nsamples, X = fish1[ii,1:3]))

# posterior distribution of niche size by species
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})

# point estimate and standard error
rbind(est = colMeans(fish.size),
      se = apply(fish.size, 2, sd))

# boxplots
#clrs <- c("blue","red", "orange") # colors for each species
boxplot(fish.size, col = clrs, pch = 16, cex = .5,
        ylab = "Niche Size", xlab = "Species")



