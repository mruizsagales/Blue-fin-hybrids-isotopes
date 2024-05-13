

#--------------------------------------------------------------------------------
# FIN-BLUE WHALE HYBRIDS ISOTOPIC NICHES
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
library(ellipse)
library(ggtext)
library(here)
library(nicheROVER) 
library(purrr)
library(stringr)
library(tidyr)
library(plyr)
library(lme4)
library(climwin)
library(readr)
library(mgcv)
library(lmerTest)
library(car)
library(nortest)
library(effects)
library(nlme)
library(lvmisc)
library(ggeffects)



# 2. Import data
merge <- read_excel("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/0_point_alignments_29_juny_total_-4_amb_sofre_blue.xlsx")
unique(merge$Species) #fin hybrid blue

which(is.na(merge$Data_capt))

merge[c(1611:1671),]$Data_capt <- as.Date(c(rep("2010/08/23", 61))) 
merge[c(1672:1701),]$Data_capt <- as.Date(c(rep("1987/08/23", 30))) #segona data inventada

# per no tenir valors NA fem aproxx en al carboni
which(is.na(merge$dC))

merge$dC <- na.approx(merge$dC )






#cabviar els dies en que no tenim la data per el dies on si que la tenim

#merge <-merge[complete.cases(merge$dC), ] # remove the rows with no d13C values in order to be able to correct for the Suess effect
#merge <-merge[complete.cases(merge$dN), ]

#merge1 <- merge %>% 
  dplyr::mutate(Species = ifelse(Whale %in% c("F13129", "F18022", "F18098"), "hybrid", "fin")) 

#merge2 <- dplyr::filter(merge1, Species == "fin") # remove blue-fin hybrids

length(unique(merge1$Whale)) # we have 30 fin whale individuals (being non-hybrids)

#merge3 <- dplyr::filter(merge2, Whale != "F18025") # select all the whales that are not F18025 (non well-sampled)

# Create a logical index for rows with non-missing values in all three columns
#index <- !is.na(merge1$dC) & !is.na(merge1$dC) & !is.na(merge1$dS)

# Select rows with non-missing values in all three columns
#selected_rows <- merge1[index, ]

selected_rows <- dplyr::filter(merge, Whale %in% c("F18098", "F18022", "F13065", "F13129" ,"F13066", "F13068", "F13073", "F13076", "BMUS_1990", "BMUS_2010"))

# 2. Time-estimation


unique(selected_rows$Whale) # (n=7; "F18098" "F18022" "F13065" "F13066" "F13068" "F13073" "F13076")

point_size <- 2
combined_df_fin_blue_hybrids <- NULL;

for (i in 1:length(unique(selected_rows$Whale))) { 
  df_i <- dplyr::filter(selected_rows, Whale == unique(selected_rows$Whale)[i]) # We select one fin whale individual from 2013
  growth_rate <- 16  # Baleen plate growth rate (centimeters per year)
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
                                           name = expression(paste(delta^{15}, "N (\u2030)"))),limits = c(-30, -10)) + labs(title= unique(combined_df_fin_blue_hybrids$Whale)[i])
  new_df <- df_i
  combined_df_fin_blue_hybrids <- rbind(combined_df_fin_blue_hybrids, new_df)
  print(isop)
}
combined_df_fin_blue_hybrids

unique(combined_df_fin_blue_hybrids$Species)

Nitrogen <- ggplot(combined_df_fin_blue_hybrids, aes(year_rev, dN, color=Whale)) + geom_line() + xlab("Year") + ylab(expression(paste(delta^{15}, "N (\u2030)"))) + theme_bw() + theme(aspect.ratio = 3/4) 
Carboni <- ggplot(combined_df_fin_blue_hybrids, aes(year_rev, dC, color=Whale)) + geom_line() + xlab("Year") + ylab(expression(paste(delta^{13}, "C(\u2030)"))) + theme_bw() + theme(aspect.ratio = 3/4)
Sofre <- ggplot(combined_df_fin_blue_hybrids, aes(year_rev, dS, color=Whale)) + geom_line() + xlab("Year") + ylab(expression(paste(delta^{34}, "S (\u2030)"))) + theme_bw() + theme(aspect.ratio = 3/4)

plot_fin <- Nitrogen/Carboni/Sofre  + 
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "a", 
                  tag_suffix = ")")
ggsave("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/png/Time_series_N_C_S_fin_and_hybrids.png", plot_fin, 
       device = png(width = 450, height = 800))

# 3. Suess and Laws effect correction

merge3 <-combined_df_fin_blue_hybrids

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
names(df)
names(df) <- c("id","Cm","d15n","dC", "dS", "Whale_1","Whale","Sex","Status","Talla_fetus_(cm)","Sexe_fetus","Length","Talla_(peus)",          
"Data_capt","Lat","Long","Edat","Year","species","days","rev_days","sample_date","year.x","year_rev","Year_from_sample_date","Year_month","Month_from_sample_date","region",                
"year.y","d13c.uncor","Laws.cor","Suess.cor","net.cor","d13c")





# 2. Statistics

library(lmerTest)
library(lme4)
library(effects)

#Nitrogen
#merge %>% group_by(Species)  %>% dplyr::summarise(sp_mean = mean(dN, na.rm=T), sd= sd(dN, na.rm=T), max= max(dN, na.rm=T), min= min(dN, na.rm=T))
#ggplot(merge, aes(Species, dN)) + geom_boxplot()

#df <- merge %>% dplyr::filter(Whale %in% c("F18098", "F18022", "F13065", "F13066", "F13068", "F13073", "F13076","F13129", "BMUS_2010", "BMUS_1990"))
unique(df$Whale)
data_without_N_Nas <- df %>% drop_na(d15n)
lmer_dN <- lmer(d15n ~ 1 + species + (1|Whale), data=data_without_N_Nas,REML = TRUE,control = lmerControl(optimizer ="Nelder_Mead"))
summary(lmer_dN)
plot(allEffects(lmer_dN))

data_without_C_Nas <- df %>% drop_na(d13c)
lmer_dN <- lmer(d13c ~ 1 + species + (1|Whale), data=data_without_C_Nas,REML = TRUE,control = lmerControl(optimizer ="Nelder_Mead"))
summary(lmer_dN)
plot(allEffects(lmer_dN))

data_without_S_Nas <- df %>% drop_na(dS)
lmer_dN <- lmer(dS ~ 1 + species + (1|Whale), data=data_without_S_Nas,REML = TRUE,control = lmerControl(optimizer ="Nelder_Mead"))
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
boost <- lme4::bootMer(lmer_dC, f, nsim=200, type="semiparametric", use.u=TRUE)
plot(boost, index= 1)

# Carbon
lmer_dC <- lmer(d13C ~ 1 + species + (1|Whale), data=df)
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
lmer_dS <- lmer(dS ~ 1 + species + (1|Whale), data=df)
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


###################################################################################
#fer-ho amb nicheROVER
###################################################################################

library(nicheROVER)
unique(df$species)
data_per_niche_rover <- df[,c(3,5,34,19)]
aggregate(data_per_niche_rover[c(1,2,3)], data_per_niche_rover[4], mean, na.rm = TRUE) # isotope means calculated for each speci

data_per_niche_rover_NA_out <- na.omit(data_per_niche_rover)
data_per_niche_rover_NA_out <- data.frame(species= data_per_niche_rover_NA_out$species,
                                          D15N= data_per_niche_rover_NA_out$d15n, 
                                          D13C= data_per_niche_rover_NA_out$d13c, 
                                          D34S = data_per_niche_rover_NA_out$dS
                                          )
# generate parameter draws from the "default" posteriors of each fish
nsamples <- 1e3
system.time({
  blue_fin_hybrid.par <- tapply(1:nrow(data_per_niche_rover_NA_out), data_per_niche_rover_NA_out$species,
                     function(ii) niw.post(nsamples = nsamples, X = data_per_niche_rover_NA_out[ii,c(2:4)]))
})

# various parameter plots
clrs <- c("red", "blue", "green") # colors for each species

# mu1 (del15N), mu2 (del13C), and Sigma12
par(mar = c(4, 4, .5, .1)+.1, mfrow = c(1,3))
niche.par.plot(blue_fin_hybrid.par, col = clrs, plot.index = 1)
niche.par.plot(blue_fin_hybrid.par, col = clrs, plot.index = 2)
niche.par.plot(blue_fin_hybrid.par, col = clrs, plot.index = 1:2)
legend("topleft", legend = names(blue_fin_hybrid.par), fill = clrs)

# all mu (del15N, del13C, del34S)
niche.par.plot(blue_fin_hybrid.par, col = clrs, plot.mu = TRUE, plot.Sigma = FALSE)
legend("topleft", legend = names(blue_fin_hybrid.par), fill = clrs)

# all mu and Sigma
par(mar = c(4.2, 4.2, 2, 1)+.1)
niche.par.plot(blue_fin_hybrid.par, col = clrs, plot.mu = TRUE, plot.Sigma = TRUE)
legend("topright", legend = names(blue_fin_hybrid.par), fill = clrs)

# 2-d projections of 10 niche regions
clrs <- c("red", "blue", "green") # colors for each species
nsamples <- 10
blue_fin_hybrid.par <- tapply(1:nrow(data_per_niche_rover_NA_out), data_per_niche_rover_NA_out$species,
                   function(ii) niw.post(nsamples = nsamples, X = data_per_niche_rover_NA_out[ii,2:4]))

# format data for plotting function
blue_fin_hybrid.data <- tapply(1:nrow(data_per_niche_rover_NA_out), data_per_niche_rover_NA_out$species, function(ii) X = data_per_niche_rover_NA_out[ii,2:4])

niche.plot(niche.par = blue_fin_hybrid.par, niche.data = blue_fin_hybrid.data, pfrac = .05,
           iso.names = expression(delta^{15}*N, delta^{13}*C, delta^{34}*S),
           col = clrs, xlab = expression("Isotope Ratio (per mil)"))

# niche overlap plots for 95% niche region sizes
nsamples <- 1000
fish.par <- tapply(1:nrow(data_per_niche_rover_NA_out), data_per_niche_rover_NA_out$species,
                   function(ii) niw.post(nsamples = nsamples, X = data_per_niche_rover_NA_out[ii,2:4]))

# Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher accuracy.
# the variable over.stat can be supplied directly to the overlap.plot function

over.stat <- overlap(blue_fin_hybrid.par, nreps = nsamples, nprob = 1e3, alpha = c(.95, 0.99))

#The mean overlap metrics calculated across iteratations for both niche 
#region sizes (alpha = .95 and alpha = .99) can be calculated and displayed in an array.
over.mean <- apply(over.stat, c(1:2,4), mean)*100
round(over.mean, 2)

over.cred <- apply(over.stat*100, c(1:2, 4), quantile, prob = c(.025, .975), na.rm = TRUE)
round(over.cred[,,,1]) # display alpha = .95 niche region

# Overlap plot.Before you run this, make sure that you have chosen your 
#alpha level.
clrs <- c("red", "blue", "green") # colors for each species
over.stat <- overlap(blue_fin_hybrid.par, nreps = nsamples, nprob = 1e3, alpha = .95)
overlap.plot(over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")

# posterior distribution of (mu, Sigma) for each species
nsamples <- 1000
fish.par <- tapply(1:nrow(data_per_niche_rover_NA_out), data_per_niche_rover_NA_out$species,
                   function(ii) niw.post(nsamples = nsamples, X = data_per_niche_rover_NA_out[ii,2:4]))

# posterior distribution of niche size by species
blue_fin_hybrid.size <- sapply(blue_fin_hybrid.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})

# point estimate and standard error
rbind(est = colMeans(blue_fin_hybrid.size),
      se = apply(blue_fin_hybrid.size, 2, sd))

# boxplots

clrs <- c("black", "red", "blue", "green") # colors for each species
boxplot(blue_fin_hybrid.size, col = clrs, pch = 16, cex = .5,
        ylab = "Niche Size", xlab = "Species")


###################################################################################
#fer-ho amb ggplot
###################################################################################
data_per_niche_rover <- df[,c(3,5,34,19,7)]

data_per_niche_rover_NA_out <- data.frame(species= data_per_niche_rover$species,
                                          D15N= data_per_niche_rover$d15n, 
                                          D13C= data_per_niche_rover$d13c
                                          #D34S = data_per_niche_rover$dS
                                          
)


library(SIBER, quietly = TRUE,
        verbose = FALSE,
        logical.return = FALSE)

library(viridis)
palette(viridis(4))

per_C_N <- na.omit(data_per_niche_rover_NA_out)
nrow(per_C_N)
max(per_C_N$D15N)
min(per_C_N$D15N)

max(per_C_N$D13C)
min(per_C_N$D13C)

#N amb C
# create the siber object
data_per_niche_rover_NA_out_SIBER <- data.frame(
  iso1= per_C_N$D13C, 
  iso2= per_C_N$D15N,
  group= per_C_N$species,
  community = rep("a",467))


####
#make ellipses overlap 
library(SIBER)
siber.example_N_overlap_b <- createSiberObject(data_per_niche_rover_NA_out_SIBER)
par(mfrow=c(1,1))
plotSiberObject(siber.example_N_overlap_b,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille')
)

# In this example, I will calculate the overlap between ellipses for groups 2
# and 3 in community 1 (i.e. the green and yellow open circles of data).

# The first ellipse is referenced using a character string representation where 
# in "x.y", "x" is the community, and "y" is the group within that community.
# So in this example: community 1, group 2
ellipse1 <- "a.hybrid" 

# Ellipse two is similarly defined: community 1, group3
ellipse2 <- "a.blue"

ellipse3 <- "a.fin"

# the overlap betweeen the corresponding 95% prediction ellipses is given by:
ellipse95.overlap <- maxLikOverlap(ellipse1, ellipse2, siber.example_N_overlap_b, 
                                   p.interval = 0.95, n = 100)
ellipse95.overlap

# so in this case, the overlap as a proportion of the non-overlapping area of 
# the two ellipses, would be
prop.95.over <- ellipse95.overlap[3] / (ellipse95.overlap[2] + 
                                          ellipse95.overlap[1] -
                                          ellipse95.overlap[3])
prop.95.over


#Overlap Nitrogen Carboni hybrid-blue = 6.237677 (0.6841612 proportion)
##Overlap Nitrogen Carboni hybrid-fin = 4.845884 (0.4489594 proportion)
##Overlap Nitrogen Carboni blue-fin = 4.834358 (0.4143555 proportion)


####
demo_data <- per_C_N %>% mutate(group = per_C_N$species, 
                                        community = rep("a",467),
                                        d13C = per_C_N$D13C, 
                                        d15N = per_C_N$D15N,
                                        .keep = "unused") 
dev.off()
first.plot <- ggplot(data = demo_data, 
                     aes(x = d13C, 
                         y = d15N)) + 
  geom_point(aes(color = group), 
             #shape = ifelse(demo_data$group == "hybrid", 16, 1),
             shape=1,
             size = 2, 
             alpha=0.75) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  #theme(text = element_text(size=16)) + 
  scale_color_manual(values = c("blue","darkgrey","red")) + theme_article(base_size = 20) + theme(aspect.ratio = 1/1) + xlim(c(-21.5,-17.5)) + ylim(c(6,13.5))

# And print our plot to screen
print(first.plot)

classic.first.plot <- first.plot + theme_classic() + 
  theme(text = element_text(size=18))

first.plot.dens_b <- ggplot(data = demo_data) + 
  geom_density(aes(x = D15N, 
                   fill = group), 
               alpha = 0.35, 
               linewidth = 0.8) +
  theme(text = element_text(size=16)) + 
  scale_fill_manual(name = "", labels = c("Blue","Fin", "Hybrid"), values = c("blue","red","green")) + xlim(6,13.5) + xlab(expression(paste(delta^{15}, "N (\u2030)"))) + ylab("Density") + theme_article(base_size = 20) + theme(aspect.ratio = 1/1, legend.position = c(0.2,0.9) )

classic.first.plot <- first.plot + theme_classic() + 
  theme(text = element_text(size=18))

first.plot.dens2_b <- ggplot(data = demo_data) + 
  geom_density(aes(x = D13C, 
                   fill = group), 
               alpha = 0.35, 
               linewidth = 0.8) +
  theme(text = element_text(size=16)) + 
  scale_fill_manual(name = "", labels = c("Blue","Fin", "Hybrid"), values = c("blue","red","green")) + xlim(-21.5,-17.5) + xlab(expression(paste(delta^{13}, "C (\u2030)"))) + ylab("Density") + theme_article(base_size = 20) + theme(aspect.ratio = 1/1, legend.position = "none" )

classic.first.plot <- first.plot + theme_classic() + 
  theme(text = element_text(size=18))

# and print to screen
print(classic.first.plot)

# use our ellipse function to generate the ellipses for plotting

# decide how big an ellipse you want to draw
p.ell <- 0.40
p.ell2 <- 0.95
# create our plot based on first.plot above adding the stat_ellipse() geometry.
# We specify thee ellipse to be plotted using the polygon geom, with fill and
# edge colour defined by our column "group", using the normal distribution and
# with a quite high level of transparency on the fill so we can see the points
# underneath. In order to get different ellipses plotted by both columns "group"
# and "community" we have to use the interaction() function to ensure both are
# considered in the aes(group = XYZ) specification. Note also the need to
# specify the scale_fill_viridis_d() function as the mapping of colors for
# points and lines is separate to filled objects and we want them to match.
ellipse.plot1 <- first.plot + 
  stat_ellipse(aes(group = interaction(group, community),
                   color = group), 
               alpha = 0, 
               linewidth=1.25,
               level = p.ell,
               type = "norm",
               #linetype = "dashed",
               geom = "polygon") + 
  stat_ellipse(aes(group = interaction(group, community),
                   color = group), 
               alpha = 0, 
               linewidth=0.75,
               linetype = "dotted",
               level = p.ell2,
               type = "norm",
               geom = "polygon") + 
  scale_color_manual(name = "", labels = c("Blue", "Fin", "Hybrid"), values = c("blue","red","green"))  +
  guides(color = guide_legend(override.aes = list(shape = NA))) + theme(legend.position = c(0.2,0.2))

print(ellipse.plot1)

#C amb S
data_per_niche_rover_NA_out <- data.frame(species= data_per_niche_rover$species,
                                          #D15N= data_per_niche_rover$d15n, 
                                          D13C= data_per_niche_rover$d13c,
                                          D34S = data_per_niche_rover$dS
                                          
)

per_C_N <- na.omit(data_per_niche_rover_NA_out)
nrow(per_C_N)

max(per_C_N$D34S)
min(per_C_N$D34S)

max(per_C_N$D13C)
min(per_C_N$D13C)

# create the siber object
data_per_niche_rover_NA_out_SIBER <- data.frame(
  iso1= per_C_N$D13C, 
  iso2= per_C_N$D34S,
  group= per_C_N$species,
  community = rep("a",232))

####
#make ellipses overlap 
library(SIBER)
siber.example_N_overlap_b <- createSiberObject(data_per_niche_rover_NA_out_SIBER)
par(mfrow=c(1,1))
plotSiberObject(siber.example_N_overlap_b,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille')
)

# In this example, I will calculate the overlap between ellipses for groups 2
# and 3 in community 1 (i.e. the green and yellow open circles of data).

# The first ellipse is referenced using a character string representation where 
# in "x.y", "x" is the community, and "y" is the group within that community.
# So in this example: community 1, group 2
ellipse1 <- "a.hybrid" 

# Ellipse two is similarly defined: community 1, group3
ellipse2 <- "a.blue"

ellipse3 <- "a.fin"

# the overlap betweeen the corresponding 95% prediction ellipses is given by:
ellipse95.overlap <- maxLikOverlap(ellipse1, ellipse3, siber.example_N_overlap_b, 
                                   p.interval = 0.95, n = 100)
ellipse95.overlap

# so in this case, the overlap as a proportion of the non-overlapping area of 
# the two ellipses, would be
prop.95.over <- ellipse95.overlap[3] / (ellipse95.overlap[2] + 
                                          ellipse95.overlap[1] -
                                          ellipse95.overlap[3])
prop.95.over


##Overlap Carboni Sofre hybrid-fin = 1.908452 (0.3445563 proportion)



####

demo_data <- per_C_N %>% mutate(group = per_C_N$species, 
                                community = rep("a",232),
                                d13C = per_C_N$D13C, 
                                d34S = per_C_N$D34S,
                                .keep = "unused") 

first.plot <- ggplot(data = demo_data, 
                     aes(x = d34S, 
                         y = d13C)) + 
  geom_point(aes(color = group), size = 2, 
             #shape = ifelse(demo_data$group == "hybrid", 16, 1),
             shape=1,
             alpha=0.75) +
  xlab(expression(paste(delta^{34}, "S (\u2030)"))) +
  ylab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme(text = element_text(size=16)) + 
  scale_color_manual(values = c("red","green")) + theme_article(base_size = 20) + theme(aspect.ratio = 1/1) + ylim(c(-21.5,-17.5)) + xlim(c(15,20))

# And print our plot to screen
print(first.plot)

first.plot.dens3_b <- ggplot(data = demo_data) + 
  geom_density(aes(x = D34S, 
                   fill = group), 
               alpha = 0.35, 
               linewidth = 0.8) +
  theme(text = element_text(size=16)) + 
  scale_fill_manual(name = "", labels = c("Fin", "Hybrid"), values = c("red","green")) + xlim(15,20) + xlab(expression(paste(delta^{34}, "S (\u2030)"))) + ylab("Density") + theme_article(base_size = 20) + theme(aspect.ratio = 1/1, legend.position = "none")

classic.first.plot <- first.plot + theme_classic() + 
  theme(text = element_text(size=18))

classic.first.plot <- first.plot + theme_classic() + 
  theme(text = element_text(size=18))

# and print to screen
print(classic.first.plot)

# use our ellipse function to generate the ellipses for plotting

# decide how big an ellipse you want to draw
p.ell <- 0.40
p.ell2 <- 0.95
# create our plot based on first.plot above adding the stat_ellipse() geometry.
# We specify thee ellipse to be plotted using the polygon geom, with fill and
# edge colour defined by our column "group", using the normal distribution and
# with a quite high level of transparency on the fill so we can see the points
# underneath. In order to get different ellipses plotted by both columns "group"
# and "community" we have to use the interaction() function to ensure both are
# considered in the aes(group = XYZ) specification. Note also the need to
# specify the scale_fill_viridis_d() function as the mapping of colors for
# points and lines is separate to filled objects and we want them to match.
ellipse.plot2 <- first.plot + 
  stat_ellipse(aes(group = interaction(group, community),
                   color = group), 
               alpha = 0, 
               linewidth=1.25,
               level = p.ell,
               type = "norm",
               #linetype = "dashed",
               geom = "polygon") + 
  stat_ellipse(aes(group = interaction(group, community),
                   color = group), 
               alpha = 0, 
               linewidth=0.75,
               linetype = "dotted",
               level = p.ell2,
               type = "norm",
               geom = "polygon") + 
  scale_color_manual(name = "", labels = c("Fin", "Hybrid"), values = c("red","green"))  +
  guides(color = guide_legend(override.aes = list(shape = NA))) + theme(legend.position = "none") 

print(ellipse.plot2)


#N amb S
data_per_niche_rover_NA_out <- data.frame(species= data_per_niche_rover$species,
                                          D15N= data_per_niche_rover$d15n, 
                                          #D13C= data_per_niche_rover$d13c,
                                          D34S = data_per_niche_rover$dS
                                          
)
per_C_N <- na.omit(data_per_niche_rover_NA_out)
nrow(per_C_N)
# create the siber object
data_per_niche_rover_NA_out_SIBER <- data.frame(
  iso1= per_C_N$D15N, 
  iso2= per_C_N$D34S,
  group= per_C_N$species,
  community = rep("a",231))

####
#make ellipses overlap 
library(SIBER)
siber.example_N_overlap_b <- createSiberObject(data_per_niche_rover_NA_out_SIBER)
par(mfrow=c(1,1))
plotSiberObject(siber.example_N_overlap_b,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille')
)

# In this example, I will calculate the overlap between ellipses for groups 2
# and 3 in community 1 (i.e. the green and yellow open circles of data).

# The first ellipse is referenced using a character string representation where 
# in "x.y", "x" is the community, and "y" is the group within that community.
# So in this example: community 1, group 2
ellipse1 <- "a.hybrid" 

# Ellipse two is similarly defined: community 1, group3
ellipse2 <- "a.blue"

ellipse3 <- "a.fin"

# the overlap betweeen the corresponding 95% prediction ellipses is given by:
ellipse95.overlap <- maxLikOverlap(ellipse1, ellipse3, siber.example_N_overlap_b, 
                                   p.interval = 0.95, n = 100)
ellipse95.overlap

# so in this case, the overlap as a proportion of the non-overlapping area of 
# the two ellipses, would be
prop.95.over <- ellipse95.overlap[3] / (ellipse95.overlap[2] + 
                                          ellipse95.overlap[1] -
                                          ellipse95.overlap[3])
prop.95.over


##Overlap Nitrogen Sofre hybrid-fin = 2.673845 (0.209836 proportion)



####

demo_data <- per_C_N %>% mutate(group = per_C_N$species, 
                                community = rep("a",231),
                                D15N = per_C_N$D15N, 
                                d34S = per_C_N$D34S,
                                .keep = "unused") 

first.plot <- ggplot(data = demo_data, 
                     aes(x = d34S, 
                         y = D15N)) + 
  geom_point(aes(color = group), size = 2, 
             #shape = ifelse(demo_data$group == "hybrid", 16, 1),
             shape=1,
             alpha=0.75) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{34}, "S (\u2030)"))) + 
  theme(text = element_text(size=16)) + 
  scale_color_manual(values = c("red","green")) + theme_article(base_size = 20) + theme(aspect.ratio = 1/1)  + ylim(c(6,13.5)) + xlim(c(15,20))

# And print our plot to screen
print(first.plot)

classic.first.plot <- first.plot + theme_classic() + 
  theme(text = element_text(size=18))

# and print to screen
print(classic.first.plot)

# use our ellipse function to generate the ellipses for plotting

# decide how big an ellipse you want to draw
p.ell <- 0.40
p.ell2 <- 0.95
# create our plot based on first.plot above adding the stat_ellipse() geometry.
# We specify thee ellipse to be plotted using the polygon geom, with fill and
# edge colour defined by our column "group", using the normal distribution and
# with a quite high level of transparency on the fill so we can see the points
# underneath. In order to get different ellipses plotted by both columns "group"
# and "community" we have to use the interaction() function to ensure both are
# considered in the aes(group = XYZ) specification. Note also the need to
# specify the scale_fill_viridis_d() function as the mapping of colors for
# points and lines is separate to filled objects and we want them to match.
ellipse.plot3 <- first.plot + 
  stat_ellipse(aes(group = interaction(group, community),
                   color = group), 
               alpha = 0, 
               linewidth=1.25,
               level = p.ell,
               type = "norm",
               #linetype = "dashed",
               geom = "polygon") + 
  stat_ellipse(aes(group = interaction(group, community),
                   color = group), 
               alpha = 0, 
               linewidth=0.75,
               linetype = "dotted",
               level = p.ell2,
               type = "norm",
               geom = "polygon") + 
  scale_color_manual(name = "", labels = c("Fin", "Hybrid"), values = c("red","green"))  +
  guides(color = guide_legend(override.aes = list(shape = NA))) + theme(legend.position = "none")

print(ellipse.plot3)



ellipse.plot1 + ellipse.plot2 + ellipse.plot3

#pells

Pells_N_C_S <- read_excel("~/Desktop/Pells_N_C_S.xls", 
                          sheet = "Pells")



merge3 <-Pells_N_C_S

merge3$id <-seq.int(nrow(merge3)) # add a sequence number to each row

merge3 <- merge3 %>% add_column(region = "Subpolar North Atlantic") # add the Subpolar North Atlantic region to each row

df1 <- merge3 # rename merge3 to df1

subset <- df1[c("id", "dC","Any","region")] # select the merge3 column names (for RSuess) and name it subset

names(subset) <- c("id", "d13c","year","region") # rename subset columnames

subset$year <- as.numeric(subset$year) # define year as a numeric variable

subset <- as.data.frame(subset) # define subset as a dataframe

df2 <- SuessR(data=subset, correct.to = 2022) # correct the Suess and Laws effect to the year 2022

fw.data_d13cor <- merge(df1,df2,by="id") # merge df1 and df2

#N amb C
data_per_niche_rover_NA_out <- data.frame(species= fw.data_d13cor$Species,
                                          D15N= fw.data_d13cor$dN, 
                                          D13C= fw.data_d13cor$d13c.cor
                                          #D34S = fw.data_d13cor$dS
                                          
)
per_C_N <- na.omit(data_per_niche_rover_NA_out)
nrow(per_C_N)
# create the siber object
data_per_niche_rover_NA_out_SIBER <- data.frame(
  iso1= per_C_N$D15N, 
  iso2= per_C_N$D13C,
  group= per_C_N$species,
  community = rep("a",43))

####
#make ellipses overlap 
library(SIBER)
data_per_niche_rover_NA_out_SIBER_no_hyb <- dplyr::filter(data_per_niche_rover_NA_out_SIBER, !group == "hybrid")
siber.example_N_overlap_b <- createSiberObject(data_per_niche_rover_NA_out_SIBER_no_hyb)
par(mfrow=c(1,1))
plotSiberObject(siber.example_N_overlap_b,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille')
)

# In this example, I will calculate the overlap between ellipses for groups 2
# and 3 in community 1 (i.e. the green and yellow open circles of data).

# The first ellipse is referenced using a character string representation where 
# in "x.y", "x" is the community, and "y" is the group within that community.
# So in this example: community 1, group 2

# Ellipse two is similarly defined: community 1, group3
ellipse1 <- "a.blue"

ellipse2 <- "a.fin"

# the overlap betweeen the corresponding 95% prediction ellipses is given by:
ellipse95.overlap <- maxLikOverlap(ellipse1, ellipse2, siber.example_N_overlap_b, 
                                   p.interval = 0.95, n = 100)
ellipse95.overlap

# so in this case, the overlap as a proportion of the non-overlapping area of 
# the two ellipses, would be
prop.95.over <- ellipse95.overlap[3] / (ellipse95.overlap[2] + 
                                          ellipse95.overlap[1] -
                                          ellipse95.overlap[3])
prop.95.over


##Overlap Nitrogen Carbon hybrid-fin = 2.003857 (0.3695364 proportion)



####

demo_data <- per_C_N %>% mutate(group = per_C_N$species, 
                                community = rep("a",43),
                                D15N = per_C_N$D15N, 
                                D13C = per_C_N$D13C,
                                .keep = "unused") 

first.plot <- ggplot(data = demo_data, 
                     aes(x = D13C, 
                         y = D15N)) + 
  geom_point(aes(color = group), size = 2, 
             shape = ifelse(demo_data$group == "hybrid", 16, 1),
             #shape=1,
             alpha=0.75) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme(text = element_text(size=16)) + 
  scale_color_manual(values = c("blue","red","green")) + theme_article(base_size = 20) + theme(aspect.ratio = 1/1)  + xlim(c(-21.5,-17.5)) + ylim(c(6,13.5))

# And print our plot to screen
print(first.plot)


first.plot.dens <- ggplot(data = demo_data) + 
  geom_density(aes(x = D15N, 
                   fill = group), 
               alpha = 0.35, 
               linewidth = 0.8) +
  theme(text = element_text(size=16)) + 
  scale_fill_manual(name = "", labels = c("Blue","Fin", "Hybrid"), values = c("blue","red","green")) + xlim(6,13.5) + xlab(expression(paste(delta^{15}, "N (\u2030)"))) + ylab("Density") + theme_article(base_size = 20) + theme(aspect.ratio = 1/1, legend.position = c(0.2,0.9) )

classic.first.plot <- first.plot + theme_classic() + 
  theme(text = element_text(size=18))

first.plot.dens2 <- ggplot(data = demo_data) + 
  geom_density(aes(x = D13C, 
                   fill = group), 
               alpha = 0.35, 
               linewidth = 0.8) +
  theme(text = element_text(size=16)) + 
  scale_fill_manual(name = "", labels = c("Blue","Fin", "Hybrid"), values = c("blue","red","green")) + xlim(-21.5,-17.5) + xlab(expression(paste(delta^{13}, "C (\u2030)"))) + ylab("Density") + theme_article(base_size = 20) + theme(aspect.ratio = 1/1, legend.position = "none")

classic.first.plot <- first.plot + theme_classic() + 
  theme(text = element_text(size=18))

# and print to screen
print(classic.first.plot)

# use our ellipse function to generate the ellipses for plotting

# decide how big an ellipse you want to draw
p.ell <- 0.40
p.ell2 <- 0.95
# create our plot based on first.plot above adding the stat_ellipse() geometry.
# We specify thee ellipse to be plotted using the polygon geom, with fill and
# edge colour defined by our column "group", using the normal distribution and
# with a quite high level of transparency on the fill so we can see the points
# underneath. In order to get different ellipses plotted by both columns "group"
# and "community" we have to use the interaction() function to ensure both are
# considered in the aes(group = XYZ) specification. Note also the need to
# specify the scale_fill_viridis_d() function as the mapping of colors for
# points and lines is separate to filled objects and we want them to match.
ellipse.plot1_p <- first.plot + 
  stat_ellipse(aes(group = interaction(group, community),
                   color = group), 
               alpha = 0, 
               linewidth=1.25,
               level = p.ell,
               type = "norm",
               #linetype = "dashed",
               geom = "polygon") + 
  stat_ellipse(aes(group = interaction(group, community),
                   color = group), 
               alpha = 0, 
               linewidth=0.75,
               linetype = "dotted",
               level = p.ell2,
               type = "norm",
               geom = "polygon") + 
  scale_color_manual(name = "", labels = c("Blue","Fin", "Hybrid"), values = c("blue","red","green"))  +
  guides(color = guide_legend(override.aes = list(shape = NA))) + theme(legend.position = c(0.2,0.2))

print(ellipse.plot1_p)

#N amb S
data_per_niche_rover_NA_out <- data.frame(species= fw.data_d13cor$Species,
                                          D15N= fw.data_d13cor$dN, 
                                          #D13C= data_per_niche_rover$d13c,
                                          D34S = fw.data_d13cor$dS
                                          
)



per_C_N <- na.omit(data_per_niche_rover_NA_out)
nrow(per_C_N)
# create the siber object
data_per_niche_rover_NA_out_SIBER <- data.frame(
  iso1= per_C_N$D15N, 
  iso2= per_C_N$D34S,
  group= per_C_N$species,
  community = rep("a",39))

####
#make ellipses overlap 
library(SIBER)
data_per_niche_rover_NA_out_SIBER_no_hyb <- dplyr::filter(data_per_niche_rover_NA_out_SIBER, !group == "hybrid")
siber.example_N_overlap_b <- createSiberObject(data_per_niche_rover_NA_out_SIBER_no_hyb)
par(mfrow=c(1,1))
plotSiberObject(siber.example_N_overlap_b,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille')
)

# In this example, I will calculate the overlap between ellipses for groups 2
# and 3 in community 1 (i.e. the green and yellow open circles of data).

# The first ellipse is referenced using a character string representation where 
# in "x.y", "x" is the community, and "y" is the group within that community.
# So in this example: community 1, group 2

# Ellipse two is similarly defined: community 1, group3
ellipse1 <- "a.blue"

ellipse2 <- "a.fin"

# the overlap betweeen the corresponding 95% prediction ellipses is given by:
ellipse95.overlap <- maxLikOverlap(ellipse1, ellipse2, siber.example_N_overlap_b, 
                                   p.interval = 0.95, n = 100)
ellipse95.overlap

# so in this case, the overlap as a proportion of the non-overlapping area of 
# the two ellipses, would be
prop.95.over <- ellipse95.overlap[3] / (ellipse95.overlap[2] + 
                                          ellipse95.overlap[1] -
                                          ellipse95.overlap[3])
prop.95.over


##Overlap Nitrogen Sofre hybrid-fin = 2.356818 (0.4773415 proportion)



####

demo_data <- per_C_N %>% mutate(group = per_C_N$species, 
                                community = rep("a",39),
                                D15N = per_C_N$D15N, 
                                d34S = per_C_N$D34S,
                                .keep = "unused") 

first.plot <- ggplot(data = demo_data, 
                     aes(x = d34S, 
                         y = D15N)) + 
  geom_point(aes(color = group), size = 2, 
             shape = ifelse(demo_data$group == "hybrid", 16, 1),
             #shape=1,
             alpha=0.75) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{34}, "S (\u2030)"))) + 
  theme(text = element_text(size=16)) + 
  scale_color_manual(values = c("blue","red","green")) + theme_article(base_size = 20) + theme(aspect.ratio = 1/1)  + ylim(c(6,13.5)) + xlim(c(15,20))

first.plot.dens3 <- ggplot(data = demo_data) + 
  geom_density(aes(x = D34S, 
                   fill = group), 
               alpha = 0.35, 
               linewidth = 0.8) +
  theme(text = element_text(size=16)) + 
  scale_fill_manual(name = "", labels = c("Blue","Fin", "Hybrid"), values = c("blue","red","green")) + xlim(17.5,20.5) + xlab(expression(paste(delta^{34}, "S (\u2030)"))) + ylab("Density") + theme_article(base_size = 20) + theme(aspect.ratio = 1/1, legend.position = "none" )

classic.first.plot <- first.plot + theme_classic() + 
  theme(text = element_text(size=18))

# And print our plot to screen
print(first.plot)

classic.first.plot <- first.plot + theme_classic() + 
  theme(text = element_text(size=18))

# and print to screen
print(classic.first.plot)

# use our ellipse function to generate the ellipses for plotting

# decide how big an ellipse you want to draw
p.ell <- 0.40
p.ell2 <- 0.95
# create our plot based on first.plot above adding the stat_ellipse() geometry.
# We specify thee ellipse to be plotted using the polygon geom, with fill and
# edge colour defined by our column "group", using the normal distribution and
# with a quite high level of transparency on the fill so we can see the points
# underneath. In order to get different ellipses plotted by both columns "group"
# and "community" we have to use the interaction() function to ensure both are
# considered in the aes(group = XYZ) specification. Note also the need to
# specify the scale_fill_viridis_d() function as the mapping of colors for
# points and lines is separate to filled objects and we want them to match.
ellipse.plot2_p <- first.plot + 
  stat_ellipse(aes(group = interaction(group, community),
                   color = group), 
               alpha = 0, 
               linewidth=1.25,
               level = p.ell,
               type = "norm",
               #linetype = "dashed",
               geom = "polygon") + 
  stat_ellipse(aes(group = interaction(group, community),
                   color = group), 
               alpha = 0, 
               linewidth=0.75,
               linetype = "dotted",
               level = p.ell2,
               type = "norm",
               geom = "polygon") + 
  scale_color_manual(name = "", labels = c("Blue","Fin", "Hybrid"), values = c("blue","red","green"))  +
  guides(color = guide_legend(override.aes = list(shape = NA))) + theme(legend.position = "none")

print(ellipse.plot2_p)

#C amb S
data_per_niche_rover_NA_out <- data.frame(species= fw.data_d13cor$Species,
                                          #D15N= fw.data_d13cor$dN, 
                                          D13C= fw.data_d13cor$d13c.cor,
                                          D34S = fw.data_d13cor$dS
                                          
)
per_C_N <- na.omit(data_per_niche_rover_NA_out)
nrow(per_C_N)
# create the siber object
data_per_niche_rover_NA_out_SIBER <- data.frame(
  iso1= per_C_N$D13C, 
  iso2= per_C_N$D34S,
  group= per_C_N$species,
  community = rep("a",39))

####
#make ellipses overlap 
library(SIBER)
data_per_niche_rover_NA_out_SIBER_no_hyb <- dplyr::filter(data_per_niche_rover_NA_out_SIBER, !group == "hybrid")
siber.example_N_overlap_b <- createSiberObject(data_per_niche_rover_NA_out_SIBER_no_hyb)
par(mfrow=c(1,1))
plotSiberObject(siber.example_N_overlap_b,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille')
)

# In this example, I will calculate the overlap between ellipses for groups 2
# and 3 in community 1 (i.e. the green and yellow open circles of data).

# The first ellipse is referenced using a character string representation where 
# in "x.y", "x" is the community, and "y" is the group within that community.
# So in this example: community 1, group 2

# Ellipse two is similarly defined: community 1, group3
ellipse1 <- "a.blue"

ellipse2 <- "a.fin"

# the overlap betweeen the corresponding 95% prediction ellipses is given by:
ellipse95.overlap <- maxLikOverlap(ellipse1, ellipse2, siber.example_N_overlap_b, 
                                   p.interval = 0.95, n = 100)
ellipse95.overlap

# so in this case, the overlap as a proportion of the non-overlapping area of 
# the two ellipses, would be
prop.95.over <- ellipse95.overlap[3] / (ellipse95.overlap[2] + 
                                          ellipse95.overlap[1] -
                                          ellipse95.overlap[3])
prop.95.over


##Overlap Carboni Sofre hybrid-fin = 1.235348 (0.3872907 proportion)



####

demo_data <- per_C_N %>% mutate(group = per_C_N$species, 
                                community = rep("a",39),
                                D13C = per_C_N$D13C, 
                                D34S = per_C_N$D34S,
                                .keep = "unused") 

first.plot <- ggplot(data = demo_data, 
                     aes(x = D34S, 
                         y = D13C)) + 
  geom_point(aes(color = group), size = 2, 
             shape = ifelse(demo_data$group == "hybrid", 16, 1),
             #shape=1,
             alpha=0.75) +
  ylab(expression(paste(delta^{13}, "C (\u2030)"))) +
  xlab(expression(paste(delta^{34}, "S (\u2030)"))) + 
  theme(text = element_text(size=16)) + 
  scale_color_manual(values = c("blue","red","green")) + theme_article(base_size = 20) + theme(aspect.ratio = 1/1)  + ylim(c(-21.5,-17.5)) + xlim(c(15,20))

# And print our plot to screen
print(first.plot)

classic.first.plot <- first.plot + theme_classic() + 
  theme(text = element_text(size=18))

# and print to screen
print(classic.first.plot)

# use our ellipse function to generate the ellipses for plotting

# decide how big an ellipse you want to draw
p.ell <- 0.40
p.ell2 <- 0.95
# create our plot based on first.plot above adding the stat_ellipse() geometry.
# We specify thee ellipse to be plotted using the polygon geom, with fill and
# edge colour defined by our column "group", using the normal distribution and
# with a quite high level of transparency on the fill so we can see the points
# underneath. In order to get different ellipses plotted by both columns "group"
# and "community" we have to use the interaction() function to ensure both are
# considered in the aes(group = XYZ) specification. Note also the need to
# specify the scale_fill_viridis_d() function as the mapping of colors for
# points and lines is separate to filled objects and we want them to match.
ellipse.plot3_p <- first.plot + 
  stat_ellipse(aes(group = interaction(group, community),
                   color = group), 
               alpha = 0, 
               linewidth=1.25,
               level = p.ell,
               type = "norm",
               #linetype = "dashed",
               geom = "polygon") + 
  stat_ellipse(aes(group = interaction(group, community),
                   color = group), 
               alpha = 0, 
               linewidth=0.75,
               linetype = "dotted",
               level = p.ell2,
               type = "norm",
               geom = "polygon") + 
  scale_color_manual(name = "", labels = c("Blue","Fin", "Hybrid"), values = c("blue","red","green"))  +
  guides(color = guide_legend(override.aes = list(shape = NA))) + theme(legend.position = "none")

print(ellipse.plot3_p)

ellipse.plot1/ellipse.plot1_p 
ellipse.plot2/ellipse.plot2_p
ellipse.plot3/ellipse.plot3_p

ellipse.plot1 +ellipse.plot3+ellipse.plot2
ellipse.plot1_p +ellipse.plot2_p+ellipse.plot3_p

first.plot.dens_b + first.plot.dens2_b + first.plot.dens3_b

first.plot.dens + first.plot.dens2 + first.plot.dens3



library(SIBER)
library(dplyr)
library(ggplot2)


siber.example <- createSiberObject(data_per_niche_rover_NA_out_SIBER)

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

# plot the raw data
par(mfrow=c(1,1))
plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = F, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "o",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'(\u2030)'),
                ylab = expression({delta}^15*N~'(\u2030)')
)


# You can add more ellipses by directly calling plot.group.ellipses()
# Add an additional p.interval % prediction ellilpse
plotGroupEllipses(siber.example, n = 100, p.interval = 0.40,
                  lty = 1, lwd = 4)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
#plotGroupEllipses(siber.example, n = 100, p.interval = 0.95,
                  ci.mean = T, lty = 1, lwd = 2)

# Calculate sumamry statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siber.example)
print(group.ML)

# add a legend
legend("topright", colnames(group.ML), 
       pch = c(1,1,1,2,2,2), col = c(1:3, 1:3), lty = 1)


## ----fit-bayes----------------------------------------------------------------

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.example, parms, priors)


## ----plot-data, fig.width = 10, fig.height = 6--------------------------------
# 
# ----------------------------------------------------------------
# Plot out some of the data and results
# ----------------------------------------------------------------

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('(\u2030)' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(SEA.B.credibles)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

print(SEA.B.modes)

## ----prob-diff-g12------------------------------------------------------------
Pg1.1_lt_g1.2 <- sum( SEA.B[,1] < SEA.B[,2] ) / nrow(SEA.B)
print(Pg1.1_lt_g1.2)

## ----prob-diff-g13------------------------------------------------------------
Pg1.1_lt_g1.3 <- sum( SEA.B[,1] < SEA.B[,3] ) / nrow(SEA.B)
print(Pg1.1_lt_g1.3)

## ----prob-diff-all------------------------------------------------------------
Pg1.1_lt_g2.1 <- sum( SEA.B[,1] < SEA.B[,4] ) / nrow(SEA.B)
print(Pg1.1_lt_g2.1)

Pg1.2_lt_g1.3 <- sum( SEA.B[,2] < SEA.B[,3] ) / nrow(SEA.B)
print(Pg1.2_lt_g1.3)

Pg1.3_lt_g2.1 <- sum( SEA.B[,3] < SEA.B[,4] ) / nrow(SEA.B)
print(Pg1.3_lt_g2.1)

Pg2.2_lt_g2.3 <- sum( SEA.B[,5] < SEA.B[,6] ) / nrow(SEA.B)
print(Pg2.2_lt_g2.3)


## ----ML-overlap---------------------------------------------------------------

overlap.G1.2.G1.3 <- maxLikOverlap("1.2", "1.3", siber.example, p = 0.95, n =)


## ----ML-overlap-proportions---------------------------------------------------
prop.of.first <- as.numeric(overlap.G1.2.G1.3["overlap"] / overlap.G1.2.G1.3["area.1"])
print(prop.of.first)

prop.of.second <- as.numeric(overlap.G1.2.G1.3["overlap"] / overlap.G1.2.G1.3["area.2"])
print(prop.of.second)

prop.of.both <- as.numeric(overlap.G1.2.G1.3["overlap"] / (overlap.G1.2.G1.3["area.1"] + overlap.G1.2.G1.3["area.2"]))
print(prop.of.both)

## ----bayesian-overlap---------------------------------------------------------
bayes.overlap.G2.G3 <- bayesianOverlap("1.2", "1.3", ellipses.posterior, 
                                       draws = 10, p.interval = 0.95,
                                       n = 360)
print(bayes.overlap.G2.G3)



##
# set the number of posterior samples that will be taken 
nsample <- 1000


# ---- Use map from purrr and niw.post from nicheROVER to estimate niches ----- 
fish_par <- data_per_niche_rover_NA_out %>% 
  split(.$species) %>% 
  purrr::map(~ dplyr::select(., D15N, D13C, D34S)) %>% 
  purrr::map(~niw.post(nsample = nsample, X = .))

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
                     levels = c("fin", "hybrid"))
  ) %>% 
  dplyr::group_by(species) %>% 
  dplyr::mutate(
    sample_number = 1:1000
  ) %>% 
  dplyr::ungroup()




df_mu_long <- df_mu %>% 
  pivot_longer(cols = -c(metric, species, sample_number), 
               names_to = "isotope", 
               values_to = "mu_est") %>% 
  dplyr::mutate(
    element =  dplyr::case_when(
      isotope == "D15N" ~ "N",
      isotope == "D13C" ~ "C",
      isotope == "D34S" ~ "S",
    ), 
    neutron =  dplyr::case_when(
      isotope == "D15N" ~ 15,
      isotope == "D13C" ~ 13,
      isotope == "D34S" ~ 34,
    ) 
  )

# ---- Extract sigma from nicheROVER object ----- 

df_sigma <- map(fish_par, pluck, 2) %>%
  imap(~ as_tibble(.x) %>%
         dplyr::mutate(
           metric = "sigma",
           id = c("D15N", "D13C","D34S"),
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

posterior_plots_pl <- posterior_plots$D15N  +
  theme(legend.position = c(0.18, 0.84)) + 
  posterior_plots$D13C + posterior_plots$D34S 

ggsave("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/png/Posterior_plots.png", posterior_plots_pl, 
       device = png(width = 800, height = 300))

# ---- Prepare sigma dataframes for plotting ----

df_sigma_cn <- df_sigma_cn %>%
  mutate(
    element_id = case_when(
      id == "D15N" ~ "N",
      id == "D13C" ~ "C",
      id == "D34S" ~ "S",
    ),
    neutron_id = case_when(
      id == "D15N" ~ 15,
      id == "D13C" ~ 13,
      id == "D34S" ~ 34,
    ),
    element_iso = case_when(
      isotopes == "D15N" ~ "N",
      isotopes == "D13C" ~ "C",
      isotopes == "D34S" ~ "S",
    ),
    neutron_iso = case_when(
      isotopes == "D15N" ~ 15,
      isotopes == "D13C" ~ 13,
      isotopes == "D34S" ~ 34,
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
p.ell <- 0.95

# ---- create blank list to dump results and category object to loop for 
species_name <- unique(df_sigma_wide$species)



all_ellipses_NC <- list()


# ---- Loop to create ellipse ------ 
for (i in 1:length(species_name)) {
  
  sigma_species <- df_sigma_wide %>% 
    filter(species %in% species_name[i])
  
  sigma_species <- filter(sigma_species, isotopes != "D34S")
  
  mu_species <- df_mu %>% 
    filter(species %in% species_name[i])
  
  ell <- NULL
  post.id <- NULL
  
  for(j in 1:length(unique(sigma_species$sample_number))) {
    
    sigma_ind <- sigma_species %>%
      filter(sample_number %in% sample_number[j]) %>% 
      dplyr::select(D15N, D13C) 
    
    Sigma <- as.matrix(sigma_ind, 2, 2)
    row.names(Sigma) <- c("D15N", "D13C")
    
    mu <- mu_species %>%
      filter(sample_number %in% sample_number[j]) %>% 
      dplyr::select(sample_number, D15N, D13C) %>% 
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
        id == "2" ~ "hybrid"
      ), level = c("fin", "hybrid")
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
               mapping = aes(x = D13C, y = D15N,
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
  
  sigma_species <- filter(sigma_species, isotopes != "D13C")
  
  mu_species <- df_mu %>% 
    filter(species %in% species_name[i])
  
  ell <- NULL
  post.id <- NULL
  
  for(j in 1:length(unique(sigma_species$sample_number))) {
    
    sigma_ind <- sigma_species %>%
      filter(sample_number %in% sample_number[j]) %>% 
      dplyr::select(D15N, D34S) 
    
    Sigma <- as.matrix(sigma_ind, 2, 2)
    row.names(Sigma) <- c("D15N", "D34S")
    
    mu <- mu_species %>%
      filter(sample_number %in% sample_number[j]) %>% 
      dplyr::select(sample_number, D15N, D34S) %>% 
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
        id == "2" ~ "hybrid"
      ), level = c("fin", "hybrid")
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
               mapping = aes(x = D34S, y = D15N,
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
  
  sigma_species <- filter(sigma_species, isotopes != "D15N")
  
  mu_species <- df_mu %>% 
    filter(species %in% species_name[i])
  
  ell <- NULL
  post.id <- NULL
  
  for(j in 1:length(unique(sigma_species$sample_number))) {
    
    sigma_ind <- sigma_species %>%
      filter(sample_number %in% sample_number[j]) %>% 
      dplyr::select(D13C, D34S) 
    
    Sigma <- as.matrix(sigma_ind, 2, 2)
    row.names(Sigma) <- c("D13C", "D34S")
    
    mu <- mu_species %>%
      filter(sample_number %in% sample_number[j]) %>% 
      dplyr::select(sample_number, D13C, D34S) %>% 
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
        id == "2" ~ "hybrid"
      ), level = c("fin", "hybrid")
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
               mapping = aes(x = D34S, y = D13C,
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


iso_long <- data_per_niche_rover_NA_out %>% pivot_longer(cols = -species,
               names_to = "isotope", 
               values_to = "value") %>% 
  mutate(
    element = case_when(
      isotope == "D15N" ~ "N",
      isotope == "D13C" ~ "C",
      isotope == "D34S" ~ "S",
    ), 
    neutron = case_when(
      isotope == "D15N" ~ 15,
      isotope == "D13C" ~ 13,
      isotope == "D34S" ~ 34,
    )
  )

iso_density1 <- iso_long %>% 
  group_split(isotope) 
iso_density1 <- ggplot(data = iso_density1[[1]]) + 
  geom_density(aes(x = value, 
                   fill = species), 
               alpha = 0.35, 
               linewidth = 0.8) +
  scale_fill_brewer(palette = "Set1", name = "Species") +
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
                   fill = species), 
               alpha = 0.35, 
               linewidth = 0.8) +
  scale_fill_brewer(palette = "Set1", name = "Species") +
  #scale_fill_viridis_d(begin = 0.25, end = 0.75, option = "D", name = "Species") +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.15, 0.75), 
        legend.title.align = 0.5,
        legend.background = element_blank(), 
        axis.title.x = element_markdown(family = "sans")) + 
  labs(x =  paste("\U03B4",
                  "<sup>", unique(iso_density2[[2]]$neutron), "</sup>",unique(iso_density2[[2]]$element), " (\u2030)",
                  sep = ""), 
       y = "Density")


iso_density3 <- iso_long %>% 
  group_split(isotope) 
iso_density3 <- ggplot(data = iso_density3[[3]]) + 
      geom_density(aes(x = value, 
                       fill = species), 
                   alpha = 0.35, 
                   linewidth = 0.8) +
  scale_fill_brewer(palette = "Set1", name = "Species") +
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
  geom_point(data = data_per_niche_rover_NA_out, aes(x = D15N, y = D13C,
                            fill = species),
             shape = 21, colour = "black", 
             stroke = 0.8,
             size = 2, alpha = 0.70) +
  scale_color_brewer(palette = "Set1", name = "Species") +
  scale_fill_brewer(palette = "Set1", name = "Species") +
  #scale_color_viridis_d(begin = 0.25, end = 0.75, option = "D", name = "Species") +
  geom_smooth(data = data_per_niche_rover_NA_out, aes(x = D15N, y = D13C, fill=species, col=species), method="lm", alpha=0.5) +
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
  geom_point(data = data_per_niche_rover_NA_out, aes(x = D15N, y = D34S,
                                                     fill = species),
             shape = 21, colour = "black", 
             stroke = 0.8,
             size = 2, alpha = 0.70) +
  scale_color_brewer(palette = "Set1", name = "Species") +
  scale_fill_brewer(palette = "Set1", name = "Species") +
  #scale_color_viridis_d(begin = 0.25, end = 0.75, option = "D", name = "Species") +
  geom_smooth(data = data_per_niche_rover_NA_out, aes(x = D15N, y = D34S, fill=species, col=species), method="lm", alpha=0.5) +
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
  geom_point(data = data_per_niche_rover_NA_out, aes(x = D13C, y = D34S,
                                                     fill = species),
             shape = 21, colour = "black", 
             stroke = 0.8,
             size = 2, alpha = 0.70) +
  scale_color_brewer(palette = "Set1", name = "Species") +
  scale_fill_brewer(palette = "Set1", name = "Species") +
  #scale_color_viridis_d(begin = 0.25, end = 0.75, option = "D", name = "Species") +
  geom_smooth(data = data_per_niche_rover_NA_out, aes(x = D13C, y = D34S, fill=species, col=species), method="lm", alpha=0.5) +
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
                       level = c("fin", "hybrid"))
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
                        level = c("fin", "hybrid")
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
                     level = c("fin", "hybrid"))
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
        legend.position = c(0.9, 0.1), aspect.ratio = 4/4) + 
  labs(x = "", 
       y = paste("Niche Size (\u2030)"))

ggsave("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/png/Niche_size.png", Niche_size, 
       device = png(width = 800, height = 600))

