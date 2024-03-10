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

# 2. Import data
merge <- read_excel("~/Desktop/0_point_alignments_29_juny_total_-4_amb_sofre.xlsx")

#merge <-merge[complete.cases(merge$dC), ] # remove the rows with no d13C values in order to be able to correct for the Suess effect
#merge <-merge[complete.cases(merge$dN), ]

merge1 <- merge %>% 
  dplyr::mutate(Species = ifelse(Whale %in% c("F13129", "F18022", "F18098"), "hybrid", "fin")) 

#merge2 <- dplyr::filter(merge1, Species == "fin") # remove blue-fin hybrids

length(unique(merge1$Whale)) # we have 30 fin whale individuals (being non-hybrids)

#merge3 <- dplyr::filter(merge2, Whale != "F18025") # select all the whales that are not F18025 (non well-sampled)

# Create a logical index for rows with non-missing values in all three columns
index <- !is.na(merge1$dC) & !is.na(merge1$dC) & !is.na(merge1$dS)

# Select rows with non-missing values in all three columns
selected_rows <- merge1[index, ]

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

ggplot(combined_df_fin_blue_hybrids, aes(year_rev, dN, color=Whale)) + geom_line()
ggplot(combined_df_fin_blue_hybrids, aes(year_rev, dC, color=Whale)) + geom_line()
ggplot(combined_df_fin_blue_hybrids, aes(year_rev, dS, color=Whale)) + geom_line() #notar l'increment de so


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
names(df) <- c("id","Cm","d15n","dC", "dS", "Whale_1","Whale","Sex","Status","Talla_fetus_(cm)","Sexe_fetus","Length","Talla_(peus)",          
"Data_capt","Lat","Long","Edat","Year","species","days","rev_days","sample_date","year.x","year_rev","Year_from_sample_date","Year_month","Month_from_sample_date","region",                
"year.y","d13c.uncor","Laws.cor","Suess.cor","net.cor","d13c")


# 2. Statistics
df_hybrids <- filter(df, species == "hybrid")
df_fin <- filter(df, species == "fin")

#Normality
shapiro.test(df_hybrids$d15n)
shapiro.test(df_hybrids$d13c)
shapiro.test(df_hybrids$dS)

library("car")
dev.off()
qqPlot(df_hybrids$d13c) #As all the points fall approximately along this reference line, we can assume normality.

# Homoscedasticity 
bartlett.test(d15n ~ species, df) #different
bartlett.test(d13c ~ species, df)
bartlett.test(dS ~ species, df)

require(ggstatsplot)
library(ggsignif)
df$species <- as.factor(df$species)
df$d15n <- as.double(df$d15n)
ggstatsplot::ggbetweenstats(data= df, x= species, y= d13c)
ggstatsplot::ggbetweenstats(data= df, x= species, y= d15n)
ggstatsplot::ggbetweenstats(data= df, x= species, y= dS)
wilcox.test(df_hybrids$d15n,df_fin$d15n)

result <- t.test(df_hybrids$d15n,df_fin$d15n, var.equal = FALSE)
result

result <- t.test(df_hybrids$d13c,df_fin$d13c, var.equal = TRUE)
result

result <- t.test(df_hybrids$dS,df_fin$dS, var.equal = TRUE)
result

# Calculate Hedges' g effect size
library(effsize)
hedges_g <- effsize::cohen.d(d = df$dS, f= df$species, method = "hedges")

###################################################################################
#fer-ho amb nicheROVER
###################################################################################

library(nicheROVER)

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
clrs <- c("red", "blue", "orange") # colors for each species

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
clrs <- c("red", "blue", "orange") # colors for each species
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
clrs <- c("red", "blue", "orange") # colors for each species
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
clrs <- c("black", "red", "blue", "orange") # colors for each species
boxplot(blue_fin_hybrid.size, col = clrs, pch = 16, cex = .5,
        ylab = "Niche Size", xlab = "Species")


###################################################################################
#fer-ho amb ggplot
###################################################################################


# set the number of posterior samples that will be taken 
nsample <- 1000


# ---- Use map from purrr and niw.post from nicheROVER to estimate niches ----- 
fish_par <- data_per_niche_rover_NA_out %>% 
  split(.$species) %>% 
  map(~ dplyr::select(., D15N, D13C, D34S)) %>% 
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
posterior_plots <- df_mu_long %>%
  split(.$isotope) %>%
  imap(
    ~ ggplot(data = ., aes(x = mu_est)) +
      geom_density(aes(fill = species), alpha = 0.5) +
      scale_fill_viridis_d(begin = 0.25, end = 0.75,
                           option = "D", name = "Species") +
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

posterior_plots$D15N  +
  theme(legend.position = c(0.18, 0.84)) + 
  posterior_plots$D13C + posterior_plots$D34S 

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
      scale_fill_viridis_d(begin = 0.25, end = 0.75,
                           option = "D", name = "Species") +
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

sigma_plots[[1]] +sigma_plots[[2]] +sigma_plots[[3]] +sigma_plots[[4]] +sigma_plots[[5]] +sigma_plots[[6]] +
  theme(legend.position = c(0.1, 0.82))


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
  
  scale_colour_viridis_d(begin = 0.25, end = 0.75, 
                         option = "D", name = "species",
  ) + 
  #scale_x_continuous(breaks = rev(seq(-20, -40, -2))) +
  #scale_y_continuous(breaks = seq(6, 16, 2)) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = expression(paste(delta ^ 15, "N")))


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
  
  scale_colour_viridis_d(begin = 0.25, end = 0.75, 
                         option = "D", name = "species",
  ) + 
  #scale_x_continuous(breaks = rev(seq(-20, -40, -2))) +
  #scale_y_continuous(breaks = seq(6, 16, 2)) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 34, "S")), 
       y = expression(paste(delta ^ 15, "N")))

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
  
  scale_colour_viridis_d(begin = 0.25, end = 0.75, 
                         option = "D", name = "species",
  ) + 
  #scale_x_continuous(breaks = rev(seq(-20, -40, -2))) +
  #scale_y_continuous(breaks = seq(6, 16, 2)) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 34, "S")), 
       y = expression(paste(delta ^ 13, "C")))


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



iso_density <- iso_long %>% 
  group_split(isotope) %>% 
  imap(
    ~ ggplot(data = .) + 
      geom_density(aes(x = value, 
                       fill = species), 
                   alpha = 0.35, 
                   linewidth = 0.8) +
      scale_fill_viridis_d(begin = 0.25, end = 0.75,
                           option = "D", name = "Species") +
      theme_bw(base_size = 10) +
      theme(axis.text = element_text(colour = "black"),
            panel.grid = element_blank(), 
            legend.position = c(0.15, 0.65), 
            legend.title.align = 0.5,
            legend.background = element_blank(), 
            axis.title.x = element_markdown(family = "sans")) + 
      labs(x =  paste("\U03B4",
                      "<sup>", unique(.$neutron), "</sup>",unique(.$element), 
                      sep = ""), 
           y = "Density")
  )

d34s_density <- iso_density[[3]] 

d15n_density <- iso_density[[2]] 

d13c_density <- iso_density[[1]]


iso_biplot_1 <- ggplot() + 
  geom_point(data = data_per_niche_rover_NA_out, aes(x = D15N, y = D13C,
                            fill = species),
             shape = 21, colour = "black", 
             stroke = 0.8,
             size = 2, alpha = 0.70) +
  scale_fill_viridis_d(begin = 0.25, end = 0.75,
                       option = "D", name = "species") +
  scale_color_viridis_d(begin = 0.25, end = 0.75,
                       option = "D", name = "species") +
  geom_smooth(data = data_per_niche_rover_NA_out, aes(x = D15N, y = D13C, fill=species, col=species), method="lm", alpha=0.5) +
  #scale_x_continuous(breaks = rev(seq(-20, -39, -1))) +
  #scale_y_continuous(breaks = seq(5, 17, 1)) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(y = expression(paste(delta ^ 13, "C")), 
       x = expression(paste(delta ^ 15, "N")))

iso_biplot_2 <- ggplot() + 
  geom_point(data = data_per_niche_rover_NA_out, aes(x = D15N, y = D34S,
                                                     fill = species),
             shape = 21, colour = "black", 
             stroke = 0.8,
             size = 2, alpha = 0.70) +
  scale_fill_viridis_d(begin = 0.25, end = 0.75,
                       option = "D", name = "species") +
  scale_color_viridis_d(begin = 0.25, end = 0.75,
                        option = "D", name = "species") +
  geom_smooth(data = data_per_niche_rover_NA_out, aes(x = D15N, y = D34S, fill=species, col=species), method="lm", alpha=0.5) +
  #scale_x_continuous(breaks = rev(seq(-20, -39, -1))) +
  #scale_y_continuous(breaks = seq(5, 17, 1)) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(y = expression(paste(delta ^ 34, "S")), 
       x = expression(paste(delta ^ 15, "N")))

iso_biplot_3 <- ggplot() + 
  geom_point(data = data_per_niche_rover_NA_out, aes(x = D13C, y = D34S,
                                                     fill = species),
             shape = 21, colour = "black", 
             stroke = 0.8,
             size = 2, alpha = 0.70) +
  scale_fill_viridis_d(begin = 0.25, end = 0.75,
                       option = "D", name = "species") +
  scale_color_viridis_d(begin = 0.25, end = 0.75,
                        option = "D", name = "species") +
  geom_smooth(data = data_per_niche_rover_NA_out, aes(x = D13C, y = D34S, fill=species, col=species), method="lm", alpha=0.5) +
  #scale_x_continuous(breaks = rev(seq(-20, -39, -1))) +
  #scale_y_continuous(breaks = seq(5, 17, 1)) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(y = expression(paste(delta ^ 34, "S")), 
       x = expression(paste(delta ^ 13, "C")))


plot <- d15n_density + ellipse_plots_1 + ellipse_plots_2  + iso_biplot_1 + d13c_density + ellipse_plots_3 + iso_biplot_2 + iso_biplot_3 +  d34s_density +
  plot_annotation(tag_levels = "a", 
                  tag_suffix = ")")


print(plot)
ggsave("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/Figure_NicheRover_ggplot.png", plot, 
       device = png(width = 600, height = 450))
#####

# ---- Estimate niche similarities with nicheROVER ----- 
over_stat <- overlap(fish_par, nreps = nsample, nprob = 1000, 
                     alpha = 0.95)


# ---- Extract niche similarities and convert into a dataframe ---- 
over_stat_df <- over_stat %>% 
  as_tibble(rownames = "species_a") %>% 
  mutate(
    id = 1:nrow(.), 
    species_a = factor(species_a, 
                       level = c("fin", "hybrid"))
  ) %>% 
  pivot_longer(cols = -c(id, species_a), 
               names_to = "species_b", 
               values_to = "mc_nr")  %>% 
  separate(species_b, into = c("species_c", "sample_number"), 
           sep = "\\.") %>% 
  select(-id) %>% 
  rename(species_b = species_c) %>% 
  mutate(
    species_b =  factor(species_b, 
                        level = c("fin", "hybrid")
    ), 
    mc_nr_perc = mc_nr * 100
  )



over_sum <- over_stat_df %>% 
  group_by(species_a, species_b) %>% 
  summarise(
    mean_mc_nr = round(mean(mc_nr_perc), digits = 2),
    qual_2.5 = round(quantile(mc_nr_perc, probs = 0.025, na.rm = TRUE), digits = 2), 
    qual_97.5 = round(quantile(mc_nr_perc, probs = 0.975, na.rm = TRUE), digits = 2)
  ) %>% 
  ungroup() %>% 
  pivot_longer(cols = -c(species_a, species_b, mean_mc_nr), 
               names_to = "percentage", 
               values_to = "mc_nr_qual") %>% 
  mutate(
    percentage = as.numeric(str_remove(percentage, "qual_"))
  ) 

# ---- plot niche similarities ----- 
#we calculate the mean overlap metric between each species. Remember that the overlap metric is directional, 
#such that it represents the probability that an individual from Species A will be found in the niche of Species B
#Species A is along the rows and Species B is along columns. 
#The plots represent the posterior probability that an individual from the species indicated by the row will be found within the niche of the species indicated by the column header

ggplot(data = over_stat_df, aes(x = mc_nr_perc)) + 
  geom_density(aes(fill = species_a)) + 
  geom_vline(data = over_sum, aes(xintercept = mean_mc_nr), 
             colour = "black", linewidth = 1) +
  geom_vline(data = over_sum, aes(xintercept = mc_nr_qual), 
             colour = "black", linewidth = 0.5, linetype = 2) +
  scale_fill_viridis_d(begin = 0.25, end = 0.75,
                       option = "D", name = "Species", 
                       alpha = 0.35) + 
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
  group_by(species) %>% 
  summarise(
    mean_niche = round(mean(niche_size), digits = 2), 
    sd_niche = round(sd(niche_size), digits = 2), 
    sem_niche = round(sd(niche_size) / sqrt(n()), digits = 2)
  ) %>% 
  ungroup()


# ---- Plot niche sizes as violin plots ------ 
ggplot(data = niche_size_df) + 
  geom_violin(
    aes(x = species, y = niche_size),
    width = 0.2) + 
  geom_point(data = niche_size_mean, aes(x = species, y = mean_niche)) +
  geom_errorbar(data = niche_size_mean, aes(x = species, 
                                            ymin = mean_niche  - sem_niche, 
                                            ymax = mean_niche  + sem_niche), 
                width = 0.05) + 
  theme_bw(base_size = 15) + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black")) + 
  labs(x = "", 
       y = "Niche Size") 

