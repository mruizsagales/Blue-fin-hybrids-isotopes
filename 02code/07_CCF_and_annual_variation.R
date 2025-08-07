#--------------------------------------------------------------------------------
# Baleen SIA in fin and blue whale hybrids (Marc Ruiz-Sagalés, 6 de juny del 2025)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# 07 - Cross correlation functions and variation along the year between isotopic oscillations
#--------------------------------------------------------------------------------

library(RColorBrewer)

setwd("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes")
source("03_Suess_Laws_effect_correction.R")
df

whale <- df[,c(2,3,34,5,7,19)]
whale1 <-whale[complete.cases(whale), ]

# CCF

ccf_results_N <- list()
ccf_results_C <- list()
ccf_results_S <- list()

whale_names <- unique(whale1$Whale) # name of individuals
num_whales <- length(whale_names) # number of individuals

# Empty matrix to save maximum correlations between the isotopic oscillations of individuals
correlation_matrix_N <- matrix(NA, nrow = num_whales, ncol = num_whales, dimnames = list(whale_names, whale_names))
correlation_matrix_C <- matrix(NA, nrow = num_whales, ncol = num_whales, dimnames = list(whale_names, whale_names))
correlation_matrix_S <- matrix(NA, nrow = num_whales, ncol = num_whales, dimnames = list(whale_names, whale_names))

# loop for each pair of individuals
for (i in 1:num_whales) {
  for (j in 1:num_whales) {
    if (i != j) {
      # dNcorr
      dN_i <- whale1[whale1$Whale == whale_names[i], "d15n"]
      dN_j <- whale1[whale1$Whale == whale_names[j], "d15n"]
      
      ccf_N <- ccf(dN_i, dN_j, lag.max = 5, plot = FALSE)
      correlation_matrix_N[i, j] <- max(abs(ccf_N$acf))  
      
      # dCcorr
      dC_i <- whale1[whale1$Whale == whale_names[i], "d13c"]
      dC_j <- whale1[whale1$Whale == whale_names[j], "d13c"]
      
      ccf_C <- ccf(dC_i, dC_j, lag.max = 5, plot = FALSE)
      correlation_matrix_C[i, j] <- max(abs(ccf_C$acf))
      
      # dScorr
      dS_i <- whale1[whale1$Whale == whale_names[i], "d34s"]
      dS_j <- whale1[whale1$Whale == whale_names[j], "d34s"]
      
      ccf_S <- ccf(dS_i, dS_j, lag.max = 5, plot = FALSE)
      correlation_matrix_S[i, j] <- max(abs(ccf_S$acf))
    }
  }
}

# Substitute the NA's for 1's 
correlation_matrix_N[is.na(correlation_matrix_N)] <- 1
correlation_matrix_C[is.na(correlation_matrix_C)] <- 1
correlation_matrix_S[is.na(correlation_matrix_S)] <- 1

# Corrplots
corrplot(correlation_matrix_N, method = "color", tl.col = "black", type = "lower", 
         outline = TRUE, title = "Cross-correlation (dNcorr)", mar=c(0,0,2,0),
         addCoef.col = "black", number.cex = 0.45, # mida del text
         cl.pos = "r", number.digits = 2)

corrplot(correlation_matrix_C, method = "color", tl.col = "black", type = "lower", 
         outline = TRUE, title = "Cross-correlation (dCcorr)", mar=c(0,0,2,0),
         addCoef.col = "black", number.cex = 0.45, # mida del text
         cl.pos = "r", number.digits = 2)

corrplot(correlation_matrix_S, method = "color", tl.col = "black", type = "lower", 
         outline = TRUE, title = "Cross-correlation (dScorr)", mar=c(0,0,2,0),
         addCoef.col = "black", number.cex = 0.45, # mida del text
         cl.pos = "r", number.digits = 2)

a<- ggcorrplot::ggcorrplot(correlation_matrix_N,method = c("square"),type = c("lower"),ggtheme = theme_article(base_size = 10, base_family = "Optima")) + scale_fill_viridis_c(option = "A", direction = -1, begin=0, end=1) + theme(legend.position = "none")
b<-ggcorrplot::ggcorrplot(correlation_matrix_C,method = c("square"),type = c("lower"),ggtheme = theme_article(base_size = 10, base_family = "Optima")) + scale_fill_viridis_c(option = "A", direction = -1, begin=0, end=1) + theme(legend.position = "none")
c<-ggcorrplot::ggcorrplot(correlation_matrix_S,method = c("square"),type = c("lower"),ggtheme = theme_article(base_size = 10, base_family = "Optima")) + scale_fill_viridis_c(option = "A", direction = -1, begin=0, end=1) 
tot <- a+b+c

# ggsave("/Users/marcruizisagales/Dropbox/00_papers/2023_hybrids_SIA_baleen/Final_plots/CCF.svg", tot, 
#        dpi = 300,  width = 12, height = 8, units = "cm")

# Intraannual variation (d15n, d13c and d34s per yearday and whale)

df$Whale <- factor(df$Whale, levels = c("BMUS_2010", "BMUS_1990", "F13129", "F18022", "F18098","F13065", "F13066", "F13068", "F13073", "F13076"))
df$species <- factor(df$species, levels = c("blue", "hybrid", "fin"))
whale_species <- df %>% distinct(Whale, species) %>% arrange(species)

# create a function to generate light variations per species
generate_variants <- function(base_color, n) {
  lighten(base_color, amount = seq(0, 0.4, length.out = n))  # up to 40% lighter
}

# assign colors per Whale by varying the species base
individual_colors <- whale_species %>%
  group_by(species) %>%
  mutate(color = generate_variants(species_colors[species[1]], n())) %>%
  ungroup() %>%
  select(Whale, color) %>%
  deframe()

a1<-ggplot(df, aes(yday(year_rev), d15n, color=Whale)) + geom_smooth(se=F) + theme_article(base_size = 10, base_family = "Optima") + theme(aspect.ratio = 1, legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()) + facet_wrap(vars(species)) + ylab(expression(paste(delta^{15}, "N (\u2030)"))) + scale_color_brewer(palette = "RdYlBu", direction = -1) + scale_color_manual(values = individual_colors)
b1<-ggplot(df, aes(yday(year_rev), d13c, color=Whale)) + geom_smooth(se=F) + theme_article(base_size = 10, base_family = "Optima") + theme(aspect.ratio = 1, axis.title.x = element_blank(), axis.text.x = element_blank()) + facet_wrap(vars(species)) + ylab(expression(paste(delta^{13}, "C (\u2030)"))) + scale_color_brewer(palette = "RdYlBu", direction = -1) + scale_color_manual(values = individual_colors)
c1<-ggplot(df, aes(yday(year_rev), d34s, color=Whale)) + geom_smooth(se=F) + theme_article(base_size = 10, base_family = "Optima") + theme(aspect.ratio = 1,legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(vars(species)) + xlab("Yearday") + ylab(expression(paste(delta^{34}, "S (\u2030)"))) + scale_color_brewer(palette = "RdYlBu", direction = -1) + scale_color_manual(values = individual_colors)
tot1 <- a1/b1/c1

# ggsave("/Users/marcruizisagales/Dropbox/00_papers/2023_hybrids_SIA_baleen/Final_plots/intraannual_variation.svg", tot1, 
#        dpi = 300,  width = 12, height = 8, units = "cm")
