#--------------------------------------------------------------------------------
# Baleen SIA in fin and blue whale hybrids (Marc Ruiz-Sagalés, 6 de juny del 2025)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# 04 - Basic statistics and plots of the baleen isotopes per individual
#--------------------------------------------------------------------------------
setwd("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes")
source("03_Suess_Laws_effect_correction.R")

df

unique(df$Whale) # check individuals

# mean + SD (range) for d15n, d13c and d34s per Whale
df %>% dplyr::group_by(Whale) %>% dplyr::reframe(mean(d15n,na.rm=T),sd(d15n,na.rm=T),range(d15n,na.rm=T),
                                                              mean(d13c,na.rm=T),sd(d13c,na.rm=T),range(d13c,na.rm=T),
                                                              mean(d34s,na.rm=T),sd(d34s,na.rm=T),range(d34s,na.rm=T))

# mean + SD (range) for d15n, d13c and d34s per Class
df %>% dplyr::group_by(species) %>% dplyr::reframe(mean(d15n,na.rm=T),sd(d15n,na.rm=T),range(d15n,na.rm=T),
                                                 mean(d13c,na.rm=T),sd(d13c,na.rm=T),range(d13c,na.rm=T),
                                                 mean(d34s,na.rm=T),sd(d34s,na.rm=T),range(d34s,na.rm=T))

# model fitting, statistical extraction, and prediction for each isotope and each individual whale

isotope= c("d15n","d13c","d34s")
species_colors <- c("blue" = "#313695","hybrid" = "#F46D43","fin" = "#A50026")
species_order <- c("blue", "hybrid", "fin")
individuals <- unique(df$Whale)

stat_results <- list() # storage for all statistics
pred_df <- data.frame()  # storage for all predictions

for(iso in isotope) { # loops over each isotope 
  for (ind in individuals) { # loops over each ind. 
    
    df_ind <- df %>% dplyr::filter(Whale == ind, is.finite(.data[[iso]]), is.finite(Cm)) # filters the dataset for one individual whale and the current isotope, keeping only finite values
    
    # AIC optimization
    k_values <- 3:29
    aic_values <- sapply(k_values, function(k) { # fits a GAM for each candidate smoothness parameter k (between 3 and 29)
      model <- try(gam(reformulate("s(Cm, k = k, bs = 'cs')", response = iso), 
                       data = df_ind, family = gaussian(), gamma = 1), silent = TRUE)
      if (!inherits(model, "try-error")) AIC(model) else NA
    })
    
    if (all(is.na(aic_values))) next
    best_k <- k_values[which.min(aic_values)] # selects the value of k that minimizes AIC
    best_k <- min(best_k, floor(nrow(df_ind)/2))  # ensures k is not too large (limited to half the number of observations)
    
    # Fit final model
    final_model <- try(gam(reformulate("s(Cm, k = best_k, bs = 'cs')", response = iso), data = df_ind, family = gaussian()), silent = TRUE) # fits the final GAM model with the best k
    
    # --- Store stats ---
    sm <- summary(final_model)
    stat_results[[paste(ind, iso, sep = "_")]] <- data.frame( # extracts key statistics from the GAM (EDF, F, p-value, AIC, and chosen k)
      Individu = ind,
      Isotope = iso,
      EDF = sum(sm$s.table[, "edf"]),
      Ref_df = sum(sm$s.table[, "Ref.df"]),
      F = sum(sm$s.table[, "F"]),
      p_value = min(sm$s.table[, "p-value"]),
      AIC = AIC(final_model),
      k_used = best_k
    )
    
    
    # Predictions
    pred <- ggpredict(final_model, terms = "Cm") # compute predicted values along the Cm gradient
    pred$species <- unique(df_ind$species)
    pred$Whale <- ind  # keep track of the individual
    pred$Isotope <- iso  # keep track of the individual
    pred_df <- rbind(pred_df, pred)  # append to master dataframe
    
  }
}

pred_df<- as.data.frame(pred_df)

# plot of the baleen isotopes per ind.

# a) d15n
pred_df_dN <- pred_df %>% filter(Isotope == "d15n")
p <- ggplot(data=pred_df_dN,aes(x = x, y = predicted)) +
  geom_point(data = df, aes(x = Cm, y = d15n, fill = species),
             shape = 21, color = "black", alpha = 0.2) +
  geom_line(aes(color = species), size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = species), alpha = 0.3) +
  labs(x = "Distance from gingiva (cm)", y = expression(delta^15*N~"(\u2030)")) +
  coord_cartesian(xlim = c(-5, 35), ylim = c(7, 12.5)) +
  scale_fill_manual(values = species_colors) +
  scale_color_manual(values = species_colors) +
  theme_article(base_size = 10, base_family = "Optima") +
  theme(aspect.ratio = 3 / 4, legend.position = "none",strip.text = element_text(size = 10), axis.title.x = element_blank(),
        axis.text.x = element_blank()) + facet_wrap(vars(Whale), ncol=5)
p

# b) d13c
pred_df_dC <- pred_df %>% filter(Isotope == "d13c")
p1 <- ggplot(data=pred_df_dC,aes(x = x, y = predicted)) +
  geom_point(data = df, aes(x = Cm, y = d13c, fill = species),
             shape = 21, color = "black", alpha = 0.2) +
  geom_line(aes(color = species), size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = species), alpha = 0.3) +
  labs(x = "Distance from gingiva (cm)", y = expression(delta^13*C~"(\u2030)")) +
  coord_cartesian(xlim = c(-5, 35), ylim = c(-21, -18)) +
  scale_fill_manual(values = species_colors) +
  scale_color_manual(values = species_colors) +
  theme_article(base_size = 10, base_family = "Optima") +
  theme(aspect.ratio = 3 / 4, legend.position = "none",strip.text = element_text(size = 10),axis.title.x = element_blank(),
        axis.text.x = element_blank()) + facet_wrap(vars(Whale), ncol=5)
p1

# c) d34s
pred_df_dS <- pred_df %>% filter(Isotope == "d34s")
p2 <- ggplot(data=pred_df_dS,aes(x = x, y = predicted)) +
  geom_point(data = df, aes(x = Cm, y = d34s, fill = species),
             shape = 21, color = "black", alpha = 0.2) +
  geom_line(aes(color = species), size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = species), alpha = 0.3) +
  labs(x = "Distance from gingiva (cm)", y = expression(delta^34*S~"(\u2030)")) +
  coord_cartesian(xlim = c(-5, 35), ylim = c(16, 20)) +
  scale_fill_manual(values = species_colors) +
  scale_color_manual(values = species_colors) +
  theme_article(base_size = 10, base_family = "Optima") +
  theme(aspect.ratio = 3 / 4, legend.position = "none",strip.text = element_text(size = 10)) + facet_wrap(vars(Whale), ncol=5)
p2

plot_ind_N_C_S = p/p1/p2 + plot_annotation(tag_levels = 'A') 
plot_ind_N_C_S

#ggsave("/Users/marcruizisagales/Dropbox/00_papers/2023_hybrids_SIA_baleen/Final_plots/plot_ind_N_C_S.svg", plot = plot_ind_N_C_S, dpi = 300, width = 15, height = 20, units = "cm", device = "svg")
