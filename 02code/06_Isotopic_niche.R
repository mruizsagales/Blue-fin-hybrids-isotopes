#--------------------------------------------------------------------------------
# Baleen SIA in fin and blue whale hybrids 
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# 06 - Isotopic niche
#--------------------------------------------------------------------------------

library(patchwork)

setwd("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/code2")
source("03_Suess_Laws_effect_correction.R")
df

whale <- df[,c(3,34,5,7,19)] # select only necessary columns
whale1 <-whale[complete.cases(whale), ] # keep rows where all columns have non-missing values
whale1$species<- factor(whale1$species, levels = c('fin','hybrid','blue'))

# Define species-specific colors
species_colors <- c("blue" = "#313695","hybrid" = "#F46D43","fin" = "#A50026")

# Compute group means and standard deviations
summary_df <- whale1 %>%
  group_by(species, Whale) %>%
  summarise(
    mean_d13c = mean(d13c, na.rm = TRUE),
    mean_d15n = mean(d15n, na.rm = TRUE),
    mean_d34s = mean(d34s, na.rm = TRUE),
    sd_d13c = sd(d13c, na.rm = TRUE),
    sd_d15n = sd(d15n, na.rm = TRUE),
    sd_d34s = sd(d34s, na.rm = TRUE),
    .groups = 'drop'
  )

# a) Isoplots d15n vs. d13c
aaa <- ggplot() +
  # # Convex hulls
  # geom_polygon(data = hull_df, aes(x = d13c, y = d15n, fill = species, group = species), alpha = 0.2, color = NA) +
  # 
  # Error bars
  geom_errorbar(data = summary_df, aes(x = mean_d13c, ymin = mean_d15n - sd_d15n, ymax = mean_d15n + sd_d15n, color = species), width = 0.1, alpha=0.75) +
  geom_errorbarh(data = summary_df, aes(y = mean_d15n, xmin = mean_d13c - sd_d13c, xmax = mean_d13c + sd_d13c, color = species), height = 0.1, alpha=0.75) +
  
  # Mean points
  geom_point(data = summary_df, aes(x = mean_d13c, y = mean_d15n, color = species), size = 1, alpha=0.95) +
  # Individual points with labels
  # geom_label_repel(
  #   data = summary_df,
  #   aes(x = mean_d13c, y = mean_d15n, label = Whale),
  #   size = 1.5,box.padding = 0.5, max.overlaps = Inf
  # ) +
  # # Labels
  # geom_text(data = summary_df, aes(x = mean_d13c, y = mean_d15n, label = species, color = species), hjust = -0.2, vjust = -0.2, size = 3) +
  
  # Aesthetics
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  #scale_shape_manual(values = location_shapes) +
  theme_article(base_size = 10, base_family = "Optima") +
  labs(
    x = expression(delta^13*C~"(\u2030)"),
    y = expression(delta^15*N~"(\u2030)"),
    color = "Species",
    fill = "Species"
  ) +
  theme(legend.position = "none", aspect.ratio = 1) + ylim(7.2,12.4) + xlim(-20.9,-18.1)

aaa

# b) Isoplots Isoplots d15n vs. d34s
bbb<-ggplot() +
  # # Convex hulls
  # geom_polygon(data = hull_df, aes(x = d13c, y = d15n, fill = species, group = species), alpha = 0.2, color = NA) +
  # 
  # Error bars
  geom_errorbar(data = summary_df, aes(x = mean_d34s, ymin = mean_d15n - sd_d15n, ymax = mean_d15n + sd_d15n, color = species), width = 0.1, alpha=0.75) +
  geom_errorbarh(data = summary_df, aes(y = mean_d15n, xmin = mean_d34s - sd_d34s, xmax = mean_d34s + sd_d34s, color = species), height = 0.1, alpha=0.75) +
  
  # Mean points
  geom_point(data = summary_df, aes(x = mean_d34s, y = mean_d15n, color = species), size = 1, alpha=0.95) +
  # Individual points with labels
  # geom_label_repel(
  #   data = summary_df,
  #   aes(x = mean_d34s, y = mean_d15n, label = Whale),
  #   size = 1.5,box.padding = 0.5, max.overlaps = Inf
  # ) +
  # # Labels
  # geom_text(data = summary_df, aes(x = mean_d13c, y = mean_d15n, label = species, color = species), hjust = -0.2, vjust = -0.2, size = 3) +
  
  # Aesthetics
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  #scale_shape_manual(values = location_shapes) +
  theme_article(base_size = 10, base_family = "Optima") +
  labs(
    x = expression(delta^34*S~"(\u2030)"),
    y = expression(delta^15*N~"(\u2030)"),
    color = "Species",
    fill = "Species"
  ) +
  theme(legend.position = "none", aspect.ratio = 1) + ylim(7.2,12.4) + xlim(16.1,19.9)

bbb

# c) Isoplots d13c vs. d34s
ccc<-ggplot() +
  # # Convex hulls
  # geom_polygon(data = hull_df, aes(x = d13c, y = d15n, fill = species, group = species), alpha = 0.2, color = NA) +
  # 
  # Error bars
  geom_errorbarh(data = summary_df, aes(y = mean_d13c, xmin = mean_d34s - sd_d34s, xmax = mean_d34s + sd_d34s, color = species), height = 0.1, alpha=0.75) +
  geom_errorbar(data = summary_df, aes(x = mean_d34s, ymin = mean_d13c - sd_d13c, ymax = mean_d13c + sd_d13c, color = species), width = 0.1, alpha=0.75) +
  
  # Mean points
  geom_point(data = summary_df, aes(y = mean_d13c, x = mean_d34s, color = species), size = 1, alpha=0.95) +
  # Individual points with labels
  # geom_label_repel(
  #   data = summary_df,
  #   aes(y = mean_d13c, x = mean_d34s, label = Whale),
  #   size = 1.5,box.padding = 0.5, max.overlaps = Inf
  # ) +
  # # Labels
  # geom_text(data = summary_df, aes(x = mean_d13c, y = mean_d15n, label = species, color = species), hjust = -0.2, vjust = -0.2, size = 3) +
  
  # Aesthetics
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  #scale_shape_manual(values = location_shapes) +
  theme_article(base_size = 10, base_family = "Optima") +
  labs(
    y = expression(delta^13*C~"(\u2030)"),
    x = expression(delta^34*S~"(\u2030)"),
    color = "Species",
    fill = "Species"
  ) +
  theme(legend.position = "none", aspect.ratio = 1) + xlim(16.1,19.9) + ylim(-20.9,-18.1)

ccc

tot <- aaa/bbb/ccc
tot

# posterior distribution of (mu, Sigma) for each species
nsamples <- 1000
whale_par <- whale1 %>% 
  split(.$species) %>% 
  map(~ dplyr::select(., d15n,d13c, d34s)) %>% 
  map(~ niw.post(nsample = nsamples, X = .))
whale_par_N_S <- whale1 %>% 
  split(.$species) %>% 
  map(~ dplyr::select(., d15n, d34s)) %>% 
  map(~ niw.post(nsample = nsamples, X = .))
whale_par_S_C <- whale1 %>% 
  split(.$species) %>% 
  map(~ dplyr::select(., d13c, d34s)) %>% 
  map(~ niw.post(nsample = nsamples, X = .))

# flatten posterior samples into a tidy tibble
posterior_df <- imap_dfr(whale_par, ~ {
  as_tibble(.x$mu) |> 
    set_names(c("d15n", "d13c", "d34s")) |> 
    mutate(species = .y)
})

posterior_df <- posterior_df %>% mutate(species = factor(species, levels = c("blue", "hybrid", "fin")))
posterior_df$species <- factor(posterior_df$species, levels = c("fin","hybrid","blue"))

pN1 <- ggplot(posterior_df, aes(y = d15n, x = fct_rev(species), fill = species)) +
  stat_eye(
    .width = c(0.66, 0.95),
    interval_color = "black",
    #adjust = 0.5,
    slab_alpha = 0.8,
    point_interval = mean_qi,
    #height=1.5,
    point_size =1,
    outline_bars = T,
    slab_color = "black",
    slab_linewidth = 0.5
  ) +
  theme_article(base_size = 10, base_family = "Optima") +
  labs(
    y = expression(paste("posterior - ", delta^15*N~"\u2030")),
    x = "",
    fill = "Species"
  ) + scale_fill_manual(values = species_colors) + ylim(7.2,12.4) + theme(aspect.ratio = 1, legend.position = "none",  
                                                                          panel.grid.major = element_line(color = "grey85", size = 0.1),
                                                                          panel.grid.minor = element_blank(),
                                                                          panel.background = element_blank(),
                                                                          axis.text.x = element_blank())

pC1 <- ggplot(posterior_df, aes(y = d13c, x = fct_rev(species), fill = species)) +
  stat_eye(
    .width = c(0.66, 0.95),
    interval_color = "black",
    #adjust = 0.5,
    slab_alpha = 0.8,
    point_interval = mean_qi,
    #height=1.5,
    point_size =1,
    outline_bars = T,
    slab_color = "black",
    slab_linewidth = 0.5
  ) +
  theme_article(base_size = 10, base_family = "Optima") +
  labs(
    y = expression(paste("posterior - ", delta^13*C~"\u2030")),
    x = "",
    fill = "Species"
  ) + scale_fill_manual(values = species_colors) + ylim(-20.9,-18.1) + theme(aspect.ratio = 1, legend.position = "none",  
                                                                             panel.grid.major = element_line(color = "grey85", size = 0.1),
                                                                             panel.grid.minor = element_blank(),
                                                                             panel.background = element_blank(),
                                                                             axis.text.x = element_blank())
pS1 <- ggplot(posterior_df, aes(y = d34s, x = fct_rev(species), fill = species)) +
  stat_eye(
    .width = c(0.66, 0.95),
    interval_color = "black",
    #adjust = 0.5,
    slab_alpha = 0.8,
    point_interval = mean_qi,
    #height=1.5,
    point_size =1,
    outline_bars = T,
    slab_color = "black",
    slab_linewidth = 0.5
  ) +
  theme_article(base_size = 10, base_family = "Optima") +
  scale_x_discrete(labels = c("fin" = "Fin", "hybrid" = "Hybrid", "blue" = "Blue")) + 
  labs(
    y = expression(paste("posterior - ", delta^34*S~"\u2030")),
    x = "Class",
    fill = "Species"
  ) + scale_fill_manual(values = species_colors) + ylim(16.1,19.9) + theme(aspect.ratio = 1, legend.position = "none",  
                                                                           panel.grid.major = element_line(color = "grey85", size = 0.1),
                                                                           panel.grid.minor = element_blank(),
                                                                           panel.background = element_blank())
pN1/pC1/pS1

# extract μ values
df_mu <- nichetools::extract_mu(whale_par)
df_mu_N_S <- extract_mu(whale_par_N_S)
df_mu_S_C <- extract_mu(whale_par_S_C)

df_mu <- df_mu %>%
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

# extract Σ values
df_sigma <- extract_sigma(whale_par,  isotope_n = 3)
df_sigma_N_S <- extract_sigma(whale_par_N_S)
df_sigma_S_C <- extract_sigma(whale_par_S_C)
df_sigma_cn <- extract_sigma(whale_par, data_format = "long",  isotope_n = 3) %>% filter(id != isotope)

# plot posterior distrubtion of μ and Σ
posterior_plots <- df_mu %>%
  split(.$isotope) %>%
  imap(
    ~ ggplot(data = ., aes(x = mu_est)) +
      geom_density(aes(fill = sample_name), alpha = 0.5) +
      scale_fill_viridis_d(begin = 0.25, end = 0.75,
                           option = "D", name = "Species") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.title.x =  element_markdown(),
            axis.title.y =  element_markdown(),
            legend.position = "none",
            legend.background = element_blank()
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

posterior_plots$d15n +
  theme(legend.position = c(0.18, 0.82)) + 
  posterior_plots$d13c + posterior_plots$d34s + theme(aspect.ratio = 1)

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
      isotope == "d15n" ~ "N",
      isotope == "d13c" ~ "C",
      isotope == "d34s" ~ "S",
    ),
    neutron_iso = case_when(
      isotope == "d15n" ~ 15,
      isotope == "d13c" ~ 13,
      isotope == "d34s" ~ 34,
    )
  )

sigma_plots <- df_sigma_cn %>%
  group_split(id, isotope) %>%
  imap(
    ~ ggplot(data = ., aes(x = post_sample)) +
      geom_density(aes(fill = sample_name), alpha = 0.5) +
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

(sigma_plots[[1]] + 
    theme(legend.position = c(0.1, 0.82)) + sigma_plots[[2]] + sigma_plots[[3]]) / (sigma_plots[[4]] + sigma_plots[[5]] + sigma_plots[[6]])

# estimate 40% niche ellipse

# d15n vs. d13c
ellipse_df <- niche_ellipse(dat_mu = df_mu, dat_sigma = df_sigma,
                            set_seed = 4,p_ell=0.4)

ellipse_plots <- ggplot() + 
  geom_polygon(data = ellipse_df,
               mapping = aes(x = d15n, y = d13c,
                             group = interaction(sample_number, sample_name),
                             color = sample_name),
               fill = NA,
               linewidth = 0.2) + 
  
  # scale_colour_viridis_d(begin = 0.25, end = 0.75, 
  #                        option = "D", name = "species",
  # ) + 
  # scale_x_continuous(breaks = rev(seq(-20, -40, -2))) +
  # scale_y_continuous(breaks = seq(6, 16, 2)) +
  theme_article(base_size = 10, base_family = "Optima") +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title = element_text(hjust = 0.5),
        legend.background = element_blank()) + 
  labs(x = expression(delta^13*C~"(\u2030)"), 
       y = expression(delta^15*N~"(\u2030)")) + scale_colour_manual(values=species_colors) + theme(aspect.ratio = 1) + ylim(7.2,12.4) + xlim(-20.9,-18.1)

ellipse_plots

# d15n vs. d34s
ellipse_df_N_S <- niche_ellipse(dat_mu = df_mu_N_S, dat_sigma = df_sigma_N_S,
                                set_seed = 4, p_ell = 0.4)

ellipse_plots1 <- ggplot() + 
  geom_polygon(data = ellipse_df_N_S,
               mapping = aes(x = d15n, y = d13c,
                             group = interaction(sample_number, sample_name),
                             color = sample_name),
               fill = NA,
               linewidth = 0.2) + 
  
  # scale_colour_viridis_d(begin = 0.25, end = 0.75, 
  #                        option = "D", name = "species",
  # ) + 
  # scale_x_continuous(breaks = rev(seq(-20, -40, -2))) +
  # scale_y_continuous(breaks = seq(6, 16, 2)) +
  theme_article(base_size = 10, base_family = "Optima") +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title = element_text(hjust = 0.5),
        legend.background = element_blank()) + 
  labs(x = expression(delta^34*S~"(\u2030)"), 
       y = expression(delta^15*N~"(\u2030)")) + scale_colour_manual(values=species_colors) + theme(aspect.ratio = 1)  + xlim(16.1,19.9) + ylim(7.2,12.4)

ellipse_plots1

# d34s vs. d13d
ellipse_df_S_C <- niche_ellipse(dat_mu = df_mu_S_C, dat_sigma = df_sigma_S_C,
                                set_seed = 4, p_ell=0.4)

ellipse_plots2 <- ggplot() + 
  geom_polygon(data = ellipse_df_S_C,
               mapping = aes(x = d15n, y = d13c,
                             group = interaction(sample_number, sample_name),
                             color = sample_name),
               fill = NA,
               linewidth = 0.2) + 
  
  scale_colour_viridis_d(begin = 0.25, end = 0.75, 
                         option = "D", name = "species",
  ) + 
  # scale_x_continuous(breaks = rev(seq(-20, -40, -2))) +
  # scale_y_continuous(breaks = seq(6, 16, 2)) +
  theme_article(base_size = 10, base_family = "Optima") +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title = element_text(hjust = 0.5),
        legend.background = element_blank()) + 
  labs(x = expression(delta^34*S~"(\u2030)"), 
       y = expression(delta^13*C~"(\u2030)")) + scale_colour_manual(values=species_colors) + theme(aspect.ratio = 1) + xlim(16.1,19.9) + ylim(-20.9,-18.1)

ellipse_plots2

plottogheter1 <- (
  (aaa + theme(plot.margin = margin(1, 1, 1, 0))) / 
    (bbb + theme(plot.margin = margin(1, 1, 1, 0))) / 
    (ccc + theme(plot.margin = margin(1, 1, 1, 0)))
) |
  (
    (pN1 + theme(plot.margin = margin(1, 0, 1, 0))) / 
      (pC1 + theme(plot.margin = margin(1, 0, 1, 0))) / 
      (pS1 + theme(plot.margin = margin(1, 0, 1, 0)))
  ) |
  (
    (ellipse_plots + theme(plot.margin = margin(1, 0, 1, 1))) / 
      (ellipse_plots1 + theme(plot.margin = margin(1, 0, 1, 1))) / 
      (ellipse_plots2 + theme(plot.margin = margin(1, 0, 1, 1))) 
  )

plottogheter1 <- plottogheter1 +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 10)) 

plottogheter1

#ggsave("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/03figures/plot_ellipses.svg", plot = plottogheter1, dpi = 300, width = 16, height = 15, units = "cm", device = "svg")

# estimate also 95% ellipses

# d15n vs. d13c
ellipse_df_95 <- niche_ellipse(dat_mu = df_mu, dat_sigma = df_sigma,
                            set_seed = 4,p_ell=0.95)

ellipse_plots_95 <- ggplot() + 
  geom_polygon(data = ellipse_df_95,
               mapping = aes(x = d15n, y = d13c,
                             group = interaction(sample_number, sample_name),
                             color = sample_name),
               fill = NA,
               linewidth = 0.2) + 
  
  # scale_colour_viridis_d(begin = 0.25, end = 0.75, 
  #                        option = "D", name = "species",
  # ) + 
  # scale_x_continuous(breaks = rev(seq(-20, -40, -2))) +
  # scale_y_continuous(breaks = seq(6, 16, 2)) +
  theme_article(base_size = 10, base_family = "Optima") +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title = element_text(hjust = 0.5),
        legend.background = element_blank()) + 
  labs(x = expression(delta^13*C~"(\u2030)"), 
       y = expression(delta^15*N~"(\u2030)")) + scale_colour_manual(values=species_colors) + theme(aspect.ratio = 1) #+ ylim(7.2,12.4) + xlim(-20.9,-18.1)

ellipse_plots_95

# d15n vs. d34s
ellipse_df_N_S_95 <- niche_ellipse(dat_mu = df_mu_N_S, dat_sigma = df_sigma_N_S,
                                set_seed = 4, p_ell = 0.95)

ellipse_plots1_95 <- ggplot() + 
  geom_polygon(data = ellipse_df_N_S_95,
               mapping = aes(x = d15n, y = d13c,
                             group = interaction(sample_number, sample_name),
                             color = sample_name),
               fill = NA,
               linewidth = 0.2) + 
  
  # scale_colour_viridis_d(begin = 0.25, end = 0.75, 
  #                        option = "D", name = "species",
  # ) + 
  # scale_x_continuous(breaks = rev(seq(-20, -40, -2))) +
  # scale_y_continuous(breaks = seq(6, 16, 2)) +
  theme_article(base_size = 10, base_family = "Optima") +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title = element_text(hjust = 0.5),
        legend.background = element_blank()) + 
  labs(x = expression(delta^34*S~"(\u2030)"), 
       y = expression(delta^15*N~"(\u2030)")) + scale_colour_manual(values=species_colors) + theme(aspect.ratio = 1)  #+ xlim(16.1,19.9) + ylim(7.2,12.4)

ellipse_plots1_95

# d34s vs. d13d
ellipse_df_S_C_95 <- niche_ellipse(dat_mu = df_mu_S_C, dat_sigma = df_sigma_S_C,
                                set_seed = 4, p_ell=0.95)

ellipse_plots2_95 <- ggplot() + 
  geom_polygon(data = ellipse_df_S_C_95,
               mapping = aes(x = d15n, y = d13c,
                             group = interaction(sample_number, sample_name),
                             color = sample_name),
               fill = NA,
               linewidth = 0.2) + 
  
  scale_colour_viridis_d(begin = 0.25, end = 0.75, 
                         option = "D", name = "species",
  ) + 
  # scale_x_continuous(breaks = rev(seq(-20, -40, -2))) +
  # scale_y_continuous(breaks = seq(6, 16, 2)) +
  theme_article(base_size = 10, base_family = "Optima") +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title = element_text(hjust = 0.5),
        legend.background = element_blank()) + 
  labs(x = expression(delta^34*S~"(\u2030)"), 
       y = expression(delta^13*C~"(\u2030)")) + scale_colour_manual(values=species_colors) + theme(aspect.ratio = 1) #+ xlim(16.1,19.9) + ylim(-20.9,-18.1)

ellipse_plots2_95

ellipse_plots_95 + ellipse_plots1_95 + ellipse_plots2_95

#ggsave("/Users/marcruizisagales/Dropbox/00_papers/2023_hybrids_SIA_baleen/Final_plots/plot_ellipses_95.svg", plot = last_plot(), dpi = 300, width = 16, height = 15, units = "cm", device = "svg")

# Determine the 95% niche similarties for each species
over_stat <- nicheROVER::overlap(whale_par, nreps = nsamples, nprob = 1000, alpha = 0.95)
over_stat_df <- extract_overlap(data = over_stat)
over_sum <- over_stat_df %>% 
  group_by(sample_name_a, sample_name_b) %>% 
  summarise(
    median_niche_overlap = round(median(niche_overlap_perc), digits = 2),
    mean_niche_overlap = round(mean(niche_overlap_perc, na.rm=T), digits = 2),
    qual_2.5 = round(quantile(niche_overlap_perc, 
                              probs = 0.025, na.rm = TRUE), digits = 2), 
    qual_97.5 = round(quantile(niche_overlap_perc, 
                               probs = 0.975, na.rm = TRUE), digits = 2)
  ) %>% 
  ungroup() %>% 
  pivot_longer(cols = -c(sample_name_a, sample_name_b, mean_niche_overlap,median_niche_overlap), 
               names_to = "percentage", 
               values_to = "niche_overlap_qual") %>% 
  mutate(
    percentage = as.numeric(str_remove(percentage, "qual_"))
  ) 

# Clean and pivot to extract summary per pair
niche_overlap_summary <- over_sum %>%
  filter(!is.na(mean_niche_overlap)) %>%
  group_by(sample_name_a, sample_name_b) %>%
  summarise(
    mean_overlap = unique(mean_niche_overlap),
    lower_95 = min(niche_overlap_qual),
    upper_95 = max(niche_overlap_qual),
    .groups = "drop"
  ) %>%
  mutate(
    direction = paste(sample_name_a, "→", sample_name_b),
    label = paste0(round(mean_overlap, 1), "% [", 
                   round(lower_95, 1), "–", round(upper_95, 1), "%]")
  )


ov <- ggplot(data = over_stat_df, aes(x = sample_name_a, 
                                      y = niche_overlap_perc, 
                                      fill = sample_name_b)) + stat_eye(position = position_dodge(width = 0.9),width = c(0.95), alpha=0.75,
                                                                        interval_color = "black",
                                                                        #adjust = 0.5,
                                                                        slab_alpha = 0.75,
                                                                        point_interval = mean_qi,
                                                                        #height=1.5,
                                                                        point_size =1,
                                                                        outline_bars = T,
                                                                        slab_color = "black",
                                                                        slab_linewidth = 0.5) +
  geom_vline(xintercept = seq(1.5, 3, 1), 
             linetype = 1, color= "lightgrey") + theme_article(base_size = 10, base_family = "Optima") +
  theme(
    panel.grid = element_blank(), 
    axis.text = element_text(colour = "black"), 
    legend.background = element_blank(),
    strip.background = element_blank(),
    legend.position = "none",
    aspect.ratio = 3/4
  ) +
  labs(x = paste("Overlap Probability (%)", "\u2013", 
                 "Niche Region Size: 95%"), 
       y = "p(Percent Overlap | X)") +
  scale_fill_manual(values = species_colors,
                    labels = c("Blue", "Hybrid", "Fin")) + 
  scale_x_discrete(labels = c("blue" = "Blue", "hybrid" = "Hybrid", "fin" = "Fin")) + ylim(0,100)

ov

# Estimate overall niche size (95%)
niche_size <- extract_niche_size(whale_par, prob=0.95)
niche_size$species <- factor(niche_size$sample_name, levels = c("blue", "hybrid", "fin"))
niche_summary <- niche_size %>%
  group_by(sample_name) %>%
  summarise(
    mean_size = mean(niche_size),
    median_size = median(niche_size),
    sd_size = sd(niche_size),
    lower_95 = quantile(niche_size, 0.025),
    upper_95 = quantile(niche_size, 0.975)
  )

ns <- ggplot(data = niche_size, 
             aes(x = niche_size, y = species)) + 
  tidybayes::stat_halfeye(aes(fill=species),width = 0.5, alpha=0.75,
                          interval_color = "black",
                          #adjust = 0.5,
                          slab_alpha = 0.75,
                          point_interval = mean_qi,
                          #height=1.5,
                          point_size =1,
                          outline_bars = T,
                          slab_color = "black",
                          slab_linewidth = 0.5) + 
  stat_summary(fun.y = median, geom = "point", 
               size = 3, 
               position = position_dodge(width = 0.9)) + 
  theme_article(base_size = 10, base_family = "Optima") + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black"), legend.title = element_blank(), aspect.ratio = 1, legend.position = "none") + 
  labs(x = "Niche Size", 
       y = "Class") + scale_fill_manual(values= species_colors) + 
  scale_fill_manual(values = species_colors,
                    labels = c("Blue", "Hybrid", "Fin")) + 
  scale_y_discrete(labels = c("blue" = "Blue", "hybrid" = "Hybrid", "fin" = "Fin"))

niche_size_and_overlap_plot <- ov + ns + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 10))

niche_size_and_overlap_plot

#ggsave("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/03figures/plot_niche_size_overlap.svg", plot = niche_size_and_overlap_plot, dpi = 300, width = 16, height = 15, units = "cm", device = "svg")


