#--------------------------------------------------------------------------------
# Baleen SIA in fin and blue whale hybrids (Marc Ruiz-Sagalés, 6 de juny del 2025)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# 09 - Violation of independence check (adapted from Buss et al., 2022; https://doi.org/10.1007/s00227-022-04131-x)
#--------------------------------------------------------------------------------

library(tidyverse)
library(lmtest)
library(tidybayes)

setwd("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes")
source("03_Suess_Laws_effect_correction.R")
df

iso_vars <- c("d15n", "d13c", "d34s")
species_col <- "species"
species_colors <- c("blue" = "#313695", "hybrid" = "#F46D43", "fin" = "#A50026")

whale <- df[,c(2,3,34,5,7,19)]
whale1 <-whale[complete.cases(whale), ] # keep rows where all columns have non-missing values

# Normality and homogeneity tests by species and isotope
normality_results <- whale1 %>%
  pivot_longer(cols = all_of(iso_vars), names_to = "Isotope", values_to = "Value") %>%
  dplyr::group_by(species, Isotope) %>%
  dplyr::summarise(Shapiro_W_p = shapiro.test(Value)$p.value, .groups = "drop") %>% 
  dplyr::mutate(
    Significance = case_when(
      Shapiro_W_p > 0.05 ~ "ns",
      Shapiro_W_p <= 0.05 & Shapiro_W_p > 0.01 ~ "*",
      Shapiro_W_p <= 0.01 & Shapiro_W_p > 0.001 ~ "**",
      Shapiro_W_p <= 0.001 ~ "***"
    )
  )

normality_results # les significatives no segueixen distribució normal

homogeneity_results <- whale1 %>%
  pivot_longer(cols = all_of(iso_vars), names_to = "Isotope", values_to = "Value") %>%
  dplyr::group_by(Isotope) %>%
  dplyr::summarise(LeveneTest_p = car::leveneTest(Value ~ as.factor(species))$"Pr(>F)"[1],.groups = "drop") %>% 
  dplyr::mutate(
    Significance = case_when(
      LeveneTest_p > 0.05 ~ "ns",
      LeveneTest_p <= 0.05 & LeveneTest_p > 0.01 ~ "*",
      LeveneTest_p <= 0.01 & LeveneTest_p > 0.001 ~ "**",
      LeveneTest_p <= 0.001 ~ "***"
    )
  )

homogeneity_results # unequal variances for d13c and d34s

# Variables
time_var <- "Cm"
id_var <- "Whale"

# Function to subsample 10 data points per individual
subsample_individuals <- function(df, n_per_id = 10) {
  df %>%
    dplyr::group_by(across(all_of(id_var))) %>%
    dplyr::slice_sample(n = n_per_id) %>%
    dplyr::ungroup()
}


# Function to check autocorrelation per isotope per individual
check_autocorrelation <- function(df) {
  df %>%
    dplyr::group_by(across(all_of(id_var))) %>%
    dplyr::summarise(
      ac_pass = all(sapply(iso_vars, function(iso) {
        model <- lm(reformulate(time_var, iso), data = cur_data())
        test  <- lmtest::bgtest(model, order = 1)
        test$p.value > 0.05
      })),
      .groups = "drop"
    ) %>%
    pull(ac_pass) %>%
    all()
}

# Subsample 200 datasets
set.seed(123)
all_subsets <- replicate(200, subsample_individuals(whale1, n_per_id = 10), simplify = F)

# Filter clean ones (no autocorrelation)
clean_subsets <- all_subsets[sapply(all_subsets, check_autocorrelation)]

# Results
length(all_subsets)     # 200 total
length(clean_subsets)   # En surten 46 de no autocorrelacionats

# Posterior overlap for clean subsets
get_overlap_matrix <- function(df) {
  par_list <- tapply(1:nrow(df), df$species, function(ii) {
    niw.post(nsamples = 1000, X = df[ii, iso_vars])
  })
  overlap(par_list, nreps = 1000, nprob = 1000, alpha = 0.95)
}

subset_overlaps <- lapply(clean_subsets, get_overlap_matrix)

# Full dataset posterior overlap
full_par <- tapply(1:nrow(whale1), whale1$species, function(ii) {
  niw.post(nsamples = 1000, X = whale1[ii, iso_vars])
})
full_overlap <- overlap(full_par, nreps = 1000, nprob = 1000, alpha = 0.95)

# Calculate 95% CI of posterior distributions from subset overlaps
subset_ci_list <- lapply(seq_along(subset_overlaps), function(i) {
  mat <- subset_overlaps[[i]]
  if (!is.null(mat) && is.array(mat)) {
    tryCatch({
      ci_array <- apply(mat * 100, c(1, 2), quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      data.frame(
        subset = i,
        expand.grid(Species.A = dimnames(ci_array)[[1]], Species.B = dimnames(ci_array)[[2]]),
        CI_low = as.vector(ci_array[1,,]),
        CI_high = as.vector(ci_array[2,,])
      )
    }, error = function(e) NULL)
  } else {
    NULL
  }
})
subset_ci_df <- do.call(rbind, subset_ci_list)

# Compute % difference between subset and full overlap (mean values)
get_overlap_diff <- function(overlap_mat, full_mat) {
  diffs <- apply(overlap_mat, c(1, 2), mean) - apply(full_mat, c(1, 2), mean)
  as.data.frame(as.table(diffs))
}

overlap_diffs <- lapply(subset_overlaps, get_overlap_diff, full_mat = full_overlap)
diff_df <- do.call(rbind, Map(cbind, subset = seq_along(overlap_diffs), overlap_diffs))
names(diff_df) <- c("Subset", "Species.A", "Species.B", "DeltaOverlap") # Clean up names
p1 <- ggplot() + 
  geom_density(data=diff_df, aes(x = DeltaOverlap)) +
  geom_point(shape = 21, fill = "black", alpha = 0.5, size = 2) +
  facet_wrap(~ paste0(Species.A, ":", Species.B), ncol = 3, scales="free_y") +
  geom_vline(xintercept = 0, linetype = "dashed", color="grey") +
  theme_minimal() +
  labs(x = "Increase of niche overlap (%) in subsampled datasets",
       y = "Subset") + theme_article(base_size = 10, base_family = "Helvetica") + theme(aspect.ratio = 3/4, axis.text.y = element_blank())

p1

# ggsave("/Users/marcruizisagales/Desktop/hybrids/overlap_subset.pdf", p1, 
#        dpi = 300,  width = 12, height = 8, units = "cm")

# Posterior niche sizes for each subset
get_niche_sizes <- function(df) {
  par_list <- tapply(1:nrow(df), df$species, function(ii) {
    niw.post(nsamples = 1000, X = df[ii, iso_vars])
  })
  sapply(par_list, function(par) apply(par$Sigma, 3, niche.size, alpha = 0.95))
}

niche_sizes_list <- lapply(clean_subsets, get_niche_sizes)
niche_size_df <- do.call(rbind, lapply(seq_along(niche_sizes_list), function(i) {
  df <- as.data.frame(niche_sizes_list[[i]])
  df$subset <- i
  df
}))
niche_size_long <- pivot_longer(niche_size_df, -subset, names_to = "Species", values_to = "SEA")

# Full dataset niche size
full_size <- sapply(full_par, function(par) apply(par$Sigma, 3, niche.size, alpha = 0.95))
colnames(full_size) <- names(full_par)  # Add column names from species
full_summary <- data.frame(Species = colnames(full_size), MeanSEA = colMeans(full_size))

# Plot Figure B style
p2 <- ggplot(niche_size_long, aes(x = SEA, y = as.character(subset), color = Species)) +
  tidybayes::stat_pointinterval(point_interval = mean_qi, alpha = 1, size=0.5) +
  geom_vline(data = full_summary, aes(xintercept = MeanSEA), linetype = "dashed", color="black") +
  facet_wrap(vars(Species)) +
  theme_article(base_size = 10, base_family = "Helvetica") +
  labs(x = "Standard Ellipse Area (‰²)", y = "Dataset") + theme(axis.text.y = element_blank(), aspect.ratio = 4/1, legend.position = "none") + ggplot2::scale_color_manual(values = species_colors)
p2

# ggsave("/Users/marcruizisagales/Desktop/hybrids/SEA_subset.pdf", p2, 
#        dpi = 300,  width = 10, height = 10, units = "cm")



clean_subsets
subset_mu_df <- purrr::map_dfr(seq_along(clean_subsets), function(i) {
  subset_data <- clean_subsets[[i]]
  par_list <- split(1:nrow(subset_data), subset_data$species) %>%
    purrr::map(~ nicheROVER::niw.post(nsamples = 1000, X = subset_data[.x, iso_vars]))
  purrr::map_dfr(names(par_list), function(sp) {
    as.data.frame(par_list[[sp]]$mu) %>%
      setNames(iso_vars) %>%
      dplyr::mutate(Species = sp, Subset = i)
  })
})

mu_long_subsets <- tidyr::pivot_longer(subset_mu_df, cols = all_of(iso_vars), names_to = "Isotope", values_to = "Value")

full_mu <- split(1:nrow(fish), fish$species) %>%
  purrr::map(~ nicheROVER::niw.post(nsamples = 1000, X = fish[.x, iso_vars]))

# Posterior means per species and isotope
mu_summary <- purrr::map_dfr(names(full_mu), function(sp) {
  mu_mat <- as.data.frame(full_mu[[sp]]$mu)  # Columns = isotopes?
  colnames(mu_mat) <- iso_vars  # Ensure correct names
  
  mu_mat %>%
    dplyr::summarise(across(everything(), mean)) %>%
    tidyr::pivot_longer(cols = everything(), names_to = "Isotope", values_to = "mu_mean") %>%
    dplyr::mutate(Species = sp)
})



# 1. δ15N plot
pN <- ggplot2::ggplot(
  mu_long_subsets %>% dplyr::filter(Isotope == "d15n"),
  ggplot2::aes(x = Value, y = as.character(Subset), color = Species)
) +
  tidybayes::stat_pointinterval(point_interval = mean_qi, alpha = 0.6, size=0.5) +
  ggplot2::geom_vline(
    data = mu_summary %>% dplyr::filter(Isotope == "d15n"),
    ggplot2::aes(xintercept = mu_mean, linetype = Species),
    color = "black"
  ) +
  theme_article(base_size = 10) +
  ggplot2::labs(x = expression(paste(delta^15, "N (", "\u2030", ")")), y = "Dataset") + scale_x_continuous(breaks = seq(7.5, 11.5, by = 1)) + ggplot2::scale_color_manual(values = species_colors) + theme(legend.position = "none", axis.text.y = element_blank(), aspect.ratio = 4/1) 

# 2. δ13C plot
pC <- ggplot2::ggplot(
  mu_long_subsets %>% dplyr::filter(Isotope == "d13c"),
  ggplot2::aes(x = Value, y = as.character(Subset), color = Species)
) +
  tidybayes::stat_pointinterval(point_interval = mean_qi, alpha = 0.6, size=0.5) +
  ggplot2::geom_vline(
    data = mu_summary %>% dplyr::filter(Isotope == "d13c"),
    ggplot2::aes(xintercept = mu_mean, linetype = Species),
    color = "black"
  ) +
  theme_article(base_size = 10) +
  ggplot2::labs(x = expression(paste(delta^13, "C (", "\u2030", ")")), y = "") + scale_x_continuous(breaks = seq(-20.5, -18.5, by = 0.25)) + ggplot2::scale_color_manual(values = species_colors) + theme(legend.position = "none", axis.text.y = element_blank()) 

# 3. δ34S plot
pS <- ggplot2::ggplot(
  mu_long_subsets %>% dplyr::filter(Isotope == "d34s"),
  ggplot2::aes(x = Value, y = as.character(Subset), color = Species)
) +
  tidybayes::stat_pointinterval(point_interval = mean_qi, alpha = 0.6, size=0.5) +
  ggplot2::geom_vline(
    data = mu_summary %>% dplyr::filter(Isotope == "d34s"),
    ggplot2::aes(xintercept = mu_mean, linetype = Species),
    color = "black"
  ) +
  theme_article(base_size = 10) +
  ggplot2::labs(x = expression(paste(delta^34, "S (", "\u2030", ")")),y = "") + 
  ggplot2::guides(color = ggplot2::guide_legend(title = "Subset Posterior"),
                  linetype = ggplot2::guide_legend(title = "Full Dataset Mean",legend.direction = "horizontal")) + 
  ggplot2::scale_color_manual(values = species_colors) + 
  theme(axis.text.y = element_blank()) + 
  scale_x_continuous(breaks = seq(17, 19, by = 0.5))


p3 <- pN + theme(aspect.ratio = 4/1) +pC+ theme(aspect.ratio = 4/1) + pS+ theme(aspect.ratio = 4/1)
p3

# ggsave("/Users/marcruizisagales/Desktop/hybrids/Posterior_mean_subset.pdf", p3, 
#        dpi = 300,  width = 30, height = 10, units = "cm")