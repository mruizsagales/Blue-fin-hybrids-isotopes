#--------------------------------------------------------------------------------
# Baleen SIA in fin and blue whale hybrids 
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# 08 - Dissimilarities in isotopic trajectory among individual whales (adapted from Sturbois et al., 2021;  https://doi.org/10.1002/ecm.1501)
#--------------------------------------------------------------------------------

library(pheatmap)
library(viridis)
library(gridExtra)

setwd("/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/02code")
source("03_Suess_Laws_effect_correction.R")
df

whale <- df[,c(3,34,5,7,19)] # select only necessary columns
whale1 <-whale[complete.cases(whale), ] # keep rows where all columns have non-missing values
whale1$species<- factor(whale1$species, levels = c('fin','hybrid','blue'))
whale1$id_cm <- rep(1:30, 10)

D <- dist(whale1[,c("d13c","d15n","d34s")])
sum(is.na(as.matrix(dist(whale1[, c("d13c", "d15n", "d34s")]))))

whales_x <- defineTrajectories(D, whale1$Whale)
Ds<-trajectoryDistances(whales_x, distance.type = "DSPD", symmetrization = "mean", add = TRUE)

pt<-c(16,16,16,16)
hsxy <- hclust(Ds, "ward.D2")
plot(hsxy,hang = -1, main="distance hybrids", cex=.6)

dist_mat <- as.matrix(Ds)
rownames(dist_mat) <- colnames(dist_mat) <- unique(whale1$Whale)

ann_df <- whale1 %>%
  distinct(Whale, species) %>%
  filter(Whale %in% rownames(dist_mat)) %>%
  arrange(match(Whale, rownames(dist_mat))) %>%
  column_to_rownames("Whale")
names(ann_df) <- "Class"
species_colors <- c("blue" = "#313695", "hybrid" = "#F46D43", "fin" = "#A50026")
ann_colors <- list(Class = species_colors)
alpha_colors <- adjustcolor(viridis::viridis(100), alpha.f = 0.9)

d<- pheatmap(dist_mat,
         clustering_distance_rows = Ds,
         clustering_distance_cols = Ds,
         cellwidth = 15,     # adjust width
         cellheight = 15,    # same as width → forces square cells
         clustering_method = "ward.D2",
         annotation_row = ann_df,
         annotation_col = ann_df,
         annotation_colors = ann_colors,
         annotation_legend = T,  # ← disables the species legend
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 8,
         color = alpha_colors
         #filename = "/Users/marcruizisagales/Documents/GitHub/Blue-fin-hybrids-isotopes/03figures/similarity_heatmap.pdf"
         ) 
d
