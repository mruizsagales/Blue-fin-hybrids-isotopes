#--------------------------------------------------------------------------------
# Baleen SIA in fin and blue whale hybrids (Marc Ruiz-Sagalés, 6 de juny del 2025)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# 10 - Temporal and spatial overlap between fin and blue whale catches
#--------------------------------------------------------------------------------

library(readxl)
library(tidyverse)
library(egg)
library(adehabitatHR)
library(sp)
library(dplyr)
library(ggplot2)
library(maps)
library(sf)
library(rphylopic)
library(ggimage)

data <- read_excel("/Users/marcruizisagales/Downloads/NA BPHY i BMUS only.xlsx") #Download data (from the International Whaling Comission dataset).
polygon_West_iceland<- st_read("/Users/marcruizisagales/Dropbox/00_papers/2022_blue_and_fin_whale_catches_NA/Blue_whale_logbooks_Na/Shapefiles/West_Iceland1.shp") # shapefile of the W Iceland feeding ground

data <- data %>% # recode
  dplyr::mutate(Sp = dplyr::recode(Sp, "10" = "B. brydei", 
                                   "11" =  "E. australis", 
                                   "12" = "E. robustus", 
                                   "15" =  "B.m. brevicauda",
                                   "20" = "B. bonaerensis", 
                                   "4" = "B. musculus", 
                                   "5" =  "B. physalus",
                                   "6" = "P. macrocephalus", 
                                   "7" = "M. novaeangliae", 
                                   "8" =  "B. borealis",
                                   "9" =  "B. acutorostrata"))


#Convert longitude and latitude from degrees, minutes and seconds to decimal units with ´lonsex2dec´ and ´latsex2dec´functions.

lonsex2dec <- function(degree, minute, direction) {
  declon = if_else(
    condition = (direction == "E"),
    true =  (degree + (minute / 60)),
    false = -(degree + (minute / 60))
  )
  return(declon)
}

longdecimal <- data %>% 
  mutate(
    longdecimal = lonsex2dec(
      data$Lon,data$Mn...20,data$...21
    )  
  )

latsex2dec <- function(degree, minute, direction) {
  declon = if_else(
    condition = (direction == "N"),
    true =  (degree + (minute / 60)),
    false = -(degree + (minute / 60))
  )
  return(declon)
}

latdecimal <- longdecimal %>% 
  mutate(
    latdecimal = latsex2dec(
      Lat,Mn...16,...17
    )  
  )

latdecimal
latdecimal_filtered <- latdecimal[complete.cases(latdecimal[, c("longdecimal", "latdecimal")]), ]

# check plot
world = map_data("world")
ggplot(latdecimal_filtered, aes(x= longdecimal, y= latdecimal, color=Sp, group=1)) + geom_point(size=0.25) + theme_article()+ theme(legend.position = "none", aspect.ratio = 1) + facet_wrap(vars(Sp)) + geom_map(data=world, map=world, aes(x=long, y=lat, map_id=region), color="white", fill="black", size=0.05, alpha=1/4) + xlim(-75,30) + ylim(25,85) + xlab("Longitude") + ylab("Latitude")

# a) Space 
points_sf_total <- st_as_sf(latdecimal_filtered, coords = c("longdecimal", "latdecimal"), crs = 4326)
points_within_polygon_West_iceland <- st_intersection(points_sf_total, polygon_West_iceland)
my_coords_West_iceland <- st_coordinates(points_within_polygon_West_iceland)
Total_West_iceland <- cbind(my_coords_West_iceland, points_within_polygon_West_iceland)
names(Total_West_iceland) <-  c("Longitude", "Latitude", "CBt", "Day", "Mon", "Year", "Sp", "Len", "L.u", "Sx", "NoF", "F1.L", "F1.S", "F2.L", "F2.S", "F.u", "Lat", "Mn...16", "...17", "Ac...18", "Lon", "Mn...20", "...21", "Ac...22", "Exp", "Sum.Ex", "Nt", "SCo", "S.q", "S.d", "S.s", "L.D", "Inf.", "Fem", "Mat", "By", "IT", "Mk", "Sa", "Txt", "Notes.Extracted.on..22.Dec.20...........14.43.16","id", "geometry")

# map

ggplot(Total_West_iceland, aes(x= Longitude, y= Latitude, color=Sp, group=1)) + 
  geom_point(size=0.25) +
  theme_article()+
  theme(legend.position = "none") + 
  facet_wrap(vars(Sp)) + 
  geom_map(data=world, map=world,
           aes(x=long, y=lat, map_id=region),
           color="white", fill="black", size=0.05, alpha=1/4) + 
  xlim(-30,-10) + 
  ylim(62,67.5) + xlab("Longitude") + ylab("Latitude") +  theme_article(base_size = 15, base_family = "Optima") + theme(aspect.ratio = 3/4, legend.position = "none")  + scale_fill_manual(values = species_colors) +
  scale_color_manual(values = species_colors) 

# create coordinates
Total_West_iceland_df <- as.data.frame(Total_West_iceland)
coords <- Total_West_iceland_df[, c("Longitude", "Latitude")]

# convert to SpatialPointsDataFrame
spdf <- sp::SpatialPointsDataFrame(
  coords = coords,
  data = as.data.frame(Total_West_iceland),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)

# Calculate kernel UD for each species
kud <- kernelUD(spdf[, "Sp"], h = "href", grid = 500)  # you can adjust bandwidth or grid
# Get 95% contour polygons
hr <- getverticeshr(kud, percent = 90)
# Convert to sf for ggplot
hr_sf <- st_as_sf(hr)
species_colors <- c("B. musculus" = "#313695", "hybrid" = "#F46D43", "B. physalus" = "#A50026")

kk <- ggplot() +
  geom_map(data = world, map = world,
           aes(x = long, y = lat, map_id = region),
           color = "white", fill = "black", size = 0.05, alpha = 0.25) +
  geom_sf(data = hr_sf, aes(fill = id, col=id), alpha = 0.1, linewidth=0.5) +
  # geom_point(data = Total_West_iceland,
  #            aes(x = Longitude, y = Latitude, color = Sp),
  #            size = 0.3) +
  scale_fill_manual(values = species_colors) +
  scale_color_manual(values = species_colors) +
  #facet_wrap(vars(Sp)) +
  xlim(-30, -10) +
  ylim(62, 68) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_article(base_size = 10, base_family = "Optima") +
  theme(
    aspect.ratio = 3 / 4,
    legend.title = element_blank(),
    legend.position = "none"
  )

kk

# convert sf to data frame
Total_West_iceland_df <- as.data.frame(Total_West_iceland)
Total_West_iceland_df$geometry <- NULL

# "B. musculus" "B. physalus" (and so on)
par(mfrow = c(1, 2))  # 1 row, 2 columns
plot(kud[["B. musculus"]], main = "B. musculus UD")
plot(kud[["B. physalus"]], main = "B. physalus UD")
par(mfrow = c(1, 1))  # reset
kerneloverlaphr(kud, method = c("HR"),
                percent = 90)

# calculate pairwise overlap
overlap <- kerneloverlaphr(kud, method = "BA", percent = 90)  # or "UDOI"
print(overlap)

##

df <- Total_West_iceland %>%
  dplyr::filter(Mon > 0)

# summary: mean and n
df_summary <- df %>%
  group_by(Sp) %>%
  summarise(
    mean_yday = mean(yday, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# Create custom y positions
df_summary <- df_summary %>%
  arrange(desc(n)) %>%  # optional: sort species to define top/bottom
  mutate(
    y_pos = seq(0.014, 0.013, length.out = n())  # Adjust as needed
  )

# UUIDs
uuid <- get_uuid(name = "Balaenoptera musculus", n = 2)
img <- get_phylopic(uuid = uuid[2])
uuid1 <- get_uuid(name = "Balaenoptera physalus", n = 2)
img1 <- get_phylopic(uuid = uuid1[2])
# Merge colors
df$yday <- yday(make_date(year = df$Year, month = df$Mon, day = df$Day))
range(df$Year)
# Plot: The area under curve sums to 1 (for each group).
species_colors <- c("B. musculus" = "#313695", "hybrid" = "#F46D43", "B. physalus" = "#A50026")
jj <- ggplot(data=df, aes(x = yday, col = Sp, fill = Sp)) +
  geom_density(aes(y = after_stat(density)), alpha = 0.1, adjust = 1.5, linewidth = 0.5) +
  #geom_vline(data = df_summary, aes(xintercept = mean_yday, col = Sp),linetype = "dotted", linewidth = 0.5) +
  geom_text(data = df_summary,
            aes(x = 360, y = y_pos, label = paste("n =", n), col = Sp),
            hjust = 1.1, vjust = 0, size = 3, show.legend = FALSE,
            inherit.aes = FALSE) +
  add_phylopic(x = 240, y = 0.0142, img = img1, alpha = 0.8, width = 40, fill = "#A50026") +
  add_phylopic(x = 260, y = 0.0132, img = img, alpha = 0.8, width = 50, fill = "#313695") +
  ylab("Density estimate") +
  ylim(0,0.015) +
  scale_x_continuous(name = "Yearday", limits = c(1, 366), breaks = seq(0, 360, 60)) +
  #labs(title = "Seasonal Distribution by Species – West Iceland", y = "Density (area = 1)") +
  theme_article(base_size = 10, base_family = "Optima") +
  theme(aspect.ratio = 3/4,
        # plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        # axis.title = element_text(size = 13),
        # axis.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.position = "none"
  ) + scale_fill_manual(values = species_colors) +
  scale_color_manual(values = species_colors) 

kk + jj + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 10))

# ggsave("/Users/marcruizisagales/Dropbox/00_papers/2023_hybrids_SIA_baleen/Final_plots/catches_denisty.svg", last_plot(), 
#        dpi = 300,  width = 17, height = 12, units = "cm")

# b) Time

monthly_catches <- df %>%
  dplyr::group_by(Sp, Mon) %>%
  dplyr::summarise(TotalCatch = n(), .groups = "drop")

# Example: get vectors for fin and sei whales
fin_catch <- monthly_catches %>%
  dplyr::filter(Sp == "B. physalus") %>%
  dplyr::arrange(Mon) %>%
  dplyr::pull(TotalCatch)

blue_catch <- monthly_catches %>%
  dplyr::filter(Sp == "B. musculus") %>%
  dplyr::arrange(Mon) %>%
  dplyr::pull(TotalCatch)

ccf_result <- ccf(fin_catch, blue_catch, plot = TRUE, main = "Cross-correlation: Fin vs Blue catches")

# Peak at lag = 0 ⇒ synchronous peaks. (posar mitjana també)

n_months <- length(fin_catch)
crit_val <- 1.96 / sqrt(n_months)
abline(h = crit_val, col = "red", lty = 2)
abline(h = -crit_val, col = "red", lty = 2)

# Years with both species
years_both <- df %>%
  dplyr::group_by(Year, Sp) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::count(Year) %>%       # count how many species per year
  dplyr::filter(n == 2) %>%    # keep only years with both
  dplyr::pull(Year)
df_filtered <- df %>%
  dplyr::filter(Year %in% years_both)
df_summary <- df_filtered %>%
  dplyr::group_by(Year, Sp) %>%
  dplyr::summarise(n = n(), mean_yday = mean(yday, na.rm = TRUE), .groups = "drop")


ggplot(df_filtered, aes(x = yday, col = Sp)) +
  geom_density(aes(y = after_stat(density)), adjust = 1, linewidth = 1) +
  geom_vline(data = df_summary,
             aes(xintercept = mean_yday, col = Sp),
             linetype = "dashed", linewidth = 1) +
  geom_text(data = df_summary,
            aes(x = 100, y = Inf, label = paste("n =", n), col = Sp),
            hjust = 1.1, vjust = 1.5, size = 4,
            inherit.aes = FALSE) +
  facet_wrap(~Year, scales = "free_y") +
  scale_x_continuous(name = "Day of year", limits = c(1, 365), breaks = seq(0, 360, 60)) +
  labs(title = "Density per year (only years with both species)", y = "Density") +
  theme_article() +
  theme(plot.title = element_text(hjust = 0.5))
