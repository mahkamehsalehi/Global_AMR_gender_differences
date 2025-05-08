# Load required libraries
library(mia)
library(ggplot2)
library(dplyr)
library(Cairo)

# ---------------------------
# 1) Load TSE object and extract PCoA coordinates
# ---------------------------

tse <- readRDS("../DATA/TSE_gender_age.rds")
pcoa <- as.data.frame(reducedDim(tse, "PCoA")) 
sample_data <- tse %>% colData(tse) %>% as.data.frame()

# Rename columns to PC1, PC2
colnames(pcoa) <- paste0("PC", 1:ncol(pcoa))

# ---------------------------
# 2) Annotate PCoA data with metadata
# ---------------------------

pcoa$age_group <- sample_data$age_category
pcoa$geo <- sample_data$geo_loc_name_country_continent_calc
pcoa$gender <- sample_data$sex_combined

# Combine geo and age_group into a single factor for faceting
pcoa$facet_group <- paste(pcoa$geo, pcoa$age_group, sep = " - ")

# ---------------------------
# 3) Plot PCoA faceted by region-age group combinations
# ---------------------------

pcoa_filtered <- pcoa %>%
  group_by(facet_group) %>%
  filter(n() >= 20) %>%         # Filter out facets with fewer than 20 samples
  ungroup()

ggplot(pcoa_filtered, aes(x = PC1, y = PC2, color = gender)) +
  geom_point(alpha = 0.7, size = 0.5) +
  facet_wrap(~ facet_group, scales = "free") +
  labs(x = "Principal Coordinate 1", y = "Principal Coordinate 2",
       color = "Gender") + 
  theme_minimal() + 
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(1, "lines")
  )

# ---------------------------
# 4) Create Figure 4: PCoA with centroids for selected regions
# ---------------------------

# Subset data for North America, Europe, and Africa
pcoa_filtered <- pcoa %>%
  filter(geo %in% c("North America", "Europe", "Africa"))

# Calculate centroids for each age group within each region
centroids <- pcoa_filtered %>%
  group_by(geo, age_group, gender) %>%
  summarize(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = 'drop'
  )

# Define text size for plot elements
s <- 20

# Construct the final Figure 4
fig4 <- ggplot(pcoa_filtered, aes(x = PC1, y = PC2, color = age_group, shape = gender)) +
  geom_point(alpha = 0.4, size = 2) + 
  geom_point(data = centroids, aes(x = PC1, y = PC2, color = age_group, fill =age_group, shape = gender), 
             size = 6, ) + 
  facet_wrap(~ geo, scales = "free") +  
  scale_color_manual(values = c("lightblue", "orange","green4", "purple")) +
  labs(x = "Principal Coordinate 1", y = "Principal Coordinate 2",
       color = "Age Group", shape = "Gender") + 
  theme_minimal() +
  guides(fill = "none") +
  theme(strip.text.x.top  = element_text(size=s),
        legend.position = "bottom", 
        legend.text = element_text(size=s),
        legend.title = element_text(size=s),
        axis.title = element_text(size= s),
        axis.text = element_text(size=s),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 10, face = "bold")
  )

# ---------------------------
# 5) Save Figure 4
# ---------------------------
CairoJPEG("../RESULTS/Fig2_PCoA.png", width = 800, height = 420, quality = 100)
print(fig4)
dev.off()
