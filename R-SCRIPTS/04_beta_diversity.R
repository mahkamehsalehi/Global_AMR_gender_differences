library(mia)
library(ggplot2)
library(dplyr)
# Creates Figure 4 
tse <- readRDS("data/TSE_gender_age.rds")
pcoa <- as.data.frame(reducedDim(tse, "PCoA")) 

sample_data <- tse %>% colData(tse) %>% as.data.frame()
# Rename columns to PC1, PC2
colnames(pcoa) <- paste0("PC", 1:ncol(pcoa))
# Add age group and region information to the PCoA data frame
pcoa$age_group <- sample_data$age_category
pcoa$geo <- sample_data$geo_loc_name_country_continent_calc
pcoa$gender <- sample_data$sex_combined

# Combine geo and age_group into a single factor for faceting
pcoa$facet_group <- paste(pcoa$geo, pcoa$age_group, sep = " - ")

# Filter out facets with fewer than 20 samples
pcoa_filtered <- pcoa %>%
  group_by(facet_group) %>%
  filter(n() >= 20) %>%
  ungroup()

# Create the PCoA plot with separate facets for each region-age group combination
ggplot(pcoa_filtered, aes(x = PC1, y = PC2, color = gender)) +
  geom_point(alpha = 0.7, size = 0.5) +            # Points for PCoA coordinates
  facet_wrap(~ facet_group, scales = "free") +    # Separate panels by region-age combination
  labs(x = "Principal Coordinate 1", y = "Principal Coordinate 2",
       color = "Gender") +                        # Axis and legend labels
  theme_minimal() +                               # Minimal theme for cleaner visuals
  theme(
    strip.text = element_text(size = 10, face = "bold"), # Font style for facet labels
    legend.position = "bottom",                   # Position the legend at the bottom
    panel.spacing = unit(1, "lines")              # Space between panels
  )



# Filter data for North America, Europe, and Africa
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

# Create the PCoA plot with facets for each region, colored by age group, and centroids
fig4 <- ggplot(pcoa_filtered, aes(x = PC1, y = PC2, color = age_group, shape = gender)) +
  geom_point(alpha = 0.4, size = 2) +               # Points for PCoA coordinates
  geom_point(data = centroids, aes(x = PC1, y = PC2, color = age_group, fill =age_group, shape = gender), 
            size = 6, ) +                # Colored triangles for centroids by age group
  facet_wrap(~ geo, scales = "free") +  
  scale_color_manual(values = c("lightblue", "orange","green4", "purple")) +
  # Custom colors for age groups
  labs(x = "Principal Coordinate 1", y = "Principal Coordinate 2",
       color = "Age Group", shape = "Gender") +                       # Axis and legend labels
  theme_minimal() +
  guides(fill = "none") +
  theme(strip.text.x.top  = element_text(size=15),
    legend.position = "bottom", 
    legend.text = element_text(size=15),
    legend.title = element_text(size=15),
    axis.title = element_text(size= 15),
    axis.text = element_text(size=15),
    panel.spacing = unit(1, "lines"), # Space between panels
    strip.text = element_text(size = 10, face = "bold") # Font style for facet labels
  )
png("RESULTS/PCoA_figure.png", width=800, height=400)
print(fig4)
dev.off()
