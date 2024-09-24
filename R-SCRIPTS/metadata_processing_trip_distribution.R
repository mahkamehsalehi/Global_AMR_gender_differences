library(readxl)
Trip_distribution_matrix <-
  read_excel("DATA/Trip_distribution_matrix.xlsx")
View(Trip_distribution_matrix)

clean_trips <-
  Trip_distribution_matrix[c(1, 5:245), c(2, 5:245)] %>% as.data.frame
rownames(clean_trips) <- clean_trips[, 1] %>% str_to_title()
colnames(clean_trips) <- clean_trips[1, ] %>% str_to_title()
clean_trips <- (clean_trips[2:242, 2:242]) 

TSE <- readRDS("DATA/TSE.rds")

tse_countries <-
  colData(TSE) %>% data.frame %>% 
  dplyr::select(matches("geo_loc_name_country_calc")) %>% 
  unique()

setdiff(names(clean_trips), (tse_countries[,]))
setdiff(tse_countries[,], names(clean_trips))

# Recode the rownames outside of mutate
rownames(clean_trips) <- recode(rownames(clean_trips),
                                "Congo (Republic)" = "Democratic Republic of the Congo",
                                "Hongkong" = "Hong Kong",
                                "Turkey Turkiye" = "Turkey",
                                "United States" = "USA")

# Recode the rownames outside of mutate
colnames(clean_trips) <- recode(colnames(clean_trips),
                                "Congo (Republic)" = "Democratic Republic of the Congo",
                                "Hongkong" = "Hong Kong",
                                "Turkey Turkiye" = "Turkey",
                                "United States" = "USA")

intersect(tse_countries[,], names(clean_trips))


# Find the intersection of tse_countries with both row and column names
common_rows <- intersect(tse_countries[,], rownames(clean_trips))
common_cols <- intersect(tse_countries[,], names(clean_trips))

# Subset clean_trips based on the common rows and columns
subset_clean_trips <- clean_trips[common_rows, common_cols]

