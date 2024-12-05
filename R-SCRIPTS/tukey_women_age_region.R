library(tidyverse)
library(vegan)
library(scater)
library(ggpubr)
library(rstatix)
library(viridis)
library(ggsignif)
library(Matrix)
library(RColorBrewer)
library(caret)
library(randomForest)
library(xgboost)
library(brms)
library(glmnet)
library(jtools)
library(multcomp)

## Helpers ************************************* ####

# Set number <0.001 to "<0.001"
to_neat <- function(x, threshold = 0.001) {
  ifelse(x > threshold, as.character(x), "<0.001")
}

# Get number of stars for P.adj significant column
get_stars <- function(x) {
  
  ifelse(x<0.0001, "****", 
         ifelse(x<0.001, "***", 
                ifelse(x<0.01, "**", 
                       ifelse(x<0.05, "*", "ns"))))
  
}

## Data **************************************** ####

TSE <- readRDS("DATA/TSE_filtered.rds")

df <- TSE %>% colData(TSE) %>% as.data.frame()

tmp <- df %>% dplyr::select(matches(
  c(
    "acc",
    "ARG_load",
    "sex_combined",
    "age_category_new",
    "World_Bank_Income_Group",
    "Usage",
    "geo_loc_name_country_continent_calc"
  )
)) %>% na.omit()

# Subset
tmp$HIC <- tmp$World_Bank_Income_Group == "High income"
tmp$High_use <- tmp$Usage > 10
tmp$geo_loc_name_country_continent_calc <-
  tmp$geo_loc_name_country_continent_calc %>% as.factor

tmp <- full_join(tmp, data.frame(adult_metadata)[, c("acc", "shannon_diversity")], by = "acc")



## ARG LOAD ***************************************************** ####
# HICs ******************* ####

data_HIC <- tmp[tmp$HIC == "TRUE" & tmp$sex_combined == "female", ]

model_HIC <- glm(
  log(ARG_load) ~ High_use + age_category_new + geo_loc_name_country_continent_calc,
  data = data_HIC
)

# Post hoc test for age_category_new
age_tukey_female_HIC <- glht(model_HIC,
                             linfct = mcp(age_category_new = "Tukey"))
# Post hoc test for region
region_tukey_female_HIC <- glht(model_HIC,
                                linfct = mcp(geo_loc_name_country_continent_calc = "Tukey"))



# LMICs ****************** ####

data_LMIC <- tmp[tmp$HIC == "FALSE" & tmp$sex_combined == "female", ]

model_LMIC <- glm(
  log(ARG_load) ~ High_use + age_category_new +
    geo_loc_name_country_continent_calc,
  data = data_LMIC
)

# Post hoc test for age_category_new
age_tukey_female_LMIC <-
  glht(model_LMIC, linfct = mcp(age_category_new = "Tukey"))

# Post hoc test for region
region_tukey_female_LMIC <-
  glht(model_LMIC, linfct = mcp(geo_loc_name_country_continent_calc = "Tukey"))

## Compute N's
GPD_age_Ns <- tmp %>%
  filter(sex_combined == "female") %>%
  dplyr::select(HIC, age_category_new) %>% 
  table %>%
  data.frame()

GPD_region_Ns <- tmp %>%
  filter(sex_combined == "female") %>%
  dplyr::select(HIC, geo_loc_name_country_continent_calc) %>% 
  table %>%
  data.frame()

## Combine LM/HIC summaries ******************* ####

age_tukey_GPD <- list(HIC = age_tukey_female_HIC, LMIC = age_tukey_female_LMIC)

age_tukey_GDP_table <- lapply(names(age_tukey_GPD), function(f) {

  X <- age_tukey_GPD[[f]]

  tukey_summary <- summary(X)

  # Confidence intervals
  # conf_intervals <- confint(X)

  # Create a tidy data frame
  results <- data.frame(
    Comparison = names(tukey_summary$test$coefficients),
    Estimate = tukey_summary$test$coefficients,
    # SE = tukey_summary$test$sigma,
    # z_value = tukey_summary$test$tstat,
    p_value = tukey_summary$test$pvalues,
    # lower2.5 = as.data.frame(conf_intervals$confint)$lwr,
    # upper97.5 = as.data.frame(conf_intervals$confint)$upr,
    `Income Group` = f
  )

  results$`Effect Direction` <- ifelse(results$Estimate > 0,
                                       gsub("-", ">",results$Comparison), 
                                       gsub("-", "<",results$Comparison))
  
  
  results <- results %>% dplyr::select(-Estimate)
  
  results <- results %>% 
    separate(Comparison, into = c("Group 1", "Group 2"), sep = " - ")
  
  # Add N1
  results <- full_join(GPD_age_Ns %>% 
                         filter(HIC == (f == "HIC")) %>% 
                         dplyr::select(-HIC) %>%
                         magrittr::set_colnames(c("Group 1", "N1")), 
                       results, 
                       multiple = "all")
  
  # Add N2
  results <- full_join(GPD_age_Ns %>% 
                         filter(HIC == (f == "HIC")) %>% 
                         dplyr::select(-HIC) %>%
                         magrittr::set_colnames(c("Group 2", "N2")), 
                       results, 
                       multiple = "all")
  
  rownames(results) <- NULL
  
  results <- results %>% drop_na()

  results$p_value <- p.adjust(results$p_value, method = "BH")
  
  results$p_value <- signif(results$p_value, digits = 3)
  
  # results$`P.adj Significant` <- ifelse(as.numeric(results$p_value) < 0.05, 
  #                                       TRUE, FALSE)

  
  results$`P.adj Significant` <- get_stars(results$p_value)
  
  results$p_value <- to_neat(results$p_value)
  
  results$`Income Group` <- f
  

  results <- results %>% dplyr::select(`Group 1`, `Group 2`, 
                                       `N1`, `N2`, `Income Group`,
                                       P.adj = p_value, `P.adj Significant`,
                                       `Effect Direction`)
  
  results$`Effect Direction` <- gsub("Middle>Aged|Middle<Aged",
                                     "Middle-Aged",
                                     results$`Effect Direction`)
  
  
  return(results)
  
}) %>% 
  do.call(rbind, .)
write.csv(age_tukey_GDP_table, "RESULTS/age_tukey_GDP_table_ARG_load.csv")


region_tukey_GPD <- list(HIC = region_tukey_female_HIC,
                         LMIC = region_tukey_female_LMIC)

region_tukey_GDP_table <- lapply(names(region_tukey_GPD), function(f) {

  X <- region_tukey_GPD[[f]]

  tukey_summary <- summary(X)

  # Confidence intervals
  # conf_intervals <- confint(X)

  # Create a tidy data frame
  results <- data.frame(
    Comparison = names(tukey_summary$test$coefficients),
    Estimate = tukey_summary$test$coefficients,
    # SE = tukey_summary$test$sigma,
    # z_value = tukey_summary$test$tstat,
    p_value = tukey_summary$test$pvalues,
    # lower2.5 = as.data.frame(conf_intervals$confint)$lwr,
    # upper97.5 = as.data.frame(conf_intervals$confint)$upr,
    `Income Group` = f
  )

  results$`Effect Direction` <- ifelse(results$Estimate > 0,
                                       gsub("-", ">",results$Comparison), 
                                       gsub("-", "<",results$Comparison))
  
  
  results <- results %>% dplyr::select(-Estimate)
  
  results <- results %>% 
    separate(Comparison, into = c("Group 1", "Group 2"), sep = " - ")
  
  # Add N1
  results <- full_join(GPD_region_Ns %>% 
                         filter(HIC == (f == "HIC")) %>% 
                         dplyr::select(-HIC) %>%
                         magrittr::set_colnames(c("Group 1", "N1")), 
                       results, 
                       multiple = "all")
  
  # Add N2
  results <- full_join(GPD_region_Ns %>% 
                         filter(HIC == (f == "HIC")) %>% 
                         dplyr::select(-HIC) %>%
                         magrittr::set_colnames(c("Group 2", "N2")), 
                       results, 
                       multiple = "all")
  
  rownames(results) <- NULL
  
  results <- results %>% drop_na()
  
  results$p_value <- p.adjust(results$p_value, method = "BH")
  
  results$p_value <- signif(results$p_value, digits = 3)
  
  # results$`P.adj Significant` <- ifelse(as.numeric(results$p_value) < 0.05, 
  #                                       TRUE, FALSE)
  
  results$`P.adj Significant` <- get_stars(results$p_value)
  
  results$p_value <- to_neat(results$p_value)
  
  results$`Income Group` <- f
  
  
  results <- results %>% dplyr::select(`Group 1`, `Group 2`, 
                                       `N1`, `N2`, `Income Group`,
                                       P.adj = p_value, `P.adj Significant`,
                                       `Effect Direction`)
  
  return(results)
  
}) %>% 
  do.call(rbind, .)
write.csv(region_tukey_GDP_table, "RESULTS/region_tukey_GDP_table_ARG_load.csv")



## Diversity **************************************************** ####
# HICs ****************** ####

data_HIC <- tmp[tmp$HIC == "TRUE" & tmp$sex_combined == "female", ]

model_HIC <- glm(
  shannon_diversity ~ High_use + age_category_new + geo_loc_name_country_continent_calc,
  data = data_HIC
)

# Post hoc test for age_category_new
age_tukey_female_HIC <-
  glht(model_HIC, linfct = mcp(age_category_new = "Tukey"))

# Post hoc test for region
region_tukey_female_HIC <- glht(model_HIC,
                                linfct = mcp(geo_loc_name_country_continent_calc = "Tukey"))


# LMICs ***************** ####

data_LMIC <- tmp[tmp$HIC == "FALSE" & tmp$sex_combined == "female", ]

model_LMIC <- glm(
  shannon_diversity ~ High_use + age_category_new +
    geo_loc_name_country_continent_calc,
  data = data_LMIC
)

# Post hoc test for age_category_new
age_tukey_female_LMIC <-
  glht(model_LMIC, linfct = mcp(age_category_new = "Tukey"))

# Post hoc test for region
region_tukey_female_LMIC <-
  glht(model_LMIC, linfct = mcp(geo_loc_name_country_continent_calc = "Tukey"))



## Combine LM/HIC summaries ******************* ####

age_tukey_GPD <- list(HIC = age_tukey_female_HIC, LMIC = age_tukey_female_LMIC)

age_tukey_GDP_table <- lapply(names(age_tukey_GPD), function(f) {
  
  X <- age_tukey_GPD[[f]]
  
  tukey_summary <- summary(X)
  
  # Confidence intervals
  # conf_intervals <- confint(X)
  
  # Create a tidy data frame
  results <- data.frame(
    Comparison = names(tukey_summary$test$coefficients),
    Estimate = tukey_summary$test$coefficients,
    # SE = tukey_summary$test$sigma,
    # z_value = tukey_summary$test$tstat,
    p_value = tukey_summary$test$pvalues,
    # lower2.5 = as.data.frame(conf_intervals$confint)$lwr,
    # upper97.5 = as.data.frame(conf_intervals$confint)$upr,
    `Income Group` = f
  )
  
  results$`Effect Direction` <- ifelse(results$Estimate > 0,
                                       gsub("-", ">",results$Comparison), 
                                       gsub("-", "<",results$Comparison))
  
  
  results <- results %>% dplyr::select(-Estimate)
  
  results <- results %>% 
    separate(Comparison, into = c("Group 1", "Group 2"), sep = " - ")
  
  # Add N1
  results <- full_join(GPD_age_Ns %>% 
                         filter(HIC == (f == "HIC")) %>% 
                         dplyr::select(-HIC) %>%
                         magrittr::set_colnames(c("Group 1", "N1")), 
                       results, 
                       multiple = "all")
  
  # Add N2
  results <- full_join(GPD_age_Ns %>% 
                         filter(HIC == (f == "HIC")) %>% 
                         dplyr::select(-HIC) %>%
                         magrittr::set_colnames(c("Group 2", "N2")), 
                       results, 
                       multiple = "all")
  
  rownames(results) <- NULL
  
  results <- results %>% drop_na()
  
  results$p_value <- p.adjust(results$p_value, method = "BH")
  
  results$p_value <- signif(results$p_value, digits = 3)
  
  results$`P.adj Significant` <- get_stars(results$p_value)
  
  results$p_value <- to_neat(results$p_value)
  
  results$`Income Group` <- f
  
  
  results <- results %>% dplyr::select(`Group 1`, `Group 2`, 
                                       `N1`, `N2`, `Income Group`,
                                       P.adj = p_value, `P.adj Significant`,
                                       `Effect Direction`)
  
  results$`Effect Direction` <- gsub("Middle>Aged|Middle<Aged",
                                     "Middle-Aged",
                                     results$`Effect Direction`)
  
  
  return(results)
  
}) %>% 
  do.call(rbind, .)
write.csv(age_tukey_GDP_table, "RESULTS/age_tukey_GDP_table_shannon.csv")


region_tukey_GPD <- list(HIC = region_tukey_female_HIC,
                         LMIC = region_tukey_female_LMIC)

region_tukey_GDP_table <- lapply(names(region_tukey_GPD), function(f) {
  
  X <- region_tukey_GPD[[f]]
  
  tukey_summary <- summary(X)
  
  # Confidence intervals
  # conf_intervals <- confint(X)
  
  # Create a tidy data frame
  results <- data.frame(
    Comparison = names(tukey_summary$test$coefficients),
    Estimate = tukey_summary$test$coefficients,
    # SE = tukey_summary$test$sigma,
    # z_value = tukey_summary$test$tstat,
    p_value = tukey_summary$test$pvalues,
    # lower2.5 = as.data.frame(conf_intervals$confint)$lwr,
    # upper97.5 = as.data.frame(conf_intervals$confint)$upr,
    `Income Group` = f
  )
  
  results$`Effect Direction` <- ifelse(results$Estimate > 0,
                                       gsub("-", ">",results$Comparison), 
                                       gsub("-", "<",results$Comparison))
  
  
  results <- results %>% dplyr::select(-Estimate)
  
  results <- results %>% 
    separate(Comparison, into = c("Group 1", "Group 2"), sep = " - ")
  
  # Add N1
  results <- full_join(GPD_region_Ns %>% 
                         filter(HIC == (f == "HIC")) %>% 
                         dplyr::select(-HIC) %>%
                         magrittr::set_colnames(c("Group 1", "N1")), 
                       results, 
                       multiple = "all")
  
  # Add N2
  results <- full_join(GPD_region_Ns %>% 
                         filter(HIC == (f == "HIC")) %>% 
                         dplyr::select(-HIC) %>%
                         magrittr::set_colnames(c("Group 2", "N2")), 
                       results, 
                       multiple = "all")
  
  rownames(results) <- NULL
  
  results <- results %>% drop_na()
  
  results$p_value <- p.adjust(results$p_value, method = "BH")
  
  results$p_value <- signif(results$p_value, digits = 3)
  
  results$`P.adj Significant` <- get_stars(results$p_value)
  
  results$p_value <- to_neat(results$p_value)
  
  results$`Income Group` <- f
  
  
  results <- results %>% dplyr::select(`Group 1`, `Group 2`, 
                                       `N1`, `N2`, `Income Group`,
                                       P.adj = p_value, `P.adj Significant`,
                                       `Effect Direction`)
  
  return(results)
  
}) %>% 
  do.call(rbind, .)
write.csv(region_tukey_GDP_table, "RESULTS/region_tukey_GDP_table_shannon.csv")




