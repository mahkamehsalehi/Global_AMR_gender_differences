## Dummy-encoded linear models

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
library(tidyverse)
library(vegan)
library(TreeSummarizedExperiment)


## Helpers *************************

# Set number <0.001 to "<0.001"
to_neat <- function(x, threshold = 0.001) {
  ifelse(x > threshold, as.character(x), "<0.001")
}

# Get number of stars for P.adj significant column
get_stars_CI_log <- function(lower, upper) {
  
  ifelse(lower > 1 | upper < 1, "*", "ns")
  
}
get_stars_CI <- function(lower, upper) {
  
  ifelse(lower > 0 | upper < 0, "*", "ns")
  
}


# Get data in adult_metadata object
source("R-SCRIPTS/prepare_data_for_lm_analysis.R")

# ARG_load and Shannon as response
# Log-Normal for ARG_load, Normal for Shannon
# Include covariates from vector dummy_var_names
# Fit separate models for HIC/LMIC

responses <- c("ARG_load", "shannon_diversity")
incomes <- 0:1

# Compare different variants: 
# 1. original: no bioproject or readcount
# 2. bioproject as random intercept
# 3. readcount as covariate
# 4. bioproject as random intercept, readcount as covariate

# Generate manuscript Figure 5 based on model 1. 

set.seed(123123)

## Model 1 ************************ ####

fit_list1 <- list()

for(r in responses) {     # Loop over responses
  for(ic in incomes) {    # Loop over income levels
    
    # Filter data
    my_data <- adult_metadata %>%
      as.data.frame() %>%
      filter(income_group_HIC == ic)
    
    if(ic == 0) {
      
      # Remove Africa, North America and Oceania since no observation for LMIC
      my_formula <- paste0(r, " ~ ",
                           paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|HIC",
                                                         x = dummy_var_names)],
                                  collapse = "+"))
      
    }
    
    if(ic == 1) {
      # Remove Africa, South America, set Europe to base level
      my_formula <- paste0(r, " ~ ",
                           paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Europe|HIC",
                                                         x = dummy_var_names)],
                                  collapse = "+"))
    }
    
    my_formula <- my_formula %>% as.formula()
    
    if(r == "shannon_diversity") {
      print(r)
      my_fit <- brm(formula = my_formula,
                    data = my_data,
                    family = gaussian(),
                    prior = c(
                      set_prior("normal(0, 1)", class = "b"),
                      set_prior("normal(0, 1)", class = "Intercept")
                    ),
                    chains = 2,
                    iter = 5000,
                    control = list(adapt_delta = 0.99, max_treedepth = 11),
                    cores = parallel::detectCores())
      
    } else if(r == "ARG_load") {
      print(r)
      my_fit <- brm(formula = my_formula,
                    data = my_data,
                    family = lognormal(),
                    prior = c(
                      set_prior("normal(0, 1)", class = "b"),
                      set_prior("normal(0, 1)", class = "Intercept")
                    ),
                    chains = 2,
                    iter = 5000,
                    control = list(adapt_delta = 0.99, max_treedepth = 11),
                    cores = parallel::detectCores())
      
    }

    
    fit_list1[[r]][[ic+1]] <- my_fit
    
  }
}

saveRDS(object = fit_list1, file = "RESULTS/FITS/dummy_lm_fit_M1.RDS")
# fit_list1 <- readRDS(file = "dummy_lm_fit_M1.RDS")

## Results
full_summary1 <- lapply(responses, function(r) {
  lapply(incomes, function(ic) {
    
    my_summary <- posterior_summary(fit_list1[[r]][[ic+1]])
    # Extract effect only
    my_summary <- my_summary[grep("^b_", rownames(my_summary)), ]
    
    my_summary <- my_summary %>% 
      data.frame() %>% 
      rownames_to_column(var = "Feature") %>% 
      mutate(Feature = gsub("b_", "", Feature), 
             Response = r) %>% 
      dplyr::select(-Est.Error) %>% 
      mutate(Significant = ifelse(Q2.5 > 0 | Q97.5 < 0, TRUE, FALSE), 
             # Sex = sex,
             income_group = ifelse(ic == 1, "HIC", "LMIC"))
    
    return(my_summary)
    
  }) %>% 
    do.call(rbind,. )
}) %>% 
  do.call(rbind,. ) %>% 
  mutate(model = 1)

## Model 2 ************************ ####

fit_list2 <- list()

for(r in responses) {     # Loop over responses
  for(ic in incomes) {    # Loop over income levels
    
    # Filter data
    my_data <- adult_metadata %>%
      as.data.frame() %>%
      filter(income_group_HIC == ic)
    
    if(ic == 0) {
      
      # Remove Africa, North America and Oceania since no observation for LMIC
      my_formula <- paste0(r, " ~ ",
                           paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|HIC",
                                                         x = dummy_var_names)],
                                  collapse = "+"))
      
    }
    
    if(ic == 1) {
      # Remove Africa, South America, set Europe to base level
      my_formula <- paste0(r, " ~ ",
                           paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Europe|HIC",
                                                         x = dummy_var_names)],
                                  collapse = "+"))
    }
    
    # Random intercept for bioproject
    my_formula <- paste0(my_formula, " + (1 | bioproject)")
    
    my_formula <- my_formula %>% as.formula()
    
    if(r == "shannon_diversity") {
      print(r)
      my_fit <- brm(formula = my_formula,
                    data = my_data,
                    family = gaussian(),
                    prior = c(
                      set_prior("normal(0, 1)", class = "b"),
                      set_prior("normal(0, 1)", class = "Intercept"), 
                      set_prior("normal(0, 1)", class = "sd", group = "bioproject")
                    ),
                    chains = 2,
                    iter = 5000,
                    control = list(adapt_delta = 0.99,
                                   max_treedepth = 11),
                    cores = parallel::detectCores())
      
    } else if(r == "ARG_load") {
      print(r)
      my_fit <- brm(formula = my_formula,
                    data = my_data,
                    family = lognormal(),
                    prior = c(
                      set_prior("normal(0, 1)", class = "b"),
                      set_prior("normal(0, 1)", class = "Intercept"), 
                      set_prior("normal(0, 1)", class = "sd", group = "bioproject")
                    ),
                    chains = 2,
                    iter = 5000,
                    control = list(adapt_delta = 0.99,
                                   max_treedepth = 11),
                    cores = parallel::detectCores())
      
    }
    
    
    
    fit_list2[[r]][[ic+1]] <- my_fit
    
  }
}

saveRDS(object = fit_list2, file = "RESULTS/FITS/dummy_lm_fit_M2.RDS")
# fit_list2 <- readRDS(file = "dummy_lm_fit_M2.RDS")

## Results
full_summary2 <- lapply(responses, function(r) {
  lapply(incomes, function(ic) {
    
    my_summary <- posterior_summary(fit_list2[[r]][[ic+1]])
    # Extract effect only
    my_summary <- my_summary[grep("^b_", rownames(my_summary)), ]
    
    my_summary <- my_summary %>% 
      data.frame() %>% 
      rownames_to_column(var = "Feature") %>% 
      mutate(Feature = gsub("b_", "", Feature), 
             Response = r) %>% 
      dplyr::select(-Est.Error) %>% 
      mutate(Significant = ifelse(Q2.5 > 0 | Q97.5 < 0, TRUE, FALSE), 
             # Sex = sex,
             income_group = ifelse(ic == 1, "HIC", "LMIC"))
    
    return(my_summary)
    
  }) %>% 
    do.call(rbind,. )
}) %>% 
  do.call(rbind,. ) %>% 
  mutate(model = 2)


## Model 3 ************************ ####

fit_list3 <- list()

for(r in responses) {     # Loop over responses
  for(ic in incomes) {    # Loop over income levels
    
    # Filter data
    my_data <- adult_metadata %>%
      as.data.frame() %>%
      filter(income_group_HIC == ic)
    
    if(ic == 0) {
      
      # Remove Africa, North America and Oceania since no observation for LMIC
      my_formula <- paste0(r, " ~ ",
                           paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|HIC",
                                                         x = dummy_var_names)],
                                  collapse = "+"))
      
    }
    
    if(ic == 1) {
      # Remove Africa, South America, set Europe to base level
      my_formula <- paste0(r, " ~ ",
                           paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Europe|HIC",
                                                         x = dummy_var_names)],
                                  collapse = "+"))
    }
    
    # Add read count as covariate
    my_formula <- paste0(my_formula, " + log(readcount)")
    
    my_formula <- my_formula %>% as.formula()
    
    
    if(r == "shannon_diversity") {
      print(r)
      my_fit <- brm(formula = my_formula,
                    data = my_data,
                    family = gaussian(),
                    prior = c(
                      set_prior("normal(0, 1)", class = "b"),
                      set_prior("normal(0, 1)", class = "Intercept")
                    ),
                    chains = 2,
                    iter = 5000,
                    control = list(adapt_delta = 0.99, max_treedepth = 11),
                    cores = parallel::detectCores())
      
    } else if(r == "ARG_load") {
      print(r)
      my_fit <- brm(formula = my_formula,
                    data = my_data,
                    family = lognormal(),
                    prior = c(
                      set_prior("normal(0, 1)", class = "b"),
                      set_prior("normal(0, 1)", class = "Intercept")
                    ),
                    chains = 2,
                    iter = 5000,
                    control = list(adapt_delta = 0.99, max_treedepth = 11),
                    cores = parallel::detectCores())
      
    }
    
    fit_list3[[r]][[ic+1]] <- my_fit
    
  }
}

saveRDS(object = fit_list3, file = "RESULTS/FITS/dummy_lm_fit_M3.RDS")
# fit_list3 <- readRDS(file = "dummy_lm_fit_M3.RDS")

## Results
full_summary3 <- lapply(responses, function(r) {
  lapply(incomes, function(ic) {
    
    my_summary <- posterior_summary(fit_list3[[r]][[ic+1]])
    # Extract effect only
    my_summary <- my_summary[grep("^b_", rownames(my_summary)), ]
    
    my_summary <- my_summary %>% 
      data.frame() %>% 
      data.frame() %>% 
      rownames_to_column(var = "Feature") %>% 
      mutate(Feature = gsub("b_", "", Feature), 
             Response = r) %>% 
      dplyr::select(-Est.Error) %>% 
      mutate(Significant = ifelse(Q2.5 > 0 | Q97.5 < 0, TRUE, FALSE), 
             # Sex = sex,
             income_group = ifelse(ic == 1, "HIC", "LMIC"))
    
    return(my_summary)
    
  }) %>% 
    do.call(rbind,. )
}) %>% 
  do.call(rbind,. ) %>% 
  mutate(model = 3)

## Model 4 ************************ ####

fit_list4 <- list()

for(r in responses) {     # Loop over responses
  for(ic in incomes) {    # Loop over income levels
    
    # Filter data
    my_data <- adult_metadata %>%
      as.data.frame() %>%
      filter(income_group_HIC == ic)
    
    if(ic == 0) {
      
      # Remove Africa, North America and Oceania since no observation for LMIC
      my_formula <- paste0(r, " ~ ",
                           paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|HIC",
                                                         x = dummy_var_names)],
                                  collapse = "+"))
      
    }
    
    if(ic == 1) {
      # Remove Africa, South America, set Europe to base level
      my_formula <- paste0(r, " ~ ",
                           paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Europe|HIC",
                                                         x = dummy_var_names)],
                                  collapse = "+"))
    }
    
    # Add readcount as covariate
    my_formula <- paste0(my_formula, " + log(readcount)")
    # Random intercept for bioproject
    my_formula <- paste0(my_formula, " + (1 | bioproject)")
    
    my_formula <- my_formula %>% as.formula()
    
    
    if(r == "shannon_diversity") {
      print(r)
      my_fit <- brm(formula = my_formula,
                    data = my_data,
                    family = gaussian(),
                    prior = c(
                      set_prior("normal(0, 1)", class = "b"),
                      set_prior("normal(0, 1)", class = "Intercept"), 
                      set_prior("normal(0, 1)", class = "sd", group = "bioproject")
                    ),
                    chains = 2,
                    iter = 5000,
                    control = list(adapt_delta = 0.99, max_treedepth = 11),
                    cores = parallel::detectCores())
      
    } else if(r == "ARG_load") {
      print(r)
      my_fit <- brm(formula = my_formula,
                    data = my_data,
                    family = lognormal(),
                    prior = c(
                      set_prior("normal(0, 1)", class = "b"),
                      set_prior("normal(0, 1)", class = "Intercept"), 
                      set_prior("normal(0, 1)", class = "sd", group = "bioproject")
                    ),
                    chains = 2,
                    iter = 5000,
                    control = list(adapt_delta = 0.99, max_treedepth = 11),
                    cores = parallel::detectCores())
      
    }
    
    fit_list4[[r]][[ic+1]] <- my_fit
    
  }
}

saveRDS(object = fit_list4, file = "RESULTS/FITS/dummy_lm_fit_M4.RDS")
# fit_list4 <- readRDS(file = "dummy_lm_fit_M4.RDS")

## Results
full_summary4 <- lapply(responses, function(r) {
  lapply(incomes, function(ic) {
    
    my_summary <- posterior_summary(fit_list4[[r]][[ic+1]])
    # Extract effect only
    my_summary <- my_summary[grep("^b_", rownames(my_summary)), ]
    
    my_summary <- my_summary %>% 
      data.frame() %>% 
      rownames_to_column(var = "Feature") %>% 
      mutate(Feature = gsub("b_", "", Feature), 
             Response = r) %>% 
      dplyr::select(-Est.Error) %>% 
      mutate(Significant = ifelse(Q2.5 > 0 | Q97.5 < 0, TRUE, FALSE), 
             # Sex = sex,
             income_group = ifelse(ic == 1, "HIC", "LMIC"))
    
    return(my_summary)
    
  }) %>% 
    do.call(rbind,. )
}) %>% 
  do.call(rbind,. ) %>% 
  mutate(model = 4)
## Combine results **************** ####

## Full summary **************
full_summary <- list(full_summary1, 
                     full_summary2, 
                     full_summary3, 
                     full_summary4) %>% 
  do.call(rbind, .)

# Recode model identifier
full_summary <- full_summary %>% 
  mutate(model = as.character(model)) %>% 
  mutate(model = recode(model,
                        "1" = "original", 
                        "2" = "+(1|bioproject)", 
                        "3" = "+log(readcount)", 
                        "4" = "+(1|bioproject)+log(readcount)"))

full_summary$model <- factor(full_summary$model, 
                             levels = c("original", 
                                        "+(1|bioproject)", 
                                        "+log(readcount)", 
                                        "+(1|bioproject)+log(readcount)"))

# Make neat
full_summary$Feature <- gsub("region", "", full_summary$Feature)
full_summary$Feature <- gsub("age_category_new", "", full_summary$Feature)
full_summary$Feature <- gsub("Usage_high", "High Antibiotic Use", full_summary$Feature)
full_summary$Feature <- gsub("sex_num_Men", "Woman", full_summary$Feature)
full_summary$Feature <- gsub("_", " ", full_summary$Feature)
full_summary$Response <- gsub("ARG_load", "ARG load", full_summary$Response)
full_summary$Response <- gsub("shannon_diversity", "Shannon", full_summary$Response)

full_summary <- full_summary %>% 
  select(Predictor = Feature,
         Estimate, Q2.5, Q97.5,
         Response,
         'Income Group' = income_group, 
         model)

full_summary <- full_summary %>% 
  mutate('exp(Estimate)' = exp(Estimate),
         'exp(Q2.5)' = exp(Q2.5),
         'exp(Q97.5)' = exp(Q97.5),)



full_summary$Predictor <- factor(full_summary$Predictor,
                                 levels = c("South America", "Oceania", "North America", 
                                            "Europe", "Asia",
                                            "Oldest Adult", "Older Adult",
                                            "Young Adult", "Teenager",
                                            "Children", "Toddler", "Infant",
                                            "High Antibiotic Use", 
                                            "Woman", 
                                            "Intercept"))


full_summary <- full_summary %>% rename(Model = model)

# Rename model 1
full_summary <- full_summary %>% 
  mutate(Model = recode(full_summary$Model, "original" = "Reference"))


# Round
# full_summary_round <- full_summary
# full_summary_round[, 2:7] <- signif(full_summary[, 2:7], 3)

## Log ARG load
# full_summary_round %>% 
#   filter(Response == "log(ARG load)") %>%
#   select(-Response) %>% 
#   write.csv(., "RESULTS/ARG_dummy_lm.csv")

## Shannon
# full_summary_round %>% 
#   filter(Response == "Shannon") %>% 
#   select(-Response) %>% 
#   write.csv(., "RESULTS/shannon_dummy_lm.csv")

## Plots ************************** ####

## Figure 5 ************ ####

p <- full_summary %>% 
  filter(Model == levels(full_summary$Model)[1]) %>%
  mutate(lower = (exp(Q2.5) - 1)*100,
         upper = (exp(Q97.5) - 1)*100) %>%
  filter(Predictor != "Intercept") %>% 
  ggplot() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  # geom_errorbar(aes(x = Predictor, ymin = `exp(Q2.5)`, ymax = `exp(Q97.5)`, color = Response),
  #               position = "dodge", width = 0.2, linewidth = 1) +
  geom_errorbar(aes(x = Predictor, ymin = lower, ymax = upper, color = Response),
                position = "dodge", width = 0.4, linewidth = 1.5) +
  facet_wrap(~`Income Group`) +
  coord_flip() +
  labs(
    # title = "95% CIs",
    y = "Effect Size (%)", x = "") +
  theme_bw(25) +
  # scale_y_continuous(breaks = c(-50, -25, 0,  25, 50, 100)) +
  theme(
    # panel.grid.major.x = element_line(color = "grey"),  # Show major gridlines for x
    # panel.grid.minor.x = element_blank(),              # Hide minor gridlines for x
    # panel.grid.major.y = element_blank(),              # Hide major gridlines for y
    # panel.grid.minor.y = element_blank()               # Hide minor gridlines for y
  ) + 
  scale_color_manual(values = c("#1f77b4", "#ff7f0e")) +
  facet_wrap(~`Income Group`, ncol = 2)



png("RESULTS/FIGURES/Fig5_lm.png",
    units = "in",
    res = 500,
    height = 7,
    width = 16)
print(p)
dev.off()


## Comparison plot ***** ####
## Supplementary figure

comparison_p <- full_summary %>% 
  # filter(Model %in% levels(full_summary$Model)[1:2]) %>% 
  mutate(lower = (exp(Q2.5) - 1)*100,
         upper = (exp(Q97.5) - 1)*100) %>%
  filter(Predictor != "Intercept") %>% 
  ggplot() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  # geom_errorbar(aes(x = Predictor, ymin = `exp(Q2.5)`, ymax = `exp(Q97.5)`, color = Response),
  #               position = "dodge", width = 0.2, linewidth = 1) +
  geom_errorbar(aes(x = Predictor, ymin = lower, ymax = upper, color = Response),
                position = "dodge", width = 0.5, linewidth = 1.5) +
  # facet_wrap(~`Income Group`) +
  coord_flip() +
  labs(
    # title = "95% CIs",
    y = "Effect Size (%)", x = "") +
  theme_bw(25) +
  # scale_y_continuous(breaks = c(-50, -25, 0,  25, 50, 100)) +
  theme(
    # panel.grid.major.x = element_line(color = "grey"),  # Show major gridlines for x
    # panel.grid.minor.x = element_blank(),              # Hide minor gridlines for x
    # panel.grid.major.y = element_blank(),              # Hide major gridlines for y
    # panel.grid.minor.y = element_blank()               # Hide minor gridlines for y
  ) + 
  scale_color_manual(values = c("#1f77b4", "#ff7f0e")) +
  facet_wrap(Model ~ `Income Group`, ncol = 2)



png("RESULTS/FIGURES/SFig_lm_comparison.png",
    units = "in",
    res = 500,
    height = 20,
    width = 20)
print(comparison_p)
dev.off()

## Print tables ******************* ####

## Comparison ********** ####
comparison_table <- full_summary

comparison_table$`exp(Estimate)` <- round(comparison_table$`exp(Estimate)`, 3)
comparison_table$`exp(Q2.5)` <- round(comparison_table$`exp(Q2.5)`, 3)
comparison_table$`exp(Q97.5)` <- round(comparison_table$`exp(Q97.5)`, 3)

comparison_table$Significant <- get_stars_CI_log(lower = comparison_table$`exp(Q2.5)`,
                                            upper = comparison_table$`exp(Q97.5)`)

comparison_table$`Estimate` <- round(comparison_table$`Estimate`, 3)
comparison_table$`Q2.5` <- round(comparison_table$`Q2.5`, 3)
comparison_table$`Q97.5` <- round(comparison_table$`Q97.5`, 3)

mono_table <- comparison_table

comparison_table <- comparison_table %>% 
  arrange(desc(Predictor), desc(Response), `Income Group`)

comparison_table <- comparison_table %>% 
  select(colnames(comparison_table)[c(1, 5, 6, 7, 2:4, 8:11)])

write.csv(comparison_table, "RESULTS/TABLES/lm_comparison.csv")
  
## Shannon ************* ####

shannon_table <- mono_table %>% 
  filter(Response == "Shannon", 
         Model == "Reference") %>% 
  select(1:4, 8:10, 6) %>% 
  arrange(desc(`Income Group`))


write.csv(shannon_table, "RESULTS/TABLES/lm_shannon.csv")



## ARG load ************ ####

ARG_table <- mono_table %>% 
  filter(Response == "ARG load", 
         Model == "Reference") %>% 
  select(1:4, 8:10, 6) %>% 
  arrange(desc(`Income Group`))


write.csv(ARG_table, "RESULTS/TABLES/lm_ARG.csv")





