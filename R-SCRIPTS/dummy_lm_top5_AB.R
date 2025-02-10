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


# Get data in adult_metadata object
source("R-SCRIPTS/prepare_data_for_lm_analysis.R")

## Dummy-encoded linear models

set.seed(123123)

# Top 5 antibiotics classes as response
# Include covariates from vector dummy_var_names
# Fit separate models for HIC/LMIC

responses <- adult_metadata %>% 
  as.data.frame() %>% 
  select(matches(substr(top_5_AB_classes, 1, 4)) & !matches("log_")) %>%
  # select(top_5_AB_classes) %>%
  colnames
incomes <- 0:1

fit_list <- list()

for(r in responses) {     # Loop over responses
  for(ic in incomes) {    # Loop over income levels
    
    # Filter data
    my_data <- adult_metadata %>%
      as.data.frame() %>%
      # filter(readcount > quantile(adult_metadata$readcount, prob = 0.1)) %>%      ## RARIFY ##
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
    
    # Add offet
    # my_formula <- paste0(my_formula, " + offset(log(readcount))")
    my_formula <- my_formula %>% as.formula()
    
    
    my_fit <- brm(formula = my_formula,
                  data = my_data,
                  family = lognormal(),
                  prior = c(
                    set_prior("normal(0, 1)", class = "b"),
                    set_prior("normal(0, 1)", class = "Intercept")
                  ),
                  chains = 2,
                  iter = 5000,
                  cores = parallel::detectCores())
    
    fit_list[[r]][[ic+1]] <- my_fit
    
  }
}

saveRDS(object = fit_list, file = "dummy_lm_fit_list_TOP5_AB.RDS")


## Results
full_summary <- lapply(responses, function(r) {
  lapply(incomes, function(ic) {
    
    my_summary <- posterior_summary(fit_list[[r]][[ic+1]], pars = "b") %>% 
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
  do.call(rbind,. )



# full_summary %>%  
#   filter(Feature != "Intercept") %>% 
#   ggplot() + 
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_errorbar(aes(x = Feature, ymin = Q2.5, ymax = Q97.5,
#                     color = Response), 
#                 position = "dodge", width = 0.2) +
#   facet_wrap(~income_group) +
#   coord_flip() +
#   labs(title = "95% CIs")



## Full summary **************
full_summary <- full_summary %>% 
  mutate(exp_Estimate = exp(Estimate), 
         exp_Q2.5 = exp(Q2.5), 
         exp_Q97.5 = exp(Q97.5))

# Make neat
full_summary$Feature <- gsub("region", "", full_summary$Feature)
full_summary$Feature <- gsub("age_category_new", "", full_summary$Feature)
full_summary$Feature <- gsub("Usage_high", "High Antibiotic Use", full_summary$Feature)
full_summary$Feature <- gsub("sex_num_Men", "Woman", full_summary$Feature)
full_summary$Feature <- gsub("_", " ", full_summary$Feature)
# full_summary$Response <- gsub("log_ARG_load", "log(ARG load)", full_summary$Response)
# full_summary$Response <- gsub("shannon_diversity", "Shannon", full_summary$Response)

full_summary <- full_summary %>% 
  select(Predictor = Feature,
         Estimate, Q2.5, Q97.5,
         'exp(Estimate)' = exp_Estimate, 
         'exp(Q2.5)' = exp_Q2.5, 
         'exp(Q97.5)' = exp_Q97.5, 
         Response,
         'Income Group' = income_group)


full_summary$Predictor <- factor(full_summary$Predictor,
                                 levels = c("South America", "Oceania", "North America", 
                                            "Europe", "Asia",
                                            "Oldest Adult", "Older Adult",
                                            "Young Adult", "Teenager",
                                            "Children", "Toddler", "Infant",
                                            "High Antibiotic Use", 
                                            "Woman", 
                                            "Intercept"))

# Round
full_summary_round <- full_summary
full_summary_round[, 2:7] <- signif(full_summary[, 2:7], 3)


## Full table
# write.csv(full_summary_round, "ARG_shannon_dummy_lm.csv")

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

## Plot

p <- full_summary %>% 
  mutate(Response = substr(Response, 1, 15)) %>% 
  mutate(lower = (exp(Q2.5) - 1)*100,
         upper = (exp(Q97.5) - 1)*100) %>%
  filter(Predictor != "Intercept") %>% 
  ggplot() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  # geom_errorbar(aes(x = Predictor, ymin = `exp(Q2.5)`, ymax = `exp(Q97.5)`, color = Response),
  #               position = "dodge", width = 0.2, linewidth = 1) +
  geom_errorbar(aes(x = Predictor, ymin = lower, ymax = upper, color = Response),
                position = position_dodge(0.5), width = 0.2, linewidth = 1) +
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
  ) 
  # scale_color_manual(values = c("#AF46B4",  "#4BB446"))

p

# png("RESULTS/FIGURES/Fig5_offset_readcount_rarify_0.1.png",
#     units = "in",
#     res = 500,
#     height = 7,
#     width = 16)
# print(p)
# dev.off()







