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

set.seed(123123)

# Get data in adult_metadata object
source("R-SCRIPTS/prepare_data_for_lm_analysis.R")

## Dummy-encoded linear models ******************************************************* ####
## Pooled model ************************ ####

# log10_ARG_formula <- paste0("log10_ARG_load ~ ",
#                             paste0(dummy_var_names[!grepl(pattern = "Africa|Asia|Infant", x = dummy_var_names)],
#                                    collapse = "+")) %>% as.formula()
# 
# brm_log10_ARG_load_dummy <- brm(formula = log10_ARG_formula,
#                                 data = adult_metadata,
#                                 family = gaussian(), 
#                                 prior = c(
#                                   set_prior("normal(0, 1)", class = "b"),
#                                   set_prior("normal(0, 1)", class = "Intercept")
#                                 ), 
#                                 chains = 2, 
#                                 iter = 5000,
#                                 cores = parallel::detectCores(),
# )
# 
# 
# shannon_formula <- paste0("shannon_diversity ~ ",
#                           paste0(dummy_var_names[!grepl(pattern = "Africa|Asia|Infant", x = dummy_var_names)],
#                                  collapse = "+")) %>% as.formula()
# 
# brm_shannon_dummy <- brm(formula = shannon_formula,
#                          data = adult_metadata,
#                          family = gaussian(), 
#                          prior = c(
#                            set_prior("normal(0, 1)", class = "b"),
#                            set_prior("normal(0, 1)", class = "Intercept")
#                          ), 
#                          chains = 2, 
#                          iter = 5000,
#                          cores = parallel::detectCores(),
# )
# 
# 
# ## Get summaries
# log10_ARG_summary <- posterior_summary(brm_log10_ARG_load_dummy, pars = "b") %>% 
#   data.frame() %>% 
#   rownames_to_column(var = "Feature") %>% 
#   mutate(Feature = gsub("b_", "", Feature), 
#          Response = "log10_ARG") %>% 
#   dplyr::select(-Est.Error) %>% 
#   mutate(Significant = ifelse(Q2.5 > 0 | Q97.5 < 0, TRUE, FALSE))
# 
# 
# shannon_summary <- posterior_summary(brm_shannon_dummy, pars = "b") %>% 
#   data.frame() %>% 
#   rownames_to_column(var = "Feature") %>% 
#   mutate(Feature = gsub("b_", "", Feature), 
#          Response = "Shannon") %>% 
#   dplyr::select(-Est.Error) %>% 
#   mutate(Significant = ifelse(Q2.5 > 0 | Q97.5 < 0, TRUE, FALSE))
# 
# 
# full_summary <- rbind(log10_ARG_summary, shannon_summary)
# 
# full_summary %>%  
#   filter(Feature != "Intercept") %>% 
#   ggplot() + 
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_errorbar(aes(x = Feature, ymin = Q2.5, ymax = Q97.5, color = Response), 
#                 position = "dodge") +
#   # facet_wrap(~Response) +
#   coord_flip() +
#   labs(title = "95% CIs")




## Separate models ********************* ####

## Do separate analyses for each HIC/LMIC, MEN/WOMEN pair ####

# ## log10_ARG_load ************
# fit_list_log10_ARG_load <- list()
# for(ic in 0:1) {
#   for(sex in c("Men", "Women")) {
#     
#     my_data <- adult_metadata %>% 
#       as.data.frame() %>% 
#       filter(sex_combined == sex, 
#              income_group_HIC == ic)
#     
#     
#     if(ic == 0) {
#       
#       # Remove North America and Oceania since no observation for LMIC
#       my_formula <- paste0("log10_ARG_load ~ ",
#                            paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|Infant|HIC|sex",
#                                                          x = dummy_var_names)],
#                                   collapse = "+")) %>% as.formula()
#     }
#     
#     if(ic == 1) {
#       # Remove South America
#       my_formula <- paste0("log10_ARG_load ~ ",
#                            paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Asia|Infant|HIC|sex",
#                                                          x = dummy_var_names)],
#                                   collapse = "+")) %>% as.formula()
#     }
# 
#     
#     my_fit <- brm(formula = my_formula,
#                   data = my_data,
#                   family = gaussian(), 
#                   prior = c(
#                     set_prior("normal(0, 1)", class = "b"),
#                     set_prior("normal(0, 1)", class = "Intercept")
#                     ), 
#                   chains = 2, 
#                   iter = 5000,
#                   cores = parallel::detectCores())
#     
#     fit_list[[sex]][[ic+1]] <- my_fit
#   }
# }
# 
# 
# ## Results
# 
# log10_ARG_full_summary <- lapply(c("Men", "Women"), function(sex) {
#   lapply(0:1, function(ic) {
#     
#     my_summary <- posterior_summary(fit_list_log10_ARG_load[[sex]][[ic+1]], pars = "b") %>% 
#       data.frame() %>% 
#       rownames_to_column(var = "Feature") %>% 
#       mutate(Feature = gsub("b_", "", Feature), 
#              Response = "log10_ARG") %>% 
#       dplyr::select(-Est.Error) %>% 
#       mutate(Significant = ifelse(Q2.5 > 0 | Q97.5 < 0, TRUE, FALSE), 
#              Sex = sex, income_group = ifelse(ic == 1, "HIC", "LMIC"))
#     
#     return(my_summary)
#     
#   }) %>% do.call(rbind,. )
# }) %>% 
#   do.call(rbind,. )
# 
# log10_ARG_full_summary %>%  
#   filter(Feature != "Intercept") %>% 
#   ggplot() + 
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_errorbar(aes(x = Feature, ymin = Q2.5, ymax = Q97.5, color = Sex), 
#                 position = "dodge", width = 0.2) +
#   facet_wrap(~income_group) +
#   coord_flip() +
#   labs(title = "95% CIs")
#     
# 
# ## Shannon Diversity ************
# # Do separate analyses for each HIC/LMIC, MEN/WOMEN pair
# fit_list_shannon <- list()
# for(ic in 0:1) {
#   for(sex in c("Men", "Women")) {
#     
#     my_data <- adult_metadata %>% 
#       as.data.frame() %>% 
#       filter(sex_combined == sex, 
#              income_group_HIC == ic)
#     
#     
#     if(ic == 0) {
#       
#       # Remove North America and Oceania since no observation for LMIC
#       my_formula <- paste0("shannon_diversity ~ ",
#                            paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|Infant|HIC|sex",
#                                                          x = dummy_var_names)],
#                                   collapse = "+")) %>% as.formula()
#     }
#     
#     if(ic == 1) {
#       # Remove South America
#       my_formula <- paste0("shannon_diversity ~ ",
#                            paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Asia|Infant|HIC|sex",
#                                                          x = dummy_var_names)],
#                                   collapse = "+")) %>% as.formula()
#     }
#     
#     
#     my_fit <- brm(formula = my_formula,
#                   data = my_data,
#                   family = gaussian(), 
#                   prior = c(
#                     set_prior("normal(0, 1)", class = "b"),
#                     set_prior("normal(0, 1)", class = "Intercept")
#                   ), 
#                   chains = 2, 
#                   iter = 5000,
#                   cores = parallel::detectCores())
#     
#     fit_list_shannon[[sex]][[ic+1]] <- my_fit
#   }
# }
# 
# 
# ## Results
# 
# shannon_full_summary <- lapply(c("Men", "Women"), function(sex) {
#   lapply(0:1, function(ic) {
#     
#     my_summary <- posterior_summary(fit_list_shannon[[sex]][[ic+1]], pars = "b") %>% 
#       data.frame() %>% 
#       rownames_to_column(var = "Feature") %>% 
#       mutate(Feature = gsub("b_", "", Feature), 
#              Response = "shannon") %>% 
#       dplyr::select(-Est.Error) %>% 
#       mutate(Significant = ifelse(Q2.5 > 0 | Q97.5 < 0, TRUE, FALSE), 
#              Sex = sex, income_group = ifelse(ic == 1, "HIC", "LMIC"))
#     
#     return(my_summary)
#     
#   }) %>% do.call(rbind,. )
# }) %>% do.call(rbind,. )
# 
# shannon_full_summary %>%  
#   filter(Feature != "Intercept") %>% 
#   ggplot() + 
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_errorbar(aes(x = Feature, ymin = Q2.5, ymax = Q97.5, color = Sex), 
#                 position = "dodge", width = 0.2) +
#   facet_wrap(~income_group) +
#   coord_flip() +
#   labs(title = "95% CIs")
# 
# 
#     
#  
# 
# 
# 
# 

## log_ARG_load ************
## Do separate analyses for each HIC/LMIC pair ####

## ARG_load ******************** ####

fit_list_log_ARG_load <- list()
for(ic in 0:1) {
  # for(sex in c("Men", "Women")) {

  my_data <- adult_metadata %>%
    as.data.frame() %>%
    filter(
      # sex_combined == sex,
      income_group_HIC == ic
    )


  if(ic == 0) {

    # Remove Africa, North America and Oceania since no observation for LMIC
    my_formula <- paste0("log_ARG_load ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+")) %>% as.formula()
  }

  if(ic == 1) {
    # Remove Africa, South America, set Europe to base level
    my_formula <- paste0("log_ARG_load ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Europe|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+")) %>% as.formula()
  }


  my_fit <- brm(formula = my_formula,
                data = my_data,
                family = gaussian(),
                prior = c(
                  set_prior("normal(0, 1)", class = "b"),
                  set_prior("normal(0, 1)", class = "Intercept")
                ),
                chains = 2,
                iter = 5000,
                cores = parallel::detectCores())

  fit_list_log_ARG_load[[ic+1]] <- my_fit
  # }
}

# saveRDS(object = fit_list_log_ARG_load, file = "fit_list_log_ARG_load.RDS")
# fit_list_log_ARG_load <- readRDS(file = "fit_list_log_ARG_load.RDS")

## Results
log_ARG_full_summary <-lapply(0:1, function(ic) {
  
  my_summary <- posterior_summary(fit_list_log_ARG_load[[ic+1]], pars = "b") %>% 
    data.frame() %>% 
    rownames_to_column(var = "Feature") %>% 
    mutate(Feature = gsub("b_", "", Feature), 
           Response = "log_ARG") %>% 
    dplyr::select(-Est.Error) %>% 
    mutate(Significant = ifelse(Q2.5 > 0 | Q97.5 < 0, TRUE, FALSE), 
           # Sex = sex,
           income_group = ifelse(ic == 1, "HIC", "LMIC"))
  
  return(my_summary)
  
}) %>% 
  do.call(rbind,. )


log_ARG_full_summary %>%  
  filter(Feature != "Intercept") %>% 
  ggplot() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(x = Feature, ymin = Q2.5, ymax = Q97.5), 
                position = "dodge", width = 0.2) +
  facet_wrap(~income_group) +
  coord_flip() +
  labs(title = "95% CIs")


## Shannon Diversity *********** ####
fit_list_shannon <- list()
for(ic in 0:1) {
  # for(sex in c("Men", "Women")) {

  my_data <- adult_metadata %>%
    as.data.frame() %>%
    filter(
      # sex_combined == sex,
      income_group_HIC == ic
    )


  if(ic == 0) {

    # Remove Africa, North America and Oceania since no observation for LMIC
    my_formula <- paste0("shannon_diversity ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+")) %>% as.formula()
  }

  if(ic == 1) {
    # Remove Africa, South America. Set Europe to base level
    my_formula <- paste0("shannon_diversity ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Europe|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+")) %>% as.formula()
  }


  my_fit <- brm(formula = my_formula,
                data = my_data,
                family = gaussian(),
                prior = c(
                  set_prior("normal(0, 1)", class = "b"),
                  set_prior("normal(0, 1)", class = "Intercept")
                ),
                chains = 2,
                iter = 5000,
                cores = parallel::detectCores())

  fit_list_shannon[[ic+1]] <- my_fit
  # }
}


# saveRDS(object = fit_list_shannon, file = "fit_list_shannon.RDS")
# fit_list_shannon <- readRDS(file = "fit_list_shannon.RDS")

## Results
shannon_full_summary <-lapply(0:1, function(ic) {
  
  my_summary <- posterior_summary(fit_list_shannon[[ic+1]], pars = "b") %>% 
    data.frame() %>% 
    rownames_to_column(var = "Feature") %>% 
    mutate(Feature = gsub("b_", "", Feature), 
           Response = "Shannon") %>% 
    dplyr::select(-Est.Error) %>% 
    mutate(Significant = ifelse(Q2.5 > 0 | Q97.5 < 0, TRUE, FALSE), 
           # Sex = sex,
           income_group = ifelse(ic == 1, "HIC", "LMIC"))
  
  return(my_summary)
  
}) %>% do.call(rbind,. )


## Full summary **************
full_summary <- rbind(shannon_full_summary, log_ARG_full_summary)

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
full_summary$Response <- gsub("log_ARG", "log(ARG load)", full_summary$Response)

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
full_summary_round %>% 
  filter(Response == "log(ARG load)") %>%
  select(-Response) %>% 
  write.csv(., "RESULTS/ARG_dummy_lm.csv")

## Shannon
full_summary_round %>% 
  filter(Response == "Shannon") %>% 
  select(-Response) %>% 
  write.csv(., "RESULTS/shannon_dummy_lm.csv")

## Plot

p <- full_summary %>% 
  mutate(lower = (exp(Q2.5) - 1)*100,
         upper = (exp(Q97.5) - 1)*100) %>%
  filter(Predictor != "Intercept") %>% 
  ggplot() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  # geom_errorbar(aes(x = Predictor, ymin = `exp(Q2.5)`, ymax = `exp(Q97.5)`, color = Response),
  #               position = "dodge", width = 0.2, linewidth = 1) +
  geom_errorbar(aes(x = Predictor, ymin = lower, ymax = upper, color = Response),
                position = "dodge", width = 0.2, linewidth = 1) +
  facet_wrap(~`Income Group`) +
  coord_flip() +
  labs(
    # title = "95% CIs",
    y = "Effect Size (%)", x = "") +
  theme_bw(25) +
  scale_y_continuous(breaks = c(-50, -25, 0,  25, 50, 100)) +
  theme(
    # panel.grid.major.x = element_line(color = "grey"),  # Show major gridlines for x
    # panel.grid.minor.x = element_blank(),              # Hide minor gridlines for x
    # panel.grid.major.y = element_blank(),              # Hide major gridlines for y
    # panel.grid.minor.y = element_blank()               # Hide minor gridlines for y
  ) + 
  scale_color_manual(values = c("black", "darkgrey"))

p
 
png("RESULTS/FIGURES/dummy_lm_CIs2.png",
    units = "in",
    res = 500,
    height = 7,
    width = 16)
print(p)
dev.off()


# library(Cairo)
# CairoJPEG("RESULTS/FIGURES/Fig5.jpg", width=1000, height=480, quality=100)
# print(p)
# dev.off()

# library(Cairo)
# # CairoJPEG("dummy_lm_CIs.jpeg",
# #           width=1000,
# #           height=480
#           # quality=100,
#           # res = 600
#           # )
# 
# CairoJPEG(filename = "dummy_lm_CIs.jpeg", width = 2000, height = 960, 
#           res = 400
#           # pointsize = 12, quality = 1000, bg = "white", res = NA
#           )
# 
# print(p)
# dev.off()











