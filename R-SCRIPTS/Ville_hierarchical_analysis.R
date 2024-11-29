log10_ARG_formula <- paste0("log10_ARG_load ~ ",
                            paste0(dummy_var_names[!grepl(pattern = "Africa|Asia|Infant", x = dummy_var_names)],
                                   collapse = "+")) %>% as.formula()

my_vars <- dummy_var_names
my_vars <- my_vars[!grepl("Asia|Africa|Infant|sex|income_group", my_vars)]

## Hierarchical analysis
library(rstan)
hier_lm_model <- stan_model("R-SCRIPTS/hierarchical_lm.stan")


N <- adult_metadata %>% nrow
G <- length(my_vars)
y <- adult_metadata$log10_ARG_load
x <- adult_metadata[, my_vars]
ind <- ifelse(adult_metadata$sex_combined == "Men" & adult_metadata$income_group_HIC == 0, 
              1, ifelse(adult_metadata$sex_combined == "Men" & adult_metadata$income_group_HIC == 1, 
                        2, ifelse(adult_metadata$sex_combined == "Women" & adult_metadata$income_group_HIC == 0, 
                                  3, 4)))

adult_metadata$ind <- ind
n_groups <- length(unique(ind))

my_data <- list(N = N, 
                G = G, 
                y = y, 
                x = x %>% as.data.frame(), 
                n_groups = n_groups, 
                ind = ind)

hier_lm_fit <- sampling(hier_lm_model, 
                        my_data, 
                        chains = 2)


## Results
coef_summary <- summary(hier_lm_fit, "coef")$summary %>%
  data.frame() %>% 
  rownames_to_column(var = "dummy") %>% 
  mutate(dummy = gsub("coef\\[", "", dummy)) %>% 
  mutate(dummy = gsub("\\]", "", dummy))%>% 
  separate(dummy, into = c("feature", "ind"), sep = ",") %>% 
  mutate(feature = as.numeric(feature), 
         ind = as.numeric(ind))

coef_summary <- full_join(coef_summary, 
                          data.frame(sex_num_Men = c("Men", "Men", "Women", "Women"), 
                                     income_group_HIC = c(0, 1, 0, 1), 
                                     ind = 1:4), 
                          by = "ind")

for(i in 1:nrow(coef_summary)) {
  coef_summary[i, "feature"] <- my_vars[as.numeric(coef_summary[i, "feature"])]
}


coef_summary %>% 
  ggplot() + 
  geom_errorbar(aes(x = feature, ymin = X2.5., ymax = X97.5., color = sex_num_Men), 
                position = "dodge") + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~income_group_HIC) +
  coord_flip()




