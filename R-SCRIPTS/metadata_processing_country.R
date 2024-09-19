library(dplyr)
library(lme4)
library(nlme)
library(readr)
# Read tse
# TSE <- readRDS("TSE.rds")
df <- TSE %>% colData() %>% as.data.frame() 
df$log10_ARG_load <- log10(df$ARG_load)

country_level_data <- read_csv("DATA/data.csv") %>%
  .[, -c(20:24)] %>%
  rename(geo_loc_name_country_calc = Country)

country_level_data <- country_level_data %>%
  mutate(geo_loc_name_country_calc = recode(geo_loc_name_country_calc,
                                            "Kazakistan" = "Kazakhstan",
                                            "Hongkong" = "Hong Kong",
                                            "Russian Federation" = "Russia",
                                            "United States of America" = "USA"))

# # Infant projects
# infant_prj<- c("PRJNA971895",
# "SRR24577893",
# "PRJNA872399",
# "PRJNA870588",
# "PRJNA630999",
# "PRJNA872399",
# "PRJNA698986",
# "PRJNA379120",
# "PRJNA396794",
# "PRJNA294605",
# "PRJNA301903",
# "PRJNA521878",
# "PRJEB51728",
# "PRJNA774819",
# "PRJNA613054",
# "PRJNA524703",
# "PRJNA806984",
# "PRJNA376566",
# "PRJNA489090",
# "PRJNA497734",
# "PRJEB55713",
# "PRJNA726032",
# "PRJEB9169",
# "PRJEB29052",
# "PRJEB24006",
# "PRJNA701480",
# "PRJNA648487",
# "PRJNA633576",
# "PRJEB52748",
# "PRJNA613032",
# "PRJNA688164",
# "PRJNA695570",
# "PRJNA327106",
# "PRJEB32135",
# "PRJNA290380",
# "PRJNA595749",
# "PRJEB52147")
# 
# 
# 
# # Sick people 
# sick_prj <- c("PRJNA728680", # MDR infections
# "PRJNA447983", # Colorectal cancer
# "PRJNA674880", # Clostridioides difficile
# "PRJNA786578", # Gastroenteritis
# "PRJNA928758", # Liver transplant infections
# "PRJNA670323", # Critically ill, ICU, infection
# "PRJNA764525", # Blood stream infection
# "PRJNA707487", # Hematopoietic cell transplantation
# "PRJNA531921", # MDR infections
# "PRJNA674880", # Clostridioides difficile
# "PRJNA650244", # COVID
# "PRJNA561398",# ICU patient
# "PRJNA754518", # Colorectal cancer
# "PRJNA637785", # FMT obesity
# "PRJNA698223", # Traveller's diarrhea
# "PRJNA746613", # Hematopoietic stem cell transplantation recipients
# "PRJNA394877", # Marrow transplant patient
# "PRJNA477326", # Blood stream infection
# "PRJNA321753", # Salmonella infection
# "PRJNA655185", # VRE and CRO infection
# "PRJNA689961", # COVID
# "PRJNA668607", # Vibrio cholerae infection
# "PRJNA297252", # Hospitalized
# "PRJEB63351", # Neuroblastoma cancer
# "PRJNA525982", # Leukemia
# "PRJNA787952", # Hematopoietic cell transplantation
# "PRJNA777429", # VRE infection
# "PRJEB35738", # Clostridioides difficile
# "PRJNA564397", # Diarrhea
# "PRJNA431482", # Bechet's disease
# "PRJNA356225", # Vogt-Koyanagi-Harada disease
# "PRJNA701982", # Blood stream infection
# "PRJNA398089", # IBD
# "PRJNA589866", # Antibiotic treatment
# "PRJNA523592", # Blood stream infection
# "PRJNA486009", # Diarrhea
# "PRJNA673918", # FMT
# "PRJNA317435", # Leukemia
# "PRJNA608678", # Vibrio cholerae infection
# "PRJNA545312", # Hematopoietic cell transplantation
# "PRJNA400628", # Urinary infection
# "PRJNA872116", # Immuno deficient, infection?
# "PRJNA283642", # Cell transplantation
# "PRJNA715586", # Leukemia
# "PRJEB9150", # Diarrhea
# "PRJNA660443", # Campylobacterial infection
# "PRJNA731589", # Colorectal cancer
# "PRJNA453965", # Breast cancer
# "PRJNA565566", # Salmonella
# "PRJNA763023", # Colorectal cancer
# "PRJNA766067", # Chirrosis with antibiotics
# "PRJDB13214", # COVID
# "PRJNA397906", # Melanoma
# "PRJNA394687", # Diarrhea
# "PRJNA420371", # C. diff
# "PRJEB55366", # Colorectal cancer
# "PRJNA389280", # IBD
# "PRJNA685914", # C. diff
# "PRJNA438847", # Cystic fibrosis
# "PRJNA641975", # FMT
# "PRJNA862908", # Gasctric infection
# "PRJNA629392", # CHAPLE disease
# "PRJEB15371", # Crohn's disease
# "PRJNA339012", # C. diff
# "PRJDB9649", # FMT
# "PRJNA311208", # Short bowel syndrome
# "PRJNA628604", # FMT, infection
# "PRJEB36140", # FMT
# "PRJNA961974", # Not sick but plant extract fermentation
# "PRJNA615162", # Intestinal systemic schlerosis
# "PRJEB1775", # Bacterial infection
# "PRJNA908301", # COVID
# "PRJNA513350", # CRE infection
# "PRJNA692325", # Colorectal cancer,
# "PRJDB10484", # Hepatic encelopathy,
# "PRJDB10829", # CRC
# "PRJDB4176", # CRC
# "PRJEB12357", # FMT
# "PRJEB12449", # CRC
# "PRJEB27928", # CRC
# "PRJEB38625", # Pancreatic cancer
# "PRJEB39023", # Clostridium difficile
# "PRJEB42013", # Pancreatic cancer
# "PRJEB43119", # Melanoma
# "PRJEB47061", # FMT colitis
# "PRJEB50260", # Ulcerative colitis
# "PRJEB60097", # Colorectal cancer
# "PRJEB65297", # end stage renal disease
# "PRJEB8094", # Antibiotic intervention
# "PRJNA1004651", # IBD
# "PRJNA240307" # Clostridium difficile FMT
# 
# 
# )  
# 
# 


# Find the geo_loc_name_country_calc values in country_level_data but not in df
missing_countries <- df %>%
  select(geo_loc_name_country_calc) %>%
  anti_join(country_level_data %>% select(geo_loc_name_country_calc), by = "geo_loc_name_country_calc")

# Print or inspect the result
print(missing_countries$geo_loc_name_country_calc %>% sort %>% unique)

# Find the geo_loc_name_country_calc values in country_level_data but not in df
missing_countries_2 <- country_level_data %>%
  select(geo_loc_name_country_calc) %>%
  anti_join(df %>% select(geo_loc_name_country_calc), by = "geo_loc_name_country_calc")


# Print or inspect the result
print(missing_countries$geo_loc_name_country_calc %>% sort %>% unique)
print(missing_countries_2$geo_loc_name_country_calc %>% sort %>% unique)

merged_data <- df %>%
  dplyr::left_join(country_level_data %>%
              unique(), by = "geo_loc_name_country_calc") %>%
  mutate(log10_ARG_load = log10(ARG_load)) # %>%
 # filter(geo_loc_name_country_calc != "" & geo_loc_name_country_calc != "uncalculated")

# Create the mapping of countries to regions
country_to_region <- c(
  "Cameroon" = "Sub-Saharan Africa",
  "Democratic Republic of the Congo" = "Sub-Saharan Africa",
  "Liberia" = "Sub-Saharan Africa",
  "Madagascar" = "Sub-Saharan Africa",
  "Mali" = "Sub-Saharan Africa",
  "Mozambique" = "Sub-Saharan Africa",
  "Romania" = "Europe & Central Asia",
  "Kazakhstan" = "Europe & Central Asia",
  "USA" = "North America",
  "Hong Kong" = "East Asia & Pacific",
  "Russia" = "Europe & Central Asia",
  "Taiwan" = "East Asia & Pacific",
  "Tanzania" = "Sub-Saharan Africa"
)

# Assign regions to the corresponding countries in merged_data
merged_data <- merged_data %>%
  mutate(Region = case_when(
    geo_loc_name_country_calc %in% names(country_to_region) ~ country_to_region[geo_loc_name_country_calc],
    TRUE ~ Region  # keep existing region if the country is not in the mapping
  ))

# Factorize and clean names
merged_data$`Koppen Climate Classification of Capital City` <- as.factor(merged_data$`Koppen Climate Classification of Capital City`)
merged_data$`World Bank Income Group` <- as.factor(merged_data$`World Bank Income Group`)
merged_data$geo_loc_name_country_calc <- as.factor(merged_data$geo_loc_name_country_calc)
merged_data$Region <- as.factor(merged_data$Region)
names(merged_data) <- gsub(" ", "_", names(merged_data))
names(merged_data) <- gsub("&", "and", names(merged_data))

# Add to colData in TSE
rownames(merged_data) <- merged_data$acc
colData(TSE) <- DataFrame(merged_data)

# Filter merged_data based on the specified conditions to remove remaining 16S, virome or isolate libraries
# merged_data <- merged_data %>%
#  filter(readcount > 0.5, ARG_load > 50, ARG_obs > 1)

# # Remove sick cohorts
# merged_data_noinf <- merged_data[!merged_data$bioproject%in%c(sick_prj, infant_prj),]
# merged_data_noinf$bioproject %>% unique %>% sort
# 
#  merged_data_noinf[,c(1,2,6,17:24)] %>% View()


# fit1<- lm(log10_ARG_load ~    Usage + `Infrastructure Index`, data=merged_data[merged_data$bioproject%in%c(sick_prj),])
# summary(fit1)

# data <- na.omit(merged_data[, c("acc",
#                                 "Region",
#                                 "Koppen Climate Classification of Capital City", 
#                                 "Usage", 
#                                 "Education Index", 
#                                 "World Bank Income Group", 
#                                 "Infrastructure Index", 
#                                 "GDP per head", 
#                                 "geo_loc_name_country_calc", 
#                                 "Climate Index",
#                                 "Corruption & Governance Index",
#                                 "geo_loc_name_country_continent_calc",
#                                 "log10_ARG_load",
#                                 "Health Spend Index",
#                                 "Logit of Average All Resistance Unstandardized",
#                                 "Logit of Average Ecoli Resistance unstandardized",
#                                 "Logit of META-DATA E.coli per cent resistance FQs*",
#                                 "Logit of META-DATA E.coli per cent resistance 3rd gen ceph*",
#                                 "Average All Resistance Unstandardized",
#                                 "Average Ecoli Resistance to 3GCeph & FQ unstandardized",
#                                 "META-DATA E.coli per cent resistance FQs*",
#                                 "META-DATA E.coli per cent resistance 3rd gen ceph*")])
# 
# # df <- data
# # 
# # # Factorize and clean names
# # df$`Koppen Climate Classification of Capital City` <- as.factor(df$`Koppen Climate Classification of Capital City`)
# # df$`World Bank Income Group` <- as.factor(df$`World Bank Income Group`)
# # df$geo_loc_name_country_calc <- as.factor(df$geo_loc_name_country_calc)
# # df$Region <- as.factor(df$Region)
# # names(df) <- gsub(" ", "_", names(df))
# # names(df) <- gsub("&", "and", names(df))
# # # 
# # # Fit the model
# # fit1 <- lm(log10_ARG_load ~ Average_All_Resistance_Unstandardized, data=df)
# # summary(fit1)
# # # Calculate R^2 value
# # r_squared <- summary(fit1)$r.squared
# # 
# # # Plotting
# # plot_1 <- ggplot(df, aes(x =Average_All_Resistance_Unstandardized, y = log10_ARG_load)) +
# #   geom_jitter(width = 1, height = 0.0, color = "black", alpha = 0.1, size = 0.1) + # Jittered points
# #   geom_smooth(method = "lm", se = FALSE, color = "red") +               # Fitted line
# #   labs(title = "Average clinical isolate resistance",
# #        x = " E. coli, Klebsiella spp, and Staphylococcus aureus\n Resistance (unstandardized)",
# #        y = "log10 ARG Load") +
# #   theme_minimal() +
# #   annotate("text", x = Inf, y = Inf, label = paste("R² =", round(r_squared, 3)),
# #            hjust = 1.1, vjust = 1.1, size = 5, color = "black")
# # 
# # # Print the plot
# # print(plot_1)
# # 
# # # Fit the model
# # fit2 <- lm(log10_ARG_load ~ `Average_Ecoli_Resistance_to_3GCeph_and_FQ_unstandardized`, data=df)
# # summary(fit2)
# # # Calculate R^2 value
# # r_squared <- summary(fit2)$r.squared
# # 
# # # Plotting
# # plot_2 <- ggplot(df, aes(x = `Average_Ecoli_Resistance_to_3GCeph_and_FQ_unstandardized`, y = log10_ARG_load)) +
# #   geom_jitter(width = 1, height = 0.0, color = "black", alpha = 0.1, size = 0.1) + # Jittered points
# #   geom_smooth(method = "lm", se = FALSE, color = "red") +               # Fitted line
# #   labs(title = "Average clinical E. coli isolate resistance",
# #        x = "Average E. coli Resistance (unstandardized)",
# #        y = "log10 ARG Load") +
# #   theme_minimal() +
# #   annotate("text", x = Inf, y = Inf, label = paste("R² =", round(r_squared, 3)),
# #            hjust = 1.1, vjust = 1.1, size = 5, color = "black")
# # 
# # # Print the plot
# # print(plot_2)
# # 
# # # Fit the model
# # fit3 <- lm(log10_ARG_load ~ Usage, data=df)
# # summary(fit3)
# # # Calculate R^2 value
# # r_squared <- summary(fit3)$r.squared
# # 
# # 
# plot_3 <- ggplot(df, aes(x = Usage , y = log10_ARG_load)) +
#   geom_jitter(width = 1, height = 0.0, color = "black", alpha = 0.1, size = 0.1) + # Jittered points
#   geom_smooth(method = "lm", se = FALSE, color = "red") +               # Fitted line
#   labs(title = "Average antibiotic use",
#        x = "Average antibiotic use",
#        y = "log10 ARG Load") +
#   theme_minimal() +
#   annotate("text", x = Inf, y = Inf, label = paste("R² =", round(r_squared, 3)),
#            hjust = 1.1, vjust = 1.1, size = 5, color = "black")
# 
# # Print the plot
# print(plot_3)
# 
# # Fit the model
# fit4 <- lm(log10_ARG_load ~ Infrastructure_Index, data=df)
# summary(fit4)
# # Calculate R^2 value
# r_squared <- summary(fit4)$r.squared
# 
# plot_4 <- ggplot(df, aes(x = Infrastructure_Index , y = log10_ARG_load)) +
#   geom_jitter(width = 0.1, height = 0.0, color = "black", alpha = 0.1, size = 0.1) + # Jittered points
#   geom_smooth(method = "lm", se = FALSE, color = "red") +               # Fitted line
#   labs(title = "Average infrastructure index",
#        x = "Average infrastructure index",
#        y = "log10 ARG Load") +
#   theme_minimal() +
#   annotate("text", x = Inf, y = Inf, label = paste("R² =", round(r_squared, 3)),
#            hjust = 1.1, vjust = 1.1, size = 5, color = "black")
# 
# # Print the plot
# print(plot_4)
# 
# # Fit the model
# fit5 <- lm(log10_ARG_load ~ `Corruption_and_Governance_Index`, data=df)
# summary(fit5)
# # Calculate R^2 value
# r_squared <- summary(fit5)$r.squared
# 
# plot_5 <- ggplot(df, aes(x = `Corruption_and_Governance_Index` , y = log10_ARG_load)) +
#   geom_jitter(width = 0.1, height = 0.0, color = "black", alpha = 0.1, size = 0.1) + # Jittered points
#   geom_smooth(method = "lm", se = FALSE, color = "red") +               # Fitted line
#   labs(title = "Average governance index",
#        x = "Average goverance index",
#        y = "log10 ARG Load") +
#   theme_minimal() +
#   annotate("text", x = Inf, y = Inf, label = paste("R² =", round(r_squared, 3)),
#            hjust = 1.1, vjust = 1.1, size = 5, color = "black")
# 
# # Print the plot
# print(plot_5)
# 
# # mixed effects
# model1 <- lm(log10_ARG_load ~ ., 
#               data = df %>% 
#                dplyr::select(-matches(c("resistance", 
#                                         "acc", "geo_loc_name_country_calc", 
#                                         "geo_loc_name_country_continent_calc", "World", "Region", "Koppen")))) 
# 
# summary(model1)
# 
# lm(log10_ARG_load~Usage, data=df) %>% summary
# 
# 
# # Reorder
# df$geo_loc_name_country_calc <- reorder(df$geo_loc_name_country_calc,
#                                         df$log10_ARG_load, median)
# 
# plot_6 <- ggplot(df, aes(x = geo_loc_name_country_calc, y = log10_ARG_load, fill = Region)) +
#      geom_boxplot() +                                            # Box plot
#      labs(title = "Distribution of log10 ARG Load by Country",
#                   x = "Country",
#                    y = "log10 ARG Load",
#                    fill = "Region") +                                  # Legend title for the fill aesthetic
#      theme_minimal() +                                           # Minimal theme
#      theme(axis.text.x = element_text(angle = 45, hjust = 1))    # Rotate x-axis labels
#  # Print the plot
#    print(plot_6)
# 
# # Reorder
# df$geo_loc_name_country_calc <- reorder(df$geo_loc_name_country_calc, df$Logit_of_Average_Ecoli_Resistance_unstandardized, median)
# 
# # Create the plot
# plot_7 <- ggplot(df, aes(x = geo_loc_name_country_calc, y = Logit_of_Average_Ecoli_Resistance_unstandardized, color = Region)) +
#   geom_point() +                                            # Box plot
#   labs(title = "Distribution Logit of Average E. coli Resistance",
#        x = "Country",
#        y = "log10 ARG Load",
#        fill = "Region") +                                  # Legend title for the fill aesthetic
#   theme_minimal() +                                           # Minimal theme
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))    # Rotate x-axis labels
# 
# # Print the plot
# print(plot_7)
# 
# 
# library(cowplot)
# # Print the plot
# cowplot::plot_grid(plot_6, plot_7, align = "hv")
# 
# library(ggplot2)
# library(dplyr)
# 
# # Reorder the countries by usage
# merged_data <- merged_data %>%
#   arrange(Usage) %>%
#   mutate(geo_loc_name_country_calc = factor(geo_loc_name_country_calc, levels = unique(geo_loc_name_country_calc)))
# 
# # Create the box plot faceted by geo_loc_name_country_continent_calc
# plot_6 <- ggplot(merged_data, aes(x = geo_loc_name_country_calc, y = log10_ARG_load)) +
#   geom_boxplot() +                                            # Box plot
#   labs(title = "Distribution of log10 ARG Load by Country and Continent",
#        x = "Country",
#        y = "log10 ARG Load") +
#   theme_minimal() +                                           # Minimal theme
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
#   facet_wrap(~ geo_loc_name_country_continent_calc, scales = "free_x")  # Facet by Continent
# 
# # Print the plot
# print(plot_6)
# 
# 
