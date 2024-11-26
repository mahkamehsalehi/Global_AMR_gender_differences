TSE <- read_rds("DATA/TSE.rds")

df <- colData(TSE) %>% as.data.frame()

# Plot age gender by bioproject category
plot_df <- df %>%
  dplyr::select(matches(c("sex_combined", "log10_ARG_load", "host_age_years", "category", "World_Bank_Income_Group")))
# Create the plot, removing NA values from host_age_years
ggplot(na.omit(plot_df), aes(x = host_age_years, y = log10_ARG_load, color = sex_combined)) +
  geom_jitter(size = 0.1, width = 0.1) +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = c(0,1,3,5,10,15,20,25,30,40,50,60,70,80))+
  theme_classic() +
  facet_grid(cols=vars(World_Bank_Income_Group))




ggplot(na.omit(plot_df), aes(x = host_age_years, y = log10_ARG_load, color = sex_combined)) +
  geom_jitter(size = 0.1, width = 0.1) +
  geom_smooth(method = "loess") +
  scale_x_continuous(breaks = c(0,1,3,5,10,15,20,25,30,40,50,60,70,80))+
  theme_classic()


df <- colData(TSE_gender) %>% as.data.frame()
tmp<-df[df$category!="Infant Study",]
tmp <- tmp[tmp$geo_loc_name_country_calc!="Zimbabwe",]
ggplot(tmp, aes(x = Usage, y = log10_ARG_load, color = sex_combined)) +
  geom_jitter(size = 0.1, width = 0.1) +
  geom_smooth(method = "loess") +
  theme_classic() 

library(jtools)
fit <- aov(log10_ARG_load~category, data=df) 

tmp<-df[df$category!="Infant Study",]
tmp <- df[df$geo_loc_name_country_calc!="Zimbawe",]
glm(ARG_div ~ sex_combined, data = df) %>% 
  summ(exp = TRUE, digits = 3)

