library(mia)
library(vegan)
library(ggpubr)
library(ggpubr)
library(patchwork)
library(tidyverse)
set.seed(4522)

devtools::load_all("~/Rpackages/microbiome/miaverse/mia")
devtools::load_all("~/Rpackages/microbiome/miaverse/miaViz")

# ---------------------------
# Data Preparation
# ---------------------------

tse <- readRDS("../DATA/TSE_filtered.rds")

# Prevalent features only (to speed up the analyses)
altExp(tse, "Prevalent") <- subsetByPrevalent(tse, assay.type="relabundance", detection=0, prevalence=1/100)

# Create age groups (different from the one in the data set) and clean the data
colData(tse) <- colData(tse) %>%
  as.data.frame() %>%
  mutate(
    age_group = case_when(
      host_age_years >= 0 & host_age_years <= 1 ~ "Infant",
      host_age_years > 1 & host_age_years <= 3 ~ "Toddler",
      host_age_years > 3 & host_age_years <= 5 ~ "Child",
      host_age_years >= 18 & host_age_years <= 100 ~ "Adult",
      TRUE ~ NA_character_
    ),
    # Order the age categories for consistent plotting
    age_group = factor(age_group, levels = c("Infant", "Toddler", "Child", "Adult"))
  ) %>%
  mutate(region = factor(gsub(" ","_", geo_loc_name_country_continent_calc))) %>%
  DataFrame()

# Drop NA cases
selected <- !is.na(tse$region) & !is.na(tse$gender) & !is.na(tse$age_group) & !is.na(tse$GDP_per_head) & !is.na(tse$Usage)
tse <- tse[, selected]

# -----------------------------------------------------------------------------------------

# There are some samples with no prevalent ARGs; we must remove these from dbRDA calculation
# because dissimilarity calculation fails otherwise (state this in the ms)
tse.sel <- altExp(tse, "Prevalent", withColData=TRUE)
empty <- apply(assay(tse.sel, "relabundance")==0, 2, all) # These have no ARGs, problem for dissimilarity calculation
tse.sel <- tse.sel[, !empty]
#tse.sel <- tse.sel[, sample(ncol(tse.sel), 500)] # test first with smaller sample size
tse.sel <- tse.sel[, sample(ncol(tse.sel))] # test first with smaller sample size

# Center and scale the numeric variables
tse.sel$GDP_per_head <- scale(tse.sel$GDP_per_head)
tse.sel$Usage <- scale(tse.sel$Usage)

# Run RDA and store results into TreeSE
tse.sel <- addRDA(
    tse.sel,
    formula = assay ~ region + gender + age_group + GDP_per_head + Usage,
    FUN = getDissimilarity,
    assay.type="relabundance",
    distance = "bray",
    na.action = na.exclude
    )
# all(rownames(reducedDim(tse.sel, "RDA"))==colnames(tse.sel))

# Create RDA plot coloured by variable
p <- plotRDA(tse.sel, "RDA", colour.by = "gender") +
       scale_color_manual(values=c("red", "blue"))
print(p)

# Publications usually request jpg
jpeg("../RESULTS/FIGURES/dbrda.jpg", width=600, height=480, quality=100)
print(p)
dev.off()

# Store results of PERMANOVA test
rda_info <- attr(reducedDim(tse.sel, "RDA"), "significance")
tab <- rda_info$permanova |> knitr::kable()
write.table(tab, file="../RESULTS/FIGURES/dbrda_table.tab", row.names=FALSE, quote=FALSE)

jpeg("../RESULTS/FIGURES/dbrda.jpg", width=600, height=480, quality=100)
print(p)
dev.off()

p2 <- plotLoadings(tse.sel, "RDA", ncomponents = 2, n = 20)

jpeg("../RESULTS/FIGURES/dbrda_with_loadings.jpg", width=600, height=1000, quality=100)
print(p/p2)
dev.off()


