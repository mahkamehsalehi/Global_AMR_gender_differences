# README: ARG Analysis Pipeline

This repository contains scripts for analyzing antibiotic resistance gene (ARG) data, including preprocessing, visualization, statistical modeling, and diversity analysis. The workflow includes data preparation, exploratory analysis, and statistical modeling to investigate ARG load across demographic and socio-economic factors.

---

## Pipeline Overview

### 1. Data Preprocessing

- **`00_functions.R`** → Contains utility functions for data cleaning.
- **`01_data_preprocessing.R`** → Processes raw data, filters metadata, normalizes gene expression values, and creates the **TreeSummarizedExperiment (TSE)** object.
- **`metadata_processing_sex.R`** → Extracts **gender, and age information** from metadata and standardizes these fields.
- **`metadata_processing_antibiotic_use.R`** → Merges antibiotic usage data (2016 estimates) with the TSE metadata.
- **`metadata_processing_bioproject.R`** → Adds bioproject titles and abstracts to the TSE metadata.
- **`metadata_processing_country.R`** → Integrates **country-level metadata** and categorizes **countries by region**.
- **`metadata_processing_JSON.R`** → Parses **JSON metadata fields** to extract structured information.
- **`metadata_processing_trip_distribution.R`** → Loads and processes a **trip distribution matrix** for country-level travel data.

---

### 2. Exploratory Analysis

- **`02_ARG_load_analysis.R`** → Visualizes **ARG load** across gender, age, and socio-economic groups.
  - Scatterplots of **log10(ARG load)** vs. **age** (colored by gender).
  - Linear and LOESS regression fits for ARG load by sex.
  - ANOVA and GLM analysis for category-based effects on ARG load.

- **`04_beta_diversity.R`** → Performs **PCoA analysis** (Principal Coordinate Analysis) to visualize **beta diversity** in ARG composition.
---

### 3. Statistical Modeling

- **`03_modeling.R`** → (Currently empty, reserved for future modeling scripts).
- **`dummy_lm.R`** → Performs linear modeling for **Figure 5**, generates **Figure S8**, and creates tables **S4, S5, S12, S13**.
- **`dummy_lm_top5_AB.R`** → Compares **top 5 antibiotic usage patterns** and generates **Figure S7**.
- **`tukey_women_age_region.R`** → Performs **Tukey's tests** to compare ARG load across age, gender, and region (outputs **Tables S6–S9**).
- **`permanova.R`** → Conducts **PERMANOVA** tests to assess group-level differences in ARG diversity.

---

| **Figure**  | **Table**  | **Script** |
|------------|------------|------------|
| **Figure 1** | - | `Figure_1_DataSummary.R` |
| **Figure 2** | - | `Figure_2_PCoA_top5_ARG.R` |
| **Figure 3** | - | `Figure_3_income_gender.R` |
| **Figure 4** | - | `Figure_4_hic_region_age.R` |
| **Figure S1** | - | `Figure_S1_S2_Usage.R` |
| **Figure S2** | - | `Figure_S1_S2_Usage.R` |
| **Figure S3** | - | `Figure_S3_income_women.R` |
| **Figure S4** | - | `Figure_S4_ARG_load.R` |
| **Figure S5** | - | `Figure_S5_ARG_diversity.R` |
| **Figure S6** | - | `Figure_S6_PCoA.R` |
| **Figure S7** | - | `dummy_lm_top5_AB.R` |
| **Figure S8** | - | `dummy_lm.R` |
| **Figure S9** | - | `Figure_subgroup_histograms.R` |
| - | **Table S1** | `PERMANOVA.R` |
| - | **Table S2** | `Table_S2_S3.R` |
| - | **Table S3** | `Table_S2_S3.R` |
| - | **Table S4** | `dummy_lm.R` |
| - | **Table S5** | `dummy_lm.R` |
| - | **Table S6** | `tukey_women_age_region.R` |
| - | **Table S7** | `tukey_women_age_region.R` |
| - | **Table S8** | `tukey_women_age_region.R` |
| - | **Table S9** | `tukey_women_age_region.R` |
| - | **Table S12** | `dummy_lm.R` |
| - | **Table S13** | `dummy_lm.R` |
---

## Usage Instructions

1. **Run data preprocessing**
   - Start with `01_data_preprocessing.R` to create the **TSE object**.
     - `metadata_processing_antibiotic_use.R`
     - `metadata_processing_bioproject.R`
     - `metadata_processing_country.R`
     - `metadata_processing_sex.R`
     - `metadata_processing_trip_distribution.R`

2. **Perform exploratory analysis**
   - Run `02_ARG_load_analysis.R` to visualize **ARG load trends**.
   - Run `04_beta_diversity.R` for **beta diversity analysis**.

3. **Perform modeling and statistical tests**
   - Use `dummy_lm.R` and `dummy_lm_top5_AB.R` for **linear modeling**.
   - Run `tukey_women_age_region.R` for **post-hoc comparisons**.
   - Run `permanova.R` to test **group-level diversity differences**.

---

