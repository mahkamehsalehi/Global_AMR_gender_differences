# Code for Global AMR project

# Directory structure

# AMR Gender Analysis Repository Structure

AMR-Gender-Analysis/  
├── **DATA/**                         # Store raw and processed data (avoid storing sensitive or large data)  
│   ├── **RAW/**                      # Raw data files (e.g., TSE Data, mapping data, SRA metadata)  
│   └── **PROCESSED/**                # Processed data for analysis (e.g., filtered datasets)  
│  
├── **R-SCRIPTS/**                    # R scripts for statistical analysis and data visualization  
│   ├── **carpentry.R**               # Script for creating TSE object  
│   ├── **01_data_preprocessing.R**   # Script for data cleaning and filtering  
│   ├── **02_ARG_load_analysis.R**    # Script for ARG load analysis (box plots, scatter plots)  
│   ├── **03_modeling.R**             # Script for modeling ARG load with gender and age  
│   └── **04_beta_diversity.R**       # Script for PCoA analysis and beta-diversity calculations  
│  
├── **BIOINFORMATIC-SCRIPTS/**        # Scripts related to bioinformatic analysis  
│   ├──  
│   └──  
│  
├── **NOTEBOOKS/**                    # Jupyter or Rmarkdown notebooks for exploratory analysis or reports  
│   └── **exploratory_analysis.Rmd**  # Exploratory analysis and data visualization  
│  
├── **OUTPUT/**                       # Store the generated results such as figures, tables, and models  
│   ├── **FIGURES/**                  # Figures such as box plots, scatter plots, and PCoA plots  
│   ├── **TABLES/**                   # Output tables (e.g., statistical results, p-values)  
│   └── **MODELS/**                   # Saved model objects (e.g., trained models for ARG load)  
│  
├── **DOCS/**                         # Documentation related to the project  
│   ├── **methodology.md**            # Description of methods and data processing steps  
│   ├── **analysis_plan.md**          # Detailed plan for analysis based on research questions  
│   └── **references.md**             # List of references and sources  
│  
├── **README.md**                     # Overview of the repository, analysis, and how to run the code  
├── **LICENSE**                       # License file for open-source usage  
├── **.gitignore**                    # List of files and directories to be ignored by git (e.g., large datasets)  
└── **requirements.txt**              # List of R package dependencies or `renv` configuration

# AMR Gender Paper

## Research Questions:
1. Is there a difference in ARG load variation between men and women across socio-economic groups?
2. Does the difference in ARG load between genders vary significantly across different age groups?
3. What are the differences in beta-diversity for genders?

## Data Sources:
- **Global AMR Data**: If the data is insufficient, supplement with Finrisk data. Due to data access issues and computational resource changes, Kata will handle Finrisk analyses.
- **Subset TSE Data**: Filter to include only samples with gender information.

## Analyses:
### ARG Load:
- **Box Plot Analysis**:  
  Plot ARG load by country/region/income group, colored by gender.  
  - **Y-axis**: ARG load  
  - **X-axis**: Country/Region/Income group  
  Include p-values for gender differences within regions and between regions.

- **Scatter Plot Analysis**:  
  Plot age (x-axis) vs. ARG load (y-axis), colored by gender.  
  Use `geom_smooth` in `ggplot2` for trend visualization.

- **Modeling**:  
  Model ARG load with gender and age included as variables, following the methodology outlined in Mahkameh's thesis.

### Beta-Diversity:
- **PCoA Analysis**:  
  Perform PCoA with Bray-Curtis distance on ARGs (without grouping by gene or class).  
  Color points by gender.

- **Differential Abundance Analysis**:  
  Group ARGs by gene.  
  Compare abundance between genders while controlling for age and region.  
  Use facet grid in `ggplot` for visualization.  
  Apply non-parametric tests (e.g., Wilcoxon test).

## References:
- **Differential Abundance Methods for ARGs**:  
  Jonsson, V., Österlund, T., Nerman, O. et al. (2016). Statistical evaluation of methods for identification of differentially abundant genes in comparative metagenomics. *BMC Genomics*, 17, 78. [Link](https://doi.org/10.1186/s12864-016-2426-8)
  
- **Socio-Demographic Stratification in AMR Research**:  
  World Health Organization (2018). *Tackling antimicrobial resistance (AMR) together: Working paper 5.0: Enhancing the focus on gender and equity*. [Link](https://www.who.int)

## Important Dates:
- **Npj Biofilms and Microbiomes**:  
  Women and Their Microbes Collection Special Issue Submission Deadline: *December 6, 2024*

- **Draft Review**:  
  Get draft to be reviewed by John Beggs, Peter Collignon, and Johan Bengtsson-Palme by *November 6*.

