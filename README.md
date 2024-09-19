# Code for Global AMR project

# TODO:
Clean up the codes to make TSE object including gender and age and country level socio-economic and antibiotic data. Main script Carpentry.R that collects smaller scripts.

Make a clean structure, following this

"DATA" -> Big files, bioinformatic results > 500 MB
"R" -> R-sripts, smaller data files (undersub dir), 
"RESULTS" figs and tables, saved as csv of pdf
"SCRIPTS" Bioinformatic scripts
"Archived" all old code not used



Title: AMR Gender Paper
Research Questions:
Is there a difference in ARG load variation between men and women across socio-economic groups?
Does the difference in ARG load between genders vary significantly across different age groups?
What are the differences in beta-diversity for genders?
Data Sources:
Global AMR Data: If the data is insufficient, supplement with Finrisk data. Due to data access issues and computational resource changes, Kata will handle Finrisk analyses.
Subset TSE Data: Filter to include only samples with gender information.
Analyses:
ARG Load:
Box Plot Analysis:
Plot ARG load by country/region/income group, colored by gender.
Y-axis: ARG load
X-axis: Country/Region/Income group
Include p-values for gender differences within regions and between regions.
Scatter Plot Analysis:
Plot age (x-axis) vs. ARG load (y-axis), colored by gender.
Use geom_smooth in ggplot2 for trend visualization.
Modeling:
Model ARG load with gender and age included as variables, following the methodology outlined in Mahkameh's thesis.
Beta-Diversity:
PcoA Analysis:
Perform PCoA with Bray-Curtis distance on ARGs (without grouping by gene or class).
Color points by gender.
Differential Abundance Analysis:
Group ARGs by gene.
Compare abundance between genders while controlling for age and region.
Use facet grid in ggplot for visualization.
Apply non-parametric tests (e.g., Wilcoxon test).
References:
Differential Abundance Methods for ARGs:
Jonsson, V., Österlund, T., Nerman, O. et al. (2016). Statistical evaluation of methods for identification of differentially abundant genes in comparative metagenomics. BMC Genomics, 17, 78. Link
Socio-Demographic Stratification in AMR Research:
World Health Organization (2018). Tackling antimicrobial resistance (AMR) together: Working paper 5.0: Enhancing the focus on gender and equity. Link
Important Dates:
Npj Biofilms and Microbiomes:
Women and Their Microbes Collection Special Issue Submission Deadline: December 6, 2024
Draft Review:
Get draft to be reviewed by John Beggs, Peter Collignon, and Johan Bengtsson-Palme by November 6.


