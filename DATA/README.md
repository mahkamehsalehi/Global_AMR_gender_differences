### This README file describes the different data sources in DATA
Subdirectories PROCESSED and RAW and their contents are in .gitignore. These contain large output files.

### Sra_metadata_jun11.txt
Sra metadata fetched based on fastq file accessions which were downloaded for the project.  
Created using bigquery_metadata_query.sh

### Classified_Sra_meta_data_jun11.tct
Sra metadata with classfications of projects based on bioproject.
Created using fecth_abstracts.py.

### gene_length.csv
Resfinder ARG gene lenghts in basepairs. Calculated using bash script.

### resfinder_phenotypes.txt
Resinfinder phenotypes file downloaded with the ResFinder database v. 2.1.1

### data.csv
Country level data with antibiotic resistance, infrasctructure, antibiotic use, socio-economics, climate etc from Anthropological and socioeconomic factors contributing to global antimicrobial resistance: a univariate and multivariable analysis.
Collignon, Peter et al.
The Lancet Planetary Health, Volume 2, Issue 9, e398 - e405
https://doi.org/10.1016/S2542-5196(18)30186-4

### total_antibiotic_consumption_estimates.csv
Estimated antibiotuc consumption data from https://www.tropicalmedicine.ox.ac.uk/gram/research/visualisation-app-antibiotic-usage-and-consumption 
Using year 2016.
Global antibiotic consumption and usage in humans, 2000–18: a spatial modelling study.
Browne, Annie et al.
The Lancet Planetary Health Volume 5, Issue 12, December 2021, Pages e893-e904

### Trip_distribution_matrix.

Trip distribution data from Peter Collignon and John Beggs. See First cut at AMR burden of Travel.docs
