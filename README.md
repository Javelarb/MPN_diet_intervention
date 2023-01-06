#Mediterranean diet intervention - MPN
by: Julio Avelar-Barragan

This repository contains all of the R and shell code used for the analysis of the manuscript:  
"Characterizing the microbiome of patients with myeloproliferative neoplasms during a Mediterranean diet intervention"  

Files and their description:  
* ancom.R - A script used for running the differential abundance test ANCOM v2.1 (https://www.tandfonline.com/doi/full/10.3402/mehd.v26.27663)
* Bash_commands.txt - Contains the parameters for functions and programs within the bash command line.  
* contig_table_maker.sh - A script used to generate a table of open reading frames x sample abundances.  
* Functional_MPN_analysis.Rmd - R markdown file containing code used to analyze the functional metagenomic data of samples.  
* megahit.sh - contains the parameters for running the program megahit.sh, which is used to assemble shotgun data.  
* metadata_github.csv - contains the sample metadata used for analysis.  
* MPN_analysis.Rmd - contains the R code used to analyze the OTU table.  
* OTU_clean.csv - A filtered OTU table containing the relative abundances of species per sample as characterized by metaphlan3.  