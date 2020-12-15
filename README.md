# ADHD
Included are the scripts for the mixed effects modeling (in the R script) and elastic net modeling analyses (in the matlab scripts) from the article 
Multimethod Investigation of the Neurobiological Basis of ADHD Symptomatology in Children Aged 9-10: Baseline Data from the ABCD Study.

There are also some auxilary analyses in the R script following the mixed effects modeling.

A postprint of this article is available at https://osf.io/425wu/.

These scripts use data from the ABCD RDS file. The RDS file can be downloaded with permission from NDA at https://nda.nih.gov/abcd.

Scripts should be run in the same directory as the RDS file. Script run order is 1) "ADHD_Analyses_Final.R" 2) "preanalysis_driver_FINAL.m"
3) "final_accumulator_metadriver_FINAL.m".

There are a few files with participant names that are used in the data cleaning phase of the R-script that I cannot share for participant privacy reasons (e.g., "philips.txt).

Matlab scripts are currently set to run elastic net regression modeling using the standard covariate approach (i.e., all covariates
but medication status). To run with no covariates or using medication status as a covariate, change the path directories in the "preanalysis driver"
script.
