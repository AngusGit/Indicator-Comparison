source('functions.R')
source('CI_30_models_functions.R')
source('functions_inla.R')
## singularity exec -B .:/home/Project ../coral_indicators.sif Rscript 00_main.R --final_year=2023 --fresh_start=true --runStage=1:5 --rerun_baselines=true

PRIMARY_DATA_PATH="../data/primary/"  #AT
DATA_PATH="../data/"   ## AT
FINAL_YEAR=2024
OUTPUT_PATH = "../output/"
FIGS_PATH = paste0(OUTPUT_PATH, "figures")
TABS_PATH = paste0(OUTPUT_PATH, "tables")

# bring in data for reef level estimates 

 source('CI_10_get_data.R')
 source('CI_20_process_data.R') # processes the data for reef level models and brings in baseline estimates from the metadata record.

#
source('CI_31_models_MA.R')

source('CI_32_models_CC.R')

# note: there is a function CI__clean_ju_data that removes years with all zeros fo juvenile counts 
# to ensure reef level models run. this had the effect of producing NA scores in some years when no Acropora juvs were observed.
# I edited to retain individual years with zeros but still removed reefs where only one year was not zeros
 source('CI_33_models_JU.R')  

 source('CI_34_models_RPI.R')

 source('CI_35_models_CO.R')

## 
source('CI_40_collation.R')

## Take the indicator level scores and compare to Coral Index indicator scores - produce plots for report
source('Compare Grades.R')

## Compare to Coral Index scores to Framework Categories - produce plots for report

source('Compare Grades Index.R')
