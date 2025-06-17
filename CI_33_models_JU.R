source('functions.R')
source('CI_33_models_JU_functions.R')
source('CI_30_models_functions.R')


## Generate baseline models
#source("CI_27_JU_baseline_models.R") # produces parameters/JUV_baseline.RData this sourced from metadata

## Combine baseline models
CI_models_JU_get_baselines() # reads in JUV_baseline, outputs modelled/JU__baseline_posteriors.RData

## Fit models
CI_models_JU_prepare_data() # input processed/juv.df.RData and  processed/spatial_lookup.RData, output modelled/data_ju.RData
CI_models_JU_prepare_nest() # input modelled/data_ju.RData, output "mods" as modelled/JU__mods.RData'
CI_models_JU_fit_models() # fits the reef level models input, modelled/data_ju.RData and  modelled/JU__mods.RData, output (via internal use of 
## CI__fit_JU_model, modelled/JU__", reef, '_', taxa, '__draws.RData, and **__modelled.RData for each reef and Taxa
CI_models_JU_cellmeans() # input modelled/data_ju.RData") and 'modelled/JU__mods.RData, and via CI__cellmeans_JU_model the reeflevel draws and out put reeflevel posteriors 
CI_models_JU_preds() #input modelled/data_ju.RData") and 'modelled/JU__mods.RData and reef level posteriors, output modelled/JU*REEFLEVEL_summary.RData
# and object mods as modelled/JU__preds.RData
## Calculate distances to baselines
CI_models_JU_distance()
#CI_models_aggregation(Indicator = 'JU', level = 'BIOREGION.agg') 
CI_models_JU_aggregation(level = 'NRM') 
# CI_models_aggregation(Indicator = 'JU', level = 'TUMRA') 
# CI_models_aggregation(Indicator = 'JU', level = 'GBRMPA.MA') 
# CI_models_aggregation(Indicator = 'JU', level = 'GBRMP') 
# CI_models_aggregation(Indicator = 'JU', level = 'ZONE') 
CI_models_JU_aggregation(level = 'REEF.d') 
