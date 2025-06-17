source('functions.R')
source('CI_31_models_MA_functions.R')
source('CI_30_models_functions.R')
#source('functions_inla.R')

## Generate baseline models
#source("CI_25_MA_baseline_models.R")
## baselines loaded from metadata record

## Combine baseline models
CI_models_MA_get_baselines()

## Fit models
CI_models_MA_prepare_data()
CI_models_MA_prepare_nest() 

CI_models_MA_fit_models() 
CI_models_MA_diagnostics()
CI_models_MA_cellmeans()
#CI_models_cellmeans(Indicator = "MA")
CI_models_MA_preds()
#CI_models_MA_partialplots()

## Calculate distances to baselines
CI_models_MA_distance()
## CI_models_MA_varify_scores()
#CI_models_aggregation(Indicator = "MA", level = 'BIOREGION.agg') 
CI_models_aggregation(Indicator = "MA", level = 'NRM') 
# CI_models_aggregation(Indicator = "MA", level = 'TUMRA') 
# CI_models_aggregation(Indicator = "MA", level = 'GBRMPA.MA') 
# CI_models_aggregation(Indicator = "MA", level = 'GBRMP') 
#CI_models_aggregation(Indicator = "MA", level = 'ZONE') 
CI_models_aggregation(Indicator = "MA", level = 'REEF.d') 

