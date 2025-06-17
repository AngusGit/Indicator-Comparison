source('functions.R')
source('CI_35_models_CO_functions.R')
source('CI_30_models_functions.R')

if (CI_isParent()) CI_startMatter()

CURRENT_STAGE <<- 7
CI__change_setting(item = 'CURRENT_STAGE', value = CURRENT_STAGE)
CI__add_stage(stage = paste0('STAGE',CURRENT_STAGE),
              title = paste0('Stage ', CURRENT_STAGE, ' - modelling Composition data'))

if (CI$setting[['FRESH_START']]) CI_clear_models_CO_data()

## Fit models
CI_models_CO_prepare_data()
CI_models_CO_calc_mds()
CI_models_CO_calc_novelty()
CI_models_CO_calc_deviation()
CI_models_CO_combine_metrics()
CI_models_CO_calc_distance()

## Standardise the output and aggregate
CI_model_CI_standardise()
# CI_models_aggregation(Indicator = 'CO', level = 'BIOREGION.agg') 
CI_models_CO_aggregation(level = 'NRM') 
# CI_models_aggregation(Indicator = 'CO', level = 'TUMRA') 
# CI_models_aggregation(Indicator = 'CO', level = 'GBRMPA.MA') 
# CI_models_aggregation(Indicator = 'CO', level = 'GBRMP') 
# CI_models_aggregation(Indicator = 'CO', level = 'ZONE') 

## Visualisations
## CI_models_CO_trajectory()
