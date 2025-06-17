source('functions.R')
source('CI_32_models_CC_functions.R')
source('CI_30_models_functions.R')


## Generate baseline models
##source("CI_26_CC_baseline_models.R")

#Baselines loaded from the meta data recored using the script CL_20_process_data.R
# CI__retrieve_data_from_metadata("CC_baseline.RData")
# CI__retrieve_data_from_metadata("CC__baseline.RData")
# CI__retrieve_data_from_metadata("CC_baseline_posteriors.RData")

## Combine baseline models
CI_models_CC_get_baselines() # outputs data/modelled/cc__baseline_posteriors.RData

## Fit models
CI_models_CC_prepare_data() #input "processed/points.analysis.data.transect.RData" and "processed/spatial_lookup.RData"  , output "modelled/data_cc.RData"
CI_models_CC_prepare_nest() # input modelled/data_cc.RData, output
CI_models_CC_fit_models()
CI_models_CC_cellmeans()
CI_models_CC_preds()

## Calculate distances to baselines
CI_models_CC_distance()

## AT added for comparison report
baselines <- get(load(file = paste0(DATA_PATH,
                                    'modelled/CC__baseline_posteriors.RData')))
cc.bio.base<-baselines |> 
  group_by(DEPTH.f, BIOREGION.agg) |> 
  summarise(value=median(value))

write.csv(cc.bio.base, file=paste0(OUTPUT_PATH,
                                    'tables/CC_bioregion_baseline_medians.csv'))
## Note the below don't seem to work due to no common variables in the last left_join.

# CI_models_aggregation(Indicator = 'CC', level = 'BIOREGION.agg') 
# CI_models_aggregation(Indicator = 'CC', level = 'NRM') 
# CI_models_aggregation(Indicator = 'CC', level = 'TUMRA') 
# CI_models_aggregation(Indicator = 'CC', level = 'GBRMPA.MA') 
# CI_models_aggregation(Indicator = 'CC', level = 'GBRMP') 
CI_models_CC_aggregation(level = 'NRM') 
