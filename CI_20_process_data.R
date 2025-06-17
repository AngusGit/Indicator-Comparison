source('functions.R')
source('CI_20_process_data_functions.R')

## AIMS data
CI_process_add_report_year()
CI_process_transect_area()
CI_process_locations_data()
CI_process_total_points()
CI_process_hc_points()
CI_process_ma_points()
CI_process_a_points()
CI_process_combine_points()

## Obtain (pre-wrangled) external data
CI_process_external_data()

## Combine AIMS and external data
CI_process_combine_data()

## Assign regions
CI_process_spatial()
CI_process_assign_regions()

## Juvenile density
CI_process_juvenile_density()
CI_process_juvenile_offset()
CI_process_juvenile_data_for_baselines()

## Unique site locations
CI_process_unique_site_location()

## RPI (recovery performance indicator) data
## this processing has been moved to 34_RPI_models_RPI()
## since it has the CC models as a dependency

## Composition data
CI_process_composition_points()
CI_process_composition_ordi_data()
CI_process_composition_data()

## Baseline models, Note we cannot fit the RPI baseline models at this
## stage because they input data from the fitted CC models
if (CI$setting[['RERUN_BASELINES']]) {
    source('../R/CI_25_MA_baseline_models.R')
    CI_25_MA_baseline_models()
    source('../R/CI_26_CC_baseline_models.R')
    CI_26_CC_baseline_models()
    source('../R/CI_27_JU_baseline_models.R')
    CI_27_JU_baseline_models()
    source('../R/CI_28_JU_IPM_models.R')
    CI_28_JU_IPM_models()
    source('../R/CI_29_JU_consequence_models.R')
    CI_29_JU_consequence_models()
}

## Check if the baseline models are present.  If they are not, then
## retrieve them from the AIMS metadata record
CI__retrieve_data_from_metadata("MA__baseline_deep.RData")
CI__retrieve_data_from_metadata("MA__baseline_shallow.RData")
CI__retrieve_data_from_metadata("MA__baseline_offshore.RData")

CI__retrieve_data_from_metadata("CC_baseline.RData")
CI__retrieve_data_from_metadata("CC__baseline.RData")
CI__retrieve_data_from_metadata("CC_baseline_posteriors.RData")

CI__retrieve_data_from_metadata("IPM_juv.RData")
CI__retrieve_data_from_metadata("JUV_baseline.RData")
CI__retrieve_data_from_metadata("JU__baseline.RData")
CI__retrieve_data_from_metadata("JU__baseline_mod.RData")
# CI__retrieve_data_from_metadata("JuvCover_Acropora_NA_Inshore_NA_ResJuv_FALSE_Morph_NA.Rdata")
# CI__retrieve_data_from_metadata("JuvCover_Acropora_NA_Offshore_NA_ResJuv_FALSE_Morph_NA.Rdata")
# CI__retrieve_data_from_metadata("JuvCover_Mound_growth_NA_Inshore_NA_ResJuv_FALSE_Morph_Merulinidae.Rdata")
# CI__retrieve_data_from_metadata("JuvCover_Mound_growth_NA_Inshore_NA_ResJuv_FALSE_Morph_Poritidae.Rdata")
# CI__retrieve_data_from_metadata("JuvCover_Mound_growth_NA_Offshore_NA_ResJuv_FALSE_Morph_Merulinidae.Rdata")
# CI__retrieve_data_from_metadata("JuvCover_Mound_growth_NA_Offshore_NA_ResJuv_FALSE_Morph_Poritidae.Rdata")
# CI__retrieve_data_from_metadata("JuvCover_Pocilloporidae_NA_Inshore_NA_ResJuv_FALSE_Morph_NA.Rdata")
# CI__retrieve_data_from_metadata("JuvCover_Pocilloporidae_NA_Offshore_NA_ResJuv_FALSE_Morph_NA.Rdata")

CI__retrieve_data_from_metadata("newdata_grid.RData")
CI__retrieve_data_from_metadata("Synoptic_text_reef.csv")
CI__retrieve_data_from_metadata("Synoptic_text_region.csv")
CI__retrieve_data_from_metadata("taxaLookup.csv")

