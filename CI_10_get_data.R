source('../R/functions.R')
source('../R/CI_10_get_data_functions.R')


# if (CI_isParent()) CI_startMatter()

#######################
## Reef names lookup 
## AIMS_REEF_MANE and REEF_ZONE required to link to new disturbance table DISTURBANCES2
## OUTSTANDING is full mapping to potential reporting zones...
## This should most probably be done via intersection of sample waypoints with spatial layers.
#######################

CURRENT_STAGE <<- 1
CI__add_setting(item = 'CURRENT_STAGE', value = 1)
CI__add_stage(stage = paste0('STAGE',CURRENT_STAGE), title = 'Stage 1 - extract data')

if (CI$setting[['FRESH_START']])
  
  CI_clear_primary_data()

## AIMS data in database
CI_get_names_lookups()
CI_get_mmp_samples_data()
CI_get_ltmp_samples_data()
CI_combine_samples_data()
CI_get_disturbances_data()
CI_get_video_code_lookup()
CI_get_mmp_juveniles_data()
CI_get_ltmp_juveniles_data()
CI_combine_juveniles_data()
CI_get_mmp_point_data()
CI_get_ltmp_point_data()
CI_combine_point_data()

## External sources
CI_get_external_data()

## Spatial data (arcgis server)
CI_get_spatial_data()
