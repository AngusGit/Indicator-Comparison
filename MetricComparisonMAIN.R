
# framework data processing
# code relies heavily on the GittHub produced for the Indicators Framework project "https://github.com/open-AIMS/CoralReefHabitat_Indicators" which has been cloned for this project
# into the folder C:\Users\athompso\OneDrive - Australian Institute of Marine Science\MMP\Reports\Visit 21 2025\Index comparison\CoralReefHabitat_Indicators
# as the various codes and functions are used/editted these will be shifted to C:\Users\athompso\OneDrive - Australian Institute of Marine Science\MMP\Reports\Visit 21 2025\Index comparison\R
# Baseline posteriors for HC, MA and Juv indicators e.g."CC_baseline_posteriors.RData" were downloaded from the GBR_DMS
# reef/depth level posteriors for hard coral cover were sourced from the LTMP dashboard preview http://tsv-ltmp-proc.aims.gov.au:3838/dashboard/

### Miscellaneous notes.... estimate distance from baseline in CI_32_models_CC_functions.R use CI_models_CC_distance(), this calls CI_index_CC()
# to aggregate to region CI_32_models_CC.R which call various functions to aggregate to desired level
# compare scores, compare medians, or compare to entire distribution simple subtraction of Coral cover index score 
#To deal with distance to baseline cf. consequence scores, consider worst case perhaps??

library(tidyverse)

################
#Hard Coral
################

#getting annual posteriors from the dashboard -  http://tsv-ltmp-proc.aims.gov.au:3838/dashboard/
# this required downloading the individual files and reading them in as below

# note the use of strsplit() to extract the components of the file names is enalbled by Murray's clever naming protocol on the dashboard that builds names in sensible manner.
hc.post.files<-list.files(path = "C:/Users/athompso/Downloads", pattern="_HARD CORAL_.*\\posteriors.csv$")
hc.posteriors<-do.call(rbind, lapply(paste0("C:/Users/athompso/Downloads/",hc.post.files), 
                                     function(x){
                                       s<-strsplit(basename(x), "_")[[1]]
                                       dat<-read.csv(x) |> 
                                         mutate(Reef=s[3],
                                                Depth=s[6],
                                                Zone=s[5])
                                       dat
                                     }
                                     ))
#######
# load the baseline distribution avaialble from the GBR-DMS
#
# this following citation should be used
# Australian Institute of Marine Science (AIMS). (2023). 
# Indicator Framework for assessing the condition of coral reef habitats in the Great Barrier Reef. 
# https://apps.aims.gov.au/metadata/view/55f56f2c-c9fd-4af7-90f9-617643975f9e, accessed 23/01/2025.



 # load("C:/Users/athompso/Downloads/CC_baseline_posteriors.RData") # note this object is called "mod" posteriors are at the level of DEPTH.f ("shallow" or "deep" and BIOREGION.agg)
# hc.ref.baseline<-mod

load("C:/Users/athompso/OneDrive - Australian Institute of Marine Science/MMP/Reports/Visit 21 2025/Index comparison/data/CC_baseline_posteriors.RData")
hc.ref.baseline<-mod


######
# This chunk pulls the indicator scores from the DMS however, the juvenile values are off
library(arrow)
library(dplyr)
library(leaflet)
library(ggplot2)
library(kableExtra)
# 
bucket <- s3_bucket("s3://gbr-dms-data-public/aims-reef-indicators-framework/data.parquet")
df <- open_dataset(bucket)
print(df$schema)
dfREEFS <- df |> filter(Shelf=="Inshore") |> collect()
juv.total<-df |> filter(Indicator=="Juvenile.density" & Metric=="Total") |> collect()
juv.reef<-juv.total |> filter(Aggregation=="reef")


##########
# as the above seems dubious try from Murrays code.
#########
# Extract raw data from  Inshore Reefs
