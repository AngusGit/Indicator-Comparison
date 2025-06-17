
CI_process_rpi_configs <- function() {
  
  
  K.limit     <<- 90
  max_init    <<- 100.0 # thresh hold for initial condition (in hard coral cover)
  min_final   <<- 0.0 # minimum threshold for final observation before next disturbance
  min_obs     <<- 4.0 # minimum number of observations
  if (RPI_PURPOSE == "reference") min_obs <<- 2.0
  if (RPI_PURPOSE == "critical") min_obs <<- 1.0
  #Start and end dates for baselines
  start_date  <<- lubridate::ymd("1990-01-01") 
  end_date    <<- lubridate::ymd("2022-01-01")
  
  DATA_DIR                  <<- paste0(DATA_PATH, "primary/")
  MODEL_DIR                 <<- paste0("../R/externalFunctions/")
  DIST_DAT_FILE             <<- "disturbance.RData"
  SAMPLE_DAT_FILE           <<- "samples.RData"
  GRP_LEV_TRANS_DAT_FILE    <<- "groups.transect.RData"
  GROUP_NAMES               <<- "spatial_lookup.RData"
  
  if (!dir.exists(paste0(DATA_PATH, "processed/RPI_baseline/")))
    dir.create(paste0(DATA_PATH, "processed/RPI_baseline/"))
  if (!dir.exists(paste0(DATA_PATH, "processed/RPI_reference/")))
    dir.create(paste0(DATA_PATH, "processed/RPI_reference/"))
  if (!dir.exists(paste0(DATA_PATH, "processed/RPI_critical/")))
    dir.create(paste0(DATA_PATH, "processed/RPI_critical/"))
  PROC_DATA_DIR             <<- paste0(DATA_PATH, "processed/")
  ## MCMC_OUTPUT_DIR           <<- paste0(DATA_PATH, "processed/RPI_baseline/")
  MCMC_OUTPUT_DIR           <<- paste0(DATA_PATH, "processed/RPI_",RPI_PURPOSE,"/")
  
  FILTER_OUT_DATA_FMT         <<- "rec.traj.trans.init%f.final%f.obs%d"
  FILTER_REFMT_OUT_DATA_FMT   <<- "rec.traj.trans.proc.init%f.final%f.obs%d"
  
}

#########################################
# Get data functions
#########################################

## run first
CI_process_rpi_data_spatial <- function() {
  
  load(paste0(DATA_PATH, 'processed/spatial_lookup.RData'))
  
  ## There are two levels of DEPTH for REEF.d, I think they are
  ##JCU sites. Will take the first level, just need a unique
  ##REEF.d There are two levels of TUMRA for two of the REEF.ds,
  ##I've just selected 1
  
  spatial_lookup<- spatial_lookup %>%
    group_by(REEF.d) %>%
    filter(DEPTH==first(DEPTH),
           !(REEF.d=="Pandora Reef deep slope" & is.na(TUMRA)) &   ## AT corrected the filter here
             !(REEF.d=="Farquharson Reef deep slope" & TUMRA=="Mandubarra")) %>% droplevels() %>%
    dplyr::select(-DEPTH )%>%
    suppressMessages() %>%
    suppressWarnings()
  
  
  save(spatial_lookup, file=paste0(DATA_PATH, "/processed/spatial_lookup_rpi.RData"))
  
}

CI_process_rpi_get_cc_data <- function() {
  
  
  load(file=paste0(DATA_PATH, "/processed/spatial_lookup_rpi.RData"))
  
  load_cc_posteriors <- function(r) { 
    if(!file.exists(file = paste0(DATA_PATH, "modelled/CC__", r, "__posteriors.RData"))) {
      return(NULL)
    } else{
      load(file=paste0(DATA_PATH, "modelled/CC__", r, "__posteriors.RData"))
      return(posteriors)
    }
  }
  
  reef.depth.names <- spatial_lookup %>%
    dplyr::pull(REEF.d) %>%
    unique() 
  reef.posteriors.all <- reef.depth.names %>%
    map_df(.f = ~ load_cc_posteriors(.x)) %>%
    mutate(value = value*100)
  
  save(reef.posteriors.all, file=paste0(DATA_PATH, "/processed/reef.posteriors.all.RData"))
  
  
}

CI_process_groups_data <- function() {
  
  load(paste0(DATA_PATH, 'primary/points.raw.RData'))
  load(paste0(DATA_PATH, "primary/video_codes.RData"))
  
  ## AIMS data
  ##Total points transect level
  total.points <- points.raw %>%
    group_by(P_CODE, REEF,DEPTH, VISIT_NO, SITE_NO) %>%
    mutate(SITE_NO = as.character(SITE_NO)) %>% 
    summarise(total.points = sum(POINTS)) %>%
    ungroup() %>%
    suppressMessages() %>%
    suppressWarnings()
  
  groups.site <- points.raw %>%
    left_join(video_codes %>% dplyr::select(GROUP_CODE, VIDEO_CODE)) %>%
    mutate(SITE_NO = as.character(SITE_NO)) %>%
    group_by(P_CODE, REEF, DEPTH, VISIT_NO, GROUP_CODE, SITE_NO) %>%
    summarise(Points = sum(POINTS)) %>%
    ungroup %>% 
    left_join(total.points) %>%
    mutate(COVER = (Points/total.points)*100) %>%
    dplyr::select(-Points, -total.points)%>%
    pivot_wider(names_from = "GROUP_CODE", values_from = "COVER", values_fill = 0) %>%
    pivot_longer(cols = c('A', 'AB', 'HC', 'SC', 'SG', 'SP', 'OT'),
                 names_to = 'GROUP_CODE', values_to = 'COVER') %>%
    mutate(P_CODE = as.factor(P_CODE),
           REEF = as.factor(REEF),
           SITE_NO = as.integer(SITE_NO),
           GROUP_CODE = as.factor(GROUP_CODE)) %>%
    suppressMessages() %>%
    suppressWarnings()
  
  save(groups.site,
       file = paste0(DATA_PATH, 'processed/groups.site.RData'))
  
}

# need to be run internally to baseline, reference or critical sets as takes RPI_PURPOSE as an input
CI_process_rpi_get_samples_data <- function() {
  
  
  load(file=paste0(DATA_PATH, "/processed/sample.reef.report.year.RData"))
  
  samples <- sample.reef.report.year %>%
    mutate(DEPTH.f = factor(ifelse(DEPTH >= 3, "deep slope", "shallow slope")),
           REEF.d = factor(paste(REEF, DEPTH.f, sep = " "))) %>%
    dplyr::select(-DEPTH, -P_CODE) %>%
    group_by(REEF, DEPTH.f, REEF.d, VISIT_NO, REPORT_YEAR) %>%
    summarise(Latitude = min(LATITUDE),
              Longitude = min(LONGITUDE),
              Date = min(Date)) %>%
    mutate(fYEAR = factor(REPORT_YEAR)) %>%
    {if (RPI_PURPOSE == 'baseline')
      filter(., REPORT_YEAR <= 2020)
      else .
    } %>%
    ungroup() %>%
    suppressMessages() %>%
    suppressWarnings()
  
  if(RPI_PURPOSE == "baseline") {
    save(samples, file=paste0(DATA_PATH, "/processed/samples_rpi_", RPI_PURPOSE, ".RData"))
  }  else
    save(samples, file=paste0(DATA_PATH, "/processed/samples_rpi.RData"))
  
  
}

CI_process_groups_data_part2 <- function() {
  
  load(file = paste0(DATA_PATH, 'processed/groups.site.RData'))
  if (RPI_PURPOSE == "baseline") {
    load(file = paste0(DATA_PATH, "processed/samples_rpi_", RPI_PURPOSE, ".RData"))
  } else 
    load(file = paste0(DATA_PATH, "processed/samples_rpi.RData"))
  load(paste0(DATA_PATH, 'processed/spatial_lookup.RData'))
  
  get.abiotic <- groups.site %>%
    mutate(DEPTH.f = factor(ifelse(DEPTH >= 3, "deep slope", "shallow slope")),
           REEF.d = factor(paste(REEF, DEPTH.f, sep=" "))) %>%
    dplyr::select(-DEPTH) %>%
    filter(GROUP_CODE == "AB") %>%
    dplyr::select(-GROUP_CODE) %>%
    group_by(P_CODE,REEF.d, VISIT_NO) %>%
    summarise(AB = mean(COVER)) %>%
    left_join(samples %>%
                dplyr::select(REEF, DEPTH.f, REEF.d, VISIT_NO,
                              REPORT_YEAR, fYEAR)) %>%
    {if (RPI_PURPOSE == "baseline")
      filter(., !is.na(REPORT_YEAR) & REPORT_YEAR <= 2020)
      else filter(., !is.na(REPORT_YEAR))
    } %>%
    droplevels() %>%
    ungroup() %>%
    dplyr::select(-VISIT_NO, -P_CODE) %>%
    suppressMessages() %>%
    suppressWarnings()
  
  ##For HC, each posterior draw is treated as a 'transect' in the pipeline
  load(file = paste0(DATA_PATH, "processed//reef.posteriors.all.RData"))
  
  reef.posteriors<- reef.posteriors.all %>%
    rename(HC=value,
           TRANSECT_NO=`.draw`) %>%
    mutate(REPORT_YEAR=as.numeric(as.character(fYEAR))) %>%
    {if (RPI_PURPOSE == "baseline")
      filter(., REPORT_YEAR<=2020)
      else .
    } %>%
    droplevels() %>%
    suppressMessages() %>%
    suppressWarnings()
  
  ##Put Abiotic and HC cover together in one dataframe
  groups.transect <- reef.posteriors %>%
    left_join(get.abiotic) %>%
    filter(!is.na(AB)) %>%
    droplevels()  %>% ## locations with NA abiotic are mostly JCU sites
    pivot_longer(cols = c(AB, HC),
                 names_to = "GROUP_CODE", values_to = "COVER") %>% 
    ## gather(AB, HC, key = "GROUP_CODE", value = COVER)  %>%      
    left_join(spatial_lookup %>%
                dplyr::select(REEF.d, REEF, DEPTH.f) %>%
                distinct) %>%
    mutate(GROUP_CODE = factor(GROUP_CODE)) %>%
    suppressMessages() %>%
    suppressWarnings()
  
  ##This forms the primary data 'groups.transect'
  if (RPI_PURPOSE == "baseline") {
    save(groups.transect,
         file = paste0(DATA_PATH, "processed/groups.transect_", RPI_PURPOSE, ".RData"))
  } else
    save(groups.transect,
         file = paste0(DATA_PATH, 'processed/groups.transect.RData'))
  
  
  
  ##Load the disturbance table and save to the pipeline-specific primary directory
  
  ##Why are there two identical records for Farquharson Reef for each year??
  
  load(paste0(DATA_PATH, "primary/disturbances.RData"))
  
  disturbance <- disturbance.reef %>% 
    mutate(DISTURBANCE = factor(DISTURBANCE),
           DEPTH.f = factor(ifelse(DEPTH >= 3, "deep slope", "shallow slope")),
           REEF.d = factor(paste(REEF, DEPTH.f, sep = " "))) %>%
    dplyr::select(-DEPTH) %>%
    left_join(samples %>%
                dplyr::select(REEF, DEPTH.f, REEF.d, VISIT_NO,
                              fYEAR, REPORT_YEAR)) %>%
    dplyr::select(-P_CODE, -VISIT_NO) %>%
    filter(RANK %in% c(NA, 1)) %>%
    droplevels() %>%
    group_by(REEF, RANK, DISTURBANCE, STORM_NAME,
             DISTURBANCE_DATE, DEPTH.f, REEF.d, REPORT_YEAR, fYEAR) %>%
    unique() %>%
    {if (RPI_PURPOSE == "baseline")
      filter(., REPORT_YEAR <= 2020)
      else .
    } %>%
    droplevels() %>%
    ungroup() %>%
    suppressMessages() %>%
    suppressWarnings()
  
  if (RPI_PURPOSE == "baseline") {
    save(disturbance,
         file = paste0(DATA_PATH, "processed/disturbance_rpi_", RPI_PURPOSE, ".RData"))
  } else
    save(disturbance,
         file = paste0(DATA_PATH, "processed/disturbance_rpi.RData"))
  
}

CI_process_rpi_data <- function() {
  
  load(paste0(DATA_PATH, "processed/groups.transect_", RPI_PURPOSE, ".RData"))
  load(paste0(DATA_PATH, "processed/samples_rpi_", RPI_PURPOSE, ".RData"))
  load(paste0(DATA_PATH, "processed/disturbance_rpi_", RPI_PURPOSE, ".RData"))
  load(paste0(DATA_PATH, "processed/sample.reef.report.year.RData"))
  
  groups.transect.2<- groups.transect %>%
    mutate(REEF.d.year = paste(REEF.d, fYEAR, sep = " "))
  
  disturbance.2 <- disturbance %>%
    mutate(REEF.d.year = paste(REEF.d, fYEAR, sep = " ")) %>%
    filter(REEF.d.year %in% unique(groups.transect.2$REEF.d.year)) %>%
    droplevels()
  
  samples.2 <- samples %>%
    mutate(REEF.d.year = paste(REEF.d, fYEAR, sep = " ")) %>%
    filter(REEF.d.year %in% unique(groups.transect.2$REEF.d.year)) %>%
    droplevels()
  
  ## Use function to assign each reef.d/report_year to a recovery trajectory
  recovery.trajectories <- extract_recovery_trajectories_transect(disturbance.2,
                                                                  samples.2,
                                                                  groups.transect.2)
  
  save(recovery.trajectories,
       file = paste0(DATA_PATH, "processed/time.series.site.transect_",
                     RPI_PURPOSE, ".RData"))
}

CI_process_rpi_filter_trajectories <- function() {
  
  if (RPI_PURPOSE == "baseline") {
    load(file = paste0(DATA_PATH, "processed/time.series.site.transect_", RPI_PURPOSE, ".RData"))
    load(paste0(DATA_PATH, "processed/groups.transect_", RPI_PURPOSE, ".RData"))
  } else {
    load(file = paste0(DATA_PATH, "processed/time.series.site.transect.RData"))
    load(paste0(DATA_PATH, 'processed/groups.transect.RData'))
  }
  load(paste0(DATA_PATH, 'processed/spatial_lookup_rpi.RData'))
  
  ## Baseline trajectories need at least 4 observations
  ## Remove missing dates - the model needs dates
  rm.ongoing.trajectories <- recovery.trajectories %>% 
    mutate(proj.site.rpid = paste(REEF.d, RP_ID, sep = " ")) %>%
    ungroup() %>%
    group_by(proj.site.rpid) %>% 
    mutate(max.report.year = max(REPORT_YEAR),
           NUM_OBS = as.numeric(length(unique(REPORT_YEAR)))) %>%
    filter(!is.na(Date), NUM_OBS > 3) %>%
    droplevels() %>%  ##Filter to trajectories that have a least 4 observations
    ungroup() %>%
    suppressMessages() %>%
    suppressWarnings()
  
  ## This was to prevent missing data, which causes issues further down the pipeline
  groups.transect.traj.reefds <- groups.transect %>%
    filter(REEF.d %in% unique(rm.ongoing.trajectories$REEF.d)) %>%
    droplevels()
  
  filt.rec.traj <- rm.ongoing.trajectories %>%
    left_join(groups.transect.traj.reefds) %>%
    left_join(spatial_lookup %>% 
                dplyr::select(ZONE, Shelf, NRM, TUMRA, BIOREGION,
                              REEF, DEPTH.f, REEF.d) %>%
                distinct()) %>% 
    mutate(BIOREGION = factor(BIOREGION),
           ZONE = factor(ZONE),
           Shelf = factor(Shelf),
           NRM = factor(NRM),
           RP_ID = factor(RP_ID),
           zone.shelf = paste(ZONE, Shelf, sep = " ")) %>%
    ungroup() %>%
    suppressMessages() %>%
    suppressWarnings()
  
  ## save the trajectory data
  if (RPI_PURPOSE == "baseline") {
    save(filt.rec.traj,
         file = paste(DATA_PATH, "processed/",
                      sprintf(FILTER_OUT_DATA_FMT,max_init,min_final,min_obs)
                      ,"_", RPI_PURPOSE, ".RData",sep=""))
  } else {
    save(filt.rec.traj,
         file = paste(DATA_PATH, "processed/",
                      sprintf(FILTER_OUT_DATA_FMT,max_init,min_final,min_obs)
                      ,".RData",sep=""))
  }
  
}

CI_calc_rpi_trajectories <- function() {
  
  
  load(file = paste(DATA_PATH, "processed/",
                    sprintf(FILTER_OUT_DATA_FMT,max_init,min_final,min_obs)
                    ,"_", RPI_PURPOSE,".RData",sep=""))
  
  ## In lists throughout, reefs are often grouped by BIOREGION
  filt.rec.traj.proc <- vector(mode="list",
                               length = length(unique(filt.rec.traj$BIOREGION)))
  names(filt.rec.traj.proc) <- unique(filt.rec.traj$BIOREGION)
  
  for(reef_group in unique(filt.rec.traj$BIOREGION)){
    
    df <- filt.rec.traj %>%
      filter(BIOREGION == reef_group) %>%
      droplevels()
    rpid.list <- vector(mode = "list",
                        length = length(unique(df$RP_ID)))
    names(rpid.list) <- unique(df$RP_ID)
    
    for(rp_id in unique(df$RP_ID)){
      rpid.df <- df %>%
        filter(RP_ID == rp_id) %>%
        droplevels()
      rpid.reformat <- rpid.df %>% 
        reformat_recovery_trajectories() %>% 
        mutate(proj.site.rpid = paste(REEF.d, RP_ID, sep=" "),
               NRM = unique(rpid.df$NRM),
               BIOREGION = factor(reef_group),
               ZONE = unique(rpid.df$ZONE),
               Shelf = unique(rpid.df$Shelf),
               zone.shelf = unique(rpid.df$zone.shelf),
               REEF = unique(rpid.df$REEF),
               DEPTH.f = unique(rpid.df$DEPTH.f),
               REEF.d = unique(rpid.df$REEF.d))
      
      rpid.list[[rp_id]] = rpid.reformat
      
    }
    
    ## Nest each list under given reef group 
    filt.rec.traj.proc[[reef_group]] <- rpid.list
  }
  
  ## save the trajectory data with HC and AB data variables
  save(filt.rec.traj.proc,
       file=paste(DATA_PATH, "processed/",
                  sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
                  ,".RData",sep=""))
}


## External Functions for growth models

CI_model_rpi_fitGrowthModel <- function() {
  
  source("externalFunctions/RPI_baseline_4_coralFitGrowthModel.R") 
  
}

CI_model_rpi_gatherThin <- function() {

  source("externalFunctions/RPI_baseline_6_gather_and_thin.R") 
  
}

CI_model_rpi_calc_peakDensity <- function() {
  
  source("externalFunctions/RPI_baseline_7_BaseTrajCalc_PeakDensity.R")
  source("externalFunctions/RPI_baseline_8_BaseTrajFilter_alphaTd.R")
  
}

# CI_model_rpi_calc10YearIncrease <- function() {
# 
#         source("../R/externalFunctions/RPI_baseline_8_predict10YearIncrease.R")
# 
# }

## only need to do the following once - perhaps move into RERUN_BASELINES section
## Note locations of external function requiref for CI_34_RPI_* scripts require
# source"*" statements in CI_34_models_RPI.R to be run
CI_34_RPI_baseline_models <- function() {
    
        RPI_PURPOSE <<- 'baseline'
        CI_process_rpi_get_samples_data()
        CI_process_groups_data_part2()
        CI_process_rpi_data()
        CI_process_rpi_filter_trajectories()
        ## Calculate mean, sd and se for HC and Abiotic for each recovery trajectory
        CI_calc_rpi_trajectories()
        ## Fit growth model
        CI_model_rpi_fitGrowthModel()
        ## Gather and thin model posteriors
        CI_model_rpi_gatherThin()       
        ## Remove non-convergence and calculate peak density
        CI_model_rpi_calc_peakDensity()
        ## CI_model_rpi_calc10YearIncrease()  AT confered with Kerryn March 2025 and decided that the peakDensity version is the correct one to use
}

### additional functions required by CI_34_RPI_reference_models()
CI__process_rpi_data_reference <- function(current.report.year,
                                           groups.transect,
                                           disturbance,
                                           samples) {
  print(paste0("Focal report year: ", current.report.year))
  groups.transect.2<- groups.transect %>%
    filter(REPORT_YEAR <= current.report.year) %>% droplevels() %>%
    mutate(REEF.d.year = paste(REEF.d, fYEAR, sep = " "))
  
  disturbance.2 <- disturbance %>%
    filter(REPORT_YEAR <= current.report.year) %>% droplevels() %>% 
    mutate(REEF.d.year = paste(REEF.d, fYEAR, sep = " ")) %>%
    filter(REEF.d.year %in% unique(groups.transect.2$REEF.d.year)) %>%
    droplevels()
  
  samples.2 <- samples %>%
    filter(REPORT_YEAR <= current.report.year) %>% droplevels() %>% 
    mutate(REEF.d.year = paste(REEF.d, fYEAR, sep = " ")) %>%
    filter(REEF.d.year %in% unique(groups.transect.2$REEF.d.year)) %>%
    droplevels()
  
  ## Use function to assign each reef.d/report_year to a recovery trajectory
  recovery.trajectories <- extract_recovery_trajectories_transect(disturbance.2,
                                                                  samples.2,
                                                                  groups.transect.2)
  
  save(recovery.trajectories,
       file = paste0(DATA_PATH, "processed/RPI_reference/time.series.site.transect.crp.",
                     current.report.year,".RData"))
  return(recovery.trajectories)
}

CI_process_rpi_data_reference <- function() {
  
  load(paste0(DATA_PATH, 'processed/groups.transect.RData'))
  load(paste0(DATA_PATH, "processed/samples_rpi.RData"))
  load(paste0(DATA_PATH, "processed/disturbance_rpi.RData"))
  load(paste0(DATA_PATH, "processed/sample.reef.report.year.RData"))
  
  ## Get a vector of years to iterate through
  report.years <- sample.reef.report.year %>%
    dplyr::select(REPORT_YEAR) %>%
    unique() %>%
    filter(!REPORT_YEAR==1993) %>%
    droplevels() %>%
    arrange(REPORT_YEAR) %>%
    pull()
  save(report.years,
       file = paste0(DATA_PATH, "modelled/report.years.RData"))
  
  ## The following creates a large file that we will subsequently
  ## break down into temporal series of focal years down
  recovery.trajectories.final_year <-
    #  CI__process_rpi_data_reference(current.report.year = CI$setting$FINAL_YEAR,
    CI__process_rpi_data_reference(current.report.year = 2024,
                                   groups.transect,
                                   disturbance,
                                   samples)
  save(recovery.trajectories.final_year,
       file = paste0(DATA_PATH, "modelled/recovery.trajectories.final_year.RData"))
  
  
}

CI_process_rpi_crp <- function() {
  
  load(file = paste0(DATA_PATH, "modelled/report.years.RData"))
  load(file = paste0(DATA_PATH, "modelled/recovery.trajectories.final_year.RData"))
  
  N <- length(report.years)
  rpi_data <- tibble(REPORT_YEAR = report.years) %>%
    mutate(n = 1:n()) %>%
    mutate(recovery.trajectories.crp =
             map2(.x = REPORT_YEAR,
                  .y = n,
                  .f = ~ {
                      rpi_data <- recovery.trajectories.final_year %>%
                      filter(REPORT_YEAR <= .x) %>%
                      mutate(proj.site.rpid = paste(REEF.d, RP_ID, sep = " ")) %>%
                      ungroup() %>%
                      group_by(proj.site.rpid) %>% 
                      mutate(max.report.year = max(REPORT_YEAR)) %>%
                      filter(max.report.year == .x) %>%
                      droplevels() %>%
                      ungroup()
                    nm <- paste0(DATA_PATH, "processed/RPI_reference/",
                                 "rpi_data_", RPI_PURPOSE, "_", .x, "_stage1.RData")
                    save(rpi_data, file = nm)
                    nm
                  }
             ))
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "_stage1.RData"))
}

CI_process_rpi_augment <- function() {
  
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "_stage1.RData"))
  
  rpi_data <- rpi_data %>%
    mutate(ongoing.df = map2(.x = recovery.trajectories.crp,
                             .y = n,
                             .f = ~ {
                               nm <- .x
                               .x <- get(load(file = nm))
                               rpi_data <- CI__34_RPI_augment_trajectories(.x)
                               nm <- str_replace(nm, 'stage1', 'stage2')
                               save(rpi_data, file = nm)
                               nm
                             }
    ))
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "_stage2.RData"))
}



CI__34_RPI_augment_trajectories <- function(recovery.trajectories.crp) {
  ## Load baseline dataframe
  load(file=paste0(DATA_PATH, "processed/RPI.baseline.RData"))
  ## load benthic data
  load(file=paste0(DATA_PATH, "processed/groups.transect.RData"))
  ## Load spatial info
  load(file=paste0(DATA_PATH, "processed/spatial_lookup_rpi.RData"))
  
  ## Add benthic data to the trajectories of interest
  ## Also add useful spatial variables
  ## Aggregate Bioregions, so that bioregions with no baseline trajectories can,
  ## .... use parameters from the most similar bioregion
  
  groups.transect %>%
    right_join(recovery.trajectories.crp %>%
                 dplyr::select(REEF.d, REPORT_YEAR, RP_ID, Date, max.report.year)) %>%
    left_join(spatial_lookup %>%
                dplyr::select(ZONE, Shelf, BIOREGION, NRM, TUMRA, REEF.d)) %>%
    mutate(proj.site.rpid = paste(REEF.d, RP_ID, sep = " ")) %>%
    spread(key = "GROUP_CODE", value = COVER) %>% 
    group_by(REEF.d, RP_ID, proj.site.rpid) %>%
    mutate(NUM_OBS = as.numeric(length(unique(REPORT_YEAR)))) %>%
    mutate(BIOREGION.rpi.agg = factor(case_when(BIOREGION %in% c("18", "29")~"bio.18and29",
                                                BIOREGION %in% c("9", "26")~"bio.9and26",
                                                BIOREGION %in% c("39")~"bio.39",
                                                BIOREGION %in% c("12","35","36")~"bio.12and35and36",
                                                BIOREGION %in% c("3", "4")~"bio.3and4",
                                                BIOREGION %in% c("27")~"bio.27",
                                                BIOREGION %in% c("11")~"bio.11",
                                                BIOREGION %in% c("37")~"bio.37",
                                                BIOREGION %in% c("23")~"bio.23",
                                                BIOREGION %in% c("6")~"bio.6",
                                                BIOREGION %in% c("26")~"bio.26",
                                                BIOREGION %in% c("16")~"bio.16",
                                                BIOREGION %in% c("17")~"bio.17",
                                                BIOREGION %in% c("34", "38")~"bio.34and38",
                                                BIOREGION %in% c("19", "22")~"bio.19and22"))) %>%
    suppressMessages() %>%
    suppressWarnings()
  
}

CI__34_RPI_define_model <- function() {
  source(paste0('externalFunctions/DefineTwoPhaseGeneralModelSingleSpeciesTypeII.R'))
  list(ode_func = general_logistic_twophase,       # RHS for ODE model
       ode_sol = general_logistic_twophase_analytic,
       like_sampler = like_sampler,                # simulation of data generation process (for pred. checks)
       varnames = vlab)
  
}

CI_process_rpi_define_model <- function() {
  
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "_stage2.RData"))
  rpi_data <- rpi_data %>%
    mutate(model = map(.x = REPORT_YEAR,
                       .f = ~ CI__34_RPI_define_model()
    ))
  save(rpi_data, file = paste0(DATA_PATH,
                               "processed/RPI_reference/rpi_datai_",
                               RPI_PURPOSE,
                               "_stage3.RData"))
  
}

CI_process_rpi_data_predict_last_obs <- function() {
  
  
  load(file = paste0(DATA_PATH,
                     "processed/RPI_reference/rpi_datai_",
                     RPI_PURPOSE,
                     "_stage3.RData"))
  ## load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data.RData"))
  load(file=paste0(PROC_DATA_DIR, "RPI.baseline.RData"))
  
  ## The following will be run in parallel.  Note that progressor can be
  ## defined inside a function, but cannot be defined in global env, so
  ## do not try to run that line interactively
  plan(multisession, workers = 10)
  p <- progressor(steps = 1)
  rpi_data <-
    rpi_data  %>%
    mutate(depth_loop = future_map2(.x = ongoing.df, .y = REPORT_YEAR,
                                    .f = ~ {
                                      p()
                                      nm <- .x
                                      .x <- get(load(file = nm))
                                      rpi_data <- CI__34_RPI_depthLoop(.x,
                                                                       RPI.baseline,
                                                                       .y)
                                      nm <- str_replace(nm, 'stage2', 'stage4')
                                      save(rpi_data, file = nm)
                                      nm
                                    }
    )
    )
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "_stage4.RData"))
}

CI__34_RPI_depthLoop <- function(ongoing.df, RPI.baseline, current.report.year) {
  count=0
  ## For each Reef.d, predict HC for the last observation in the
  ## trajectory of interest Each Reef.d should only have 1 trajectory
  ## The last observation of the trajectory should be the current
  ## report year
  depth.list=vector(mode="list",
                    length=length(unique(ongoing.df$DEPTH.f)))
  names(depth.list)<- unique(ongoing.df$DEPTH.f)
  
  for (d in unique(ongoing.df$DEPTH.f)) {
    print(d)
    ongoing.depth <- ongoing.df %>% filter(DEPTH.f==d) %>% droplevels()
    baseline.depth <- RPI.baseline %>% filter(DEPTH.f==d) %>% droplevels
    site.list <- vector(mode = "list", length = length(unique(ongoing.depth$REEF.d)))
    names(site.list) <- unique(ongoing.depth$REEF.d)
    for(s in unique(ongoing.depth$REEF.d)){
      count = count+1
      print(paste("site count", count, sep = " "))
      print(paste("Site", s, sep = " "))
      ## Select the corresponding benthic data for the reef's
      ## trajectory
      ongoing.site.a <- ongoing.depth %>%
        filter(REEF.d == s) %>%
        droplevels() %>%
        ungroup()
      obs <- as.numeric(ongoing.site.a$NUM_OBS[1])
      if (obs > 1) {
        rpid <- as.character(unique(ongoing.site.a$RP_ID))
        print(paste("rpid", rpid, sep = " "))
        baseline.bioregion <- as.character(unique(ongoing.site.a$BIOREGION.rpi.agg))
        print(paste("baseline bioregion", baseline.bioregion, sep = " "))
        ## Select the model parameters from the corresponding
        ## bioregion group
        baseline.BIOREGION <- baseline.depth %>%
          filter(BIOREGION.rpi.agg == baseline.bioregion) %>%
          droplevels() %>%
          ungroup()
        if (nrow(baseline.BIOREGION) == 0) next
        else {
          ## Predict HC for the last observation of the
          ## trajectory
          try({
            current.traj.preds <- predict.ongoing.random(baseline.BIOREGION, 
                                                         ongoing.site.a, 
                                                         model,1)  %>% 
              left_join(ongoing.site.a %>% 
                          dplyr::select(ZONE, Shelf, NRM, TUMRA,
                                        BIOREGION, REEF, DEPTH.f,
                                        REEF.d, RP_ID, proj.site.rpid) %>%
                          distinct) %>%
              mutate(REPORT_YEAR = current.report.year) %>%
              suppressWarnings() %>% suppressMessages()
            site.list[[s]] = current.traj.preds })
        }
      } else
        print(paste("Not enough observations to predict for", s, current.report.year, sep=" "))
    }
    site.list.df <- do.call("rbind", site.list)
    depth.list[[d]] = site.list.df
  }
  
  pred.df.crp.full <- do.call("rbind", depth.list)
  if (is.null(pred.df.crp.full)) return(NULL)
  save(pred.df.crp.full,
       file = paste0(PROC_DATA_DIR, "pred.df.crp.full.",
                     current.report.year, ".random.RData"))
  ## Check if there are many NaN predictions
  nans <- pred.df.crp.full %>%
    filter(is.nan(HC_PRED)) %>%
    droplevels() %>%
    group_by(REEF.d) %>%
    summarise(nan.ss = n())
  ## Remove NaNs
  pred.df.crp <- pred.df.crp.full %>%
    filter(!is.nan(HC_PRED)) %>%
    droplevels
  toc()
  save(pred.df.crp,
       file = paste0(PROC_DATA_DIR, "pred.df.crp.",
                     current.report.year, ".random.RData"))
  rm(RPI.baseline)
  rm(baseline.depth)
  list(pred.df.crp.full = paste0(PROC_DATA_DIR, "pred.df.crp.full.",
                                 current.report.year, ".random.RData"),
       pred.df.crp = paste0(PROC_DATA_DIR, "pred.df.crp.",
                            current.report.year, ".random.RData")
  )
}

CI_process_rpi_gather_preds <- function() {
  
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "_stage4.RData"))
  rpi_data <- rpi_data %>%
    mutate(pred.df.crp.temporal = map(.x = depth_loop,
                                      .f = ~ {
                                        nm <- .x
                                        if (file.exists(nm)) {
                                          .x <- get(load(file = nm))
                                          if(is.null(.x)) return(NULL)
                                          nm <- .x$pred.df.crp
                                          if (file.exists(nm)) {
                                            .x <- get(load(file = nm))
                                            return(.x) 
                                          } else return(NULL)
                                        } else return(NULL)
                                      }
    ))
  
  pred.df.crp.temporal <- rpi_data %>%
    dplyr::select(pred.df.crp.temporal) %>%
    unnest(pred.df.crp.temporal)
  save(pred.df.crp.temporal, file = paste0(DATA_PATH,
                                           "processed/RPI_reference/pred.df.crp.temporal.random_",
                                           RPI_PURPOSE,
                                           "_.RData"))
  
  ## load benthic data
  load(file=paste0(DATA_PATH, "processed/groups.transect.RData"))
  ## Load spatial info
  load(file=paste0(DATA_PATH, "processed/spatial_lookup_rpi.RData"))
  
  modelled.HC <- groups.transect %>%
    filter(GROUP_CODE == "HC") %>%
    droplevels() %>%
    rename(modelled.cover = COVER)
  
  modelled.v.pred.dist <- modelled.HC %>%
    left_join(pred.df.crp.temporal %>%
                rename(expected.cover = HC_PRED)) %>%
    dplyr::select(-TUMRA, -NRM, -BIOREGION, -NRM, -ZONE, -Shelf) %>%
    left_join(spatial_lookup) %>%
    suppressMessages() %>% suppressWarnings()
  
  save(modelled.v.pred.dist, file=paste0(DATA_PATH,
                                         "/modelled/modelled.v.pred.dist.random.RData"))
  
}

CI_process_rpi_calc_recovery_index <- function() {
  
  
  load(file=paste0(DATA_PATH,
                   "/modelled/modelled.v.pred.dist.random.RData"))
  calculate.RPI.score <- modelled.v.pred.dist %>%
    mutate(distance.metric = log2(modelled.cover/expected.cover),
           cap.dist.met = as.numeric(case_when(distance.metric < -2 ~ -2,
                                               distance.metric > 2 ~ 2,
                                               distance.metric > -2 &
                                                 distance.metric < 2 ~  distance.metric)),
           rescale.dist.met = scales::rescale(cap.dist.met, to = c(0,1))) %>%
    filter(!distance.metric %in% NaN) %>%
    droplevels() ##There are some negative predictions and so some NaNs are produced by the log2 calculation
  
  save(calculate.RPI.score,
       file=paste0(DATA_PATH, "/modelled/calculate.RPI.score.random.RData"))
  
  
  ##     ## ################################################################################
  ##     ## Summarise RPI score
  ## ## The following is not really used for anything...
  ##     load(file=paste0(DATA_PATH, "/modelled/calculate.RPI.score.random.RData"))
  
  ##     current.pred.summary <- calculate.RPI.score %>%
  ##         group_by(ZONE, Shelf, TUMRA, NRM,BIOREGION,
  ##                  DEPTH.f, REEF, REEF.d, REPORT_YEAR) %>%
  ##         rename(index = rescale.dist.met) %>%
  ##         tidybayes::median_hdci(modelled.cover, index, expected.cover) %>%
  ##         mutate(fYEAR = factor((REPORT_YEAR))) %>%
  ##         group_by(ZONE, Shelf, TUMRA, NRM,
  ##                  BIOREGION, REEF, DEPTH.f, REEF.d) %>%
  ##         arrange(ZONE, Shelf, TUMRA, NRM, BIOREGION,
  ##                 REEF, DEPTH.f, REEF.d, REPORT_YEAR) %>%
  ##         tidyr::fill(index, index.upper, index.lower)  #For disturbance years, the index score from the last time it was measured is carried over
  
  ##     save(current.pred.summary,
  ##          file=paste0(DATA_PATH, "/modelled/current.pred.summary.random.RData"))
  
  ##############################
  ## Low coral cover correction
  
  ## when coral cover is very low (< 6) and time since
  ## disturbance is > 3yrs, and the upper 95% interval of the
  ## index is < 0.5, we will override the index to be 0 (along
  ## with the intervals)
  
  ## Start by calculating the time since disturbance for each proj.site.rpid/REPORT_YEAR
  load(file = paste0(DATA_PATH, "/modelled/calculate.RPI.score.random.RData"))
  load(file=paste0(DATA_PATH, "/processed/samples_rpi.RData"))
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     "reference",
                     "_stage4.RData"))
  
  time_since_disturbance <- rpi_data %>%
    dplyr::select(REPORT_YEAR, recovery.trajectories.crp) %>%
    mutate(TEMP = map(.x = recovery.trajectories.crp,
                      .f = ~ {
                        x <- get(load(.x))
                        x %>% group_by(REEF.d) %>%
                          filter(REPORT_YEAR == max.report.year) %>%
                          dplyr::select(REPORT_YEAR, REEF.d, proj.site.rpid)
                      })) %>%
    dplyr::select(TEMP) %>%
    unnest(TEMP) %>%
    left_join(
      samples %>% dplyr::select(REPORT_YEAR, REEF.d, Date)
    ) %>%
    group_by(proj.site.rpid) %>%
    arrange(REPORT_YEAR) %>%
    mutate(Time_since_disturbance = as.numeric(Date - first(Date))/365.25)
  
  RPI_reference_posteriors <- calculate.RPI.score %>%
    group_by(ZONE, Shelf, TUMRA, NRM, BIOREGION, DEPTH.f,
             REEF, REEF.d, REPORT_YEAR) %>%
    rename(index = rescale.dist.met) %>%
    mutate(index.upper = HDInterval::hdi(index)[2],
           modelled.cover = median(modelled.cover)) %>%
    left_join(time_since_disturbance %>%
                dplyr::select(REPORT_YEAR, REEF.d, Time_since_disturbance)) %>%
    mutate(ref.under = ifelse(index.upper < 0.5, TRUE, FALSE)) %>%
    mutate(old.index = index) %>% 
    mutate(index = ifelse(Time_since_disturbance > 3 &
                            modelled.cover < 6 &
                            !ref.under, 0, index)) %>%
    dplyr::rename(.draw = TRANSECT_NO)
  
  save(RPI_reference_posteriors,
       file = paste0(DATA_PATH, "/modelled/RPI_reference_posteriors.RData"))
  
}



## reference models
CI_34_RPI_reference_models <- function() {
  
  RPI_PURPOSE <<- 'reference'
  
  CI_process_rpi_get_samples_data()#ok
  CI_process_groups_data_part2() #ok
  CI_process_rpi_data_reference() #ok
  CI_process_rpi_crp() #ok
  CI_process_rpi_augment() #ok
  CI_process_rpi_define_model() #ok
  
  ### seems to be missing the function to fit the growth model he as after this all goes wonky 
  CI_process_rpi_data_predict_last_obs() #ok but only just runs on my machine
  CI_process_rpi_gather_preds()  # ok 
  CI_process_rpi_calc_recovery_index() #ok
  
}

## additional functions required for critical 

CI_process_rpi_critical_filter_trajectories <- function() {
  
  load(file = paste0(DATA_PATH, "modelled/report.years.RData"))
  load(file = paste0(DATA_PATH, "modelled/recovery.trajectories.final_year.RData"))
  ## load benthic data
  load(file=paste0(DATA_PATH, "processed/groups.transect.RData"))
  ## Load spatial info
  load(file=paste0(DATA_PATH, "processed/spatial_lookup_rpi.RData"))
  
  N <- length(report.years)
  rpi_data <- tibble(REPORT_YEAR = report.years) %>%
    mutate(n = 1:n()) %>%
    mutate(filt.rec.traj.critical =
             map2(.x = REPORT_YEAR,
                  .y = n,
                  .f = ~ {
                      rm.ongoing.trajectories <- recovery.trajectories.final_year %>%
                      mutate(proj.site.rpid = paste(REEF.d, RP_ID, sep = " ")) %>%
                      ungroup() %>%
                      filter(REPORT_YEAR <= .x) %>%
                      group_by(proj.site.rpid) %>%
                      mutate(max.report.year = max(REPORT_YEAR),
                             NUM_OBS = as.numeric(length(unique(REPORT_YEAR)))) %>%
                      filter(!is.na(Date), NUM_OBS > 3,
                             max.report.year == .x) %>%
                      droplevels() %>%
                      ungroup()
                    
                    ## This was to prevent missing data, which causes issues further down the pipeline
                    groups.transect.traj.reefds <- groups.transect %>%
                      filter(REEF.d %in% unique(rm.ongoing.trajectories$REEF.d)) %>%
                      droplevels()
                    
                    filt.rec.traj.critical <- rm.ongoing.trajectories %>%
                      left_join(groups.transect.traj.reefds) %>%
                      left_join(spatial_lookup %>%
                                  dplyr::select(ZONE, Shelf, NRM, TUMRA,
                                                BIOREGION, REEF, DEPTH.f, REEF.d) %>%
                                  distinct) %>%
                      mutate(BIOREGION = factor(BIOREGION),
                             ZONE= factor(ZONE),
                             Shelf=factor(Shelf),
                             NRM=factor(NRM),
                             RP_ID=factor(RP_ID),
                             zone.shelf=paste(ZONE, Shelf, sep=" ")) %>%
                      ungroup() %>%
                      suppressMessages() %>% suppressWarnings()
                    
                    nm <- paste0(DATA_PATH, "processed/RPI_reference/",
                                 "rpi_data_", RPI_PURPOSE, "_", .x, "_stage1.RData")
                    save(filt.rec.traj.critical, file = nm)
                    nm
                  }
             ))
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "_stage1.RData"))
}

CI_process_rpi_add_benthos <- function() {
  
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "_stage1.RData"))
  N <- max(rpi_data$n)
  rpi_data <- rpi_data %>%
    mutate(filt.rec.traj.proc.critical =
             pmap(.l = list(REPORT_YEAR, n, filt.rec.traj.critical),
                  .f = ~ {
                    nm <- ..3
                    .x <- get(load(nm)) 
                    filt.rec.traj.proc.critical <- CI__add_benthic_group_data(.x)
                    nm <- str_replace(nm, "stage1", "stage2")
                    save(filt.rec.traj.proc.critical, file = nm)
                    nm
                  }
             ))
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "_stage2.RData"))
}

CI__add_benthic_group_data <- function(filt.rec.traj.critical) {
  
  filt.rec.traj.proc.critical <- vector(mode="list",
                                        length = length(unique(filt.rec.traj.critical$BIOREGION)))
  names(filt.rec.traj.proc.critical) <- unique(filt.rec.traj.critical$BIOREGION)
  
  for(reef_group in unique(filt.rec.traj.critical$BIOREGION)){
    
    df<- filt.rec.traj.critical %>% filter(BIOREGION == reef_group) %>% droplevels()
    rpid.list<- vector(mode="list", length=length(unique(df$RP_ID)))
    names(rpid.list)<-unique(df$RP_ID)
    
    for(rp_id in unique(df$RP_ID)){
      
      rpid.df<- df %>% filter(RP_ID==rp_id) %>% droplevels()
      rpid.reformat<- rpid.df %>%
        reformat_recovery_trajectories() %>%
        mutate(proj.site.rpid=paste(REEF.d, RP_ID, sep=" "),
               NRM=unique(rpid.df$NRM),
               BIOREGION=factor(reef_group),
               ZONE=unique(rpid.df$ZONE),
               Shelf=unique(rpid.df$Shelf),
               zone.shelf=unique(rpid.df$zone.shelf),
               REEF=unique(rpid.df$REEF),
               DEPTH.f=unique(rpid.df$DEPTH.f),
               REEF.d=unique(rpid.df$REEF.d))
      
      rpid.list[[rp_id]]=rpid.reformat
    }
    
    ## Nest each list under given reef group
    filt.rec.traj.proc.critical[[reef_group]] <- rpid.list
  }
  filt.rec.traj.proc.critical
}

CI_process_rpi_fit_coral_growth_model <- function() {
  
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "_stage2.RData"))
  N <- max(rpi_data$n)
  rpi_data <- rpi_data %>%
    mutate(filt.rec.traj.proc.mcmc.res.critical =
             pmap(.l = list(REPORT_YEAR, n, filt.rec.traj.proc.critical),
                  .f = ~ {
                    YR <- ..1
                    nm <- ..3
                    .x <- get(load(nm)) 
                    if (length(.x) == 0) {
                      filt.rec.traj.proc.mcmc.res.critical <- NULL
                    } else 
                      filt.rec.traj.proc.mcmc.res.critical <-
                      CI__fit_coral_growth_model(.x, YR = YR)
                    nm <- str_replace(nm, "stage2", "stage3")
                    save(filt.rec.traj.proc.mcmc.res.critical, file = nm)
                    nm
                  }
             ))
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "_stage3.RData"))
}

CI__fit_coral_growth_model <- function(filt.rec.traj.proc.critical, YR) {
  ## import model definition ###
  source(paste0('externalFunctions/DefineTwoPhaseGeneralModelSingleSpeciesTypeII.R'))
  ##*******************************************************************************
  ##MMP
  ## build MCMC sampler configuration
  conf <- list(chains = 4,        # number of independent chains to use
               iter = 8000,       # number of sampling iterations (or number of iterations between diagnostic checks)  #default 4000
               burnin = 8000,     # number of iterations per burning/warmup step #default 4000
               CPUs = 4,          # number of CPUs available for parallel chains (optimal CPUs = chains)
               nadapt = 1,        # number of adaptation steps
               initscale = 0.1,
               Rthresh = 1.1,     # stopping criteria threshold for Gelman-Rubin statistic diagnostic check
               ESSthresh = 400,  # stopping criteria threshold for Effective Sample Size diagnostic check
               maxChecks = 5,
               convcheck = TRUE,  # repeat iterations until stopping criteria are satisfied #has been added #default TRUE
               maxInits = 10000
  )
  
  ## Fit Growth Model ####
  
  ## ## load filter results
  ## load(paste(PROC_DATA_DIR, sprintf(FILTER_REFMT_OUT_DATA_FMT,
  ##                                   max_init,min_final,min_obs)
  ##           ,".critical.RData",sep=""))
  
  ## Perform model calibration to each recovery trajectory using MCMC
  filt.rec.traj.proc.mcmc.res.critical <- vector("list", length = length(filt.rec.traj.proc.critical))
  names(filt.rec.traj.proc.mcmc.res.critical)<- unique(names(filt.rec.traj.proc.critical))
  count=0
  for (i in unique(names(filt.rec.traj.proc.critical))) {
    
    reef.group.list<- filt.rec.traj.proc.critical[[i]]
    
    rpid.list<- vector(mode='list', length=length(reef.group.list))
    names(rpid.list)<-unique(names(reef.group.list))
    
    for(jj in unique(names(reef.group.list))){
      
      count=count+1
      print(paste("trajectory number", count, sep=" "))
      ## ensure data is sorted in ascending time order
      traj <- reef.group.list[[jj]] %>% arrange(REPORT_YEAR) %>%
        mutate(NUM_OBS=as.numeric(length(unique(REPORT_YEAR))))
      
      ## one ahead, remove last (N) observation
      n <- dim(traj)[1]
      traj <- traj[1:(n-1), ]
      
      ## build data object
      data <- list(nVisits = length(traj$HC[-1]),         # number of visits excluding initial visit
                   c0 = traj$HC[1],                       # cover of initial visit
                   t0 = traj$T[1]/365.0,                  # time of initial visit in years
                   ts = traj$T[-1]/365.0,                 # time of final visit in years
                   K = K.limit - traj$AB[length(traj$HC)],    # carrying capacity cover
                   C = traj$HC[-1],                       # cover value time series
                   Serr = traj$HC_sd[-1])                 # cover standard deviation time series
      
      ## update upper bound for T_d
      model$upper[4] <- traj$T[length(traj$T)]/365.0
      print(paste0("MCMC sampling for rpid ", jj, " in bioregion ", i))
      ## Store model in list
      samples <- adaptMCMC_fit_ode_model(data,model,conf)
      
      if (length(samples)!=0) {
        ## Store predictions in list
        rpid.list[[jj]]=samples
      }
      
      save(samples, file=paste0(DATA_PATH, "modelled/mcmc.samples.rpid.", jj, ".bioregion.", i,
                                "_",YR, "_.RData"))
      
    }
    
    filt.rec.traj.proc.mcmc.res.critical[[i]] <- rpid.list
    
    
  }
  
  ## save results
  filt.rec.traj.proc.mcmc.res.critical
}


CI_process_rpi_gather_posteriors <- function() {
  
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "_stage3.RData"))
  N <- max(rpi_data$n)
  rpi_data <- rpi_data %>%
    mutate(mcmc.bioregions.mpsrf.critical =
             pmap(.l = list(REPORT_YEAR, n, filt.rec.traj.proc.mcmc.res.critical,
                            filt.rec.traj.proc.critical),
                  .f = ~ {
                    YR <- ..1
                    nm <- ..3
                    .x <- get(load(nm)) 
                    nm1 <- ..4
                    .x1 <- get(load(nm1)) 
                    
                    if (length(.x) == 0) {
                      mcmc.bioregions.mpsrf.critical <- NULL
                    } else { 
                      mcmc.bioregions.mpsrf.critical <-
                        CI__gather_posteriors(.x, YR = YR, .x1)
                    }
                    nm <- str_replace(nm, "stage3", "stage4")
                    save(mcmc.bioregions.mpsrf.critical, file = nm)
                    nm
                  }
             )
    )
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "_stage4.RData"))
}

CI__gather_posteriors <- function(filt.rec.traj.proc.mcmc.res.critical, YR,
                                  filt.rec.traj.proc.critical, type = "critical") {
  current.report.year <- YR 
  
  mcmc.list=vector(mode='list', length=length(filt.rec.traj.proc.mcmc.res.critical))
  names(mcmc.list)<- unique(names(filt.rec.traj.proc.mcmc.res.critical))
  
  mpsrf.list.group=vector(mode='list', length=length(unique(names(filt.rec.traj.proc.mcmc.res.critical))))
  names(mpsrf.list.group)<- unique(names(filt.rec.traj.proc.mcmc.res.critical))
  
  for (i in unique(names(filt.rec.traj.proc.mcmc.res.critical))) {
    
    reef.group.list<- filt.rec.traj.proc.mcmc.res.critical[[i]]
    
    rpid.list<- vector(mode='list', length=length(unique(names(reef.group.list))))
    names(rpid.list)<-unique(names(reef.group.list))
    
    mpsrf.list<- vector(mode='list', length=length(unique(names(reef.group.list))))
    names(mpsrf.list)<-unique(names(reef.group.list))
    
    
    for(mm in unique(names(reef.group.list))){
      
      ## print(paste("RP_ID", mm, sep=" "))
      
      spatial.var<- filt.rec.traj.proc.critical[[i]][[mm]]
      mcmc <- reef.group.list[[mm]]
      if (type == 'critical') {
        mcmc.df<- ggs(mcmc) %>%
          mutate(RP_ID=factor(mm), BIOREGION=factor(i)) %>%
          as.data.frame
      } else {
        mcmc.df<- ggs(mcmc) %>%
          mutate(previous.RP_ID=factor(mm), BIOREGION=factor(i)) %>%
          as.data.frame
      }
      
      rpid.list[[mm]]=mcmc.df
      
      if (type == 'critical') {
        mpsrf= data.frame(mpsrf=gelman.diag(reef.group.list[[mm]])[[2]],
                          ess=effectiveSize(reef.group.list[[mm]])[[2]]) %>%
          mutate(RP_ID=factor(mm),
                 NRM=spatial.var$NRM[1],
                 BIOREGION=factor(i),
                 ZONE=spatial.var$ZONE[1],
                 Shelf=spatial.var$Shelf[1],
                 REEF=spatial.var$REEF[1],
                 DEPTH.f=spatial.var$DEPTH.f[1],
                 REEF.d=spatial.var$REEF.d[1],
                 proj.site.rpid=spatial.var$proj.site.rpid[1],
                 REPORT_YEAR=current.report.year)
      } else if (type == 'previous') {
        mpsrf= data.frame(mpsrf=gelman.diag(reef.group.list[[mm]])[[2]],
                          ess=effectiveSize(reef.group.list[[mm]])[[2]]) %>%
          mutate(previous.RP_ID=factor(mm),
                 NRM=spatial.var$NRM[1],
                 BIOREGION=factor(i),
                 ZONE=spatial.var$ZONE[1],
                 Shelf=spatial.var$Shelf[1],
                 REEF=spatial.var$REEF[1],
                 DEPTH.f=spatial.var$DEPTH.f[1],
                 REEF.d=spatial.var$REEF.d[1],
                 proj.site.rpid=spatial.var$proj.site.rpid[1],
                 REPORT_YEAR=current.report.year)
      }
      
      mpsrf.list[[mm]]=mpsrf
      
    }
    
    
    rpid.list.df=do.call('rbind', rpid.list)
    
    mpsrf.list.df=do.call('rbind', mpsrf.list)
    
    mcmc.list[[i]]=rpid.list.df
    
    mpsrf.list.group[[i]]=mpsrf.list.df
    
  }
  
  if (type == "critical") {
    mpsrf.group.df.critical<- do.call('rbind', mpsrf.list.group)
    save(mpsrf.group.df.critical, file=paste0(PROC_DATA_DIR, "mpsrf.group.df.critical_",YR,"_.RData"))
    
    mcmc.group.df.critical<- do.call('rbind', mcmc.list)
    ## add mpsrf to mcmc
    mcmc.bioregions.mpsrf.critical<- mcmc.group.df.critical %>% ungroup() %>%
      left_join(mpsrf.group.df.critical)
    
    ## save mcmc
    save(mcmc.bioregions.mpsrf.critical, file=paste0(PROC_DATA_DIR, "mcmc.bioregions.mpsrf.critical_",YR,"_.RData"))
    
    return(mcmc.bioregions.mpsrf.critical)
  }
  if (type == "previous") {
    mpsrf.group.df.previous<- do.call('rbind', mpsrf.list.group)
    save(mpsrf.group.df.previous, file=paste0(PROC_DATA_DIR, "mpsrf.group.df.previous_",YR,"_.RData"))
    
    mcmc.group.df.previous<- do.call('rbind', mcmc.list)
    ## add mpsrf to mcmc
    mcmc.bioregions.mpsrf.previous<- mcmc.group.df.previous %>% ungroup() %>%
      left_join(mpsrf.group.df.previous)
    
    ## save mcmc
    save(mcmc.bioregions.mpsrf.previous, file=paste0(PROC_DATA_DIR, "mcmc.bioregions.mpsrf.previous_",YR,"_.RData"))
    
    return(mcmc.bioregions.mpsrf.previous)
  }
}

CI_process_rpi_filter_thin <- function() {
  
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "_stage4.RData"))
  N <- max(rpi_data$n)
  rpi_data <- rpi_data %>%
    mutate(mcmc.converge.info.critical =
             pmap(.l = list(REPORT_YEAR, n, mcmc.bioregions.mpsrf.critical),
                  .f = ~ {
                    YR <- ..1
                    nm <- ..3
                    .x <- get(load(nm)) 
                    if (length(.x) == 0) {
                      mcmc.converge.info.critical <- NULL
                    } else { 
                      ## This step takes a long time......
                      mcmc.converge.info.critical <- .x %>%
                        ungroup() %>%
                        group_by(proj.site.rpid) %>%
                        mutate(iterations=max(Iteration)) %>%
                        dplyr::select(ZONE, Shelf, NRM, BIOREGION,
                                      REEF, DEPTH.f, REEF.d, RP_ID,
                                      proj.site.rpid, REPORT_YEAR,
                                      iterations, mpsrf, ess) %>%
                        unique()
                    }
                    nm <- str_replace(nm, "stage4", "stage5")
                    save(mcmc.converge.info.critical, file = nm)
                    nm
                  }
             )
    )
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "_stage5.RData"))
}

CI_process_rpi_optimise <- function() {
  
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "_stage5.RData"))
  N <- max(rpi_data$n)
  rpi_data <- rpi_data %>%
    mutate(mcmc.traj.ess.critical =
             pmap(.l = list(REPORT_YEAR, n,
                            mcmc.converge.info.critical,
                            mcmc.bioregions.mpsrf.critical),
                  .f = ~ {
                    YR <- ..1
                    nm <- ..3
                    .x <- get(load(nm)) 
                    nm1 <- ..4
                    .x1 <- get(load(nm1))
                    
                    if (length(.x) == 0) {
                      mcmc.traj.ess.critical <- NULL
                    } else { 
                      mcmc.table.critical <- .x %>%
                        mutate(ess.thinning.interval=iterations/ess)
                      mcmc.traj.ess.critical<- .x1 %>%
                        ungroup() %>%
                        left_join(mcmc.table.critical) %>%
                        mutate(thin.to.ess=round(ess.thinning.interval)) %>%
                        group_by(NRM, BIOREGION, ZONE, Shelf, REEF,
                                 DEPTH.f, REEF.d, RP_ID, proj.site.rpid,
                                 REPORT_YEAR,
                                 Chain, Parameter) %>%
                        slice(which(row_number() %% thin.to.ess == 1)) %>%
                        mutate(zone.shelf=paste(ZONE, Shelf, sep=" ")) %>%
                        suppressMessages() %>%
                        suppressWarnings()
                      mcmc.traj.ess.critical
                    }
                    nm <- str_replace(nm, "stage5", "stage6")
                    save(mcmc.traj.ess.critical, file = nm)
                    nm
                  }
             )
    ) %>%
    ## Remove non-convergence
    mutate(rm.poor.chains.critical =
             map(.x = mcmc.traj.ess.critical,
                 .f = ~ {
                   nm <- .x
                   .x <- get(load(nm)) 
                   if (length(.x) == 0) {
                     rm.poor.chains.critical <- NULL
                   } else {
                     rm.poor.chains.critical <- .x %>%
                       filter(!mpsrf>1.4) %>% droplevels() ##Highest mpsrf was 1.39,
                     ##From visual assessment of these chains and densities,
                     ##... they aren't bad enough to throw out
                     ##Still leave some filter in for future safeguarding
                   }
                   
                   nm <- str_replace(nm, "stage6", "stage6a")
                   save(rm.poor.chains.critical, file = nm)
                   nm
                 }
             )
    )
  
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "_stage6.RData"))
}

CI__predict_last_obs <- function(filt.rec.traj.critical,
                                 rm.poor.chains.critical,
                                 YR,
                                 type = "critical"
) {
  current.report.year <- YR 
  ## load benthic data
  load(file=paste0(DATA_PATH, "processed/groups.transect.RData"))
  ## Load spatial info
  load(file=paste0(DATA_PATH, "processed/spatial_lookup_rpi.RData"))
  source(paste0('externalFunctions/DefineTwoPhaseGeneralModelSingleSpeciesTypeII.R'))
  
  if (type == 'critical') {
    rpids.to.predict.for <- filt.rec.traj.critical %>%
      ungroup() %>%
      ## { if(type == "critical")
      ##       filter(., RP_ID %in% rm.poor.chains.critical$RP_ID) %>% droplevels() 
      ##       else .
      ## } %>%
      filter(RP_ID %in% rm.poor.chains.critical$RP_ID) %>% droplevels() %>%
      dplyr::select(ZONE, Shelf, NRM, TUMRA, BIOREGION, REEF, DEPTH.f, REEF.d,
                    RP_ID, proj.site.rpid, REPORT_YEAR, max.report.year,
                    NUM_OBS, Date) %>% unique()
  } else if (type == 'previous') {
    rpids.to.predict.for <- filt.rec.traj.critical %>%
      ungroup() %>%
      dplyr::select(ZONE, Shelf, NRM, TUMRA, BIOREGION, REEF, DEPTH.f, REEF.d,
                    RP_ID, proj.site.rpid, REPORT_YEAR,
                    max.report.year, NUM_OBS, Date) %>%
      unique()
  }
  print(rpids.to.predict.for)
  
  ongoing.df<- groups.transect %>% right_join(rpids.to.predict.for) %>%
    spread(key="GROUP_CODE", value=COVER)
  print(nrow(ongoing.df))
  ## define the model for prediction
  model <- list(ode_func = general_logistic_twophase,       # RHS for ODE model
                ode_sol = general_logistic_twophase_analytic,
                like_sampler = like_sampler,                # simulation of data generation process (for pred. checks)
                varnames = vlab)
  count=0
  
  site.list<- vector(mode="list", length=length(unique(ongoing.df$RP_ID)))
  names(site.list)<- unique(ongoing.df$RP_ID)
  
  print(unique(rm.poor.chains.critical$previous.RP_ID))
  for(rpid in unique(ongoing.df$RP_ID)){
    
    count=count+1
    print(paste("site count", count, sep=" "))
    
    ongoing.site.a<- ongoing.df %>% filter(RP_ID==rpid) %>% droplevels() %>%
      ungroup()
    obs<-as.numeric(ongoing.site.a$NUM_OBS[1])
    reef.d<- as.character(ongoing.site.a$REEF.d[1])
    rpid<-as.character(unique(ongoing.site.a$RP_ID))
    
    print(paste("rpid", rpid, "for", reef.d, sep=" "))
    
    print(type)
    if (type == 'critical') {
      rpid.parameters <- rm.poor.chains.critical %>%
        filter(RP_ID == rpid) %>% droplevels() %>%
        ungroup()
    } else if (type == 'previous') {
      rpid.parameters <- rm.poor.chains.critical %>%
        filter(REEF.d == reef.d) %>% droplevels() %>%
        ungroup()
    }
    
    if (nrow(rpid.parameters)==0) next
    
    else {
      
      ##use predict.ongoing.random with N=1
      try({
        if (type == 'critical') {
          current.traj.preds <- predict.ongoing.random(rpid.parameters,
                                                       ongoing.site.a,
                                                       model,1)  %>%
            left_join(ongoing.site.a %>%
                        dplyr::select(ZONE, Shelf, NRM, TUMRA, BIOREGION, REEF, DEPTH.f, REEF.d, RP_ID, proj.site.rpid) %>% distinct) %>%
            mutate(REPORT_YEAR=current.report.year)
          
          #save(current.traj.preds, file=paste0(MCMC_OUTPUT_DIR, "current.traj.preds.", s, ".depth.", d, ".", y, ".RData"))
          
        } else if (type == 'previous') {
          current.traj.preds <- predict.previous.random(rpid.parameters,
                                                        ongoing.site.a,
                                                        model,1)  %>%
            left_join(ongoing.site.a %>%
                        dplyr::select(ZONE, Shelf, NRM, TUMRA, BIOREGION, REEF, DEPTH.f, REEF.d, RP_ID, proj.site.rpid) %>% distinct) %>%
            mutate(REPORT_YEAR=current.report.year)
        }
        
        
        site.list[[rpid]]=current.traj.preds })
      
    }
    
  }
  
  pred.df.critical.obs3.full<-do.call("rbind", site.list)
  
  ## save(pred.df.critical.obs3.full, file=paste0(PROC_DATA_DIR, "pred.df.critical.obs3.full.", current.report.year, ".random.RData"))
  
  nans<- pred.df.critical.obs3.full %>%
    filter(is.nan(HC_PRED)) %>% droplevels() %>%
    group_by(REEF.d) %>%
    summarise(nan.ss=n())
  
  pred.df.critical.obs3<- pred.df.critical.obs3.full %>%
    filter(!is.nan(HC_PRED)) %>% droplevels
  
  ## save(pred.df.critical.obs3, file=paste0(PROC_DATA_DIR, "pred.df.critical.obs3.", current.report.year, ".random.RData"))
  pred.df.critical.obs3
}

CI_process_rpi_predict_last <- function() {
  
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "_stage6.RData"))
  N <- max(rpi_data$n)
  rpi_data <- rpi_data %>%
    mutate(mcmc.traj.ess.critical =
             pmap(.l = list(REPORT_YEAR, n,
                            filt.rec.traj.critical,
                            rm.poor.chains.critical),
                  .f = ~ {
                    YR <- ..1
                    nm <- ..3
                    .x <- get(load(nm)) 
                    nm1 <- ..4
                    .x1 <- get(load(nm1))
                    if (is.null(.x)) { 
                      pred.df.critical.obs3 <- NULL
                    } else if (nrow(.x) == 0) {
                      pred.df.critical.obs3 <- NULL
                    } else { 
                      pred.df.critical.obs3 <- CI__predict_last_obs(.x, .x1, YR)
                    }
                    nm <- str_replace(nm, "stage.", "stage7")
                    save(pred.df.critical.obs3, file = nm)
                    nm
                  }
             )
    )  
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "_stage7.RData"))
}

CI_process_rpi_critical_filter_trajectories_less3 <- function() {
  
  load(file = paste0(DATA_PATH, "modelled/report.years.RData"))
  load(file = paste0(DATA_PATH, "modelled/recovery.trajectories.final_year.RData"))
  ## load benthic data
  load(file=paste0(DATA_PATH, "processed/groups.transect.RData"))
  ## Load spatial info
  load(file=paste0(DATA_PATH, "processed/spatial_lookup_rpi.RData"))
  
  N <- length(report.years)
  rpi_data <- tibble(REPORT_YEAR = report.years) %>%
    mutate(n = 1:n()) %>%
    mutate(rm.ongoing.trajectories =
             map2(.x = REPORT_YEAR,
                  .y = n,
                  .f = ~ {
                    rm.ongoing.trajectories <- recovery.trajectories.final_year %>%
                      mutate(proj.site.rpid = paste(REEF.d, RP_ID, sep = " ")) %>%
                      ungroup() %>%
                      filter(REPORT_YEAR <= .x) %>%
                      group_by(proj.site.rpid) %>%
                      mutate(max.report.year = max(REPORT_YEAR),
                             NUM_OBS = as.numeric(length(unique(REPORT_YEAR)))) %>%
                      filter(!is.na(Date), (NUM_OBS> 1 & NUM_OBS<=3),
                             max.report.year == .x) %>%
                      droplevels() %>%
                      ungroup()
                    
                    nm <- paste0(DATA_PATH, "processed/RPI_reference/",
                                 "rpi_data_", RPI_PURPOSE, "_less3a_", .x, "_stage1.RData")
                    save(rm.ongoing.trajectories, file = nm)
                    nm
                  }
             )) %>%
    mutate(filt.rec.traj.critical.lessthan3 =
             pmap(.l = list(REPORT_YEAR, n, rm.ongoing.trajectories),
                  .f = ~ {
                    .x <- ..1
                    rm.ongoing.trajectories <- get(load(..3))
                    filt.rec.traj.critical.lessthan3 <- rm.ongoing.trajectories %>%
                      left_join(groups.transect) %>%
                      left_join(spatial_lookup %>%
                                  dplyr::select(ZONE, Shelf, NRM, TUMRA,
                                                BIOREGION, REEF, DEPTH.f, REEF.d) %>%
                                  distinct) %>%
                      mutate(BIOREGION = factor(BIOREGION),
                             ZONE= factor(ZONE),
                             Shelf=factor(Shelf),
                             NRM=factor(NRM),
                             RP_ID=factor(RP_ID),
                             zone.shelf=paste(ZONE, Shelf, sep=" ")) %>%
                      ungroup() %>%
                      suppressMessages() %>% suppressWarnings()
                    
                    nm <- paste0(DATA_PATH, "processed/RPI_reference/",
                                 "rpi_data_", RPI_PURPOSE, "_less3b_", .x, "_stage1.RData")
                    save(filt.rec.traj.critical.lessthan3, file = nm)
                    nm
                  }
             ))
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "less3_stage1.RData"))
}

CI_process_rpi_critical_find_previous_trajectory <- function() {
  
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "less3_stage1.RData"))
  load(file = paste0(DATA_PATH, "modelled/recovery.trajectories.final_year.RData"))
  recovery.trajectories.critical <- recovery.trajectories.final_year 
  
  ## load benthic data
  load(file=paste0(DATA_PATH, "processed/groups.transect.RData"))
  ## Load spatial info
  load(file=paste0(DATA_PATH, "processed/spatial_lookup_rpi.RData"))
  
  rpi_data <- rpi_data %>%
    mutate(N = max(n)) %>%
    mutate(filt.rec.traj.previous =
             pmap(.l = list(REPORT_YEAR, n,
                            rm.ongoing.trajectories),
                  .f = ~ {
                    .x <- ..1
                    .y <- ..2
                    nm <- ..3
                    rm.ongoring.trajectories <- get(load(nm))
                    
                    reef.ds.to.look.for<- rm.ongoing.trajectories %>% 
                      dplyr::select(REEF.d, DEPTH.f, RP_ID, max.report.year) %>%
                      unique()
                    
                    recovery.trajectories.critical <- recovery.trajectories.critical %>%
                      filter(REPORT_YEAR <= .x) %>%
                      droplevels()  
                    
                    possibilities<- recovery.trajectories.critical %>%
                      mutate(proj.site.rpid=paste(REEF.d, RP_ID, sep=" ")) %>%
                      ungroup() %>%
                      group_by(proj.site.rpid) %>%
                      mutate(max.report.year=max(REPORT_YEAR),
                             NUM_OBS=as.numeric(length(unique(REPORT_YEAR)))) %>%
                      filter(REEF.d %in% reef.ds.to.look.for$REEF.d,
                             NUM_OBS>=3) %>% droplevels() %>%
                      ungroup() %>%
                      dplyr::select(REEF.d, max.report.year) %>% unique() %>%
                      rename(other.max.report.year=max.report.year)
                    
                    select.closest<- reef.ds.to.look.for %>%
                      left_join(possibilities) %>%
                      filter(!other.max.report.year==max.report.year) %>%
                      mutate(year.difference=max.report.year-other.max.report.year) %>%
                      group_by(REEF.d, max.report.year) %>%
                      filter(year.difference==min(year.difference)) %>% droplevels() %>%
                      mutate(other.proj.site.max.year=paste(REEF.d, other.max.report.year))
                    
                    previous.rpids<- recovery.trajectories.critical %>%
                      mutate(proj.site.rpid=paste(REEF.d, RP_ID, sep=" ")) %>%
                      ungroup() %>%
                      group_by(proj.site.rpid) %>%
                      mutate(max.report.year=max(REPORT_YEAR),
                             NUM_OBS=as.numeric(length(unique(REPORT_YEAR)))) %>%
                      ungroup() %>%
                      filter(REEF.d %in% select.closest$REEF.d) %>%
                      mutate(proj.site.max.year=paste(REEF.d, max.report.year)) %>%
                      filter(proj.site.max.year %in% select.closest$other.proj.site.max.year) %>%
                      rename(previous.RP_ID=RP_ID)
                    
                    filt.rec.traj.previous<- previous.rpids %>%
                      left_join(groups.transect) %>%
                      left_join(spatial_lookup %>%
                                  dplyr::select(ZONE, Shelf, NRM, TUMRA, BIOREGION, REEF, DEPTH.f, REEF.d) %>% distinct) %>%
                      mutate(BIOREGION = factor(BIOREGION),
                             ZONE= factor(ZONE),
                             Shelf=factor(Shelf),
                             NRM=factor(NRM),
                             previous.RP_ID=factor(previous.RP_ID),
                             zone.shelf=paste(ZONE, Shelf, sep=" ")) %>%
                      ungroup()
                    
                    nm <- paste0(DATA_PATH, "processed/RPI_reference/",
                                 "rpi_data_", RPI_PURPOSE, "_less3_", .x, "_stage2.RData")
                    save(filt.rec.traj.previous, file = nm)
                    nm
                  }
             ))
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "less3_stage2.RData"))
}

CI_process_rpi_critical_parameterise_model <- function() {
  
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "less3_stage2.RData"))
  rpi_data <- rpi_data %>%
    mutate(filt.rec.traj.proc.previous =
             pmap(.l = list(REPORT_YEAR, n,
                            filt.rec.traj.previous),
                  .f = ~ {
                    .x <- ..1
                    .y <- ..2
                    nm <- ..3
                    filt.rec.traj.previous <- get(load(nm))
                    
                    filt.rec.traj.proc.previous <- CI__parameterise_model(filt.rec.traj.previous)
                    ## save the trajectory data with HC and AB data variables
                    nm <- paste0(DATA_PATH, "processed/RPI_reference/",
                                 "rpi_data_", RPI_PURPOSE, "_less3_", .x, "_stage3.RData")
                    save(filt.rec.traj.proc.previous, file = nm)
                    nm
                  }
             ))
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "less3_stage3.RData"))
}

CI__parameterise_model <- function(filt.rec.traj.previous){
  
  filt.rec.traj.proc.previous <- vector(mode="list",
                                        length=length(unique(filt.rec.traj.previous$BIOREGION)))
  names(filt.rec.traj.proc.previous)<-unique(filt.rec.traj.previous$BIOREGION)
  
  for(reef_group in unique(filt.rec.traj.previous$BIOREGION)){
    
    df<- filt.rec.traj.previous %>% filter(BIOREGION==reef_group) %>% droplevels()
    rpid.list<- vector(mode="list", length=length(unique(df$previous.RP_ID)))
    names(rpid.list)<-unique(df$previous.RP_ID)
    
    for(rp_id in unique(df$previous.RP_ID)){
      
      rpid.df<- df %>% filter(previous.RP_ID==rp_id) %>% droplevels()
      rpid.reformat<- rpid.df %>%
        rename(RP_ID=previous.RP_ID) %>% ##function doesn't recognise 'previous.RP_ID'
        #mutate(GROUP_CODE=factor(GROUP_CODE)) %>%
        reformat_recovery_trajectories() %>%
        rename(previous.RP_ID=RP_ID) %>% ##change it back again
        mutate(proj.site.rpid=paste(REEF.d, previous.RP_ID, sep=" "),
               NRM=unique(rpid.df$NRM),
               BIOREGION=factor(reef_group),
               ZONE=unique(rpid.df$ZONE),
               Shelf=unique(rpid.df$Shelf),
               zone.shelf=unique(rpid.df$zone.shelf),
               REEF=unique(rpid.df$REEF),
               DEPTH.f=unique(rpid.df$DEPTH.f),
               REEF.d=unique(rpid.df$REEF.d))
      
      rpid.list[[rp_id]]=rpid.reformat
      
    }
    
    ## Nest each list under given reef group
    filt.rec.traj.proc.previous[[reef_group]] <- rpid.list
  }
  filt.rec.traj.proc.previous
}

CI_process_rpi_fit_coral_growth_model_less3 <- function() {
  
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "less3_stage3.RData"))
  N <- max(rpi_data$n)
  rpi_data <- rpi_data %>%
    mutate(filt.rec.traj.proc.mcmc.res.previous =
             pmap(.l = list(REPORT_YEAR, n, filt.rec.traj.proc.previous),
                  .f = ~ {
                    YR <- ..1
                    nm <- ..3
                    .x <- get(load(nm)) 
                    if (length(.x) == 0) {
                      filt.rec.traj.proc.mcmc.res.previous <- NULL
                    } else 
                      filt.rec.traj.proc.mcmc.res.previous <-
                      CI__fit_coral_growth_model(.x, YR = YR)
                    nm <- str_replace(nm, "stage3", "stage4")
                    save(filt.rec.traj.proc.mcmc.res.previous, file = nm)
                    nm
                  }
             ))
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "less3_stage4.RData"))
}

CI_process_rpi_gather_posteriors_less3 <- function() {
  
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "less3_stage4.RData"))
  N <- max(rpi_data$n)
  rpi_data <- rpi_data %>%
    mutate(mcmc.bioregions.mpsrf.previous =
             pmap(.l = list(REPORT_YEAR, n, filt.rec.traj.proc.mcmc.res.previous,
                            filt.rec.traj.proc.previous),
                  .f = ~ {
                    YR <- ..1
                    nm <- ..3
                    .x <- get(load(nm)) 
                    nm1 <- ..4
                    .x1 <- get(load(nm1)) 
                    
                    if (length(.x) == 0) {
                      mcmc.bioregions.mpsrf.previous <- NULL
                    } else { 
                      mcmc.bioregions.mpsrf.previous <-
                        CI__gather_posteriors(.x, YR = YR, .x1, type = "previous")
                    }
                    nm <- str_replace(nm, "stage4", "stage5")
                    save(mcmc.bioregions.mpsrf.previous, file = nm)
                    nm
                  }
             )
    )
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "less3_stage5.RData"))
}

CI_process_rpi_filter_thin_less3 <- function() {
  
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "less3_stage5.RData"))
  N <- max(rpi_data$n)
  rpi_data <- rpi_data %>%
    mutate(mcmc.converge.info.previous =
             pmap(.l = list(REPORT_YEAR, n, mcmc.bioregions.mpsrf.previous),
                  .f = ~ {
                    
                    YR <- ..1
                    nm <- ..3
                    .x <- get(load(nm)) 
                    if (length(.x) == 0) {
                      mcmc.converge.info.previous <- NULL
                    } else { 
                      ## This step takes a long time......
                      mcmc.converge.info.previous <- .x %>%
                        ungroup() %>%
                        group_by(proj.site.rpid) %>%
                        mutate(iterations=max(Iteration)) %>%
                        dplyr::select(ZONE, Shelf, NRM, BIOREGION,
                                      REEF, DEPTH.f, REEF.d, previous.RP_ID,
                                      proj.site.rpid, REPORT_YEAR,
                                      iterations, mpsrf, ess) %>%
                        unique()
                    }
                    nm <- str_replace(nm, "stage5", "stage6")
                    save(mcmc.converge.info.previous, file = nm)
                    nm
                  }
             )
    )
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "less3_stage6.RData"))
}

CI_process_rpi_optimise_less3 <- function() {
  
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "less3_stage6.RData"))
  N <- max(rpi_data$n)
  rpi_data <- rpi_data %>%
    mutate(mcmc.traj.ess.previous =
             pmap(.l = list(REPORT_YEAR, n,
                            mcmc.converge.info.previous,
                            mcmc.bioregions.mpsrf.previous),
                  .f = ~ {
                    
                    YR <- ..1
                    nm <- ..3
                    .x <- get(load(nm)) 
                    nm1 <- ..4
                    .x1 <- get(load(nm1))
                    
                    if (length(.x) == 0) {
                      mcmc.converge.info.previous <- NULL
                      mcmc.traj.ess.previous <- NULL
                    } else { 
                      mcmc.table.previous <- .x %>%
                        mutate(ess.thinning.interval=iterations/ess)
                      mcmc.traj.ess.previous<- .x1 %>%
                        ungroup() %>%
                        left_join(mcmc.table.previous) %>%
                        mutate(thin.to.ess=round(ess.thinning.interval)) %>%
                        group_by(NRM, BIOREGION, ZONE, Shelf, REEF,
                                 DEPTH.f, REEF.d, previous.RP_ID, proj.site.rpid,
                                 REPORT_YEAR,
                                 Chain, Parameter) %>%
                        slice(which(row_number() %% thin.to.ess == 1)) %>%
                        mutate(zone.shelf=paste(ZONE, Shelf, sep=" ")) %>%
                        suppressMessages() %>%
                        suppressWarnings()
                      mcmc.traj.ess.previous
                    }
                    nm <- str_replace(nm, "stage6", "stage7")
                    save(mcmc.traj.ess.previous, file = nm)
                    nm
                  }
             )
    ) %>%
    ## Remove non-convergence
    mutate(rm.poor.chains.previous =
             map(.x = mcmc.traj.ess.previous,
                 .f = ~ {
                   nm <- .x
                   .x <- get(load(nm)) 
                   if (length(.x) == 0) {
                     rm.poor.chains.previous <- NULL
                   } else {
                     rm.poor.chains.previous <- .x %>%
                       filter(!mpsrf>1.4) %>% droplevels() ##Highest mpsrf was 1.39,
                     ##From visual assessment of these chains and densities,
                     ##... they aren't bad enough to throw out
                     ##Still leave some filter in for future safeguarding
                   }
                   
                   nm <- str_replace(nm, "stage7", "stage7a")
                   save(rm.poor.chains.previous, file = nm)
                   nm
                 }
             )
    )
  
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "less3_stage7.RData"))
}

CI_process_rpi_predict_last_less3 <- function() {
  
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "less3_stage7.RData"))
  N <- max(rpi_data$n)
  rpi_data <- rpi_data %>%
    mutate(mcmc.traj.ess.previous =
             pmap(.l = list(REPORT_YEAR, n,
                            filt.rec.traj.critical.lessthan3,
                            rm.poor.chains.previous),
                  .f = ~ {
                    YR <- ..1
                    nm <- ..3
                    .x <- get(load(nm))
                    nm1 <- ..4
                    .x1 <- get(load(nm1))
                    
                    if (is.null(.x) | is.null(.x1)) {
                      pred.df.previous.obs3 <- NULL
                    } else if (nrow(.x) == 0) {
                      pred.df.previous.obs3 <- NULL
                    } else { 
                      pred.df.previous.obs3 <- CI__predict_last_obs(.x, .x1, YR,
                                                                    type = "previous")
                    }
                    nm <- str_replace(nm, "stage.", "stage8")
                    save(pred.df.previous.obs3, file = nm)
                    nm
                  }
             )
    )  
  
  save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                               RPI_PURPOSE,
                               "less3_stage8.RData"))
}


CI_process_rpi_critical_gather_preds <- function() {
  
  
  ## get predictions from critical
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "_stage7.RData"))
  rpi_data <- rpi_data %>%
    mutate(pred.df.critical.obs3.temporal =
             map(.x = mcmc.traj.ess.critical, #pred.df.critical.obs3,
                 .f = ~ {
                   nm <- .x
                   if (file.exists(nm)) {
                     .x <- get(load(file = nm))
                     return(.x) 
                   } else return(NULL)
                 }))
  
  pred.df.critical.obs3.temporal <- rpi_data %>%
    dplyr::select(pred.df.critical.obs3.temporal) %>%
    unnest(pred.df.critical.obs3.temporal) %>%
    mutate(trajectory.standard = "current")
  save(pred.df.critical.obs3.temporal, file = paste0(DATA_PATH,
                                                     "processed/RPI_reference/pred.df.critical.temporal.obs3.random_",
                                                     RPI_PURPOSE,
                                                     "_.RData"))
  
  ## Got up to here
  ## get predictions from critical less than 3
  load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                     RPI_PURPOSE,
                     "less3_stage8.RData"))
  rpi_data <- rpi_data %>%
    mutate(pred.df.critical.lessthan3.temporal =
             map(.x = mcmc.traj.ess.previous, 
                 .f = ~ {
                   nm <- .x
                   if (file.exists(nm)) {
                     .x <- get(load(file = nm))
                     return(.x) 
                   } else return(NULL)
                 }))
  
  pred.df.critical.lessthan3.temporal <- rpi_data %>%
    dplyr::select(pred.df.critical.lessthan3.temporal) %>%
    unnest(pred.df.critical.lessthan3.temporal) %>%
    mutate(trajectory.standard = "previous")
  save(pred.df.critical.lessthan3.temporal, file = paste0(DATA_PATH,
                                                          "processed/RPI_reference/pred.df.critical.temporal.lessthan3.random_",
                                                          RPI_PURPOSE,
                                                          "_.RData"))
  
}

CI_process_rpi_calc_critical_recovery_index <- function() {
  
  RPI_PURPOSE <- 'critical'
  load(file = paste0(DATA_PATH,
                     "processed/RPI_reference/pred.df.critical.temporal.obs3.random_",
                     RPI_PURPOSE,
                     "_.RData"))
  load(file = paste0(DATA_PATH,
                     "processed/RPI_reference/pred.df.critical.temporal.lessthan3.random_",
                     RPI_PURPOSE,
                     "_.RData"))
  
  critical.preds.distribution <- bind_rows(pred.df.critical.lessthan3.temporal,
                                           pred.df.critical.obs3.temporal)
  
  save(critical.preds.distribution,
       file = paste0(DATA_PATH, "/processed/RPI_reference/critical.preds.distribution.RData"))
  
  load(paste0(DATA_PATH, 'processed/groups.transect.RData'))
  load(file=paste0(DATA_PATH, "/processed/spatial_lookup_rpi.RData"))
  load(file = paste0(DATA_PATH, "/processed/RPI_reference/critical.preds.distribution.RData"))
  
  modelled.HC <- groups.transect %>%
    filter(GROUP_CODE=="HC") %>%
    droplevels() %>%
    dplyr::rename(modelled.cover = COVER)
  
  modelled.v.pred.dist.critical <- modelled.HC %>%
    left_join(critical.preds.distribution %>%
                dplyr::rename(expected.cover = HC_PRED)) %>%
    dplyr::select(-TUMRA, -NRM, -BIOREGION, -NRM, -ZONE, -Shelf) %>%
    left_join(spatial_lookup)
  
  save(modelled.v.pred.dist.critical,
       file=paste0(DATA_PATH, "/modelled/modelled.v.pred.dist.critical.random.RData"))
  
  ## Calculate recovery performance index:
  
  ##There are some negative predictions and so some NaNs are
  ## produced by the log2 calculation
  calculate.critical.score <- modelled.v.pred.dist.critical %>%
    mutate(distance.metric = log2(modelled.cover/expected.cover),
           cap.dist.met = as.numeric(case_when(distance.metric < -2 ~ -2,
                                               distance.metric > 2 ~2,
                                               distance.metric > -2 &
                                                 distance.metric < 2 ~  distance.metric)),
           rescale.dist.met = scales::rescale(cap.dist.met, to = c(0,1))) %>%
    filter(!distance.metric %in% NaN) %>%
    droplevels()
  
  save(calculate.critical.score,
       file = paste0(DATA_PATH, "/modelled/calculate.critical.score.random.RData"))
  
  
  load(file = paste0(DATA_PATH, "/modelled/calculate.critical.score.random.RData"))
  ##############################
  ## Low coral cover correction
  
  ## when coral cover is very low (< 6) and time since
  ## disturbance is > 3yrs, and the upper 95% interval of the
  ## index is < 0.5, we will override the index to be 0 (along
  ## with the intervals)
  
  load(file = paste0(DATA_PATH, "modelled/recovery.trajectories.final_year.RData"))
  
  time_since_disturbance <- recovery.trajectories.final_year %>%
    mutate(proj.site.rpid = paste(REEF.d, RP_ID)) %>%
    group_by(proj.site.rpid) %>%
    arrange(REPORT_YEAR) %>%
    mutate(Time_since_disturbance = as.numeric(Date - first(Date))/365.25)
  
  RPI_critical_posteriors <- calculate.critical.score %>%
    group_by(ZONE, Shelf, TUMRA, NRM, BIOREGION, DEPTH.f,
             REEF, REEF.d, REPORT_YEAR) %>%
    rename(index = rescale.dist.met) %>%
    mutate(index.upper = HDInterval::hdi(index)[2],
           modelled.cover = median(modelled.cover)) %>%
    left_join(time_since_disturbance %>%
                dplyr::select(REPORT_YEAR, REEF.d, Time_since_disturbance)) %>%
    mutate(critical.under = ifelse(index.upper < 0.5, TRUE, FALSE)) %>%
    mutate(old.index = index) %>% 
    mutate(index = ifelse(Time_since_disturbance > 3 &
                            modelled.cover < 4 &
                            !critical.under, 0, index)) %>%
    dplyr::rename(.draw = TRANSECT_NO)
  
  save(RPI_critical_posteriors,
       file = paste0(DATA_PATH, "/modelled/RPI_critical_posteriors.RData"))
  
}

CI_34_RPI_critical_models <- function() {
  
  RPI_PURPOSE <<- 'critical'
  CI_process_rpi_get_samples_data() #ok
  CI_process_groups_data_part2() #ok
  
  ## Make the reference data
  CI_process_rpi_data_reference() #ok
  
  #######################################################
  ##Filter trajectories with 4 obs or more             ##
  ##Parameterise them using all but the last obs       ##
  ##Predict for the last obs                           ##
  #######################################################
  CI_process_rpi_critical_filter_trajectories() #ok
  CI_process_rpi_add_benthos() #ok
  CI_process_rpi_fit_coral_growth_model() #ok
  CI_process_rpi_gather_posteriors() #ok
  CI_process_rpi_filter_thin() #ok
  CI_process_rpi_optimise() #ok
  CI_process_rpi_predict_last() # ok but should check the warnings relating to RP_ID being uninitiated, which is occurring due to nulls
  #######################################################################################
  ##Filter trajectories with 2 o3 obs
  ##Find a previous trajectory for the same Reef.d
  ##Use the parameters from the previous trajectory to predict
  ## for last obs of the current trajectory
  ##If there is no previous trajectory, add the Reef.d to a list called 'needs.neighbour'
  ##Meaning that reef will need to use parameters from a nearest neighbour
  #######################################################################################
  CI_process_rpi_critical_filter_trajectories_less3() #ok
  CI_process_rpi_critical_find_previous_trajectory() #ok but warnings could be checked
  CI_process_rpi_critical_parameterise_model() #ok
  CI_process_rpi_fit_coral_growth_model_less3() #ok
  CI_process_rpi_gather_posteriors_less3() #ok
  CI_process_rpi_filter_thin_less3() #ok
  CI_process_rpi_optimise_less3() #ok
  CI_process_rpi_predict_last_less3() #ok
  
  ## Put them both together
  CI_process_rpi_critical_gather_preds() #ok
  CI_process_rpi_calc_critical_recovery_index() #ok
  
}

#######
# unused??

CI_process_rpi_crp_critical <- function() {
    
        load(file = paste0(DATA_PATH, "modelled/report.years.RData"))
        load(file = paste0(DATA_PATH, "modelled/recovery.trajectories.final_year.RData"))
        load(file = paste0(DATA_PATH, 'processed/groups.transect.RData'))
        load(paste0(DATA_PATH, 'processed/spatial_lookup_rpi.RData'))
        
        N <- length(report.years)
        rpi_data <- tibble(REPORT_YEAR = report.years) %>%
            mutate(n = 1:n()) %>%
            mutate(rm.ongoing.trajectories =
                       map2(.x = REPORT_YEAR,
                            .y = n,
                            .f = ~ {
                                   recovery.trajectories.final_year %>%
                                   filter(REPORT_YEAR <= .x) %>%
                                   mutate(proj.site.rpid = paste(REEF.d, RP_ID, sep = " ")) %>%
                                   ungroup() %>%
                                   group_by(proj.site.rpid) %>% 
                                   mutate(max.report.year = max(REPORT_YEAR),
                                          NUM_OBS = as.numeric(length(unique(REPORT_YEAR)))) %>%
                                   filter(max.report.year == .x) %>%
                                   droplevels() %>%
                                   ungroup()
                            }
                            ))
        ## The following is to prevent missing data which causes issues further down the pipeline
        rpi_data <- rpi_data %>%
            mutate(groups.transect.traj.reefds =
                       map(.x = rm.ongoing.trajectories,
                           .f = ~ groups.transect %>%
                               filter(REEF.d %in%
                                      unique(.x$REEF.d)) %>%
                               droplevels()
                           )) %>%
            mutate(filt.rec.traj.critical =
                       map2(.x = rm.ongoing.trajectories,
                            .y = groups.transect.traj.reefds,
                            .f = ~ .x %>%
                                left_join(.y) %>%
                                left_join(spatial_lookup %>%
                                          dplyr::select(ZONE, Shelf, NRM,
                                                        TUMRA, BIOREGION, REEF,
                                                        DEPTH.f, REEF.d) %>%
                                          distinct) %>%
                                mutate(BIOREGION = factor(BIOREGION),
                                       ZONE= factor(ZONE),
                                       Shelf=factor(Shelf),
                                       NRM=factor(NRM),
                                       RP_ID=factor(RP_ID),
                                       zone.shelf=paste(ZONE, Shelf, sep=" ")) %>%
                                ungroup() %>%
                                suppressMessages() %>%
                                suppressWarnings()
                            ))

        save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                                     RPI_PURPOSE,
                                    "_stage1.RData"))
        
}

CI_model_rpi_filter_thin <- function() {
    
        load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                           RPI_PURPOSE,
                           "_stage3.RData"))
        N <- rpi_data %>% pull(REPORT_YEAR) %>% unique() %>% length()
        rpi_data <- rpi_data %>%
            mutate(mcmc.converge.info.critical =
                       map2(.x = mcmc.bioregions.mpsrf.critical,
                            .y = n,
                            .f = ~ {
                                   .x %>%
                                   ungroup() %>%
                                   group_by(proj.site.rpid) %>%
                                   mutate(iterations=max(Iteration)) %>%
                                   dplyr::select(ZONE, Shelf, NRM, BIOREGION, REEF,
                                                 DEPTH.f, REEF.d, RP_ID,
                                                 proj.site.rpid, REPORT_YEAR,
                                                 iterations, mpsrf, ess) %>%
                                   unique()
                               
                           }))
        
        save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                                     RPI_PURPOSE,
                                    "_stage4.RData"))
}

CI__gather <- function(filt.rec.traj.proc.mcmc.res.critical, YR) {
  mcmc.list=vector(mode='list', length=length(filt.rec.traj.proc.mcmc.res.critical))
  names(mcmc.list)<- unique(names(filt.rec.traj.proc.mcmc.res.critical))
  
  mpsrf.list.group=vector(mode='list', length=length(unique(names(filt.rec.traj.proc.mcmc.res.critical))))
  names(mpsrf.list.group)<- unique(names(filt.rec.traj.proc.mcmc.res.critical))
  
  for (i in unique(names(filt.rec.traj.proc.mcmc.res.critical))) {
    
    print(paste("bioregion", i, sep=" "))
    
    reef.group.list<- filt.rec.traj.proc.mcmc.res.critical[[i]]
    
    rpid.list<- vector(mode='list', length=length(unique(names(reef.group.list))))
    names(rpid.list)<-unique(names(reef.group.list))
    
    mpsrf.list<- vector(mode='list', length=length(unique(names(reef.group.list))))
    names(mpsrf.list)<-unique(names(reef.group.list))
    
    
    for(mm in unique(names(reef.group.list))){
      
      print(paste("RP_ID", mm, sep=" "))
      
      spatial.var<- filt.rec.traj.proc.critical[[i]][[mm]]
      mcmc <- reef.group.list[[mm]]
      mcmc.df<- ggs(mcmc) %>% mutate(RP_ID=factor(mm), BIOREGION=factor(i)) %>% as.data.frame
      
      rpid.list[[mm]]=mcmc.df
      
      mpsrf= data.frame(mpsrf=gelman.diag(reef.group.list[[mm]])[[2]],
                        ess=effectiveSize(reef.group.list[[mm]])[[2]]) %>%
        mutate(RP_ID=factor(mm), NRM=spatial.var$NRM[1], BIOREGION=factor(i), ZONE=spatial.var$ZONE[1],
               Shelf=spatial.var$Shelf[1], REEF=spatial.var$REEF[1],
               DEPTH.f=spatial.var$DEPTH.f[1], REEF.d=spatial.var$REEF.d[1], proj.site.rpid=spatial.var$proj.site.rpid[1],
               REPORT_YEAR=current.report.year)
      
      
      
      mpsrf.list[[mm]]=mpsrf
      
    }
    
    
    rpid.list.df=do.call('rbind', rpid.list)
    
    mpsrf.list.df=do.call('rbind', mpsrf.list)
    
    mcmc.list[[i]]=rpid.list.df
    
    mpsrf.list.group[[i]]=mpsrf.list.df
    
  }
  
  
  mpsrf.group.df.critical<- do.call('rbind', mpsrf.list.group)
  save(mpsrf.group.df.critical, file=paste0(PROC_DATA_DIR, "mpsrf.group.df.critical.",YR,".RData"))
  
  mcmc.group.df.critical<- do.call('rbind', mcmc.list)
  #add mpsrf to mcmc
  if (is.null(mcmc.group.df.critical)) {
    mcmc.bioregions.mpsrf.critical <- NULL
  } else {
    mcmc.bioregions.mpsrf.critical<- mcmc.group.df.critical %>%
      ungroup() %>%
      left_join(mpsrf.group.df.critical)
  }
  #save mcmc
  save(mcmc.bioregions.mpsrf.critical, file=paste0(PROC_DATA_DIR, "mcmc.bioregions.mpsrf.critical.",YR,".RData"))
  mcmc.bioregions.mpsrf.critical
}

CI_model_rpi_gather <- function() {
   
        load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                           RPI_PURPOSE,
                           "_stage2.RData"))

        ## Unlist the posteriors and convert to dataframes
        ## Get the mpsrf statistic, which speaks to convergence.
        ## Thin the chains to Effective Sample Size

        N <- rpi_data %>% pull(REPORT_YEAR) %>% unique() %>% length()
        rpi_data <- rpi_data %>%
            mutate(mcmc.bioregions.mpsrf.critical =
                       pmap(.l = list(filt.rec.traj.proc.mcmc.res.critical,
                                      n,
                                      REPORT_YEAR),
                            .f = ~ {
                               CI__gather(..1, ..3)
                            }
                            )
                    
                   )
        
        save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                                     RPI_PURPOSE,
                                    "_stage3.RData"))
}

CI__benthic_groups <- function(filt.rec.traj.critical) {
  filt.rec.traj.proc.critical <- vector(mode="list",
                                        length = length(unique(filt.rec.traj.critical$BIOREGION)))
  names(filt.rec.traj.proc.critical) <- unique(filt.rec.traj.critical$BIOREGION)
  
  for(reef_group in unique(filt.rec.traj.critical$BIOREGION)){
    df<- filt.rec.traj.critical %>% filter(BIOREGION==reef_group) %>% droplevels()
    rpid.list<- vector(mode="list", length=length(unique(df$RP_ID)))
    names(rpid.list)<-unique(df$RP_ID)
    for(rp_id in unique(df$RP_ID)){
      rpid.df<- df %>% filter(RP_ID==rp_id) %>% droplevels()
      rpid.reformat<- rpid.df %>%
        reformat_recovery_trajectories() %>%
        mutate(proj.site.rpid=paste(REEF.d, RP_ID, sep=" "),
               NRM=unique(rpid.df$NRM),
               BIOREGION=factor(reef_group),
               ZONE=unique(rpid.df$ZONE),
               Shelf=unique(rpid.df$Shelf),
               zone.shelf=unique(rpid.df$zone.shelf),
               REEF=unique(rpid.df$REEF),
               DEPTH.f=unique(rpid.df$DEPTH.f),
               REEF.d=unique(rpid.df$REEF.d))
      rpid.list[[rp_id]]=rpid.reformat
    }
    ## Nest each list under given reef group
    filt.rec.traj.proc.critical[[reef_group]] <- rpid.list
  }
  filt.rec.traj.proc.critical
}

CI_process_rpi_benthic <- function() {
    
        load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                           RPI_PURPOSE,
                           "_stage1.RData"))
        ## Calculate mean, sd and se for HC and Abiotic for each recovery trajectory
        rpi_data <- rpi_data %>%
            mutate(filt.rec.traj.proc.critical =
                       map(.x = filt.rec.traj.critical,
                           .f = ~ CI__benthic_groups(.x)
                           )
                   )
        
        save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                                     RPI_PURPOSE,
                                    "_stage2.RData"))
        
}

CI__model_calibration_each_recovery_traj <- function(filt.rec.traj.proc.critical,
                                                     YR,
                                                     conf) {
  
  source(paste0('externalFunctions/DefineTwoPhaseGeneralModelSingleSpeciesTypeII.R'))
  ## YR <- filt.rec.traj.proc.critical[[1]] %>% pull(REPORT_YEAR) %>% max()
  filt.rec.traj.proc.mcmc.res.critical <- vector("list", length = length(filt.rec.traj.proc.critical))
  names(filt.rec.traj.proc.mcmc.res.critical) <- unique(names(filt.rec.traj.proc.critical))
  count=0
  for (i in unique(names(filt.rec.traj.proc.critical))) {
    reef.group.list<- filt.rec.traj.proc.critical[[i]]
    rpid.list<- vector(mode='list', length=length(reef.group.list))
    names(rpid.list)<-unique(names(reef.group.list))
    for(jj in unique(names(reef.group.list))){
      count=count+1
      print(paste("trajectory number", count, sep=" "))
      ## ensure data is sorted in ascending time order
      traj <- reef.group.list[[jj]] %>% arrange(REPORT_YEAR) %>%
        mutate(NUM_OBS=as.numeric(length(unique(REPORT_YEAR))))
      
      ## one ahead, remove last (N) observation
      n <- dim(traj)[1]
      traj <- traj[1:(n-1), ]
      
      ## build data object
      data <- list(nVisits = length(traj$HC[-1]),         # number of visits excluding initial visit
                   c0 = traj$HC[1],                       # cover of initial visit
                   t0 = traj$T[1]/365.0,                  # time of initial visit in years
                   ts = traj$T[-1]/365.0,                 # time of final visit in years
                   K = K.limit - traj$AB[length(traj$HC)],    # carrying capacity cover
                   C = traj$HC[-1],                       # cover value time series
                   Serr = traj$HC_sd[-1])                 # cover standard deviation time series
      
      ## update upper bound for T_d
      model$upper[4] <- traj$T[length(traj$T)]/365.0
      print(paste0("MCMC sampling for rpid ", jj, " in bioregion ", i))
      ## Store model in list
      samples <- adaptMCMC_fit_ode_model(data,model,conf)
      
      if (length(samples)!=0) {
        ## Store predictions in list
        rpid.list[[jj]]=samples
      }
      
      save(samples, file=paste0(MCMC_OUTPUT_DIR, "mcmc.samples.rpid.",
                                YR,".",jj, ".bioregion.", i, ".RData"))
      
    }
    
    filt.rec.traj.proc.mcmc.res.critical[[i]] <- rpid.list
    
    
  }
  
}
CI_model_rpi_critical <- function() {
    
        load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                           RPI_PURPOSE,
                           "_stage2.RData"))
        N <- rpi_data %>% pull(REPORT_YEAR) %>% unique() %>% length()

        ## import model definition ###
        source(paste0('externalFunctions/DefineTwoPhaseGeneralModelSingleSpeciesTypeII.R'))
        ## build MCMC sampler configuration
        conf <- list(chains = 4,        # number of independent chains to use
                     iter = 8000,       # number of sampling iterations (or number of iterations between diagnostic checks)  #default 4000
                     burnin = 8000,     # number of iterations per burning/warmup step #default 4000
                     CPUs = 4,          # number of CPUs available for parallel chains (optimal CPUs = chains)
                     nadapt = 1,        # number of adaptation steps
                     initscale = 0.1,
                     Rthresh = 1.1,     # stopping criteria threshold for Gelman-Rubin statistic diagnostic check
                     ESSthresh = 400,  # stopping criteria threshold for Effective Sample Size diagnostic check
                     maxChecks = 5,
                     convcheck = TRUE,  # repeat iterations until stopping criteria are satisfied #has been added #default TRUE
                     maxInits = 10000
                     )

        ## Fit Growth Model ####
        ## data will also be stored in <MCMC_OUTPUT_DIR>/mcmc.samples.rpid.<YR>.<jj>.bioregion.<i>.RData"
        rpi_data <- rpi_data %>%
            mutate(filt.rec.traj.proc.mcmc.res.critical =
                       pmap(.l = list(filt.rec.traj.proc.critical,
                                      n,
                                      REPORT_YEAR),
                            .f = ~ {
                               CI__model_calibration_each_recovery_traj(..1, ..3, conf)
                               }
                           )
                   )
        
        save(rpi_data, file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
                                     RPI_PURPOSE,
                                    "_stage2.RData"))

}



CI_models_RPI_aggregation <- function(level = 'NRM') {
  
  
  mods <- get(load(file = paste0(DATA_PATH, 'modelled/RPI__scores_reef_year.RData')))
  spatial_lookup <- get(load(file = paste0(DATA_PATH, 'processed/spatial_lookup.RData')))
  
  
  mods <-  mods %>%
    dplyr::select(Scores) %>% 
    unnest(Scores) %>% 
    dplyr::select(fYEAR, REEF.d, .draw, Metric, .value) %>%
    # add in infill for disturbance years
    group_by(REEF.d, .draw, Metric) |> 
    arrange(as.numeric(as.character(fYEAR))) |> 
    tidyr::fill(.value)  |> 
    ungroup()  |> 
    left_join(spatial_lookup %>%
                dplyr::select(REEF.d, MMP.Report,!!level)  |> 
                distinct()
    ) %>%
    filter(!is.na(!!sym(level))) %>%           ## exclude all reefs outside boundary
    filter(MMP.Report=="TRUE") |> 
    droplevels() |> 
    dplyr::select(-MMP.Report) |>   
    group_by(!!sym(level), Metric, fYEAR, .draw) %>%
    summarise(.value = mean(.value, na.rm=TRUE)) %>%
    ungroup() |> 
    group_by(fYEAR, !!sym(level), Metric) |> 
    summarise(Score=median(.value, na.rm=TRUE),
              q95=quantile(.value, 0.95, na.rm=TRUE)) |> 
    mutate(BelowPar=ifelse(q95<0.5, 0,1),
           Indicator="Performance") |> 
    ungroup()
  
  save(mods,
       file = paste0(DATA_PATH, 'modelled/RPI__scores_', level,'_year.RData'))

}











# test <- function() {
#     load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
#                        RPI_PURPOSE,
#                        "less3_stage3.RData"))
#     nm <- rpi_data[7,'filt.rec.traj.proc.previous'][[1]][[1]]
#     .x <- get(load(nm)) 
#     ## This must be a list
#     class(.x)
#     str(.x)
#     }
# 
# 
# 
# test <- function() {
#     load(file = paste0(DATA_PATH, "processed/RPI_reference/rpi_data_",
#                        RPI_PURPOSE,
#                        "less3_stage4.RData"))
#     nm <- rpi_data[7,'filt.rec.traj.proc.mcmc.res.previous'][[1]][[1]]
#     .x <- get(load(nm)) 
#     nm1 <- rpi_data[7,'filt.rec.traj.proc.previous'][[1]][[1]]
#     .x1 <- get(load(nm1)) 
#     .x
# 
#     }



# 
# 
# CI_models_RPI_distance <- function() {
# 
#         baselines <- get(load(file = paste0(DATA_PATH,
#                                             'modelled/CC__baseline_posteriors.RData')))
#         mods <- get(load(file = paste0(DATA_PATH, "modelled/CC__preds.RData")))
#         load(file=paste0(DATA_PATH, 'processed/site.location.RData'))
# 
#         mods <- mods %>%
#             mutate(Pred = map(.x = Pred,
#                               .f = ~ .x %>%
#                                   left_join(site.location %>%
#                                             dplyr::select(REEF, REEF.d, BIOREGION.agg,
#                                                           DEPTH.f) %>%
#                                             distinct())
#                               )) %>%
#             mutate(Scores = map(.x = Pred,
#                                 .f = ~ CI__index_CC(.x, baselines) %>%
#                                     filter(Metric %in% c('rescale.dist.metric',
#                                                          'pcb.rescale.dist.metric')) %>%
#                                     mutate(fYEAR = factor(fYEAR, levels = unique(fYEAR))) %>%
#                                     arrange(fYEAR, .draw)
#                                )) %>%
#             dplyr::select(-any_of(c("data", "newdata","Full_data", "Pred", "Summary"))) %>%
#             mutate(Summary = map(.x = Scores,
#                                  .f = ~ .x %>%
#                                      dplyr::select(-any_of(c(
#                                                 "P_CODE",
#                                                 "Model",
#                                                 "value",
#                                                 "baseline",
#                                                 "DEPTH.f"))) %>%
#                                      group_by(fYEAR, REEF, REEF.d, BIOREGION.agg, Metric) %>%
#                                      summarise_draws(median, mean, sd,
#                                                      HDInterval::hdi,
#                                                      `p<0.5` = ~ mean(.x < 0.5)
#                                                      )
#                                  ),
#                    Below = map(.x = Summary,
#                                .f = ~ .x %>%
#                                    ungroup() %>%
#                                    mutate(Below = ifelse(upper < 0.5, 1, 0),
#                                           PBelow = ifelse(`p<0.5` > 0.9, 1, 0)) %>%
#                                    dplyr::select(fYEAR, Metric, Below, PBelow) %>%
#                                    distinct()
#                                )
#                    ) %>%
#             suppressMessages() %>%
#             suppressWarnings()
# 
#         save(mods,
#               file = paste0(DATA_PATH, 'modelled/CC__scores_reef_year.RData'))
# 
# }

CI_models_RPI_distance <- function() {
    
        ## Get reference indices posteriors
        load(file = paste0(DATA_PATH, "/modelled/RPI_reference_posteriors.RData"))
        ## Get critical indices posteriors
        load(file = paste0(DATA_PATH, "/modelled/RPI_critical_posteriors.RData"))

        mods <- RPI_reference_posteriors %>%
            ungroup() %>%
            dplyr::select(fYEAR, REEF.d, .draw, REEF, BIOREGION.agg,
                          DEPTH.f, .value = index) %>%
            mutate(Metric = "reference") %>%
            bind_rows(
                RPI_critical_posteriors %>%
                ungroup() %>%
                dplyr::select(fYEAR, REEF.d, .draw, REEF, BIOREGION.agg,
                              DEPTH.f, .value = index) %>%
                mutate(Metric = "critical")
            ) %>%
            arrange(REEF.d, fYEAR, .draw, Metric) %>%
            mutate(fYEAR = factor(fYEAR)) %>%
            group_by(REEF.d) %>%
            summarise(Scores = list(cur_data_all()), .groups = "drop") %>%
            dplyr::select(-any_of(c("data", "newdata","Full_data", "Pred", "Summary"))) %>%
            mutate(Summary = map(.x = Scores,
                                 .f = ~ .x %>%
                                     dplyr::select(-any_of(c(
                                                "P_CODE",
                                                "Model",
                                                "value",
                                                "baseline",
                                                "DEPTH.f"))) %>%
                                     group_by(fYEAR, REEF, REEF.d, BIOREGION.agg, Metric) %>%
                                     tidybayes::summarise_draws(median, mean, sd,
                                                     HDInterval::hdi,
                                                     `p<0.5` = ~ mean(.x < 0.5)
                                                     )
                                 ),
                   Below = map(.x = Summary,
                               .f = ~ .x %>%
                                   ungroup() %>%
                                   mutate(Below = ifelse(upper < 0.5, 1, 0),
                                          PBelow = ifelse(`p<0.5` > 0.9, 1, 0)) %>%
                                   dplyr::select(fYEAR, Metric, Below, PBelow) %>%
                                   distinct()
                               )
                   ) %>% 
            suppressMessages() %>%
            suppressWarnings()

        save(mods,
              file = paste0(DATA_PATH, 'modelled/RPI__scores_reef_year.RData'))
}

