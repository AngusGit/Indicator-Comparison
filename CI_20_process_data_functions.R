CI_clear_processed_data <- function() {
    ## Clear the /data/processed folder
    unlink("../data/processed/*.*", recursive = TRUE)
}

CI_process_add_report_year <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'add_year',
                   label = "Add report year", status = 'pending')
    CI_tryCatch({
        load(paste0(PRIMARY_DATA_PATH, 'sample.reef.RData'))

        sample.reef.report.year<-sample.reef %>%
            mutate(REPORT_YEAR = ifelse(P_CODE == "RM", VISIT_NO+1992,
                            ifelse(P_CODE == 'IN', VISIT_NO+2004, year(Date)))) %>%
            mutate(REPORT_YEAR = ifelse(P_CODE == 'AP' & VISIT_NO %in% c('2','4'), NA, REPORT_YEAR)) %>%
            #filter(REPORT_YEAR < (CI$setting$FINAL_YEAR + 1))
          filter(REPORT_YEAR < 2025)

        save(sample.reef.report.year,
             file = paste0(DATA_PATH, 'processed/sample.reef.report.year.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'add_year',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Add report year.'), return=NULL)
}

# define juvenile transect area for each project
CI_process_transect_area <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'transect_area',
                   label = "Add transect area", status = 'pending')
    CI_tryCatch({

        juv.tran.area <- data.frame(P_CODE = c("RM","IN","GH","RR","AP"),
                                    tran.area = c(8.5,34,34,34,34)) 

        save(juv.tran.area,
             file = paste0(DATA_PATH, 'processed/juv.tran.area.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'transect_area',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Add transect area.'), return=NULL)
}

CI_process_locations_data <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'locations_data',
                   label = "Process locations data", status = 'pending')
    CI_tryCatch({

        load(paste0(DATA_PATH, "processed/sample.reef.report.year.RData"))

        locations <- sample.reef.report.year %>%
            dplyr::select(P_CODE, REEF, SITE_NO, DEPTH, LATITUDE, LONGITUDE) %>%
            mutate(name = paste(REEF, SITE_NO, sep=' ')) %>%
            group_by(REEF, SITE_NO, name) %>%
            summarise(LATITUDE = mean(LATITUDE),
                      LONGITUDE = mean(LONGITUDE))%>%
            mutate(SITE_NO = as.character(SITE_NO)) %>%
            ungroup() %>%
            suppressMessages() %>%
            suppressWarnings()

        save(locations,
             file = paste0(DATA_PATH, 'processed/locations.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'locations_data',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Process locations data.'), return=NULL)
}

CI_process_total_points <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'total_points',
                   label = "Process total points data", status = 'pending')
    CI_tryCatch({

        load(paste0(DATA_PATH, "primary/points.raw.RData"))
        load(paste0(DATA_PATH, "processed/sample.reef.report.year.RData"))

        ## Total points transect level
        total.points.aims <- points.raw %>%
           left_join(sample.reef.report.year) %>%   
           # filter(REPORT_YEAR < (CI$setting$FINAL_YEAR + 1)) %>%  
          filter(REPORT_YEAR < 2025) %>% 
            group_by(P_CODE, REEF, DEPTH, VISIT_NO,
                     REPORT_YEAR, SITE_NO, TRANSECT_NO) %>%
            mutate(SITE_NO = as.character(SITE_NO)) %>% # needed as jcu sites are named
            summarise(total.points = sum(POINTS)) %>%
            ungroup() %>%
            suppressMessages() %>%
            suppressWarnings()

        save(total.points.aims,
             file = paste0(DATA_PATH, 'processed/total.points.aims.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'total_points',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Process total points data.'), return=NULL)
}


CI_process_hc_points <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'hcc_points',
                   label = "Process HCC points data", status = 'pending')
    CI_tryCatch({

        load(paste0(DATA_PATH, "primary/video_codes.RData"))
        load(paste0(DATA_PATH, "primary/points.raw.RData"))
        load(file = paste0(DATA_PATH, 'processed/total.points.aims.RData'))

        ## hard coral points                            
        hc.aims <- points.raw %>%
            left_join(video_codes %>%
                      dplyr::select(GROUP_CODE, VIDEO_CODE)) %>%
            filter(GROUP_CODE == 'HC')%>%
            mutate(SITE_NO = as.character(SITE_NO)) %>%
            group_by(P_CODE, REEF, DEPTH, VISIT_NO,
                     SITE_NO, TRANSECT_NO) %>%
            summarise(HC = sum(POINTS)) %>%
            ungroup() %>% 
            right_join(total.points.aims %>%
                       dplyr::select(-total.points)) %>%
            mutate(HC = ifelse(is.na(HC), 0, HC)) %>%
            suppressMessages() %>%
            suppressWarnings()

        save(hc.aims,
             file = paste0(DATA_PATH, 'processed/hc.aims.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'hcc_points',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Process HCC points data.'), return=NULL)
}


CI_process_ma_points <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'ma_points',
                   label = "Process MA points data", status = 'pending')
    CI_tryCatch({

        load(paste0(DATA_PATH, "primary/video_codes.RData"))
        load(paste0(DATA_PATH, "primary/points.raw.RData"))
        load(file = paste0(DATA_PATH, 'processed/total.points.aims.RData'))

        ma.aims <- points.raw %>%
            left_join(video_codes %>%
                      dplyr::select(GROUP_CODE, BENTHOS_CODE, COMP_2021, VIDEO_CODE)) %>%
            filter(BENTHOS_CODE == 'MA')%>%
            mutate(SITE_NO = as.character(SITE_NO)) %>%
            group_by(P_CODE, REEF, DEPTH, VISIT_NO,
                     SITE_NO, TRANSECT_NO) %>%
            summarise(MA = sum(POINTS)) %>%
            ungroup() %>% 
            right_join(total.points.aims %>%
                       dplyr::select(-total.points)) %>%
            mutate(MA = ifelse(is.na(MA), 0, MA)) %>%
            suppressMessages() %>%
            suppressWarnings()

        save(ma.aims,
             file = paste0(DATA_PATH, 'processed/ma.aims.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'ma_points',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Process MA points data.'), return=NULL)
}

CI_process_a_points <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'a_points',
                   label = "Process A points data", status = 'pending')
    CI_tryCatch({

        load(paste0(DATA_PATH, "primary/video_codes.RData"))
        load(paste0(DATA_PATH, "primary/points.raw.RData"))
        load(file = paste0(DATA_PATH, 'processed/total.points.aims.RData'))

        a.aims <- points.raw %>%
            left_join(video_codes %>%
                      dplyr::select(GROUP_CODE, BENTHOS_CODE, VIDEO_CODE)) %>%
            filter(GROUP_CODE == 'A' | BENTHOS_CODE == 'ST')%>%
            mutate(SITE_NO = as.character(SITE_NO)) %>%
            group_by(P_CODE, REEF, DEPTH, VISIT_NO,
                     SITE_NO, TRANSECT_NO) %>%
            summarise(A = sum(POINTS)) %>%
            ungroup() %>% 
            right_join(total.points.aims %>%
                       dplyr::select(-total.points)) %>%
            mutate(A = ifelse(is.na(A), 0, A)) %>%
            suppressMessages() %>%
            suppressWarnings()

        save(a.aims,
             file = paste0(DATA_PATH, 'processed/a.aims.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'a_points',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Process A points data.'), return=NULL)
}

##Points response variables
CI_process_combine_points <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'combine_points',
                   label = "Combine processed points", status = 'pending')
    CI_tryCatch({

        load(paste0(DATA_PATH, "processed/hc.aims.RData"))
        load(paste0(DATA_PATH, "processed/ma.aims.RData"))
        load(paste0(DATA_PATH, "processed/a.aims.RData"))
        load(paste0(DATA_PATH, "processed/total.points.aims.RData"))
        load(paste0(DATA_PATH, "processed/locations.RData"))

        points.analysis.data.transect.aims <- hc.aims %>% 
            left_join(locations %>%
                      dplyr::select(REEF, SITE_NO, LATITUDE, LONGITUDE))%>%
            left_join(a.aims) %>%
            left_join(ma.aims) %>%
            left_join(total.points.aims) %>%
            mutate(REEF = as.factor(REEF),
                   reef.site = as.factor(paste0(REEF, '.', SITE_NO)),
                   reef.site.tran = as.factor(paste0(reef.site, '.', TRANSECT_NO)),
                   Project = as.factor(ifelse(P_CODE == 'RM', 'LTMP', 'MMP'))) %>%
            suppressMessages() %>%
            suppressWarnings()

        save(points.analysis.data.transect.aims,
             file = paste0(DATA_PATH, 'processed/points.analysis.data.transect.aims.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'combine_points',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Combine processed points.'), return=NULL)
}

CI_load_and_parse_file <- function(file) {
    df <- get(load(file))
    if(!assertthat::are_equal(colnames(df),
                          c("P_CODE", "REEF", "DEPTH", "VISIT_NO",      
                            "SITE_NO", "TRANSECT_NO", "HC", "REPORT_YEAR",
                            "LATITUDE", "LONGITUDE", "A", "MA",
                            "total.points", "reef.site", "reef.site.tran",
                            "Project"))) {
        
       # CI_log('ERROR', LOG_FILE, Category = '--Data processing--',
               # msg = paste0('Format of ', basename(file),
               #              ' file does not match expectations'))
        return(NULL)
    }
    df
}

CI_process_external_data <- function() {
    ## Will want to alter this function to read in csv versions of the data
    ## this will be easier to gaurentee that the data are compatible with
    ## the current version of R.  RData files can be created via either `save()`
    ## or `saveRDS()` and these are different formats.
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'external_data',
                   label = "Obtain external data", status = 'pending')
    CI_tryCatch({
        files <- list.files(path = paste0(DATA_PATH, 'external'),
                            pattern = 'points.analysis.data.transect.*.RData',
                            full.names = TRUE)
        points.analysis.data.transect.external <- do.call('rbind',
                                                         lapply(files, CI_load_and_parse_file))
    
        save(points.analysis.data.transect.external,
             file = paste0(DATA_PATH, 'processed/points.analysis.data.transect.external.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'external_data',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Obtain external data'), return=NULL)
}

CI_process_combine_data <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'combine_sources',
                   label = "Combine all sources", status = 'pending')
    CI_tryCatch({

        files <- list.files(path = paste0(DATA_PATH, 'processed'),
                            pattern = 'points.analysis.data.transect\\..*\\.RData',
                            full.names = TRUE)
        points.analysis.data.transect <- do.call('rbind',
                                                         lapply(files, CI_load_and_parse_file))

        save(points.analysis.data.transect,
             file = paste0(DATA_PATH, 'processed/points.analysis.data.transect.RData'))

        points.analysis.data <- points.analysis.data.transect %>%
            dplyr::select(-reef.site.tran, -TRANSECT_NO) %>%
            group_by(P_CODE, REEF, DEPTH, LATITUDE, LONGITUDE, VISIT_NO,
                     SITE_NO, REPORT_YEAR, Project, reef.site) %>%
            summarise(across(c("HC","MA","A","total.points"),sum)) %>%
            ungroup()

        save(points.analysis.data,
             file = paste0(DATA_PATH, 'processed/points.analysis.data.RData'))
        
        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'combine_sources',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Combine all sources'), return=NULL)
}


CI_process_spatial <- function() {
    ## This is for the shapefiles stored in this repo, for others
    ## set CI_get_spatial_data() in CI_10_get_data_functions.R
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'process_spatial',
                   label = "Process spatial", status = 'pending')
    CI_tryCatch({

        ## Bioregions
        bregions.sf <- read_sf(paste0(DATA_PATH, "spatial/bioregions/Marine_Bioregions_of_the_Great_Barrier_Reef__Reef_.shp")) %>%
                            st_transform(crs=4326) %>%
                            st_make_valid()
        save(bregions.sf,
             file = paste0(DATA_PATH, 'processed/bregions.sf.RData'))

        ## NRM spatial
        nrm.sf  <- read_sf(paste0(DATA_PATH, "spatial/NRM Regions/NRM_MarineRegions.shp")) %>% # this is not included in the CI_get_spatial_data function, or anywhere else I can see?
            st_transform(crs=4326)%>%
            st_make_valid()
        save(nrm.sf,
             file = paste0(DATA_PATH, 'primary/nrm.sf.RData'))
        
        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'process_spatial',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Process spatial'), return=NULL)
}


CI_process_assign_regions <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'assign_regions',
                   label = "Assign regions", status = 'pending')
    CI_tryCatch({
        
        load(paste0(DATA_PATH, 'processed/points.analysis.data.transect.RData'))
        load(paste0(DATA_PATH, 'processed/bregions.sf.RData'))
        load(paste0(DATA_PATH, 'primary/nrm.sf.RData'))   
        load(paste0(DATA_PATH, 'primary/tumra.RData'))
        load(paste0(DATA_PATH, 'primary/gbrmpa.management.RData'))
        load(paste0(DATA_PATH, 'primary/gbrmpa.RData'))
        load(paste0(DATA_PATH, 'primary/zones.RData'))

        bioregions <- points.analysis.data.transect %>%
            dplyr::select(P_CODE, REEF, SITE_NO, DEPTH, LATITUDE, LONGITUDE) %>%
            distinct() %>%
            rename(Longitude = LONGITUDE,
                   Latitude = LATITUDE) %>%
            mutate(SITE_NO = as.character(SITE_NO)) %>%
            st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326) %>%
            st_join(., bregions.sf, join = st_nearest_feature)%>%
            st_drop_geometry()  %>%
            dplyr::select(REEF, DEPTH, SITE_NO, BIOREGION) %>% 
            group_by(REEF, DEPTH) %>%
            filter(SITE_NO == first(SITE_NO)) %>%
            droplevels() %>%
            dplyr::select(-SITE_NO) %>%
            ungroup()

        ## NRM spatial
        nrm <- points.analysis.data.transect %>%
            rename(Latitude = LATITUDE, Longitude = LONGITUDE) %>%
            st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326) %>%
            st_join(., nrm.sf, join = st_within)%>%            # nrm.sf not yet created
            st_drop_geometry()%>%
            rename(NRM=NAME) %>%
            dplyr::select(REEF,DEPTH, NRM) %>%
            distinct()


        ## TUMRA
        tumra <- points.analysis.data.transect %>%
            rename(Latitude = LATITUDE, Longitude = LONGITUDE) %>%
            st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326) %>%
            st_join(., tumra, join = st_within)%>%
            st_drop_geometry()%>%
            rename(TUMRA = Name) %>%
            dplyr::select(REEF, DEPTH, TUMRA) %>% 
            distinct() %>%
            suppressMessages() %>%
            suppressWarnings()

        ## GBRMPA management areas
        gbrmpa.ma <- points.analysis.data.transect %>%
            rename(Latitude = LATITUDE, Longitude = LONGITUDE) %>%
            st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326) %>%
            ## st_join(., gbrmpa.management, join = st_within)%>%
            st_join(., gbrmpa.management, join = st_nearest_feature)%>%
            st_drop_geometry()%>%
            rename(GBRMPA.MA = Name) %>%
            dplyr::select(REEF, DEPTH, GBRMPA.MA) %>% 
            distinct() %>%
            suppressMessages() %>%
            suppressWarnings()
        
        ## Whole GBRMP 
        gbrmp <- points.analysis.data.transect %>%
            rename(Latitude = LATITUDE, Longitude = LONGITUDE) %>%
            st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326) %>%
            ## st_join(., gbrmpa, join = st_within)%>%
            st_join(., gbrmpa, join = st_nearest_feature)%>%
            st_drop_geometry()%>%
            rename(GBRMP = Name) %>%
            dplyr::select(REEF, DEPTH, GBRMP) %>% 
            distinct() %>%
            suppressMessages() %>%
            suppressWarnings()

        ## Zones (Northern, Central, Southern)
        zones <- points.analysis.data.transect %>%
            rename(Latitude = LATITUDE, Longitude = LONGITUDE) %>%
            st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326) %>%
            ## st_join(., zones, join = st_within)%>%
            st_join(., zones, join = st_nearest_feature)%>%
            st_drop_geometry()%>%
            rename(ZONE = Name) %>%
            dplyr::select(REEF, DEPTH, ZONE) %>% 
            distinct() %>%
            suppressMessages() %>%
            suppressWarnings()

        ## Latitude and Longitude
        latlong <- points.analysis.data.transect %>%
            rename(Latitude = LATITUDE, Longitude = LONGITUDE) %>%
            dplyr::select(REEF, DEPTH, Latitude, Longitude) %>%
            distinct() %>%
            group_by(REEF, DEPTH) %>%
            summarise(Latitude = mean(Latitude),
                   Longitude = mean(Longitude)) %>%
            ungroup() %>% 
            suppressWarnings() %>%
            suppressMessages()
            
        mmpReefs<- c("Snapper North","Snapper South","Low Isles",
                     "Fitzroy East","Fitzroy West","Fitzroy Island", "Franklands East","Franklands West","High East","High West",
                     "Barnards","Dunk North", "Dunk South","Bedarra",
                     "Palms East","Palms West", "Lady Elliot","Pandora Reef", "Pandora", "Havannah Island", "Havannah", "Magnetic",
                     "Daydream","Dent","Double Cone", "Hook","Pine","Seaforth","Shute Harbour","Border Island","Hayman Island","Langford and Bird Isles",
                     "North Keppel","Middle","Keppels South", "Barren","Pelican")
        ## Join REEF, BIOREGION, NRM, Lat/longs and save
        spatial_lookup <- bioregions %>% 
            left_join(nrm) %>% 
            left_join(tumra) %>%
            left_join(gbrmpa.ma) %>%
            left_join(gbrmp) %>%
            left_join(zones) %>%
            left_join(latlong) %>% 
            distinct() %>%
            mutate(BIOREGION = as.character(BIOREGION),
                   BIOREGION.agg = as.factor(case_when(BIOREGION %in% c("4", "3") ~"4:3",
                                                       BIOREGION %in% c("35", "36") ~"35:36",
                                                       !BIOREGION %in% c("4", "3", "35", "36")
                                                       ~BIOREGION)),
                   DEPTH.f = factor(case_when(DEPTH >3 ~"deep slope",
                                              DEPTH <= 3 ~"shallow slope")),
                   Shelf = ifelse(BIOREGION %in% c('16','17','18','19','22','23','29'),
                                  'Inshore','Offshore'),
                   REEF.d = factor(paste(REEF, DEPTH.f))) %>%
                distinct() %>%
          mutate(MMP.Report=ifelse(REEF %in% mmpReefs, TRUE, FALSE)) |> 
            suppressWarnings() %>%
            suppressMessages()
## adjust Lady Elliot to reflect MMP reporting as Burdekin
        spatial_lookup <- spatial_lookup |> 
          mutate(NRM=ifelse(REEF=="Lady Elliot", "Burdekin", NRM))
        
        save(spatial_lookup,
             file = paste0(DATA_PATH, 'processed/spatial_lookup.RData'))
        
        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'assign_regions',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Assign regions'), return=NULL)
}

# AT edits to impose zeros to juvenile data

CI_process_juvenile_density <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'juvenile_density',
                   label = "Juvenile density", status = 'pending')
    CI_tryCatch({
        
        load(paste0(DATA_PATH, 'primary/juveniles.RData'))
      load(paste0(DATA_PATH, 'primary/juvenile.samples.RData'))

        ## convert to density per m2 of transect
        juv.analysis.data <- juveniles %>%
          right_join(juvenile.samples) %>%
          mutate(GENUS=ifelse(is.na(GENUS), 'Acropora', GENUS),
                 ABUNDANCE=ifelse(is.na(ABUNDANCE),0,ABUNDANCE)) |>
            group_by(P_CODE, REEF, DEPTH, VISIT_NO, SITE_NO) %>%
            summarise(ABUNDANCE = sum(ABUNDANCE)) %>%
            ungroup() %>%
          mutate(GENUS='total.juv') %>%
          rbind(juveniles) %>% 
          pivot_wider(names_from = GENUS,
                      values_from = ABUNDANCE,
                      values_fill = 0) 

        save(juv.analysis.data,
             file = paste0(DATA_PATH, 'processed/juv.analysis.data.RData'))

        
        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'juvenile_density',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Juvenile density'), return=NULL)
}

CI_process_juvenile_offset <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'juvenile_offset',  # corrected spelling of juvinile
                   label = "Juvenile offset", status = 'pending')
    CI_tryCatch({
        
        load(paste0(DATA_PATH, 'processed/juv.tran.area.RData'))
        #mmp.points <- read.csv(paste0(DATA_PATH, 'primary/mmp.points.csv'), as.is=TRUE)
      load(paste0(PRIMARY_DATA_PATH, 'points.raw.RData'))
        mmp.points <- points.raw |> filter(P_CODE %in% c("IN", "AP", "GH", "RR"))
        video_codes <- read.csv(paste0(DATA_PATH, 'primary/video_codes.csv'), as.is=TRUE)
        
        ## MMP
        mmp.points.site <- mmp.points %>%
            left_join(video_codes) %>%
            filter(!BENTHOS_CODE == 'IN') %>% #remove indeterminate points
            group_by(P_CODE, REEF, DEPTH, VISIT_NO, SITE_NO) %>%
            summarise(n.point = sum(POINTS)) %>%
            ungroup() %>% 
          suppressMessages() %>%
          suppressWarnings()

        mmp.cover.algae.site <- mmp.points %>%
            left_join(video_codes) %>%
            filter(BENTHOS_CODE %in% c('AO','CA','MA','TA','ST')) %>% 
            droplevels() %>%
            group_by(P_CODE, REEF, DEPTH, VISIT_NO, SITE_NO, BENTHOS_CODE) %>%
            summarise(points = sum(POINTS)) %>%
            ungroup() %>%
            pivot_wider(names_from = BENTHOS_CODE,
                        values_from = points,
                        values_fill = list(points = 0)) %>%
            right_join(mmp.points.site) %>%
            left_join(juv.tran.area) %>%
            mutate(avail.sub = (MA+ST+TA+CA+AO)/n.point,
                   avail.area=tran.area*avail.sub) %>% 
            suppressMessages() %>%
            suppressWarnings()

        ## LTMP
       # ltmp.points <- read.csv(paste0(DATA_PATH, 'primary/ltmp.points.csv'), as.is=TRUE)
        load(paste0(PRIMARY_DATA_PATH, 'points.raw.RData'))
       ltmp.points <- points.raw |> filter(P_CODE %in% c("RM")) # note that RMRAP and RAP have been converted to RM in this dataframe
        ltmp.points.site <- ltmp.points %>%
            left_join(video_codes) %>%
            filter(!BENTHOS_CODE == 'IN' & VISIT_NO > 14) %>% #remove indeterminate points and pre juvenile samples
            group_by(P_CODE, REEF, VISIT_NO, SITE_NO) %>%
            summarise(n.point = sum(POINTS)) %>%
            ungroup() %>%
            suppressMessages() %>%
            suppressWarnings()

        ltmp.cover.algae.site <- ltmp.points %>%
            left_join(video_codes) %>%
            filter(BENTHOS_CODE %in% c('AO','CA','MA','TA','ST')) %>% 
            droplevels() %>%
            group_by(P_CODE, REEF, VISIT_NO, SITE_NO, BENTHOS_CODE) %>%
            summarise(points = sum(POINTS)) %>%
            ungroup() %>%
            pivot_wider(names_from = BENTHOS_CODE,
                        values_from = points,
                        values_fill = list(points = 0)) %>%
            right_join(ltmp.points.site) %>%
            left_join(juv.tran.area) %>%
            mutate(avail.sub = (MA+ST+TA+CA+AO)/n.point,
                   avail.area = tran.area*avail.sub,
                   DEPTH = as.integer(ifelse(REEF == 'Middle Reef', 2, 8))) %>%
            suppressMessages() %>%
            suppressWarnings()
        
        juvenile.offset <- ltmp.cover.algae.site %>% 
            rbind(mmp.cover.algae.site) %>%
            suppressMessages() %>%
            suppressWarnings()

        save(juvenile.offset,
             file = paste0(DATA_PATH, 'processed/juvenile.offset.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'juvenile_offset',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Juvenile offset'), return=NULL)
}

CI_process_juvenile_data_for_baselines <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'juvenile_baseline',
                   label = "Juvenile baseline data", status = 'pending')
    CI_tryCatch({
        
        load(paste0(DATA_PATH, 'processed/juv.analysis.data.RData'))
        load(paste0(DATA_PATH, 'processed/juvenile.offset.RData'))
        load(paste0(DATA_PATH, 'processed/points.analysis.data.RData'))
        load(paste0(DATA_PATH, 'processed/sample.reef.report.year.RData'))
        load(paste0(DATA_PATH, 'processed/spatial_lookup.RData'))

        ## Prepare data for baseline model
        juv.df <- juv.analysis.data %>%
            mutate(total.juv.excl.turb = total.juv - Turbinaria,
                   Non_Acropora = total.juv - Acropora, 
                   Non_Acr_exTurb = Non_Acropora - Turbinaria,
                   SITE_NO = as.factor(as.character(SITE_NO))) %>%
            left_join(sample.reef.report.year %>%
                      mutate(SITE_NO = as.factor(as.character(SITE_NO))) %>%
                      select(-Date),
                      by = c("P_CODE", "REEF", "DEPTH", "VISIT_NO", "SITE_NO")) %>%
            filter(!is.na(REPORT_YEAR) & REPORT_YEAR>2006) %>%
            left_join(spatial_lookup %>%
                      dplyr::select(REEF, DEPTH, DEPTH.f, REEF.d, BIOREGION, BIOREGION.agg, NRM) %>%
                      distinct()
                      ) %>%
            ## left_join(spatial_lookup %>%
            ##           dplyr::select(REEF, DEPTH, BIOREGION, NRM) %>%
            ##           distinct()) %>%
            mutate(P_CODE = as.factor(P_CODE),
                   NRM = as.factor(NRM),
                   #DEPTH.f = factor(case_when(DEPTH > 3 ~"deep slope",    # these fields are already factors in spatial_lookup
                   #               DEPTH <= 3 ~"shallow slope")),
                   #REEF.d = factor(paste(REEF, DEPTH.f)),
                   #REPORT_YEAR = as.numeric(as.character(REPORT_YEAR)), # already a numeric field in "sample.reef.report.year"
                   #BIOREGION = as.character(BIOREGION),
                   #BIOREGION.agg = as.factor(case_when(BIOREGION %in% c("4", "3") ~"4:3",
                    #                                   BIOREGION %in% c("35", "36") ~"35:36",
                     #                                  !BIOREGION %in% c("4", "3", "35", "36")
                      #                                 ~BIOREGION))
                   ) %>%
            left_join(juvenile.offset %>%
                      mutate(SITE_NO = as.factor(SITE_NO))) %>%
            filter(DEPTH < 9.1) %>% 
            mutate(fYEAR = factor(REPORT_YEAR),
                   Site = factor(paste0(REEF.d, SITE_NO))) %>%
            dplyr::select(P_CODE, NRM, BIOREGION, BIOREGION.agg, REEF, DEPTH, DEPTH.f,
                          REEF.d, VISIT_NO, REPORT_YEAR, fYEAR, SITE_NO, Site,
                          LATITUDE, LONGITUDE, total.juv, total.juv.excl.turb,
                          Turbinaria, Non_Acropora, Non_Acr_exTurb,
                          Acropora, avail.area) %>% #Select relevant variables
            filter(!is.na(avail.area)) %>%
            droplevels() %>% 
            suppressMessages() %>%
            suppressWarnings()

        save(juv.df,
             file = paste0(DATA_PATH, 'processed/juv.df.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'juvenile_baseline',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Juvenile baseline data'), return=NULL)
}



CI_process_unique_site_location <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'unique_locations',
                   label = "Unique site locations", status = 'pending')
    CI_tryCatch({
        
        ## create a unique identifier for each P_CODE,REEF,DEPTH,SITE
        ## these will need to link back to the baseline estimates at
        ## the level of P_CODE,REEF,DEPTH ultimately there will need
        ## to be a lookup table that has fields for Bioregion, nrm,
        ## other regional aggregations, Latitude, longitude, reef,
        ## site, depth, reef.site.depth, year the above is the level
        ## we intend to estimate index score distributions that are
        ## then aggregated to regional summaries.
        load(file=paste0(DATA_PATH, 'processed/points.analysis.data.RData'))
        load(paste0(DATA_PATH, 'processed/spatial_lookup.RData'))

        site.location <- points.analysis.data %>% 
            dplyr::select(P_CODE, REEF, DEPTH, SITE_NO, LATITUDE, LONGITUDE) %>%
            distinct() %>%
            left_join(spatial_lookup %>%
                      dplyr::select(REEF, DEPTH, DEPTH.f, REEF.d, BIOREGION, BIOREGION.agg, NRM) %>%
                      distinct()
                      ) %>%
            mutate(
                BIOREGION = as.character(BIOREGION),
                   BIOREGION.agg = as.factor(case_when(BIOREGION %in% c("4", "3") ~"4:3",
                                                       BIOREGION %in% c("35", "36") ~"35:36",
                                                       !BIOREGION %in% c("4", "3", "35", "36")
                                                       ~BIOREGION)),
                   DEPTH.f = factor(case_when(DEPTH >3 ~"deep slope",
                                              DEPTH <= 3 ~"shallow slope")),
                   Site = factor(paste(REEF, DEPTH.f, SITE_NO)),
                   Shelf = ifelse(BIOREGION %in% c('16','17','18','19','22','23','29'),
                                  'Inshore','Offshore'),
                   REEF.d = factor(paste(REEF, DEPTH.f))) %>%
            distinct() %>%
            suppressMessages() %>%
            suppressWarnings()

        save(site.location,
             file = paste0(DATA_PATH, 'processed/site.location.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'unique_locations',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Unique site locations'), return=NULL)
}

CI_process_composition_points <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'comp.points',
                   label = "Composition points", status = 'pending')
    CI_tryCatch({
        load(file = paste0(DATA_PATH, 'primary/points.raw.RData')) 
        codes <- read.csv(paste0(DATA_PATH, 'primary/video_codes.csv'))

        points.zeros <- points.raw%>% 
            pivot_wider(names_from = VIDEO_CODE,
                        values_from = POINTS,
                        values_fill = list(POINTS = 0)) %>%
            pivot_longer(cols = num_range(prefix = "", range = 1:997),
                         names_to = "VIDEO_CODE",
                         values_to = "n.points") %>%
            mutate(VIDEO_CODE = as.integer(VIDEO_CODE)) %>%
            left_join(codes) %>%
            filter(!GROUP_CODE == 'IN') %>%
            suppressMessages() %>%
            suppressWarnings()

        ## total points as 
        total.points.transect <- points.zeros %>%
            group_by(REEF, DEPTH, VISIT_NO, SITE_NO, TRANSECT_NO) %>%
            summarise(total.points = sum(n.points)) %>%
            ungroup() %>% 
            suppressMessages() %>%
            suppressWarnings()

        save(points.zeros, total.points.transect,
             file = paste0(DATA_PATH, 'processed/comp.points.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'comp.points',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Composition points'), return=NULL)
}


CI_process_composition_ordi_data <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'ordi.data',
                   label = "Ordination data", status = 'pending')
    CI_tryCatch({
        
        load(file = paste0(DATA_PATH, 'processed/comp.points.RData'))
        codes <- read.csv(paste0(DATA_PATH, 'primary/video_codes.csv'))
        load(file = paste0(DATA_PATH, 'processed/sample.reef.report.year.RData'))

        ## Community composition for ordination
        ordiData <- points.zeros %>% 
            left_join(codes %>%
                      dplyr::select(VIDEO_CODE, GROUP_CODE, COMP_2021)) %>%
            filter(GROUP_CODE %in% c('HC','SC')) %>%
            group_by(P_CODE, REEF, DEPTH, VISIT_NO, SITE_NO, TRANSECT_NO, COMP_2021) %>%
            summarise(n.points = sum(n.points)) %>%
            ungroup() %>%
            left_join(total.points.transect) %>%
            mutate(cover = 100*(n.points/total.points)) %>%
            group_by(P_CODE, REEF, DEPTH, VISIT_NO, SITE_NO, COMP_2021) %>%
            summarise(cover=mean(cover))%>%
            ungroup() %>%
            pivot_wider(id_cols = c(P_CODE, REEF,DEPTH, VISIT_NO, SITE_NO),
                        names_from = COMP_2021,
                        values_from = cover) %>%
            left_join(sample.reef.report.year %>%
                      dplyr::select(P_CODE, REEF, DEPTH, VISIT_NO, SITE_NO, REPORT_YEAR)) %>%
           # filter(REPORT_YEAR < (CI$setting$FINAL_YEAR + 1)) %>% 
          filter(REPORT_YEAR < (FINAL_YEAR + 1)) %>% 
            suppressMessages() %>%
            suppressWarnings()
  
        save(ordiData, file = paste0(DATA_PATH, 'processed/ordiData.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'ordi.data',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Ordination data'), return=NULL)
}

CI_process_composition_data <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'comp.data',
                   label = "Composition data", status = 'pending')
    CI_tryCatch({
        
        load(file = paste0(DATA_PATH, 'processed/comp.points.RData'))
        codes <- read.csv(paste0(DATA_PATH, 'primary/video_codes.csv'))
        load(file = paste0(DATA_PATH, 'processed/sample.reef.report.year.RData'))

        comp2021.transect.temp <- points.zeros %>% 
            left_join(codes %>%
                      dplyr::select(VIDEO_CODE, GROUP_CODE, COMP_2021)) %>%
            group_by(P_CODE, REEF, DEPTH, VISIT_NO,
                     SITE_NO, TRANSECT_NO, GROUP_CODE,COMP_2021) %>%
            summarise(n.points = sum(n.points)) %>%
            ungroup() %>% 
            left_join(total.points.transect) %>%
            mutate(cover = 100*(n.points/total.points)) %>%
            left_join(sample.reef.report.year %>%
                      dplyr::select(P_CODE, REEF, DEPTH, VISIT_NO,
                                    SITE_NO, REPORT_YEAR, Date)) %>%
            #filter(REPORT_YEAR < (CI$setting$FINAL_YEAR + 1)) %>% 
          filter(REPORT_YEAR < (FINAL_YEAR + 1)) %>% 
            suppressMessages() %>%
            suppressWarnings()

        ## remove times when COR* codes were used to excess greater
        ## than 10% of the coral cover.
        Cor <- comp2021.transect.temp  %>%
            filter(GROUP_CODE %in% c('HC','SC')) %>%
            mutate(cor = case_when(COMP_2021 %in%
                                    c('COR_CBCF','COR_CE','COR_CL','COR_CMCS') ~ "1"),
                   Cor = ifelse(is.na(cor), "other", "COR")) %>%
            group_by(REEF, DEPTH, VISIT_NO, Cor) %>%
            summarise(points = sum(n.points)) %>%
            ungroup() %>%
            pivot_wider(names_from = Cor,
                        values_from = points) %>%
            mutate(cor.prop = COR/other) %>%
            filter(cor.prop < 0.1) %>% 
            dplyr::select(REEF, DEPTH, VISIT_NO)

        comp.transect<-comp2021.transect.temp %>%
            right_join(Cor) %>%
            filter (GROUP_CODE %in% c('OT','SC','HC') &
                    !COMP_2021 %in% c('IN','OT','SC_OTH'))

        save(comp.transect, file = paste0(DATA_PATH, 'processed/comp.transect.RData'))
        
        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'comp.data',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Composition data'), return=NULL)
}


CI__retrieve_data_from_metadata <- function(file) {
    in_file <- paste0("https://data.aims.gov.au/data-download/55f56f2c-c9fd-4af7-90f9-617643975f9e/",
                      file)
    out_file <- paste0(DATA_PATH, "parameters/", file)
    if (!file.exists(file)) download.file(in_file, out_file) 
}

