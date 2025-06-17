######################################################################
## Notes made by Angus                                              ##
##                                                                  ##
## For all mmp data use mmp_site_name as REEF                       ##
## For al LTMP data use AIMS_REEF_NAME as REEF                      ##
######################################################################

CI_clear_primary_data <- function() {
    ## Clear the /data/primary folder
    unlink("../data/primary/*.csv")
    unlink("../data/primary/*.sql")
    unlink("../data/primary/*.zip")
    unlink("../data/primary/*.RData")
    unlink("../data/primary/*.geojson")
    unlink("../data/primary/*.json")
    unlink("../data/primary/GIS/*.*", recursive = TRUE)
}


###############################################################################################
## Get data functions
###############################################################################################
CI_get_names_lookups <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'names_lookup',
                   label = "Names lookup", status = 'pending')
    CI_tryCatch({

        ## MMP
        writeLines("select p_code, min(lat_dd), min(long_dd), nrm_region, shelf,catchment,fullreef_id,reef_id,reef_name,mmp_site_name as reef,AIMS_REEF_NAME,reef_zone         
          from v_in_sample
           where p_code in('IN','AP','RR','GH')
            and sample_type= 'PPOINT'
            and visit_no >0
            and mmp_site_name not like 'Cape%'
            group by p_code,nrm_region, shelf,catchment,fullreef_id,reef_id,reef_name,mmp_site_name,AIMS_REEF_NAME,reef_zone",
          paste0(PRIMARY_DATA_PATH, "name.lookup.sql"))

        system(paste0("java -jar dbExport.jar ", PRIMARY_DATA_PATH, "name.lookup.sql ",
                      PRIMARY_DATA_PATH, "mmp.name.lookup.csv reef reefmon"),
               ignore.stdout = TRUE )

        mmp_names_lookup <- read.csv(paste0(PRIMARY_DATA_PATH, "mmp.name.lookup.csv"))

        ## LTMP
        writeLines("select p_code, min(lat_dd), min(long_dd), nrm_region, shelf,catchment,fullreef_id,reef_id,reef_name,mmp_site_name as reef,AIMS_REEF_NAME,reef_zone         
          from v_in_sample
           where p_code in('RM','RAP','RMRAP')
            and sample_type= 'PPOINT'
            and visit_no >0
            and mmp_site_name not like 'Cape%'
            group by p_code,nrm_region, shelf,catchment,fullreef_id,reef_id,reef_name,mmp_site_name,AIMS_REEF_NAME,reef_zone",
          paste0(PRIMARY_DATA_PATH, "name.lookup.sql"))

        system(paste0("java -jar dbExport.jar ", PRIMARY_DATA_PATH, "name.lookup.sql ",
                      PRIMARY_DATA_PATH, "name.lookup.csv reef reefmon"),
               ignore.stdout = TRUE)
        ltmp_names_lookup<-read.csv(paste0(PRIMARY_DATA_PATH, "name.lookup.csv"))

        ## Combine
        names.lookup<-mmp_names_lookup %>% 
            dplyr::select(AIMS_REEF_NAME, REEF_ZONE,REEF,NRM_REGION) %>% #AT added NRM_REGION as not extracted elsewhere
            full_join(ltmp_names_lookup %>%
                      dplyr::select(AIMS_REEF_NAME, REEF,REEF_ZONE,NRM_REGION)) %>% 
            distinct() %>%
            suppressMessages() %>%
            suppressWarnings()

        save(names.lookup, file = paste0(DATA_PATH, "primary/names.lookup.RData"))
            
        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'names_lookup',status = 'success')
    }, logFile=LOG_FILE, Category='--Data extraction--',
    msg=paste0('names lookups.'), return=NULL)

}

CI_get_mmp_samples_data <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'mmp_samples',
                   label = "MMP samples", status = 'pending')
    CI_tryCatch({
        writeLines("select p_code, min(lat_dd), min(long_dd), mmp_site_name as reef,site_no,depth,
          visit_no, year_code as report_year, to_char(min(sample_date),'DD/MM/YYYY') as sample_date      
          from v_in_sample
           where p_code in ('IN','AP','RR','GH')
            and visit_no >0
            group by p_code,mmp_site_name,site_no,depth,visit_no, year_code",
          paste0(PRIMARY_DATA_PATH, "mmp.sample.sql"))

        system(paste0("java -jar dbExport.jar ", PRIMARY_DATA_PATH, "mmp.sample.sql ",
                      PRIMARY_DATA_PATH, "mmp.sample.csv reef reefmon"),
               ignore.stdout = TRUE)

        sample.reef.mmp  <- read.csv(paste0(PRIMARY_DATA_PATH, "mmp.sample.csv")) %>%
            mutate(LATITUDE = as.numeric(as.character(MIN.LAT_DD.)),
                   LONGITUDE = as.numeric(as.character(MIN.LONG_DD.)),
                   Date = as.Date(SAMPLE_DATE, format='%d/%m/%Y')
                   ) %>%
            dplyr::select(-MIN.LONG_DD.,-MIN.LAT_DD.,-SAMPLE_DATE)

        save(sample.reef.mmp, file = paste0(DATA_PATH, 'primary/sample.reef.mmp.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                          item = 'mmp_samples',status = 'success')
        
    }, logFile=LOG_FILE, Category='--Data extraction--',
    msg=paste0('MMP samples.'), return=NULL)
}

CI_get_ltmp_samples_data <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'ltmp_samples',
                   label = "LTMP samples", status = 'pending')
    CI_tryCatch({
        writeLines("select p_code, min(lat_dd), min(long_dd),aims_reef_name as reef, site_no,visit_no, ltmp_report_year as report_year, to_char(min(sample_date),'DD/MM/YYYY') as sample_date  
           from v_in_sample 
           where p_code in ('RM','RAP','RMRAP')
            and visit_no >0
            and sample_type in ('PPOINT', 'VPOINT','VPOINTNOB')
            group by p_code, aims_reef_name,site_no,visit_no, ltmp_report_year",
          paste0(PRIMARY_DATA_PATH, "ltmp.sample.sql"))

        system(paste0("java -jar dbExport.jar ", PRIMARY_DATA_PATH, "ltmp.sample.sql ",
                      PRIMARY_DATA_PATH, "ltmp.sample.csv reef reefmon"),
               ignore.stdout = TRUE)

        sample.reef.ltmp <- read.csv(paste0(PRIMARY_DATA_PATH, "ltmp.sample.csv")) %>%
            mutate(DEPTH = as.integer(ifelse(REEF == 'Middle Reef', 2, 8)),
                   LATITUDE = as.numeric(as.character(MIN.LAT_DD.)),
                   LONGITUDE = as.numeric(as.character(MIN.LONG_DD.)),
                   Date = as.Date(SAMPLE_DATE, format = '%d/%m/%Y'),
                   P_CODE = 'RM') %>%
            dplyr::select(-MIN.LONG_DD., -MIN.LAT_DD., -SAMPLE_DATE) 

        save(sample.reef.ltmp, file = paste0(DATA_PATH, 'primary/sample.reef.ltmp.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                          item = 'ltmp_samples',status = 'success')
        
    }, logFile=LOG_FILE, Category='--Data extraction--',
    msg=paste0('LTMP samples.'), return=NULL)
}

CI_combine_samples_data <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'combine_samples',
                   label = "Combine samples", status = 'pending')
    CI_tryCatch({
        load(file = paste0(DATA_PATH, 'primary/sample.reef.mmp.RData'))
        load(file = paste0(DATA_PATH, 'primary/sample.reef.ltmp.RData'))

        sample.reef <- sample.reef.mmp %>%
            rbind(sample.reef.ltmp) %>%
            mutate(P_CODE = factor(P_CODE),
                   REEF = factor(REEF),
                   DEPTH = as.integer(ifelse(P_CODE == 'RM', '8',
                                      ifelse(REEF == 'Middle Reef','2', DEPTH))))

        save(sample.reef, file = paste0(DATA_PATH, 'primary/sample.reef.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                          item = 'combine_samples',status = 'success')
        
    }, logFile=LOG_FILE, Category='--Data extraction--',
    msg=paste0('Combine samples.'), return=NULL)
}

###########################################################################
## Extract  disturbances data...    
## uses new disturbance table for 2021 "DISTURBANCES2"
###########################################################################

CI_get_disturbances_data <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'disturbances_samples',
                   label = "Disturbances data", status = 'pending')
    CI_tryCatch({
        writeLines("select aims_reef_name, depth, reef_zone,p_code, visit_no, rank, disturbance_type as disturbance, storm_name,to_char(distrubance_date,'DD/MM/YYYY') as disturbance_date
  from disturbances2
   where p_code in ('IN', 'RM','RMRAP','RAP','AP','GH','RR')
           and sample_type in ('PPOINT','VPOINT')",
          paste0(PRIMARY_DATA_PATH, "disturbance.sql"))

        system(paste0("java -jar dbExport.jar ", PRIMARY_DATA_PATH, "disturbance.sql ",
                      PRIMARY_DATA_PATH, "disturbance.csv reef disturbance"),
               ignore.stdout = TRUE)

        load(file = paste0(DATA_PATH, "primary/names.lookup.RData"))

        disturbance.reef <- read.csv(paste0(PRIMARY_DATA_PATH, "disturbance.csv")) %>%
            mutate(P_CODE = ifelse(P_CODE %in% c('RM', 'RMRAP', 'RAP'), 'RM', P_CODE)) %>%
            full_join(names.lookup) %>%
            mutate(REEF = ifelse(P_CODE == "RM", as.character(AIMS_REEF_NAME), REEF)) %>%
            dplyr::select(-AIMS_REEF_NAME, -REEF_ZONE) %>%
            mutate(DISTURBANCE_DATE = as.Date(DISTURBANCE_DATE, format = '%d/%m/%Y'),
                   DEPTH = as.integer(ifelse(P_CODE == "RM", '8',
                                      ifelse(REEF == 'Middle Reef', '2', DEPTH)))) %>%
            suppressMessages() %>%
            suppressWarnings()
  
                                 
        save(disturbance.reef, file = paste0(DATA_PATH, 'primary/disturbances.RData'))
        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                          item = 'disturbances_samples',status = 'success')
        
    }, logFile=LOG_FILE, Category='--Data extraction--',
    msg=paste0('Disturbances data.'), return=NULL)
}


##################################
## Extract video code lookup ##
##################################
CI_get_video_code_lookup <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'video_lookup',
                   label = "Video code lookup", status = 'pending')
    CI_tryCatch({

        writeLines("select video_code, b.benthos_code, b.benthos_desc, v.group_code, group_desc,genus_2021, family_2021, comp_2021, comp_2021_description, AI_CODE_2021
          from vpoint_codes v, benthos b
          where b.benthos_code=v.benthos_code",
          paste0(PRIMARY_DATA_PATH, "video_codes.sql"))

        system(paste0("java -jar dbExport.jar ", PRIMARY_DATA_PATH, "video_codes.sql ",
                      PRIMARY_DATA_PATH, "video_codes.csv reef reefmon"),
               ignore.stdout = TRUE)

        video_codes <- read.csv(paste0(PRIMARY_DATA_PATH, "video_codes.csv"))
                                 
        save(video_codes, file = paste0(DATA_PATH, 'primary/video_codes.RData'))
        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                          item = 'video_lookup',status = 'success')
        
    }, logFile=LOG_FILE, Category='--Data extraction--',
    msg=paste0('Video code lookup.'), return=NULL)
}


CI_get_mmp_juveniles_data <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'mmp_juv',
                   label = "MMP juveniles", status = 'pending')
    CI_tryCatch({

        writeLines("select p_code, mmp_site_name as reef, depth, visit_no, site_no,b.genus,sum(abundance) abundance
           from v_in_sample s, demog b, juvenile_taxa j
           where s.sample_id = b.sample_id
           and b.genus = j.taxa
           and p_code  in('IN','RR','AP','GH') 
           and visit_no >0 
           and size_class in ('000-002','002-005')
           and s.mmp_site_name not like 'Cape%'
           group by p_code, mmp_site_name, depth, visit_no, site_no, b.genus",
          paste0(PRIMARY_DATA_PATH, "juv.mmp.sql"))

        system(paste0("java -jar dbExport.jar ", PRIMARY_DATA_PATH, "juv.mmp.sql ",
                      PRIMARY_DATA_PATH, "juv.mmp.csv reef reefmon"),
               ignore.stdout = TRUE)

        juv.mmp <- read.csv(paste0(PRIMARY_DATA_PATH, 'juv.mmp.csv'), as.is = TRUE) %>%
            filter(SITE_NO != '3') # remove the only place with a site 3 in the mmp dataset

        save(juv.mmp, file = paste0(DATA_PATH, 'primary/juv.mmp.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                          item = 'mmp_juv',status = 'success')
        
    }, logFile=LOG_FILE, Category='--Data extraction--',
    msg=paste0('MMP juveniles.'), return=NULL)
}

CI_get_mmp_juveniles_samples <- function() {
  CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                 item = 'mmp_juv_samp',
                 label = "MMP juveniles samples", status = 'pending')
  CI_tryCatch({

writeLines("select p_code, mmp_site_name as reef, depth, visit_no, site_no
           from v_in_sample 
           where p_code in('IN','RR','AP','GH') 
           and visit_no >0 
           and sample_type = 'DEMOG'
           group by p_code, mmp_site_name, visit_no, site_no, depth",
           paste0(PRIMARY_DATA_PATH, "juv.mmp.samples.sql"))
    
    system(paste0("java -jar dbExport.jar ", PRIMARY_DATA_PATH, "juv.mmp.samples.sql ",
                  PRIMARY_DATA_PATH, "juv.mmp.samples.csv reef reefmon"),
           ignore.stdout = TRUE)
           
           juv.mmp.samples<-read.csv(paste0(PRIMARY_DATA_PATH, 'juv.mmp.samples.csv'), as.is = TRUE) %>%
             filter(SITE_NO != '3')
           
           save(juv.mmp.samples, file = paste0(DATA_PATH, 'primary/juv.mmp.samples.RData'))
           
           CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                             item = 'mmp_juv_samp',status = 'success')
           
           }, logFile=LOG_FILE, Category='--Data extraction--',
msg=paste0('MMP juveniles samples.'), return=NULL)
    }          
           
CI_get_ltmp_juveniles_data <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'ltmp_juv',
                   label = "LTMP juveniles", status = 'pending')
    CI_tryCatch({

        writeLines("select p_code,aims_reef_name as reef, visit_no, site_no, genus,sum(abundance) abundance
           from v_in_sample s, demog b 
           where s.sample_id = b.sample_id
           and p_code in ('RM','RMRAP','RAP') 
           and visit_no >0
           and transect_no<6
           group by p_code,aims_reef_name,visit_no, site_no, genus",
          paste0(PRIMARY_DATA_PATH, "juv.ltmp.sql"))

        system(paste0("java -jar dbExport.jar ", PRIMARY_DATA_PATH, "juv.ltmp.sql ",
                      PRIMARY_DATA_PATH, "juv.ltmp.csv reef reefmon"),
               ignore.stdout = TRUE)

        juv.ltmp = read.csv(paste0(PRIMARY_DATA_PATH, 'juv.ltmp.csv'), as.is=TRUE) %>%
            mutate(DEPTH = as.integer(ifelse(REEF == 'Middle Reef', 2, 8)),
                   P_CODE = ifelse(P_CODE %in% c('RM', 'RMRAP', 'RAP'), 'RM', P_CODE))

        save(juv.ltmp, file = paste0(DATA_PATH, 'primary/juv.ltmp.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                          item = 'ltmp_juv',status = 'success')
        
    }, logFile=LOG_FILE, Category='--Data extraction--',
    msg=paste0('LTMP juveniles.'), return=NULL)
}

CI_get_ltmp_juveniles_samples <- function() {
  CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                 item = 'ltmp_juv_samp',
                 label = "LTMP juveniles samples", status = 'pending')
  CI_tryCatch({
    
    writeLines("select p_code, aims_reef_name as reef, visit_no, site_no
           from v_in_sample 
           where p_code in('RM','RMRAP','RAP') 
           and visit_no >14 
           and sample_type = 'JUVCOR'
           group by p_code, aims_reef_name, visit_no, site_no",
               paste0(PRIMARY_DATA_PATH, "juv.ltmp.samples.sql"))
    
    system(paste0("java -jar dbExport.jar ", PRIMARY_DATA_PATH, "juv.ltmp.samples.sql ",
                  PRIMARY_DATA_PATH, "juv.ltmp.samples.csv reef reefmon"),
           ignore.stdout = TRUE)
    
    juv.ltmp.samples<-read.csv(paste0(PRIMARY_DATA_PATH, 'juv.ltmp.samples.csv'), as.is = TRUE) |>
          mutate(DEPTH = as.integer(ifelse(REEF == 'Middle Reef', 2, 8)),
             P_CODE = ifelse(P_CODE %in% c('RM', 'RMRAP', 'RAP'), 'RM', P_CODE))
    
    save(juv.ltmp.samples, file = paste0(DATA_PATH, 'primary/juv.ltmp.samples.RData'))
    
    CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                      item = 'ltmp_juv_samp',status = 'success')
    
  }, logFile=LOG_FILE, Category='--Data extraction--',
  msg=paste0('LTMP juveniles samples.'), return=NULL)
}          

CI_combine_juveniles_samples <- function() {
  CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                 item = 'combine_juv_samples',
                 label = "Combine juveniles samples", status = 'pending')
  CI_tryCatch({
    
    load(file = paste0(DATA_PATH, 'primary/juv.mmp.samples.RData'))
    load(file = paste0(DATA_PATH, 'primary/juv.ltmp.samples.RData'))
    
    juvenile.samples <- juv.mmp.samples %>%
      rbind(juv.ltmp.samples)
    
    save(juvenile.samples, file = paste0(PRIMARY_DATA_PATH, 'juvenile.samples.RData'))
    
    CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                      item = 'combine_juv_samles',status = 'success')
    
  }, logFile=LOG_FILE, Category='--Data extraction--',
  msg=paste0('Combine juvenile samples.'), return=NULL)
}


CI_combine_juveniles_data <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'combine_juv',
                   label = "Combine juveniles", status = 'pending')
    CI_tryCatch({

        load(file = paste0(DATA_PATH, 'primary/juv.mmp.RData'))
        load(file = paste0(DATA_PATH, 'primary/juv.ltmp.RData'))

        juveniles <- juv.mmp %>%
            rbind(juv.ltmp)

        save(juveniles, file = paste0(PRIMARY_DATA_PATH, 'juveniles.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                          item = 'combine_juv',status = 'success')
        
    }, logFile=LOG_FILE, Category='--Data extraction--',
    msg=paste0('Combine juveniles.'), return=NULL)
}


CI_get_mmp_point_data <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'mmp_points',
                   label = "MMP points", status = 'pending')
    CI_tryCatch({

        writeLines("select p_code, mmp_site_name as reef, depth, visit_no, site_no, transect_no, b.video_code, count(point_no) points
           from v_in_sample s, vpoint_codes c, rm_vpoint b
           where s.sample_id = b.sample_id
           and b.video_code=c.video_code
            and p_code in ('IN','AP','GH','RR') 
           and visit_no >0 
           and s.mmp_site_name not like 'Cape%'
           and s.mmp_site_name not like 'Connor%'
           and group_code not like 'IN'
           group by p_code, mmp_site_name, depth, visit_no,site_no,transect_no, b.video_code",
          paste0(PRIMARY_DATA_PATH, "mmp.points.sql"))

        system(paste0("java -jar dbExport.jar ", PRIMARY_DATA_PATH, "mmp.points.sql ",
                      PRIMARY_DATA_PATH, "mmp.points.csv reef reefmon"),
               ignore.stdout = TRUE)

        mmp.points <- read.csv(paste0(PRIMARY_DATA_PATH, 'mmp.points.csv')) 

        save(mmp.points, file = paste0(DATA_PATH, 'primary/mmp.points.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                          item = 'mmp_points',status = 'success')
        
    }, logFile=LOG_FILE, Category='--Data extraction--',
    msg=paste0('MMP points.'), return=NULL)
}


CI_get_ltmp_point_data <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'ltmp_points',
                   label = "LTMP points", status = 'pending')
    CI_tryCatch({

        writeLines("select p_code, aims_reef_name as reef, depth, visit_no, site_no, transect_no, b.video_code, count(point_no) points
           from v_in_sample s, vpoint_codes c, rm_vpoint b
           where s.sample_id = b.sample_id
           and b.video_code=c.video_code
           and p_code in ('RAP', 'RM','RMRAP') 
           and visit_no >0 
           and transect_no<6
           and group_code not like 'IN'
           group by p_code, aims_reef_name, depth, visit_no,site_no, transect_no,b.video_code",
          paste0(PRIMARY_DATA_PATH, "ltmp.points.sql"))

        system(paste0("java -jar dbExport.jar ", PRIMARY_DATA_PATH, "ltmp.points.sql ",
                      PRIMARY_DATA_PATH, "ltmp.points.csv reef reefmon"),
               ignore.stdout = TRUE)

        ltmp.points <- read.csv(paste0(PRIMARY_DATA_PATH, 'ltmp.points.csv')) 

        save(ltmp.points, file = paste0(DATA_PATH, 'primary/ltmp.points.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                          item = 'ltmp_points',status = 'success')
        
    }, logFile=LOG_FILE, Category='--Data extraction--',
    msg=paste0('LTMP points.'), return=NULL)
}

CI_combine_point_data <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'combine_points',
                   label = "Combine points", status = 'pending')
    CI_tryCatch({

        load(file = paste0(DATA_PATH, 'primary/mmp.points.RData'))
        load(file = paste0(DATA_PATH, 'primary/ltmp.points.RData'))

        points.raw <- ltmp.points %>% 
            mutate(DEPTH = ifelse(REEF == "Middle Reef", 2, 8) ,
                   P_CODE = ifelse(P_CODE %in% c('RM', 'RMRAP', 'RAP'), 'RM', P_CODE)) %>%
            rbind(mmp.points)

        save(points.raw, file = paste0(PRIMARY_DATA_PATH, 'points.raw.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                          item = 'combine_points',status = 'success')
        
    }, logFile=LOG_FILE, Category='--Data extraction--',
    msg=paste0('Combine points.'), return=NULL)
}

## External data are defined by two files:
## - points.data.<P_CODE>.csv
##   this is the observed spatio-temporal benthic points data
## - samples.<P_CODE>.csv
##   this provides the sample location geodata
CI_get_external_data <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'external_data',
                   label = "External data", status = 'pending')
    CI_tryCatch({

        ##JCU data total points jcu exclude the early data based on 25
        ## points per transect, and impose a visit_no beginning at the
        ## first year of 50 points per transect.
        files <- list.files(path = paste0(DATA_PATH, 'external'),
                            pattern = 'points.data.*.csv',
                            full.names = TRUE)
        ext.points <- purrr::map(files,
                                 .f = ~ read_csv(.x) %>%
                                     dplyr::filter(total.points > 25) %>% 
                                     ## this record is they only one
                                     ## in the data set that is not
                                     ## whole numbers
                                     filter(!(SITE_NO=='HY2' & TRANSECT=='5' & YEAR=='2012')) %>%  
                                     mutate(VISIT_NO=YEAR-2002,
                                            REPORT_YEAR=YEAR) %>%
                                     filter(REPORT_YEAR < (CI$setting$FINAL_YEAR + 1)) %>%
                                     rename(TRANSECT_NO=TRANSECT) %>%
                                     suppressMessages() %>%
                                     suppressWarnings())
        ## Total points (transect level)
        total.points.ext <- purrr::map(.x = ext.points,
                                       .f = ~ .x %>%
                                           dplyr::select(P_CODE, REEF,DEPTH, VISIT_NO,
                                                         REPORT_YEAR, SITE_NO, TRANSECT_NO,
                                                         total.points) %>%
                                           group_by(P_CODE, REEF, DEPTH, VISIT_NO,
                                                    REPORT_YEAR, SITE_NO, TRANSECT_NO) %>%
                                           summarise(total.points = sum(total.points)) %>%
                                           ungroup() %>% 
                                           suppressMessages() %>%
                                           suppressWarnings())

        ## hard coral points  
        hc.ext <- purrr::map(.x = ext.points,
                             .f = ~ .x %>% 
                                 dplyr::select(P_CODE, REEF, DEPTH, VISIT_NO,
                                               REPORT_YEAR, SITE_NO, TRANSECT_NO, HC) %>%
                                 group_by(P_CODE, REEF, DEPTH, VISIT_NO, REPORT_YEAR,
                                          SITE_NO, TRANSECT_NO) %>%
                                 summarise(HC = sum(HC))%>%
                                 ungroup() %>%
                                 suppressMessages() %>%
                                 suppressWarnings())

        ## MA points
        ma.ext <- purrr::map(.x = ext.points,
                             .f = ~ .x %>% 
                                 dplyr::select(P_CODE, REEF, DEPTH, VISIT_NO,
                                               REPORT_YEAR, SITE_NO, TRANSECT_NO, MA) %>%
                                 group_by(P_CODE, REEF, DEPTH, VISIT_NO, REPORT_YEAR,
                                          SITE_NO, TRANSECT_NO) %>%
                                 summarise(MA = sum(MA)) %>%
                                 ungroup() %>%
                                 suppressMessages() %>%
                                 suppressWarnings())

        ## Algae points 
        a.ext <- purrr::map(.x = ext.points,
                            .f = ~ .x %>% 
                                dplyr::select(P_CODE, REEF, DEPTH, VISIT_NO, REPORT_YEAR,
                                              SITE_NO, TRANSECT_NO, A) %>%
                                group_by(P_CODE, REEF, DEPTH, VISIT_NO, REPORT_YEAR,
                                         SITE_NO, TRANSECT_NO) %>%
                                summarise(A = sum(A)) %>%
                                ungroup() %>%
                                suppressMessages() %>%
                                suppressWarnings())
        ## samples data
        files_samples <- list.files(path = paste0(DATA_PATH, 'external'),
                                    pattern = 'sample.*.csv',
                                    full.names = TRUE)
        ext.samples <- purrr::map(files_samples,
                                  .f = ~ read_csv(.x) %>%
                                      dplyr::select(P_CODE,REEF,SITE_NO,Lat,Long) %>%
                                      unique %>%
                                      rename(LONGITUDE="Long",
                                             LATITUDE="Lat") %>%
                                      mutate(SITE_NO=as.character(SITE_NO)) %>%
                                      suppressMessages() %>%
                                      suppressWarnings())

        points.analysis.data.transect.jcu <-
            purrr::pmap(.l = list(hc.ext, ma.ext, a.ext,
                                  ext.samples, total.points.ext),
                        .f = ~ ..1 %>%
                            left_join(..2) %>%
                            left_join(..3) %>%
                            left_join(..4) %>%
                            left_join(..5) %>%
                            mutate(REEF = as.factor(REEF),
                                   reef.site = as.factor(paste0(REEF,'.',SITE_NO)),
                                   reef.site.tran = as.factor(paste0(reef.site,'.',TRANSECT_NO)),
                                   Project = as.factor(P_CODE)) %>%
                            suppressMessages() %>%
                            suppressWarnings()
                        ) %>%
            bind_rows() %>%
            dplyr::select(P_CODE, REEF, DEPTH, VISIT_NO,      
                            SITE_NO, TRANSECT_NO, HC, REPORT_YEAR,
                            LATITUDE, LONGITUDE, A, MA,
                            total.points, reef.site, reef.site.tran,
                            Project)
        
        save(points.analysis.data.transect.jcu,
             file = paste0(DATA_PATH, 'external/points.analysis.data.transect.jcu.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                          item = 'external_data',status = 'success')
        
    }, logFile=LOG_FILE, Category='--Data extraction--',
    msg=paste0('External data'), return=NULL)
}

#######
# bring in JCU point data as wrangled  by AT 
#C:\Users\athompso\OneDrive - Australian Institute of Marine Science\Coral_Indicators\R
# https://aimsgovau-my.sharepoint.com/:u:/g/personal/a_thompson_aims_gov_au/ERqAfnNP5TZOsBLt7LaUDacBVgE5dSv2uUX_Z_vKFVFapQ?e=qkSyaC
#######CI_get_external_data <- function() {
CI_get_external_data <- function() {
CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
               item = 'external_data',
               label = "External data", status = 'pending')
CI_tryCatch({
  
  ##JCU data total points jcu exclude the early data based on 25
  ## points per transect, and impose a visit_no beginning at the
  ## first year of 50 points per transect.
  load(file = paste0(DATA_PATH, 'external/jcu.points.group.RData'))
  load(file = paste0(DATA_PATH, 'external/jcu.loc.RData'))
 
  # eork through the below to makke data compatible iwth aims data
  
  
   ext.points <- purrr::map(files,
                           .f = ~ read_csv(.x) %>%
                             dplyr::filter(total.points > 25) %>% 
                             ## this record is they only one
                             ## in the data set that is not
                             ## whole numbers
                             filter(!(SITE_NO=='HY2' & TRANSECT=='5' & YEAR=='2012')) %>%  
                             mutate(VISIT_NO=YEAR-2002,
                                    REPORT_YEAR=YEAR) %>%
                             filter(REPORT_YEAR < (CI$setting$FINAL_YEAR + 1)) %>%
                             rename(TRANSECT_NO=TRANSECT) %>%
                             suppressMessages() %>%
                             suppressWarnings())
  ## Total points (transect level)
  total.points.ext <- purrr::map(.x = ext.points,
                                 .f = ~ .x %>%
                                   dplyr::select(P_CODE, REEF,DEPTH, VISIT_NO,
                                                 REPORT_YEAR, SITE_NO, TRANSECT_NO,
                                                 total.points) %>%
                                   group_by(P_CODE, REEF, DEPTH, VISIT_NO,
                                            REPORT_YEAR, SITE_NO, TRANSECT_NO) %>%
                                   summarise(total.points = sum(total.points)) %>%
                                   ungroup() %>% 
                                   suppressMessages() %>%
                                   suppressWarnings())
  
  ## hard coral points  
  hc.ext <- purrr::map(.x = ext.points,
                       .f = ~ .x %>% 
                         dplyr::select(P_CODE, REEF, DEPTH, VISIT_NO,
                                       REPORT_YEAR, SITE_NO, TRANSECT_NO, HC) %>%
                         group_by(P_CODE, REEF, DEPTH, VISIT_NO, REPORT_YEAR,
                                  SITE_NO, TRANSECT_NO) %>%
                         summarise(HC = sum(HC))%>%
                         ungroup() %>%
                         suppressMessages() %>%
                         suppressWarnings())
  
  ## MA points
  ma.ext <- purrr::map(.x = ext.points,
                       .f = ~ .x %>% 
                         dplyr::select(P_CODE, REEF, DEPTH, VISIT_NO,
                                       REPORT_YEAR, SITE_NO, TRANSECT_NO, MA) %>%
                         group_by(P_CODE, REEF, DEPTH, VISIT_NO, REPORT_YEAR,
                                  SITE_NO, TRANSECT_NO) %>%
                         summarise(MA = sum(MA)) %>%
                         ungroup() %>%
                         suppressMessages() %>%
                         suppressWarnings())
  
  ## Algae points 
  a.ext <- purrr::map(.x = ext.points,
                      .f = ~ .x %>% 
                        dplyr::select(P_CODE, REEF, DEPTH, VISIT_NO, REPORT_YEAR,
                                      SITE_NO, TRANSECT_NO, A) %>%
                        group_by(P_CODE, REEF, DEPTH, VISIT_NO, REPORT_YEAR,
                                 SITE_NO, TRANSECT_NO) %>%
                        summarise(A = sum(A)) %>%
                        ungroup() %>%
                        suppressMessages() %>%
                        suppressWarnings())
  ## samples data
  files_samples <- list.files(path = paste0(DATA_PATH, 'external'),
                              pattern = 'sample.*.csv',
                              full.names = TRUE)
  ext.samples <- purrr::map(files_samples,
                            .f = ~ read_csv(.x) %>%
                              dplyr::select(P_CODE,REEF,SITE_NO,Lat,Long) %>%
                              unique %>%
                              rename(LONGITUDE="Long",
                                     LATITUDE="Lat") %>%
                              mutate(SITE_NO=as.character(SITE_NO)) %>%
                              suppressMessages() %>%
                              suppressWarnings())
  
  points.analysis.data.transect.jcu <-
    purrr::pmap(.l = list(hc.ext, ma.ext, a.ext,
                          ext.samples, total.points.ext),
                .f = ~ ..1 %>%
                  left_join(..2) %>%
                  left_join(..3) %>%
                  left_join(..4) %>%
                  left_join(..5) %>%
                  mutate(REEF = as.factor(REEF),
                         reef.site = as.factor(paste0(REEF,'.',SITE_NO)),
                         reef.site.tran = as.factor(paste0(reef.site,'.',TRANSECT_NO)),
                         Project = as.factor(P_CODE)) %>%
                  suppressMessages() %>%
                  suppressWarnings()
    ) %>%
    bind_rows() %>%
    dplyr::select(P_CODE, REEF, DEPTH, VISIT_NO,      
                  SITE_NO, TRANSECT_NO, HC, REPORT_YEAR,
                  LATITUDE, LONGITUDE, A, MA,
                  total.points, reef.site, reef.site.tran,
                  Project)
  
  save(points.analysis.data.transect.jcu,
       file = paste0(DATA_PATH, 'external/points.analysis.data.transect.jcu.RData'))
  
  CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                    item = 'external_data',status = 'success')
  
}, logFile=LOG_FILE, Category='--Data extraction--',
msg=paste0('External data'), return=NULL)
}






CI_get_spatial_data <- function() {
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'spatial_data',
                   label = "Spatial data", status = 'pending')
    CI_tryCatch({

        ## Bioregions (reef)
        bioregions <- st_read(paste0("https://services8.arcgis.com/",
                                "ll1QQ2mI4WMXIXdm/arcgis/rest/services/",
                                "Reef_marine_bioregions_of_the_Great_Barrier_Reef/",
                                "FeatureServer/0/query?outFields=*&where=1%3D1&f=geojson"),
                         quiet = T) %>%
            select(BIOREGION, DESCRIP, geometry) %>% 
            mutate(Region = "BIOREGION") %>%
            st_make_valid()

        save(bioregions,
             file = paste0(DATA_PATH, 'primary/bioregions.RData'))
        
        
        ## TUMRA
        tumra <- st_read(paste0("https://services8.arcgis.com/",
                                "ll1QQ2mI4WMXIXdm/arcgis/rest/services/",
                                "Traditional_Use_of_Marine_Resources_Agreement_areas/",
                                "FeatureServer/0/query?outFields=*&where=1%3D1&f=geojson"),
                         quiet = T) %>%
            select(NAME, geometry) %>% 
            rename(Name = NAME) %>%
            ## group_by(Name) %>%
            ## summarise(geometry = st_combine(.))  %>%
            mutate(Region = "TUMRA") %>%
            st_make_valid()

        save(tumra,
             file = paste0(DATA_PATH, 'primary/tumra.RData'))

        ## GBRMPA management
        gbrmpa.management <- st_read(paste0("https://services8.arcgis.com/",
                                 "ll1QQ2mI4WMXIXdm/arcgis/rest/services/",
                                 "Management_Areas_of_the_Great_Barrier_Reef_Marine_Park/",
                                 "FeatureServer/0/query?outFields=*&where=1%3D1&f=geojson"),
                          quiet = T) %>%
            select(AREA_DESCR, geometry) %>%
            rename(Name = AREA_DESCR)%>%
            mutate(Region = "GBRMPA_Management")

        save(gbrmpa.management,
             file = paste0(DATA_PATH, 'primary/gbrmpa.management.RData'))

        ## GBR
        gbrmpa <- st_read(paste0("https://services8.arcgis.com/",
                                 "ll1QQ2mI4WMXIXdm/arcgis/rest/services/",
                                 "Great_Barrier_Reef_Marine_Park_Boundary/",
                                 "FeatureServer/0/query?outFields=*&where=1%3D1&f=geojson"),
                          quiet = T) %>%
            mutate(Region = "GBR", Name = "GBRMP")%>%
            select(Name, Region, geometry) 

        save(gbrmpa,
             file = paste0(DATA_PATH, 'primary/gbrmpa.RData'))

        ## Zones
        northern.bbox <- st_bbox(c(xmin = 142, xmax = 155, ymin = -15.4, ymax = 0)) %>%
            st_as_sfc() %>%
            st_sf(crs = st_crs(gbrmpa)) %>%
            mutate(Zone = 'Northern')
        central.bbox <- rbind(c(142,-20.7),
                               c(148.7,-20.7),
                               c(152,-19.6),
                               c(152,-15.4),
                               c(142,-15.4)) %>%
            st_linestring() %>%
            st_cast("POLYGON") %>%
            st_sfc(crs = st_crs(gbrmpa)) %>%
            st_sf() %>%
            mutate(Zone = 'Central')

       southern.bbox <- rbind(c(142,-20.7),
                              c(148.7,-20.7),
                              c(152,-19.6),
                              c(155,-19.6),
                              c(155,-25),
                              c(142,-25)) %>%
            st_linestring() %>%
            st_cast("POLYGON") %>%
            st_sfc(crs = st_crs(gbrmpa)) %>%
            st_sf() %>%
            mutate(Zone = 'Southern')
       zones.bbox <- rbind(northern.bbox, central.bbox, southern.bbox)
 
       zones <- gbrmpa %>%
           st_intersection(zones.bbox) %>%
           mutate(Region = 'Zones',
                  Name = Zone) %>%
           dplyr::select(-Zone) %>%
           suppressMessages() %>%
           suppressWarnings()
       
        save(zones,
             file = paste0(DATA_PATH, 'primary/zones.RData'))
        
        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'spatial_data',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('Spatial data'), return=NULL)
}
