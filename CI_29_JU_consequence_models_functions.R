CI__get_JUV_MA_site_data <- function() {
    load(file = paste0(DATA_PATH, 'processed/points.analysis.data.transect.aims.RData'))
    load(file = paste0(DATA_PATH, 'processed/juv.df.RData'))
    load(paste0(DATA_PATH, 'processed/spatial_lookup.RData'))
    load(file = paste0(DATA_PATH, 'primary/disturbances.RData'))

    juv.ma.dat.site <- points.analysis.data.transect.aims %>%
        dplyr::select(-P_CODE, -TRANSECT_NO, -reef.site.tran) %>%
        group_by(REEF, VISIT_NO, DEPTH, SITE_NO, REPORT_YEAR,
                 LATITUDE, LONGITUDE, reef.site, Project) %>%
        summarise_all(sum) %>%
        ungroup() %>%
        left_join(juv.df %>%
                  dplyr::select(REEF, VISIT_NO, DEPTH, SITE_NO,
                                total.juv, Acropora, avail.area)) %>%
        filter(!is.na(total.juv)) %>%  # these are missing samples
        mutate(reef.depth = paste0(REEF,DEPTH, sep=''),
               reef.depth.site = as.factor(paste0(reef.depth, SITE_NO, sep = '')),
               MAp = MA/A,
               Acropora.d = Acropora/avail.area,
               total.juv.d = total.juv/avail.area) %>%
        left_join(spatial_lookup %>%
                  dplyr::select(REEF, DEPTH, Shelf, DEPTH.f) %>%
                  distinct()) %>%
        mutate(Habitat =as.factor(paste(Shelf, DEPTH.f))) %>%
        left_join(disturbance.reef %>%
                  dplyr::select(-P_CODE) %>%
                  filter(RANK==1)) %>%
        suppressMessages() %>%
        suppressWarnings()

    juv.no.dist <- juv.ma.dat.site %>%
        group_by(Habitat, REEF, reef.site, DEPTH) %>%
        arrange(REEF, reef.site, DEPTH, REPORT_YEAR)%>%
        mutate(d.lag = lag(DISTURBANCE)) %>%
        ungroup() %>%
        filter(DISTURBANCE == 'n') %>% 
        suppressMessages() %>%
        suppressWarnings()

    ## apply lagsd.lag
    juv.ma.site <- juv.no.dist %>% 
        group_by(Habitat, REEF, reef.site, DEPTH) %>%
        arrange(REEF, reef.site, DEPTH, REPORT_YEAR)%>%
        mutate(l.MAp = lag(MAp),
               l.MAp2 = lag(n = 2, MAp),
               duration = REPORT_YEAR - lag(REPORT_YEAR),
               MApLag = ifelse(is.na(duration), MAp,
                        ifelse(duration == 1 &
                               !is.na(lag(duration)) &
                               lag(duration) == 1,
                               (MAp + l.MAp + l.MAp2)/3,
                               (MAp+l.MAp)/2))
               ) %>%
        ungroup() %>% 
        suppressMessages() %>%
        suppressWarnings()
    save(juv.ma.site,
         file = paste0(DATA_PATH, 'processed/juv.ma.site.RData'))
    #return(juv.ma.site)
}

