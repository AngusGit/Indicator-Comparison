
CI_models_collate_indices <- function() {
  

        ## site.location <- get(load(file = paste0(DATA_PATH, 'processed/site.location.RData')))
        spatial_lookup <- get(load(file = paste0(DATA_PATH, "processed/spatial_lookup.RData")))
        files <- list.files(path = paste0(DATA_PATH, "modelled"),
                              pattern = "..__scores_reef_year.RData",
                              full.names = TRUE)
        indices <- purrr::map_df(.x = files,
                              .f = ~ {
                                  indicator <- str_replace(.x,
                                                        #   ".*modelled/([A-Z]{2})__.*",
                                                        ".*modelled/([A-Z]{2}).*",
                                                           "\\1") 
                                  indicator <-  case_when(
                                      indicator == "CC" ~ "Coral.cover",
                                      indicator == "MA" ~ "Macroalgae",
                                      indicator == "JU" ~ "Juvenile.density",
                                      indicator == "CO" ~ "Community.composition",
                                      indicator == "RP" ~ "Recovery.performance"
                                  )
                                
                                  
                                  x <- get(load(.x)) %>%
                                      dplyr::select(Summary) %>%
                                      unnest(Summary) 
                                  x <- x %>%
                                          left_join(spatial_lookup %>%
                                                    dplyr::select(REEF.d, REEF, DEPTH, BIOREGION, DEPTH.f) %>%
                                                    distinct()) 
                                  
                                  x <- x %>%  
                                      mutate(#Level = {{level}},
                                             Indicator = {{indicator}},
                                             Reference = case_when(
                                                 Metric == 'pcb.rescale.dist.metric' ~ 'Critical',
                                                 Metric == 'rescale.dist.metric' ~ 'Baseline',
                                                 Metric == 'rescale.consequence.metric' ~ 'Critical',
                                                 Metric == 'distance.metric' ~ 'Baseline',
                                                 Metric == 'Total' ~ 'Baseline',
                                                 Metric == 'Acropora' ~ 'Critical',
                                                 Metric == 'Reference' ~ 'Baseline',
                                                 Metric == 'Critical' ~ 'Critical',
                                                 Metric == 'Reference' ~ 'reference',
                                                 Metric == 'Critical' ~ 'critical'
                                                 )) %>%
                                    dplyr::select(
                                                    Year = fYEAR,
                                                    REEF,
                                                    REEF.d
                                                    Depth = DEPTH.f,
                                                    Metric,
                                                    variable,
                                                    Median = median,
                                                    Lower = lower,
                                                    Upper = upper) %>% 
                                      suppressMessages() %>%
                                      suppressWarnings()
                                  }
                              )
            
        
        save(indices, file = paste0(DATA_PATH, 'modelled/indices.RData'))
        write_csv(indices, file = paste0(TABS_PATH, "/Indices.csv"))
        
        

 
}

## tests

tests <- function(level = 'NRM', value = 'Burdekin') {
    cat(paste0('\n##-----', level, ' (', value, ') tn.reefs from indices in 2000 ---------\n'))
    indices %>%
        filter(Level == level,
               Year == 2000,
               Indicator == 'Coral.cover',
               Reference == 'Baseline') %>%
        as.data.frame %>%
        print()

    
    cat(paste0('\n##-----all unique ', level, ' (', value, ') reefs in site.locations ---------\n'))
    if (is.na(value)) {
        SS <- site.location %>%
            filter(is.na(!!sym(level))) %>%
            dplyr::select(REEF) %>%
            distinct() %>%
            pull(REEF) %>%
            unique
    } else {
        SS <- site.location %>%
            filter(!!sym(level) == value) %>%
            dplyr::select(REEF) %>%
            distinct() %>%
            pull(REEF) %>%
            unique
    }
    print(SS)
    
    cat(paste0('\n##-----all unique ', level, ' (', value, ') reefs in points.analysis.data.transect in 2000---------\n'))
    points.analysis.data.transect %>%
        filter(REPORT_YEAR == 2000, REEF %in% SS) %>%
        dplyr::select(REEF) %>%
        unique %>%
        print()
}

## tests(level = 'NRM', value = 'Burdekin')
## tests(level = 'NRM', value = 'Burnett Mary')
## tests(level = 'NRM', value = 'Cape York')
## tests(level = 'NRM', value = 'Fitzroy')
## tests(level = 'NRM', value = 'Mackay Whitsunday')
## tests(level = 'NRM', value = 'Wet Tropics')
## tests(level = 'GBRMP', value = 'GBRMP')
