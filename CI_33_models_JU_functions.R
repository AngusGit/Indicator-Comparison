

CI_clear_models_JU_data <- function() {
    files <- list.files(path = paste0(DATA_PATH, "modelled"),
                        pattern = "JU.*|data_ju.*",
                        full.names = TRUE)
    unlink(files)
}

CI_models_JU_get_baselines <- function() {
    load(file = paste0(DATA_PATH, "parameters/JUV_baseline.RData"))

        baselines <- JUV_baseline %>%
            ## mutate(DEPTH.f = ifelse(DEPTH.f == 'deep', 'deep slope', 'shallow slope'))
            mutate(DEPTH.f = ifelse(str_detect(DEPTH.f, 'deep'), 'deep slope', 'shallow slope'))
        
        save(baselines,
             file = paste0(DATA_PATH, 'modelled/JU__baseline_posteriors.RData'))

    }

CI__clean_ju_data <- function(data) {
    df.a <- data %>%
        ungroup() %>%
        filter(!REPORT_YEAR == "NA") %>%
        droplevels %>%
        dplyr::select(-VISIT_NO)
# remove reefs with only a single observation    
    report.year.ss <- df.a %>%
        dplyr::select(REEF.d, REPORT_YEAR) %>%
        unique() %>% 
        group_by(REEF.d) %>%
        summarise(report.year.ss = n()) %>%
        ungroup()
    df.a <- df.a %>%
        left_join(report.year.ss) %>%
        filter(!report.year.ss == "1") %>%
        droplevels
    
    detect.zeros <- df.a %>% 
        group_by(REEF.d, fYEAR) %>%
        summarise(juv.transect.sum = sum(value)) %>%
        ungroup
    
    df.b<- df.a %>% 
        left_join(detect.zeros) %>% 
        filter(!juv.transect.sum==0) %>% 
        droplevels
    ## Finally, model also terminates if
    ## there is only 1 level of 'year'
    ## greater than 0 (i.e. no variance)
    no.variance <- df.b %>%
        dplyr::select(REEF.d, fYEAR, juv.transect.sum) %>% 
        filter(juv.transect.sum>0) %>%
        droplevels() %>%
        dplyr::select(REEF.d, fYEAR) %>%
        unique() %>%
        group_by(REEF.d) %>%
        summarise(ss.years.morethanzero = n()) %>%
        ungroup

    #df.b %>%  AT alteration to retain years with valid zeros when other years at the reef were ok, models still run and maybe this also because of the change to cell mean models
      df.a %>% 
        left_join(no.variance) %>% 
        filter(ss.years.morethanzero>1) %>% 
        droplevels %>%
        dplyr::select(-any_of(c("juv.transect.sum",
                                "report.year.ss")))
}

CI_models_JU_prepare_data <- function() {
    
        load(file = paste0(DATA_PATH, "processed/juv.df.RData"))
        load(file = paste0(DATA_PATH, "processed/spatial_lookup.RData"))
        
        #juv.df<-juv.df |> filter (REEF=="Daydream")
  #places where no Acropora juvs were observed that will be removed by    CI__clean_ju_data   

  #For the case study only      
        juv.df<-juv.df |> 
          left_join(spatial_lookup |> 
                      dplyr::select(REEF,Shelf) |> 
                      unique()) |> 
          filter(Shelf=="Inshore")
        
        juv <- juv.df %>% 
            group_by(REEF.d) %>%
            mutate(LONGITUDE = mean(LONGITUDE),
                   LATITUDE = mean(LATITUDE),
                   ) %>%
            ungroup() %>%
            dplyr::select(-any_of(c("total.juv.excl.turb",
                                    "Turbinaria", "Non_Acropora",
                                    "Non_Acr_exTurb"))) %>%
            dplyr::rename("Total" = total.juv) %>%
            pivot_longer(cols = c(Total, Acropora),
                         names_to = 'Taxa',
                         values_to = 'value') %>% 
            group_by(DEPTH.f, Taxa) %>%
            nest()

        data_ju <- juv %>%
            mutate(data = map(.x = data,
                              .f = ~ CI__clean_ju_data(.x) %>%
                                  suppressMessages() %>%
                                  suppressWarnings()
                              )
                   )
        
        save(data_ju,
              file = paste0(DATA_PATH, 'modelled/data_ju.RData'))
}

CI_models_JU_prepare_nest <- function() {
            ## nest the data
 load(paste0(DATA_PATH, 'modelled/data_ju.RData'))
 
        mods <- data_ju %>% 
            unnest(data) %>%
            ungroup() %>% 
            nest(data=everything(), .by=c(REEF.d, DEPTH.f,Taxa)) %>%
          mutate(n = 1:n())
        
        ## Prepare the data
        mods <- mods %>%
            mutate(newdata = map(.x = data,
                                 .f = ~ .x %>%
                                     droplevels() %>% 
                                     ## tidyr::expand(Site = Site, fYEAR = fYEAR, # replace
                                     tidyr::expand(REEF.d = REEF.d, fYEAR = fYEAR,
                                                   Site = NA, #Transect = NA, # removed as no Transect level data used 
                                                   avail.area=1,   # is this appropriate to add in ??
                                                   value = NA) %>%
                                     distinct()
                                 ),
                   Full_data = pmap(.l = list(data, newdata),
                                    .f = ~ ..1 %>% bind_rows(..2))
                   )
        save(mods,
             file = paste0(DATA_PATH, 'modelled/JU__mods.RData'))
        
  }


CI__fit_JU_model <- function(form, data, family='nbinomial', n, N) {
    
  reef <- unique(data$REEF.d)
  taxa <- na.omit(unique(data$Taxa))
        
  environment(form)<-environment()
  
        mod <- inla(formula = form,
                    E = data$avail.area,
                    data = data,
                    family = family, 
                    control.predictor = list(link = 1, compute = TRUE),
                    control.compute = list(
                        dic = TRUE, cpo = TRUE, waic = TRUE,
                        config = TRUE) 
                    )
     
        save(mod, file = paste0(DATA_PATH, "modelled/JU__", reef, '_', taxa, '__model.RData'))
        draws <- inla.posterior.sample(n=1000, mod, seed=123) %>%
            suppressWarnings() %>%
            suppressMessages()
        save(draws, file = paste0(DATA_PATH, "modelled/JU__", reef, '_', taxa, '__draws.RData'))
}


CI_models_JU_fit_models <- function() {
    
        load(file = paste0(DATA_PATH, "modelled/data_ju.RData"))
        load(file = paste0(DATA_PATH, 'modelled/JU__mods.RData'))  
       
       
        
        form <- value ~ 0+fYEAR +   ## added the 0+ to make a cell means formula
            f(Site , model='iid')
        
        ## Fit the models - output models and draws to
        purrr::pwalk(.l = list(mods$Full_data, mods$n),
                     .f = ~ CI__fit_JU_model(form = form, 
                                             data = ..1,
                                             family = 'nbinomial',
                                             n = ..2,
                                             N = nrow(mods))
                   )

}

CI__cellmeans_JU_model <- function(obs_data, Full_data, newdata, n, N) {
    
  reef <- unique(newdata$REEF.d)
  taxa <- na.omit(unique(obs_data$Taxa))
        
        draws <- get(load(file = paste0(DATA_PATH, "modelled/JU__", reef, '_', taxa, '__draws.RData')))
        cellmeans <- sapply(draws, function(x)
            x[[2]][(nrow(Full_data)-nrow(newdata)+1):nrow(Full_data)]) 
        posteriors <- newdata %>%
            dplyr::select(fYEAR, REEF.d) %>%
            cbind(exp(cellmeans)) %>%
            pivot_longer(cols = matches('[0-9]'), names_to = 'Rep') %>%
            mutate(REEF.d = reef,
                   .draw = as.integer(Rep)) %>%
            dplyr::select(-Rep) %>%
            left_join(obs_data %>%
                      dplyr::select(fYEAR, REEF.d) %>%
                      distinct()) %>% 
            suppressWarnings() %>%
            suppressMessages()
        save(posteriors, file = paste0(DATA_PATH, "modelled/JU__", reef, '_', taxa, '__posteriors.RData'))
    
}

CI_models_JU_cellmeans <- function() {
    

        load(file = paste0(DATA_PATH, "modelled/data_ju.RData"))
        load(file = paste0(DATA_PATH, 'modelled/JU__mods.RData'))
        ## Calculate cellmeans
        cellmeans <- purrr::pwalk(.l = list(mods$data, mods$Full_data, mods$newdata, mods$n, nrow(mods)),
                                 .f = ~ CI__cellmeans_JU_model(..1, ..2, ..3, ..4, ..5)) 
        
}

CI_models_JU_preds <- function() {
    
        load(file = paste0(DATA_PATH, "modelled/data_ju.RData"))
        load(file = paste0(DATA_PATH, 'modelled/JU__mods.RData'))
        
        mods <- mods %>%
            mutate(Pred = map2(.x = REEF.d, .y = Taxa,
                               .f = ~ get(load(file = paste0(DATA_PATH,
                                                             "modelled/JU__", .x, '_', .y, '__posteriors.RData'))) %>%
                                   dplyr::select(fYEAR, REEF.d, .draw, value)),
                   Summary = pmap(.l = list(data, Pred),
                                  .f = ~  ..2 %>% posterior::as_draws() %>%
                                      group_by(fYEAR, REEF.d) %>%  
                                      tidybayes::summarise_draws(mean,
                                                                 sd,
                                                                 median,
                                                                 HDInterval::hdi) %>%
                                      left_join(..1 %>%
                                                dplyr::select(fYEAR, Site) %>%
                                                distinct()) %>%
                                      suppressMessages() %>%
                                      suppressWarnings()
                              )
               )
        pwalk(.l = list(mods$REEF.d, mods$Taxa, mods$Summary),
              .f = ~ ..3 %>% save(file = paste0(DATA_PATH, "modelled/JU__", ..1, '_', ..2, '__summary.RData')))
        
        save(mods,
             file = paste0(DATA_PATH, 'modelled/JU__preds.RData'))
     
}

# this function called internally by the following function "CI_models_JU_distance()"

CI__index_JU <- function(dat, taxa, baselines) {    
    dat %>%
        left_join(baselines %>%
                  filter(Taxa == taxa) %>%
                  dplyr::rename(baseline = value)) %>%
        mutate(
            distance.met = log2(baseline/value),
            cap.dist.met = as.numeric(case_when(distance.met < -3 ~ -3,  # The capping at log2(-3 or 3) is outrageous as represents 8* the baseline to score a 1! 
                                              distance.met > 3 ~3,
                                              distance.met > -3 & distance.met < 3 ~
                                                  distance.met)),
            rescale.dist.metric = scales::rescale(cap.dist.met, from = c(-3, 3), to = c(1, 0))) %>%
        dplyr::select(-any_of(ends_with("met"))) %>%
        pivot_longer(cols = ends_with('metric'), names_to = 'Metric', values_to = '.value') %>%
        filter(!is.na(REEF.d)) %>% 
        suppressMessages() %>%
        suppressWarnings()
}


CI_models_JU_distance <- function() {
    
  
        Baselines <- get(load(file = paste0(DATA_PATH,
                                            'modelled/JU__baseline_posteriors.RData')))
        mods <- get(load(file = paste0(DATA_PATH, "modelled/JU__preds.RData")))
        load(file=paste0(DATA_PATH, 'processed/site.location.RData'))
         # add in Acropora juvenile limits based on Manu's modelling
        load(file=paste0(DATA_PATH, 'parameters/IPM_juv.RData'))
        
        .draw=tibble(.draw=seq(from=1, to=1000, by=1))
        # updated code to include the IPM_juv estimate - noting this is a little different to the value reported in the indicators tech report.
        Acr.baseline<- site.location |> 
          dplyr::select(DEPTH.f, BIOREGION.agg, Shelf) |>  
          unique() |> 
          left_join(IPM_juv |> rename(value=mean)) |> 
          mutate(Taxa='Acropora') |> 
          cross_join(.draw) |> 
          dplyr::select(-Shelf)
        # remove density of Acropora estimated from INLA models of observed values and replace with that derived from IPM model
        baselines<- Baselines |> 
          filter(Taxa=='Total') |> 
          rbind(Acr.baseline)

        mods <- mods %>% 
                    mutate(Pred = map(.x = Pred,   
                              .f = ~ .x %>%
                                  left_join(site.location %>%
                                            dplyr::select(REEF, REEF.d, BIOREGION.agg,
                                                          DEPTH.f) %>%
                                            distinct()) %>%
                                  suppressMessages() %>%
                                  suppressWarnings())) %>% 
            mutate(Scores = map2(.x = Pred, .y = Taxa,
                                .f = ~ CI__index_JU(.x, .y, baselines) %>%
                                    filter(Metric %in% c('rescale.dist.metric',
                                                         'pcb.rescale.dist.metric')) %>%   
                                    mutate(fYEAR = factor(fYEAR, levels = unique(fYEAR))) %>%
                                    arrange(fYEAR, .draw)
                               )) 
      
      
      
   mods_summary<-mods   |>  
            dplyr::select(Scores) %>%
            unnest(Scores) %>%
            dplyr::select(-Metric) %>%
            dplyr::rename(Metric = Taxa) %>%
            group_by(REEF.d, DEPTH.f) %>%
            summarise(data = list(cur_data_all()), .groups = "drop") %>%
            dplyr::rename(Scores = data) %>%
            mutate(Summary = map(.x = Scores,
                                 .f = ~ .x %>%
                                     dplyr::select(-any_of(c(
                                                "P_CODE",
                                                "Model",
                                                "value",
                                                "baseline"))) %>%
                                     group_by(fYEAR, REEF, REEF.d, DEPTH.f, BIOREGION.agg, Metric) %>%
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

        save(mods_summary,
              file = paste0(DATA_PATH, 'modelled/JU__scores_reef_year.RData'))
        

}

CI_models_JU_aggregation <- function(level = 'NRM') {
  
  
  mods <- get(load(file = paste0(DATA_PATH, 'modelled/JU__scores_reef_year.RData')))
  spatial_lookup <- get(load(file = paste0(DATA_PATH, 'processed/spatial_lookup.RData')))
  
  
  mods <-  mods %>%
    dplyr::select(Scores) %>% 
    unnest(Scores) %>% 
    dplyr::select(fYEAR, REEF.d, .draw, Metric, .value) %>%
    left_join(spatial_lookup %>%
                dplyr::select(REEF.d, MMP.Report,!!level)  |> 
                distinct()
    ) %>%
    filter(!is.na(!!sym(level))) %>%           ## exclude all reefs outside boundary
    filter(MMP.Report=="TRUE") |> 
    droplevels() |> 
    dplyr::select(-MMP.Report) |>   
    group_by(!!sym(level), Metric, fYEAR, .draw) %>%
    summarise(.value = mean(.value)) %>%
    ungroup() |> 
    group_by(fYEAR, !!sym(level), Metric) |> 
    summarise(Score=median(.value),
              q95=quantile(.value, 0.95)) |> 
    mutate(BelowPar=ifelse(q95<0.5, 0,1),
           Indicator="Juvenile") |> 
    ungroup()
  
  save(mods,
       file = paste0(DATA_PATH, 'modelled/JU__scores_', level,'_year.RData'))
}

# some handy wrangling of nested objects
# summarise mods
# apply a number to each row to facilitate navigating the rows.
# xx<-mods |> mutate(n=1:n()) |> select(n, everything())
# 
# #get the data from a row /list item as a dataframe
# xx<-mods$Full_data[[1]] |> as.data.frame()
# 
# # scores for all reefs and taxa
# xx<-mods |> select(n, Scores) |> unnest("Scores")