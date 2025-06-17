CI_clear_models_CC_data <- function() {
    files <- list.files(path = paste0(DATA_PATH, "modelled"),
                        pattern = "CC.*|data_cc.*",
                        full.names = TRUE)
    unlink(files)
}

CI_models_CC_get_baselines <- function() {
    

        load(file=paste0(DATA_PATH, 'processed/site.location.RData'))

        ## Shallow posteriors (reef) ===========================================
        load(file = paste0(DATA_PATH, "parameters/CC_baseline.RData"))

        baselines <- mod %>%
            mutate(DEPTH.f = ifelse(DEPTH.f == 'deep', 'deep slope', 'shallow slope'))
        
        save(baselines,
             file = paste0(DATA_PATH, 'modelled/CC__baseline_posteriors.RData'))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                              item = 'get_baselines',status = 'success')

    
}


CI_models_CC_prepare_data <- function() {
    
        load(file = paste0(DATA_PATH, "processed/points.analysis.data.transect.RData"))
        load(file = paste0(DATA_PATH, "processed/spatial_lookup.RData"))

        tran.data <- points.analysis.data.transect %>%
          left_join(spatial_lookup %>%
                      dplyr::select(REEF, DEPTH, DEPTH.f, REEF.d, BIOREGION, BIOREGION.agg, NRM,Shelf) %>%
                      distinct()) %>%
          mutate(P_CODE = as.factor(P_CODE),
                 NRM = as.factor(NRM),
                 REEF.d=factor(paste(REEF, DEPTH.f)),
                Site=factor(paste(REEF.d,SITE_NO)),
                   Transect=factor(paste(Site,TRANSECT_NO)),
                   fYEAR=factor(as.numeric(as.character(REPORT_YEAR))))%>%
            group_by(Transect) %>%
            mutate(Obs=factor(1:n())) %>%
            ungroup %>%
          dplyr::select(P_CODE, Shelf,NRM, BIOREGION, BIOREGION.agg, REEF, DEPTH, DEPTH.f,
                        REEF.d, VISIT_NO, REPORT_YEAR, fYEAR, SITE_NO, Site,Transect,
                        LATITUDE, LONGITUDE, HC, A, MA, total.points)       
        
        ## Model terminated when there were reefs with only 1 level of
        ## the "REPORT_YEAR" factor. Remove Could use Murray's
        ## function in place of this
        report.year.ss <- tran.data %>%
            dplyr::select(Site, REPORT_YEAR) %>%
            unique() %>% 
            group_by(Site) %>%
            summarise(report.year.ss = n())

        df.a <- tran.data %>%
            left_join(report.year.ss) %>%
            filter(!report.year.ss == "1") %>%
            droplevels

        data_cc <- df.a %>% filter(Shelf=="Inshore") %>% droplevels()
        
        save(data_cc,
              file = paste0(DATA_PATH, 'modelled/data_cc.RData'))
        
}

CI_models_CC_prepare_nest <- function() {
    
        load(file = paste0(DATA_PATH, "modelled/data_cc.RData"))
        
        ## nest the data
       
        mods <- data_cc %>% 
          nest(data=everything(), .by=REEF.d) %>%
          mutate(n = 1:n())
        
        ## Prepare the data
        mods <- mods %>%
            mutate(newdata = map(.x = data,
                                   .f = ~ .x %>%
                                       droplevels() %>% 
                                       ## tidyr::expand(Site = Site, fYEAR = fYEAR, # replace
                                       tidyr::expand(REEF.d = REEF.d, fYEAR = fYEAR,
                                                     Site=NA, Transect = NA, Obs = NA,
                                                     HC = NA, total.points = NA) %>%
                                       distinct()
                                   ),
                   Full_data = pmap(.l = list(data, newdata),
                                    .f = ~ ..1 %>% bind_rows(..2))
                   )
        save(mods,
              file = paste0(DATA_PATH, 'modelled/CC__mods.RData'))
}

CI__fit_CC_model <- function(form, data, family='binomial', n, N) {
  environment(form)<-environment()
  reef <- unique(data$REEF.d)
        ## first try betabinomial
        res <- tryCatch(
            mod <- inla(formula = form,
                        data = data,
                        Ntrials = data$total.points,
                        family = family, 
                        ## control.family=list(link='logit'),
                        control.predictor = list(link = 1, compute = TRUE),
                        control.compute = list(
                            dic = TRUE, cpo = TRUE, waic = TRUE,
                            config = TRUE) 
                        ),
            error = function(e) e
            )
        ## else try a binomial with Obs-level RE
        if (inherits(res, "error")) {
            form1 <- update(form, .~.+f(Obs, model = 'iid'))
            mod <- inla(formula = form1,
                        data = data,
                        Ntrials = data$total.points,
                        family = 'binomial', 
                        ## control.family=list(link='logit'),
                        control.predictor = list(link = 1, compute = TRUE),
                        control.compute = list(
                            dic = TRUE, cpo = TRUE, waic = TRUE,
                            config = TRUE) 
                        )
        } else {
            mod <- res
        }
        
        save(mod, file = paste0(DATA_PATH, "modelled/CC__", reef, '__model.RData'))
        draws <- inla.posterior.sample(n=1000, mod, seed=123) %>%
            suppressWarnings() %>%
            suppressMessages()
        save(draws, file = paste0(DATA_PATH, "modelled/CC__", reef, '__draws.RData'))
    
}

CI_models_CC_fit_models <- function() {
    

        load(file = paste0(DATA_PATH, "modelled/data_cc.RData"))
        load(file = paste0(DATA_PATH, 'modelled/CC__mods.RData'))

        form <- HC ~ fYEAR +
            f(Site , model='iid') +
            f(Transect , model='iid')## +
            ## f(Obs, model = "iid")
        
        ## Fit the models - output models and draws to
        ## DATA_PATH/modelled/CC__.*__.RData
        purrr::pwalk(.l = list(mods$Full_data,mods$n),
                     .f = ~ CI__fit_CC_model(form = form, data = ..1,
                                             family = 'betabinomial', n = ..2, N = nrow(mods))
                   )
}

CI__cellmeans_CC_model <- function(obs_data, Full_data, newdata, n, N) {
    
        reef <- unique(newdata$REEF.d)
        draws <- get(load(file = paste0(DATA_PATH, "modelled/CC__", reef, '__draws.RData')))
        cellmeans <- sapply(draws, function(x)
            x[[2]][(nrow(Full_data)-nrow(newdata)+1):nrow(Full_data)]) 
        posteriors <- newdata %>%
            dplyr::select(fYEAR, REEF.d) %>%
            cbind(plogis(cellmeans)) %>%
            pivot_longer(cols = matches('[0-9]'), names_to = 'Rep') %>%
            mutate(REEF.d = reef,
                   .draw = as.integer(Rep)) %>%
            dplyr::select(-Rep) %>%
            left_join(obs_data %>%
                      dplyr::select(fYEAR, REEF.d) %>%
                      distinct()) %>% 
            suppressWarnings() %>%
            suppressMessages()
        save(posteriors, file = paste0(DATA_PATH, "modelled/CC__", reef, '__posteriors.RData'))
   
}

CI_models_CC_cellmeans <- function() {
   
        load(file = paste0(DATA_PATH, "modelled/data_cc.RData"))
        load(file = paste0(DATA_PATH, 'modelled/CC__mods.RData'))
        ## Calculate cellmeans
        cellmeans <- purrr::pwalk(.l = list(mods$data, mods$Full_data, mods$newdata, mods$n, nrow(mods)),
                                 .f = ~ CI__cellmeans_CC_model(..1, ..2, ..3, ..4, ..5)) 
}

CI_models_CC_preds <- function() {
    
        load(file = paste0(DATA_PATH, "modelled/data_cc.RData"))
        load(file = paste0(DATA_PATH, 'modelled/CC__mods.RData'))
        mods <- mods %>%
            mutate(Pred = map(.x = REEF.d,
                              .f = ~ get(load(file = paste0(DATA_PATH,
                                                            "modelled/CC__", .x, '__posteriors.RData'))) %>%
                                  dplyr::select(fYEAR, REEF.d, .draw, value)),
                   Summary = pmap(.l = list(data, Pred),
                                  .f = ~ ..2 %>% posterior::as_draws() %>%
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
        pwalk(.l = list(mods$REEF.d, mods$Summary),
              .f = ~ ..2 %>% save(file = paste0(DATA_PATH, "modelled/CC__", ..1, '__summary.RData')))

        save(mods,
              file = paste0(DATA_PATH, 'modelled/CC__preds.RData'))
        
}


CI__index_CC <- function(dat, baselines) {
    dat %>%
        left_join(baselines %>%
                  dplyr::rename(baseline = value)) %>%
        mutate(
            calc.met = plogis(log2(value/baseline)),  # we need to consider an upper bound for coral cover at which point scores are 1
            rescale.dist.metric = ifelse(value >= baseline,
                                         my_rescale(calc.met,
                                                    from = list(plogis(log2(1/baseline)), 0.5),
                                                    to = c(1, 0.5)),
                                         calc.met),
            pcb.distance.met = log2(value/0.2),
            pcb.cap.dist.met = as.numeric(case_when(pcb.distance.met < -1 ~ -1,
                                                    pcb.distance.met > 1 ~ 1,
                                                    pcb.distance.met >= -1 & pcb.distance.met <= 1 ~
                                                        pcb.distance.met)),
            pcb.rescale.dist.metric = scales::rescale(pcb.cap.dist.met,
                                                      from = c(-1,1),
                                                      to = c(0, 1))) %>%
        dplyr::select(-any_of(ends_with("met"))) %>%
        pivot_longer(cols = ends_with('metric'), names_to = 'Metric', values_to = '.value') %>%
        filter(!is.na(REEF.d)) %>% 
        suppressMessages() %>%
        suppressWarnings()
    
    
}

source('CI_30_models_functions.R') # for my_rescale()
CI_models_CC_distance <- function() {
    
        baselines <- get(load(file = paste0(DATA_PATH,
                                            'modelled/CC__baseline_posteriors.RData')))
        mods <- get(load(file = paste0(DATA_PATH, "modelled/CC__preds.RData")))
        load(file=paste0(DATA_PATH, 'processed/site.location.RData'))
        ## load(file = paste0(DATA_PATH, "processed/spatial_lookup.RData"))

        mods <- mods %>%
            mutate(Pred = map(.x = Pred,
                              .f = ~ .x %>%
                                  left_join(site.location %>%
                                            dplyr::select(REEF, REEF.d, BIOREGION.agg,
                                                          DEPTH.f) %>%
                                            distinct()) )) %>%
            mutate(Scores = map(.x = Pred,
                                .f = ~ CI__index_CC(.x, baselines) %>%
                                    filter(Metric %in% c('rescale.dist.metric',      
                                                        'pcb.rescale.dist.metric')) %>%
                                    mutate(fYEAR = factor(fYEAR, levels = unique(fYEAR))) %>%
                                    arrange(fYEAR, .draw)
                               )) %>%
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
              file = paste0(DATA_PATH, 'modelled/CC__scores_reef_year.RData'))
        
}



CI_models_CC_aggregation <- function(level = 'NRM') {
  
  
  mods <- get(load(file = paste0(DATA_PATH, 'modelled/CC__scores_reef_year.RData')))
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
           Indicator="CoralCover") |> 
    ungroup()
  
  save(mods,
       file = paste0(DATA_PATH, 'modelled/CC__scores_', level,'_year.RData'))
  
  
}


## CI_models_CC_aggregation <- function(level = 'NRM') {
##     CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
##                    item = paste0('agg_', level),
##                    label = paste0("Aggregate to ", level), status = 'pending')
##     CI_tryCatch({

##         mods <- get(load(file = paste0(DATA_PATH, 'modelled/CC__scores_reef_year.RData')))
##         spatial_lookup <- get(load(file = paste0(DATA_PATH, 'processed/spatial_lookup.RData')))

##         ## Number of reefs below 0.5
##         mods.n <- mods %>%
##             dplyr::select(REEF.d, Below) %>%
##             unnest(Below) %>%
##             left_join(spatial_lookup %>%
##                       dplyr::select(REEF.d, !!level) %>%
##                       distinct()
##                       ) %>%
##             filter(!is.na(!!sym(level))) %>%           ## exclude all reefs outside boundary
##             group_by(fYEAR, Metric, !!sym(level)) %>%
##             summarise(n.below = sum(Below),
##                       n.Pbelow = sum(PBelow),
##                       tn.reefs = n()) %>%
##             ungroup() %>%
##             group_by(!!sym(level)) %>%
##             nest() %>%
##             dplyr::rename(Below = data) %>%
##             ungroup() %>%
##             suppressMessages() %>%
##             suppressWarnings()

##         ## Scores
##         mods <- mods %>%
##             dplyr::select(Scores) %>%
##             unnest(Scores) %>%
##             dplyr::select(fYEAR, REEF.d, .draw, Metric, .value) %>%
##             left_join(spatial_lookup %>%
##                       dplyr::select(REEF.d, !!level) %>%
##                       distinct()
##                       ) %>%
##             filter(!is.na(!!sym(level))) %>%           ## exclude all reefs outside boundary
##             group_by(!!sym(level)) %>%
##             summarise(data = list(cur_data_all()), .groups = "drop") %>%
##             mutate(Scores = map(.x = data,
##                                 .f = ~ .x %>%
##                                     ungroup() %>%
##                                     group_by(fYEAR, !!sym(level), .draw, Metric) %>%
##                                     summarise(.value = mean(.value))
##                                 ),
##                    Summary = map(.x = Scores,
##                                  .f = ~ .x %>%
##                                      group_by(fYEAR, !!sym(level), Metric) %>%
##                                      summarise_draws(median, mean, sd,
##                                                      HDInterval::hdi,
##                                                      `p<0.5` = ~ mean(.x < 0.5))
##                                  )
##                    ) %>%
##             suppressMessages() %>%
##             suppressWarnings()

##         ## Combine
##         mods <- mods %>%
##             left_join(mods.n) %>%
##             mutate(Summary = map2(.x = Summary, .y = Below,
##                                   .f = ~ .x %>% left_join(.y))) %>%
##             suppressMessages() %>%
##             suppressWarnings()

##         save(mods,
##               file = paste0(DATA_PATH, 'modelled/CC__scores_', level,'_year.RData'))

##         CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
##                               item = paste0('agg_', level), status = 'success')

##     }, logFile=LOG_FILE, Category='--CC models--',
##     msg=paste0('Aggregate to ', level), return=NULL)
## }
