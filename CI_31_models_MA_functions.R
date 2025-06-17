CI_clear_models_MA_data <- function() {
    files <- list.files(path = paste0(DATA_PATH, "modelled"),
                        pattern = "MA.*|data_ma.*",
                        full.names = TRUE)
    unlink(files)
}

CI_models_MA_get_baselines <- function() {
   
        load(file=paste0(DATA_PATH, 'processed/site.location.RData'))

        ## Shallow posteriors (reef) ===========================================
        #load(file = paste0(DATA_PATH, "parameters/MA__baseline_shallow.RData"))

        ## Get the set of reefs to predict for and
        ## nest according to whether they should be predicted against:
        ## 1. inshore shallow model
        ## 2. inshore deep model
        ## 3. offshore model
        newdata <- site.location %>%
            mutate(Model = ifelse(Shelf == 'Offshore',
                                  'offshore',
                           ifelse(DEPTH.f == 'shallow slope',
                                  'shallow',
                                  'deep'))) %>%
            dplyr::group_by(Model) %>%
            nest() %>%
            mutate(data = map2(.x = data, .y = Model,
                               .f = ~ {
                                   model <- ifelse(.y == 'offshore', 'BIOREGION.agg', 'REEF')
                                   .x %>% 
                                       dplyr::group_by_at(model) %>%
                                       summarise(across(c(LATITUDE, LONGITUDE), mean))
                                   }
                              )
                   )

        baselines <- newdata %>%
            mutate(Preds = map2(.x = data, .y = Model,
                                .f = ~ {
                                    load(file = paste0(DATA_PATH, "parameters/MA__baseline_",.y,".RData"))
                                    .x %>%
                                    mutate(Case = 1:n()) %>%
                                    left_join(CI__make_predictions(.x, mod, mesh) %>%
                                              as.data.frame() %>%
                                              pivot_longer(cols = everything(), names_to = 'Case') %>%
                                              mutate(Case = as.integer(Case)) %>%
                                              group_by(Case) %>%
                                              mutate(.draw = 1:n()) %>%
                                              ungroup) %>%
                                        (if(.y == "shallow") { 
                                             . %>% mutate(REEF.d = factor(paste(REEF, 'shallow slope'))) %>%
                                             dplyr::select(REEF.d, LATITUDE, LONGITUDE, .draw, value)
                                         } else if(.y == 'deep') {
                                             . %>% mutate(REEF.d = factor(paste(REEF, 'deep slope'))) %>%
                                             dplyr::select(REEF.d, LATITUDE, LONGITUDE, .draw, value)
                                         } else {
                                             . %>% dplyr::select(BIOREGION.agg, LATITUDE, LONGITUDE, .draw, value) 
                                        }) %>% 
                                        left_join(site.location %>%
                                                  dplyr::select(P_CODE, REEF, REEF.d, BIOREGION.agg) %>%
                                                  distinct()) %>%
                                        suppressMessages() %>%
                                        suppressWarnings()
                                }
                                )
                   ) %>%
            dplyr::select(Preds) %>%
            unnest(c(Preds)) %>% 
            ungroup() %>% 
            suppressMessages() %>%
            suppressWarnings()
                                    
        save(baselines,
             file = paste0(DATA_PATH, 'modelled/MA__baseline_posteriors.RData'))

        
}

CI_models_MA_prepare_data <- function() {
    
        load(file = paste0(DATA_PATH, "processed/points.analysis.data.transect.RData"))
        load(paste0(DATA_PATH, 'processed/spatial_lookup.RData'))
        
        data <- points.analysis.data.transect %>%
            mutate(DEPTH.f = factor(case_when(DEPTH > 3 ~"deep slope",
                                              DEPTH <= 3 ~"shallow slope")),
                   REEF.d = factor(paste(REEF, DEPTH.f)),
                   Site = factor(paste(REEF.d, SITE_NO)),
                   Transect = factor(paste(Site, TRANSECT_NO)),
                   fYEAR = factor(as.numeric(as.character(REPORT_YEAR)))) %>%
            group_by(Site) %>%
            mutate(Obs=factor(1:n())) %>%
            ungroup() %>%
          left_join(spatial_lookup |> 
                      dplyr::select(REEF.d, Shelf)  |> 
                      distinct()) |> 
          filter(Shelf=="Inshore") |> 
            suppressMessages() %>%
            suppressWarnings()

        ## Model terminated when there were reefs with only 1 level of
        ## the "REPORT_YEAR" factor. Remove
        report.year.ss <- data %>% 
            dplyr::select(Site, fYEAR) %>% 
            distinct() %>% 
            group_by(Site) %>% 
            summarise(report.year.ss = n()) %>%
            ungroup() %>%
            suppressMessages() %>%
            suppressWarnings()

        df.a <- data %>%
            left_join(report.year.ss) %>% 
            filter(!report.year.ss == "1") %>% 
            droplevels() %>%
            suppressMessages() %>%
            suppressWarnings()

        data = df.a

        save(data,
              file = paste0(DATA_PATH, 'modelled/data_ma.RData'))
}

# in the below site was added to .f
CI_models_MA_prepare_nest <- function() {
   
        load(file = paste0(DATA_PATH, "modelled/data_ma.RData"))
        
        mods <- data %>% 
          nest(data=everything(), .by=REEF.d) %>%
          mutate(n = 1:n())
        
        ## Prepare the data
        mods <- mods %>%
            mutate(newdata = map(.x = data,
                                   .f = ~ .x %>%
                                       droplevels() %>% 
                                       #tidyr::expand(Site = Site, fYEAR = fYEAR, # replace
                                       tidyr::expand(REEF.d = REEF.d, fYEAR = fYEAR,
                                                     Site=NA, Transect = NA, 
                                                     Obs = NA,
                                                     MA = NA, A = NA) %>%
                                       distinct()
                                   ),
                   Full_data = pmap(.l = list(data, newdata),
                                    .f = ~ ..1 %>% bind_rows(..2))
                   )
        save(mods,
              file = paste0(DATA_PATH, 'modelled/MA__mods.RData'))
        
}

CI__best_model <- function(model_type) {
    disp <- lapply(model_type, `[[`, "diags") %>% lapply(`[[`, "Disp")
    dispp <- lapply(model_type, `[[`, "diags") %>% lapply(`[[`, "Dispp")
    ks <- lapply(model_type, `[[`, "diags") %>% lapply(`[[`, "KS")
    ksp <- lapply(model_type, `[[`, "diags") %>% lapply(`[[`, "KSp")
    waic <- lapply(model_type, `[[`, "waic") %>% lapply('[[', "WAIC")
    mse <- lapply(model_type, `[[`, "mse") %>% lapply('[[', "MSE")
    ## first vet by evidence of issues with KS or dispersion
    wch <- which(unlist(ksp)>0.05 &
                 unlist(dispp)>0.05)
    ## if all fail (not withstanding that these are anti-conservative tests)
    ## retain all models
    if (length(wch) == 0) wch <- names(model_type) 
    ## Compare dispersion and WAIC (weighted most to WAIC)
    r <- rbind(
        rank(-abs(1-unlist(disp)))*0.1,
        rank(-unlist(waic[wch]))
        ## rank(-unlist(mse[wch]))
    ) %>% colSums()
    ## Select best model
    best <- names(which.max(r))
    best    
}

CI__model_MA_diagnostics <- function(mod, obsdata, type, m, reef) {
    ## Information Criterion
    waic <- list(WAIC = mod$waic$waic)
    ## DHARMa diagnostics
    preds <- posterior_predict.inla(mod, newdata = obsdata, ndraws = 250, new_random_levels = FALSE) %>% t()
    fitted_median_inla <- apply(preds, 1, mean)
    ## Check whether there are any missing values in either
    ## fitted_median_inla or obsdata.
    if (any(is.na(obsdata$MA))) { 
        wch <- which(is.na(obsdata$MA))
        preds <- preds[-wch,]
        obsdata <- obsdata[-wch,]
        fitted_median_inla <- fitted_median_inla[-wch]
    }
    mod.resids <- DHARMa::createDHARMa(
                              simulatedResponse = preds,
                              observedResponse = obsdata$MA,
                              fittedPredictedResponse = fitted_median_inla,
                              integerResponse = TRUE
                          )
    ggsave(paste0(FIGS_PATH, "/", type, "_", m, "_DHARMa_", reef, ".png"),
           patchwork::wrap_elements(~(ks <- DHARMa::testUniformity(mod.resids))) +
           patchwork::wrap_elements(~DHARMa::plotResiduals(mod.resids)) +
           patchwork::wrap_elements(~(disp <- DHARMa::testDispersion(mod.resids))),
           width = 12, height = 4)
    diags <- list(
        KS =  ks$statistic[[1]],
        KSp =  ks$p.value[[1]],
        Disp = disp$statistic[[1]],
        Dispp = disp$p.value[[1]]
    )
    mse <- list(MSE = mean((fitted_median_inla - obsdata$MA)^2))
    list(mod = mod,
         waic = waic,
         diags = diags,
         mse = mse
         )
}

## No longer need to pass in form and family - remove these arguments
CI__fit_MA_model <- function(fulldata, obsdata, n, N) {
    
        ## site <- unique(data$Site)
        type <- "MA"
        reef <- unique(fulldata$REEF.d)
        
       # environment(form)<-environment()
        model_type <- list(
            list(form = MA ~ fYEAR + f(Site , model='iid') + f(Transect , model='iid'),
                 family = "binomial"),
            list(form = MA ~ fYEAR + f(Site , model='iid') + f(Transect , model='iid') +
                     f(Obs, model='iid'),
                 family = "binomial"),
            list(form = MA ~ fYEAR + f(Site , model='iid') + f(Transect , model='iid'),
                 family = "betabinomial")
         ) %>% setNames(c("binomial","binomialURE", "betabinomial"))
        ## Fit each of the candidate models and determine which is "best"
        for (m in names(model_type)) {
           
            ## fit the model
            mod <- tryCatch({
                inla(formula = model_type[[m]]$form,
                            data = fulldata,
                            Ntrials = fulldata$A,
                            family = model_type[[m]]$family, 
                            ## control.family=list(link='logit'),
                            control.predictor = list(link = 1, compute = TRUE),
                            control.compute = list(
                                dic = TRUE, cpo = TRUE, waic = TRUE,
                                config = TRUE) 
                            )
            }, error = function(e) e)
           if (inherits(mod, "simpleError")) {
                model_type[[m]] <- NULL
            } else {
                ## Collect model diagnostics/characteristics This will
                ## also save the DHARMa plots in
                ## FIGS_PATH/<type>_<m>_DHARMa.*
                model_type[[m]] <- model_type[[m]] %>%
                    append(CI__model_diagnostics(mod, obsdata, type, m, reef))
            }
        }
        
        best_model <- CI__best_model(model_type)
        mod <- model_type[[best_model]]$mod
        
        
        save(mod, file = paste0(DATA_PATH, "modelled/MA__", reef, '__model.RData'))
        draws <- inla.posterior.sample(n=1000, mod, seed=123) %>%
            suppressWarnings() %>%
            suppressMessages()
        save(draws, file = paste0(DATA_PATH, "modelled/MA__", reef, '__draws.RData'))
    
}

## No longer need to pass in form and family - remove these arguments
CI__fit_MA_model <- function(fulldata, obsdata, n, N) {
    
        ## site <- unique(data$Site)
        type <- "MA"
        reef <- unique(fulldata$REEF.d)
        
       # environment(form)<-environment()
        model_type <- list(
            list(form = MA ~ fYEAR + f(Site , model='iid') + f(Transect , model='iid'),
                 family = "binomial"),
            list(form = MA ~ fYEAR + f(Site , model='iid') + f(Transect , model='iid') +
                     f(Obs, model='iid'),
                 family = "binomial"),
            list(form = MA ~ fYEAR + f(Site , model='iid') + f(Transect , model='iid'),
                 family = "betabinomial")
         ) %>% setNames(c("binomial","binomialURE", "betabinomial"))
        ## Fit each of the candidate models and determine which is "best"
        for (m in names(model_type)) {
           
            ## fit the model
            mod <- tryCatch({
                inla(formula = model_type[[m]]$form,
                            data = fulldata,
                            Ntrials = fulldata$A,
                            family = model_type[[m]]$family, 
                            ## control.family=list(link='logit'),
                            control.predictor = list(link = 1, compute = TRUE),
                            control.compute = list(
                                dic = TRUE, cpo = TRUE, waic = TRUE,
                                config = TRUE) 
                            )
            }, error = function(e) e)
           if (inherits(mod, "simpleError")) {
                model_type[[m]] <- NULL
            } else {
                ## Collect model diagnostics/characteristics This will
                ## also save the DHARMa plots in
                ## FIGS_PATH/<type>_<m>_DHARMa.*
                model_type[[m]] <- model_type[[m]] %>%
                    append(CI__model_diagnostics(mod, obsdata, type, m, reef))
            }
        }
        
        best_model <- CI__best_model(model_type)
        mod <- model_type[[best_model]]$mod
        
        
        save(mod, file = paste0(DATA_PATH, "modelled/MA__", reef, '__model.RData'))
        draws <- inla.posterior.sample(n=1000, mod, seed=123) %>%
            suppressWarnings() %>%
            suppressMessages()
        save(draws, file = paste0(DATA_PATH, "modelled/MA__", reef, '__draws.RData'))
    
}
## use a cell means model as per juveniles 
## this not used for framework or case study but should check with Murray to see if maybe it would be better
# CI__fit_MA_model<- function(form, data, n, N) {
#   
#   reef <- unique(data$REEF.d)
#   environment(form)<-environment()
#   
#   mod <- inla(formula = form,
#               Ntrials = fulldata$A,
#               data = fulldata,
#               family = "binomial", 
#               control.predictor = list(link = 1, compute = TRUE),
#               control.compute = list(
#                 dic = TRUE, cpo = TRUE, waic = TRUE,
#                 config = TRUE) 
#   )
#   
#   save(mod, file = paste0(DATA_PATH, "modelled/MA__", reef, '_', taxa, '__model.RData'))
#   draws <- inla.posterior.sample(n=1000, mod, seed=123) %>%
#     suppressWarnings() %>%
#     suppressMessages()
#   save(draws, file = paste0(DATA_PATH, "modelled/MA__", reef, '_', taxa, '__draws.RData'))
# }
# 



CI__cellmeans_MA_model <- function(obs_data, Full_data, newdata, n, N) {
    
        reef <- unique(newdata$REEF.d)
        
       
        draws <- get(load(file = paste0(DATA_PATH, "modelled/MA__", reef, '__draws.RData')))
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
        save(posteriors, file = paste0(DATA_PATH, "modelled/MA__", reef, '__posteriors.RData'))
    
}



CI_models_MA_fit_models <- function() {
   
    #CI_tryCatch({

        ## load(file = paste0(DATA_PATH, "modelled/data_ma.RData"))
        load(file = paste0(DATA_PATH, 'modelled/MA__mods.RData'))

        ## For each reef/depth, fit 3 alternative models (binomial,
        ## binomial with unit-level RE, and beta-binomial) and select
        ## the "best" one.  The model and draws are saved to
        ## DATA_PATH/modelled/MA__.*__.RData
        purrr::pwalk(.l = list(mods$Full_data, mods$data, mods$n),
                     .f = ~ CI__fit_MA_model(fulldata = ..1, obsdata = ..2,
                                             n = ..3, N = nrow(mods))
                   )

       
   # }, logFile=LOG_FILE, Category='--MA models--',
   # msg=paste0('Fit MA models'), return=NULL)
}

## using cell means model as per juveniles 
## this not used for framework or case study but should check with Murray to see if maybe it would be better
# CI_models_MA_fit_models<- function() {
#   
#   load(file = paste0(DATA_PATH, "modelled/data_ma.RData"))
#   load(file = paste0(DATA_PATH, 'modelled/MA__mods.RData'))  
#   
#  # form<- MA ~ fYEAR + f(Site , model='iid') + f(Transect , model='iid') +
#  #    f(Obs, model='iid')
#  #  
#   form <- MA ~ 0+fYEAR +   ## added the 0+ to make a cell means formula
#     f(Site , model='iid')+
#     + f(Transect , model='iid')
#   
#   ## Fit the models - output models and draws to
#   purrr::pwalk(.l = list(mods$Full_data, mods$n),
#                .f = ~ CI__fit_MA_model(form = form, 
#                                        data = ..1,
#                                        family = 'binomial',
#                                        n = ..3,
#                                        N = nrow(mods))
#   )
#   
# }

CI__DHARMa <- function(type, reef, preds, data, fitted_median_inla) {
    mod.resids <- DHARMa::createDHARMa(
                              simulatedResponse = preds,
                              observedResponse = data$MA,
                              fittedPredictedResponse = fitted_median_inla,
                              integerResponse = TRUE
                          )
    ggsave(paste0(FIGS_PATH, "/", type, "_DHARMa_", reef, ".png"),
           patchwork::wrap_elements(~DHARMa::testUniformity(mod.resids)) +
           patchwork::wrap_elements(~DHARMa::plotResiduals(mod.resids)) +
           patchwork::wrap_elements(~DHARMa::testDispersion(mod.resids)),
           width = 12, height = 4)
}

# 
CI__diagnostics_MA_model <- function(REEF.d, data, n, N) {
    
        type <- "MA"
        reef <- as.character(unlist(REEF.d))
      

        ## PIT
        mod <- get(load(file = paste0(DATA_PATH, "modelled/MA__", reef, '__model.RData')))
        ggsave(paste0(FIGS_PATH, "/", type, "_PIT_", reef, ".png"),
               wrap_plots(
                   pit_qq_plot(mod, i.mod = 1:nrow(data), logit_scale = TRUE),
                   pit_plot(mod, i.mod = 1:nrow(data)),
                   pit_resid_plot(mod, i.mod = 1:nrow(data)),
                   pit_hist_plot(mod, i.mod = 1:nrow(data))) +
               plot_layout(ncol = 3),
               width = 15,
               height = 6)

        ## DHARMa plots
        preds <- posterior_predict.inla(mod, newdata = data, ndraws = 250, new_random_levels = FALSE) %>% t()
        fitted_median_inla <- apply(preds, 1, mean)
        CI__DHARMa(type = type, reef = reef, preds, data, fitted_median_inla)

}

CI_models_MA_diagnostics <- function() {
  

        ## Generate diagnostics
        purrr::pwalk(.l = list(mods$REEF.d, mods$data, mods$n, nrow(mods)),
                     .f = ~ CI__diagnostics_MA_model(..1, ..2,
                                                     n = ..3,
                                                     N = ..4))

    
}





CI_models_MA_cellmeans <- function() {
  
        load(file = paste0(DATA_PATH, "modelled/data_ma.RData"))
        load(file = paste0(DATA_PATH, 'modelled/MA__mods.RData'))
        ## Calculate cellmeans
        cellmeans <- purrr::pwalk(.l = list(mods$data, mods$Full_data, mods$newdata, mods$n, nrow(mods)),
                                 .f = ~ CI__cellmeans_MA_model(..1, ..2, ..3, ..4, ..5)) 
        
      
       
   
}

CI_models_MA_preds <- function() {
   
   
        load(file = paste0(DATA_PATH, "modelled/data_ma.RData"))
        load(file = paste0(DATA_PATH, 'modelled/MA__mods.RData'))
        mods <- mods %>%
            mutate(Pred = map(.x = REEF.d,
                              .f = ~ get(load(file = paste0(DATA_PATH,
                                                            "modelled/MA__", .x, '__posteriors.RData'))) %>%
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
              .f = ~ ..2 %>% save(file = paste0(DATA_PATH, "modelled/MA__", ..1, '__summary.RData')))

        save(mods,
              file = paste0(DATA_PATH, 'modelled/MA__preds.RData'))
        
       
}

# CI_models_MA_partialplots <- function() {
#     CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
#                    item = 'partial',
#                    label = "Generate partial plots", status = 'pending')
#     CI_tryCatch({
# 
#         load(file = paste0(DATA_PATH, "modelled/data_ma.RData"))
#         load(file = paste0(DATA_PATH, 'modelled/MA__mods.RData'))
# 
#         ## Raw means
#         mods <- mods %>%
#             mutate(RawMeans = map(.x = data,
#                                   .f = ~ .x %>% group_by(fYEAR) %>%
#                                       summarise(Value = mean(MA/A, na.rm = TRUE))
#                                   )) %>%
#             mutate(ModelledMeans = map(.x = REEF.d,
#                                        .f = ~ get(load(file = paste0(DATA_PATH,
#                                                                      "modelled/MA__", .x,
#                                                                      "__summary.RData")))))
#         ##Continue from here - add partial plot with modelled trends overlayed onto raw means 
#         
#         CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
#                               item = 'partials',status = 'success')
# 
#     }, logFile=LOG_FILE, Category='--MA models--',
#     msg=paste0('Generate partial plots'), return=NULL)
# }
#                                                        
CI__index_MA <- function(dat, baselines, consequence) {
    dat %>%
        left_join(baselines %>%
                  dplyr::select(-any_of(c("LONGITUDE", "LATITUDE"))) %>%
                  dplyr::rename(baseline = value)) %>%
        mutate(
            calc.met = plogis(log2(baseline/value)),
            distance.metric = ifelse(value >= baseline,
                                     my_rescale(calc.met,
                                                from = list(plogis(log2(baseline/1)), 0.5),
                                                to = c(0, 0.5)),
                                     calc.met),
            consequence.metric = ifelse(value <= consequence,
                                        scales::rescale(value,
                                                      from = c(0, consequence), 
                                                      to = c(1, 0)),
                                        0),
            combined.metric = (distance.metric + consequence.metric)/2 ) %>%
            dplyr::select(-calc.met) |> 
        pivot_longer(cols = ends_with('metric'), names_to = 'Metric', values_to = '.value') %>%  # this pivot produces nonsense duplicated rows for Metric=='consequence.metric and variable==calc.met
        filter(!is.na(REEF)) %>% 
        suppressMessages() %>%
        suppressWarnings()
}


CI_models_MA_distance <- function() {
   
        load(file=paste0(DATA_PATH, 'processed/site.location.RData'))
        baselines <- get(load(file = paste0(DATA_PATH,
                                            'modelled/MA__baseline_posteriors.RData')))
        mods <- get(load(file = paste0(DATA_PATH, "modelled/MA__preds.RData")))
        ## At the moment ma_from_juv_consequence is just the upper CI
        ## of MA values inverse predicted from a model relating
        ## Acropora juvenile density to MA cover
      
        ma_from_juv_consequence <- get(load(file = paste0(DATA_PATH,
                                                          'modelled/ma_from_juv_consequence.RData'))) %>%
            pull(ymax) %>%
            round(2)

  # note my running of the code (21-MAr-25) returns a value of ma_from_juv_consequence of 0.28 as opposed to the value 
  #  of 0.22 reported in the final methods document....
        mods <- mods %>%
          mutate(Pred = map(.x = Pred,
                              .f = ~ .x %>%
                                  left_join(site.location %>%
                                            dplyr::select(REEF, REEF.d, BIOREGION.agg,
                                                          DEPTH.f) %>%
                                            distinct()) )) %>% 
            mutate(Scores = map(.x = Pred,
                                .f = ~ CI__index_MA(.x, baselines, ma_from_juv_consequence) %>%
                                    filter(Metric %in% c('distance.metric',
                                                         'consequence.metric')) %>%
                                    mutate(fYEAR = factor(fYEAR, levels = unique(fYEAR))) %>%
                                    arrange(fYEAR, .draw)
                               )) %>%
            dplyr::select(-data, -newdata,-Full_data, -Pred, -Summary) %>%
            mutate(Summary = map(.x = Scores,
                                 .f = ~ .x %>%
                                     dplyr::select(-any_of(c(
                                                "P_CODE",
                                                "Model",
                                                "value",
                                                "baseline",
                                                "DEPTH.f"))) %>%
                                     group_by(fYEAR, REEF, REEF.d, BIOREGION.agg, Metric) %>%
                                     tidybayes::summarise_draws( median, mean, sd,
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
              file = paste0(DATA_PATH, 'modelled/MA__scores_reef_year.RData'))
        
}


# CI_models_MA_varify_scores <- function() {
#     CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
#                    item = 'varify_scores',
#                    label = "Varify scores", status = 'pending')
#     CI_tryCatch({
# 
#         mods <- get(load(file = paste0(DATA_PATH, 'modelled/MA__scores_reef_year.RData')))
# 
#         pwalk(.l = list(mods$REEF.d, mods$Scores),
#               .f = ~ {
#                   g1 <- ..2 %>%
#                       filter(Metric == 'distance.metric') %>%
#                       droplevels() %>% 
#                       ggplot(aes(x = value)) +
#                       geom_density(aes(fill = 'Cover'), alpha = 0.5, trim = TRUE) +
#                       geom_density(aes(x = baseline, fill = 'Baseline'), alpha = 0.5, trim = TRUE) +
#                       facet_grid(fYEAR ~ ., scales = 'free') +
#                       theme_bw() +
#                       theme(strip.text.y = element_text(angle = 0)) +
#                       ggtitle(..1)
#                   g2 <- ..2 %>%
#                       filter(Metric == 'distance.metric') %>%
#                       droplevels() %>% 
#                       ggplot(aes(x = .value)) +
#                       geom_density(aes(fill = 'Distance metric')) +
#                       facet_grid(fYEAR ~ ., scales = 'free') +
#                       theme_bw() +
#                       theme(strip.text.y = element_text(angle = 0)) +
#                       ggtitle(..1)
# 
#                   ggsave(filename = paste0(FIGS_PATH, '/MA__varify_scores_',..1,'.png'),
#                          g1 + g2,
#                          width = 10, height = 10)
#               }
#               )
#         zip(zipfile = paste0(FIGS_PATH, "/MA__varify_scores.zip"),
#             files = list.files(path = paste0(FIGS_PATH),
#                                pattern = "MA__varify_scores_.*.png",
#                                full.names = TRUE),
#             flags = "-9rX")
#         CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
#                               item = 'varify_scores',status = 'success')
# 
#     }, logFile=LOG_FILE, Category='--Data processing--',
#     msg=paste0('Varify scores'), return=NULL)
# }


 CI_models_MA_aggregation <- function(level = 'NRM') {
  

         mods <- get(load(file = paste0(DATA_PATH, 'modelled/MA__scores_reef_year.RData')))
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
           filter(MMP.Report=="TRUE") |>  # this is a filter for MMP case study
           droplevels() |> 
           dplyr::select(-MMP.Report) |>   
           group_by(!!sym(level), Metric, fYEAR, .draw) %>%
           summarise(.value = mean(.value)) %>%
           ungroup() |> 
           group_by(fYEAR, !!sym(level), Metric) |> 
          summarise(Score=median(.value),
                    q95=quantile(.value, 0.95)) |> 
           mutate(BelowPar=ifelse(q95<0.5, 0,1),
                  Indicator="Macroalgae") |> 
         ungroup()

         save(mods,
               file = paste0(DATA_PATH, 'modelled/MA__scores_', level,'_year.RData'))


 }

 
# CI_template <- function() {
#     CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
#                    item = 'template',
#                    label = "Template", status = 'pending')
#     CI_tryCatch({
# 
#         load(file = paste0(DATA_PATH, "modelled/template.RData"))
#         
#         ## save(template,
#         ##       file = paste0(DATA_PATH, 'modelled/template.RData'))
#         
#         CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
#                               item = 'template',status = 'success')
# 
#     }, logFile=LOG_FILE, Category='--Data processing--',
#     msg=paste0('Template'), return=NULL)


