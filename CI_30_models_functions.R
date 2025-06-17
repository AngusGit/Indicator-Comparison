CI__make_predictions <- function(newdata, mod, mesh) {
    draws <- inla.posterior.sample(1000, result=mod, seed=123) %>%
        suppressWarnings()

    points.grid <- newdata %>% 
        dplyr::select(LONGITUDE, LATITUDE)
    ## sets up an interpolator to actual location points from our new data
    proj.grid <- inla.mesh.projector(mesh, loc = as.matrix(points.grid))

    ## Rearrange draws so that it is a matrix of cellmeans, rather than a list
    cellmeans <- sapply(draws, function(x) x[['latent']])

    ## Index the cell means for fixed effects
    i.mod <- sapply(c('APredictor','^Predictor','spatial.field','Site','Intercept',
                      'fYEAR'),
                    function(x) grep(x, draws[[1]]$latent %>% rownames, perl=TRUE))

    ## get the partial predictions for the spatial field (all nodes on the mesh)
    cellmeans.spatial <- cellmeans[i.mod[["spatial.field"]],] 

    ##Get the predictions for fixed effects (covariate).  These
    ## are the intercept and the partial effects of the covariates
    ## Generate model matrix
    Xmat <- model.matrix(~1, data = newdata)

    wch <- grep('Intercept', names(i.mod))
    ii = unlist(i.mod[wch]) 

    ## multiply the predictions by the fixed effects for the covariates
    cellmeans.full.1 <- (cellmeans[ii,]) %*% t(Xmat)

    ## inla.mesh.project uses proj.grid to convert from mesh nodes
    ## to points beyond nodes the conversion is applied to the
    ## parameters from cellmeans.spatial transpose from rows to
    ## columns
    cellmeans.spatial <- t(inla.mesh.project(proj.grid, field = cellmeans.spatial))

    ## add the fixed and spatial effects together
    cellmeans.full.2 <- cellmeans.full.1 + cellmeans.spatial
    
    ## Backtransform
    if (any(grepl("poisson|nbinomial", mod$.args$family))) {
        cellmeans.spatial.2 <- cellmeans.full.2 %>%
            as.matrix() %>%
            exp()
    } else {
        cellmeans.spatial.2 <- cellmeans.full.2 %>%
            as.matrix() %>%
            plogis()
    }

    cellmeans.spatial.2
}


CI__cellmeans_model <- function(Indicator, obs_data, Full_data, newdata, n, N) {
    
        draws <- get(load(file = paste0(DATA_PATH, "modelled/", Indicator, "__", reef, '__draws.RData')))
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
        save(posteriors, file = paste0(DATA_PATH, "modelled/", Indicator, "__", reef, '__posteriors.RData'))
    
}

# CI_models_cellmeans <- function(Indicator = 'CC') {
#     
#         load(file = paste0(DATA_PATH, "modelled/data_ma.RData"))
#         load(file = paste0(DATA_PATH, 'modelled/MA__mods.RData'))
#         ## Calculate cellmeans
#         cellmeans <- purrr::pwalk(.l = list(mods$data, mods$Full_data, mods$newdata, mods$n, nrow(mods)),
#                                  .f = ~ CI__cellmeans_model(Indicator, ..1, ..2, ..3, ..4, ..5)) 
# }        


CI_models_aggregation <- function(Indicator = 'CC', level = 'NRM') {
    
        mods <- get(load(file = paste0(DATA_PATH, 'modelled/', Indicator, '__scores_reef_year.RData')))
        spatial_lookup <- get(load(file = paste0(DATA_PATH, 'processed/spatial_lookup.RData')))

        ## Number of reefs below 0.5
        mods.n <- mods %>%
            dplyr::select(REEF.d, Below) %>%
            unnest(Below) %>%
            left_join(spatial_lookup %>%
                      dplyr::select(REEF.d, !!level, Shelf) %>%
                      distinct()
                      ) %>%
            filter(!is.na(!!sym(level))) %>%           ## exclude all reefs outside boundary 
            nest(data = everything()) %>%
            mutate(WithShelf = map(.x = data,
                                   .f = ~ .x %>% 
                                       group_by(fYEAR, Metric, !!sym(level),
                                                Shelf, across(any_of('Taxa'))) %>%
                                       summarise(n.below = sum(Below, na.rm = TRUE),
                                                 n.Pbelow = sum(PBelow, na.rm = TRUE),
                                                 tn.reefs = n()) 
                                   ),
                   NoShelf = map(.x = data,
                                 .f = ~ .x %>% 
                                     mutate(Shelf = "All") %>% 
                                     group_by(fYEAR, Metric, !!sym(level),
                                              Shelf, across(any_of('Taxa'))) %>%
                                     summarise(n.below = sum(Below, na.rm = TRUE),
                                               n.Pbelow = sum(PBelow, na.rm = TRUE),
                                               tn.reefs = n()) 
                                 ),
                   Combined = map2(.x = NoShelf, .y = WithShelf,
                                   .f = ~ .x %>% rbind(.y))
                   ) %>%
            dplyr::select(Combined) %>% 
            unnest(Combined) %>%
            group_by(!!sym(level)) %>%
            nest() %>%
            dplyr::rename(Below = data) %>% 
            ungroup() %>%
            suppressMessages() %>%
            suppressWarnings()

        ## Scores 
        mods <- b <- mods %>%
            dplyr::select(Scores) %>% 
            unnest(Scores) %>% 
            dplyr::select(fYEAR, REEF.d, .draw, any_of('Taxa'), Metric, .value) %>%
            left_join(spatial_lookup %>%
                      dplyr::select(REEF.d, !!level, Shelf) %>%
                      distinct()
                      ) %>%
            filter(!is.na(!!sym(level))) %>%           ## exclude all reefs outside boundary 
            group_by(!!sym(level)) %>%
            summarise(data = list(cur_data_all()), .groups = "drop") %>% 
            ## if the supplied posterior only contains a single draw,
            ## duplicate this 1000 times
            mutate(data = map(.x = data,
                             .f = ~ { if (max(.x$.draw == 1)) {
                                         .x %>%
                                             group_by(fYEAR, REEF.d, !!sym(level),
                                                      across(any_of('Taxa')), Shelf, Metric) %>%
                                             sample_n(size = 1000, replace = TRUE) %>%
                                             mutate(.draw = 1:n())

                                      } else .x
                                      })) %>%
            mutate(
                ## Aggregate scores marginalising over shelf
               
                Scores = map(.x = data,
                                .f = ~ .x %>%
                                    ungroup() %>% 
                                    group_by(fYEAR, !!sym(level), .draw,
                                             across(any_of('Taxa')), Metric) %>%
                                    summarise(.value = mean(sample(.value, replace = TRUE),
                                                            na.rm = TRUE)) %>%
                                    mutate(Shelf = 'All')
                                ),
                
                ## Aggregate scores to shelf level
                
                ScoresShelf = map(.x = data,
                             .f = ~ .x %>%
                                 ungroup() %>% 
                                 group_by(fYEAR, !!sym(level), .draw,
                                          across(any_of('Taxa')), Shelf, Metric) %>%
                                 summarise(.value = mean(sample(.value, replace = TRUE),
                                                         na.rm = TRUE)) 
                             ),
                
                ## Combine together
                Scores = map2(.x = Scores, .y = ScoresShelf,
                              .f = ~ .x %>% rbind(.y)),
                ## Summarise
                Summary = map(.x = Scores,
                              .f = ~ .x %>%
                                  group_by(fYEAR, !!sym(level),
                                           across(any_of('Taxa')), Shelf, Metric) %>%
                                  summarise_draws(median = ~ median(.x, na.rm = TRUE),
                                                  mean = ~ mean(.x, na.rm = TRUE),
                                                  sd = ~ sd(.x, na.rm = TRUE),
                                                  HDInterval::hdi,
                                                  `p<0.5` = ~ mean(.x < 0.5, na.rm = TRUE))
                              )
            ) %>% 
            suppressMessages() %>%
            suppressWarnings()

        ## Add the Below values to the Summary
        mods <- mods %>%
            left_join(mods.n) %>% 
            mutate(Summary = map2(.x = Summary, .y = Below,
                                  .f = ~ .x %>% left_join(.y))) %>% 
            suppressMessages() %>%
            suppressWarnings()
        
        save(mods,
              file = paste0(DATA_PATH, 'modelled/', Indicator, '__scores_', level,'_year.RData'))
        
}


my_rescale <- function(x, to = c(0,1), from = range(x, na.rm = TRUE, finite = TRUE)) {
    (x - from[[1]])/(from[[2]] - from[[1]]) * diff(to) + to[1]
}
