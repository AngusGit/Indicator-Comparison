CI_26_CC_baseline_models <- function() {
    
    CI__add_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                   item = 'cc_baseline',
                   label = "CC baseline model", status = 'pending')

    CI_tryCatch({
        # source('../R/functions.R')
        # source('../R/CI_26_CC_baseline_models_functions.R')
        # source('../R/CI_30_models_functions.R')

      source('functions.R')
      source('CI_26_CC_baseline_models_functions.R')
      source('CI_30_models_functions.R')
      
        load(paste0(DATA_PATH, 'processed/points.analysis.data.transect.RData'))
        load(paste0(DATA_PATH, 'processed/spatial_lookup.RData'))
        load(paste0(DATA_PATH, 'primary/disturbances.RData'))

################################################################################
############ PART 1 - preparation of hex grid              #####################
################################################################################
        points <- points.analysis.data.transect %>%
            ungroup() %>%
            left_join(spatial_lookup |> dplyr::select(-TUMRA) |> distinct()) %>%  # multiple TUMRA at some reefs e.g Farquharson Reef
            mutate(P_CODE = as.factor(P_CODE),
                   NRM = as.factor(NRM),
                   DEPTH.f = factor(case_when(DEPTH >3 ~ "deep",
                                              DEPTH <= 3 ~ "shallow")),
                   REEF.d = factor(paste(REEF, DEPTH.f)),
                   REPORT_YEAR = as.numeric(as.character(REPORT_YEAR)),
                   BIOREGION = as.character(BIOREGION),
                   BIOREGION.agg = as.factor(case_when(BIOREGION %in% c("4", "3") ~"4:3",
                                                       BIOREGION %in% c("35", "36") ~"35:36",
                                                       !BIOREGION %in% c("4", "3", "35", "36") ~BIOREGION))) %>%
            ungroup() %>%
            filter(DEPTH < 9.1) %>%
            droplevels() %>%
            mutate(fYEAR = factor(REPORT_YEAR),
                   Site = factor(paste0(REEF.d, SITE_NO)),
                   Transect = factor(paste0(Site, TRANSECT_NO))) %>% 
            dplyr::select(-any_of(c('Hs70', '-Hs95', '-k490', '-chl', '-chl.wet',
                                    '-total.algae', '-MA', '-MA_BROWN', '-MA_RED',
                                    '-Lobophora', '-Sargassum', '-reef.site')))

        ## Construct the polygon sampling design for bioregions
        getSource <- FALSE
        while (!getSource) {
            try({
                gbr <- st_read(paste0("https://services8.arcgis.com/",
                                      "ll1QQ2mI4WMXIXdm/arcgis/rest/services/",
                                      "Great_Barrier_Reef_Marine_Park_Boundary/",
                                      "FeatureServer/0/query?outFields=*&where=1%3D1&f=geojson"),
                               quiet = T) %>%
                    select(geometry) %>% 
                    st_make_valid() %>%
                    st_transform(3857)
                getSource <- TRUE
            },
            silent = TRUE)
        }

        ## load(paste0(DATA_PATH, "processed/bregions.sf.RData"))
        load(paste0(DATA_PATH, "primary/bioregions.RData"))
        ## bioregion <- bregions.sf %>%
        bioregion <- bioregions %>%
            st_make_valid() %>% 
            st_transform(3857)

        ## Grid up the gbr into 500x500m hexagons
        ## This takes a long time
        hex_grid <- gbr %>%
            st_make_grid(n = c(500,500), square = FALSE) %>%
            st_as_sf() %>%
            mutate(grid_ID = 1:n()) %>%
            st_intersection(gbr)

        ## Explore this for Bioregion 27
        bb <- bioregion %>%
            filter(BIOREGION == 27) %>%
            st_bbox()

        g <- ggplot() +
            geom_sf(data = hex_grid %>% st_crop(bb), fill = NA) +
            geom_sf(data = bioregion %>% st_crop(bb), color = 'blue') +
            coord_sf(xlim = bb[c(1,3)],
                     ylim = bb[c(2,4)]) +
            theme_bw() +
            coord_sf(expand = 0)
        
        ggsave(filename = paste0(FIGS_PATH, "/hex_grid__bioregion_27.png"),
               g,
               height = 6,
               width = 6)


        ## Further filter the hex grid to just those hexagons that sit over
        ## reef
        hex_grid.over_reef <- hex_grid %>%
            filter(lengths(st_intersects(., bioregion)) >0 )

        bb <- bioregion %>%
            filter(BIOREGION == 27) %>%
            st_bbox()
        g <- ggplot() +
            geom_sf(data = bioregion %>% st_crop(bb), color = 'blue') +
            geom_sf(data = hex_grid.over_reef %>% st_crop(bb), fill = NA, fill = 'red') +
            coord_sf(xlim = bb[c(1,3)],
                     ylim = bb[c(2,4)]) +
            theme_bw() +
            coord_sf(expand = 0)
        
        ggsave(filename = paste0(FIGS_PATH, "/hex_grid__bioregion_27_reef.png"),
               g,
               height = 6,
               width = 6)

        ## Simplify the bioregion object such that it just contains a
        ## multipolygon per bioregion This will prevent duplications when
        ## intersecting with the hex grid
        bioregion_sum <- bioregion %>%
            group_by(BIOREGION) %>%
            summarise()

        ## Cut the non-reef parts of the hexagons out
        hex_grid.just_reef <-
            hex_grid.over_reef %>%
            st_intersection(bioregion_sum)

        bb <- bioregion %>%
            filter(BIOREGION == 27) %>%
            st_bbox()
        g <- ggplot() +
            geom_sf(data = bioregion %>% st_crop(bb), color = 'blue') +
            geom_sf(data = hex_grid.just_reef %>% st_crop(bb), fill = NA, fill = 'red') +
            coord_sf(xlim = bb[c(1,3)],
                     ylim = bb[c(2,4)]) +
            theme_bw() +
            coord_sf(expand = 0)
        ggsave(filename = paste0(FIGS_PATH, "/hex_grid__bioregion_27_justreef.png"),
               g,
               height = 6,
               width = 6)

        ## Calculate the reef area and weight within each hexagon
        newdata_grid <- hex_grid.just_reef %>%
            mutate(Reef.area = st_area(.)) %>%
            st_centroid() %>%
            st_transform(4326) %>%
            mutate(Longitude = st_coordinates(.)[,1],
                   Latitude = st_coordinates(.)[,2]) %>%
            group_by(BIOREGION) %>%
            mutate(Weight = as.numeric(Reef.area)/sum(as.numeric(Reef.area))) %>%
            ungroup()
        save(newdata_grid,
             file = paste0(DATA_PATH, "parameters/newdata_grid.RData"))

        ## Visual check
        bb <- hex_grid.just_reef %>%
            filter(BIOREGION == 27) %>%
            slice(1:10) %>%
            st_bbox()

        g <- ggplot() +
            geom_sf(data = hex_grid.just_reef %>% st_crop(bb), fill='orange') +
            geom_sf_text(data = hex_grid.just_reef %>% st_crop(bb), aes(label = grid_ID)) +
            geom_sf(data = bioregion %>% st_crop(bb), fill = NA, colour = 'blue') +
            theme_bw() +
            coord_sf(expand = 0)

        ggsave(filename = paste0(FIGS_PATH, "/new_grid.png"),
               g,
               height = 6,
               width = 6)

################################################################################
############ PART 2 - fit the spatial models               #####################
################################################################################

                                        #The baseline model structure is set up to treat Site as the unit of
                                        #sampling. That is, points are summed across transects
        points.site <- points %>% 
            group_by(P_CODE, REEF, DEPTH, VISIT_NO, SITE_NO, REPORT_YEAR, LATITUDE, 
                     LONGITUDE, Project, BIOREGION, NRM, DEPTH.f, REEF.d, BIOREGION.agg, fYEAR, 
                     Site) %>% 
            summarise(HC = sum(HC),
                      total.points = sum(total.points))

        points.full <- points.site %>%
            filter(REPORT_YEAR < 2011) %>%
            droplevels() %>%
            group_by(REEF.d) %>%
            mutate(LONGITUDE = mean(LONGITUDE),
                   LATITUDE = mean(LATITUDE))


        bndry <- gbr %>%
            st_transform(4326) %>%
            st_cast('POINT') %>%
            mutate(Longitude = st_coordinates(.)[,1],
                   Latitude = st_coordinates(.)[,2]) %>%
            select(Longitude, Latitude) %>%
            st_drop_geometry() %>%
            as.matrix() %>%
            inla.nonconvex.hull()

        ## Model separately per depth
        mods <- points.site %>%
            group_by(DEPTH.f) %>%
            nest()

        mods <- mods %>%
            ## STEP 1 - establish the full spatio-temporal grid eventually
            ## replace with an expression that will provide unique
            ## long/lat/year combinations for all possible locations for now,
            ## I will use the locations defined in the data
            mutate(points.full.grid = map(.x = data,
                                          .f = ~ .x %>% 
                                              dplyr::select(LONGITUDE, LATITUDE) %>%
                                              distinct() %>%
                                              tidyr::crossing(fYEAR = .x$fYEAR)
                                          )
                   ) %>%
            ## STEP 2 - establish a spatio-temporal grid of the actual
            ## data.class This is all observed locations for all monitored
            ## years
            mutate(points.data.grid = map(.x = data,
                                          .f = ~ .x %>%
                                              dplyr::select(LONGITUDE, LATITUDE) %>%
                                              distinct() %>%
                                              tidyr::crossing(fYEAR = .x$fYEAR)
                                          )
                   ) %>%
            ## STEP 3 - create a matrix of locations
            mutate(point.coords = map(.x = points.data.grid,
                                      .f = ~ .x %>%
                                          dplyr::select(LONGITUDE, LATITUDE) %>%
                                          distinct() %>%
                                          as.matrix()
                                      )
                   ) %>%
            ## STEP 4 - create mesh
            mutate(mesh = map(.x = point.coords,
                              .f = ~ inla.mesh.2d(
                                      loc = .x,
                                      boundary = bndry,
                                      max.edge = c(0.5,3)*0.95,
                                      offset = c(0.95,3),
                                      cutoff = 0.95/20) ##40 secs with cutoff = 0.95/5
                              )
                   ) %>%
            ## STEP 5 - generate the SPDE
            mutate(spde = map(.x = mesh,
                              .f = ~ inla.spde2.matern(.x, alpha = 2)
                              )
                   ) %>% 
            ## define the spatial indicies for convenient access to model
            ## outputs
            mutate(i.spatial = map(.x = spde,
                                   .f = ~ inla.spde.make.index('spatial.field',
                                                               n.spde = .x$n.spde)
                                   )
                   ) %>%
            mutate(A.est = map2(.x = data, .y = mesh,
                                .f = ~ inla.spde.make.A(mesh = .y,
                                                        loc = as.matrix(.x[,c('LONGITUDE', 'LATITUDE')])
                                                        )
                                )
                   ) %>%
            ## STEP 6 - generate the data stack
            mutate(stack.est = pmap(.l = list(data, A.est, i.spatial),
                                    .f = ~inla.stack(data = list(y = ..1$HC,
                                                                 Total = ..1$total.points),
                                                     A = list(..2, 1, 1, 1),
                                                     effects = list(
                                                         c(list(Intercept = 1), ..3),
                                                         list(fYEAR = ..1$fYEAR),
                                                         list(Project = ..1$Project),
                                                         list(Site = ..1$Site)
                                                     ),
                                                     tag = 'est'
                                                     )
                                    )
                   ) %>%
            ## STEP 7 - define formula
            mutate(form = map(.x = data,
                              .f = ~ {
                                  y ~ 0 + #Intercept +
                                      f(spatial.field,
                                        model = spde
                                        ## group = spatial.field.group,
                                        ## control.group = list(model = "ar1")) +
                                        ) +
                                      f(Project, model = 'iid') +
                                      f(Site, model = 'iid') +
                                      f(fYEAR, model = 'iid')
                              }
                              )
                   ) %>%
            ## STEP 8 - remove formula terms for terms that have fewer than 2
            ## levels
            mutate(form = pmap(.l = list(form, data),
                               .f = ~ CI__clean_inla_formula(form = ..1, data = ..2) 
                               )
                   ) %>%
            ## STEP 8 - fit the model
            mutate(mod = pmap(.l = list(form, stack.est, spde),
                              .f = ~ fitModel(..1, ..2, spde = ..3) 
                              )
                   )

        save(mods,
             file=paste0(DATA_PATH, 'parameters/CC__baseline.RData'))

################################################################################
####### PART 3 - get predictions for polygon samples ###########################
################################################################################

        newdata<- newdata_grid %>%
            st_drop_geometry() %>%
            dplyr::select(BIOREGION, grid_ID, Longitude, Latitude, Weight) %>%
            mutate(grid_ID=as.factor(grid_ID),
                   BIOREGION=as.character(BIOREGION),
                   BIOREGION.agg=as.factor(case_when(BIOREGION %in% c("4", "3")~"4:3",
                                                     BIOREGION %in% c("35", "36")~"35:36",
                                                     !BIOREGION %in% c("4", "3", "35", "36")~BIOREGION))) %>%
            dplyr::rename(LONGITUDE = Longitude,
                          LATITUDE = Latitude)
        
        ## Predict to the level
        mods <- mods %>%
            mutate(Preds = pmap(.l = list(mod, mesh),
                                .f = ~ CI__make_predictions(newdata, mod = ..1, mesh = ..2)))

        ## Aggregate to BIOREGION level
        mods <- mods %>%
            mutate(bio.draws = map(.x = Preds,
                                   .f = ~ newdata %>%
                                       dplyr::select(grid_ID) %>%
                                       cbind(as.data.frame(t(.x))) %>%
                                       pivot_longer(cols=matches('[0-9]'), names_to='Rep') %>%
                                       right_join(newdata) 
                                   )
                   )
        ## sum of weighted medians for every Rep
        mods <- mods %>%
            mutate(bio.draws.mesh = map(.x = bio.draws,
                                        .f = ~ .x %>% group_by(BIOREGION.agg, Rep) %>%
                                            summarise(value = sum(value*(Weight/sum(Weight)))) %>%
                                            mutate(.draw = as.integer(str_replace(Rep, 'V',''))) %>%
                                            dplyr::select(-Rep) %>% 
                                            arrange(BIOREGION.agg, .draw)
                                        )
                   )

        mod <- mods %>%
            dplyr::select(bio.draws.mesh) %>%
            unnest(bio.draws.mesh)

        save(mod, file = paste0(DATA_PATH, "parameters/CC_baseline.RData"))
        ## load(file = paste0(DATA_PATH, "parameters/CC_baseline.RData"))
        ## save(mod, file = paste0(DATA_PATH, "parameters/CC_baseline_posteriors.RData"))
        ## save(mod, file = paste0(DATA_PATH, "parameters/CC_baseline_posteriors.RData"))

        CI__change_status(stage = paste0('STAGE',CI$setting$CURRENT_STAGE),
                          item = 'cc_baseline',status = 'success')

    }, logFile=LOG_FILE, Category='--Data processing--',
    msg=paste0('CC baseline model.'), return=NULL)
}
