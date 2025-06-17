## --------------------------- #
##
## Script name: JI_DefineBaselineThreshold.R
##
## Purpose of script: Explore historical records of hard coral
## juvenile density to produce a distribution of values over a defined
## time period and per bio region. These values are used as threshold
## regions for evaluating the juvenile density values per year (see
## indicator score)
##
## Author: Manuel Gonzalez-Rivero and Murray Logan
##
## Date Created: 17/09/2022
##
##
## --------------------------- #
##
## Notes :
##  *) Using all data up to and including REPORT_YEAR 2008
##  *) Bioregions with small sample sizes are aggregated based on similar
##     reef communities
##  *) The baseline model structure is set up to treat Site as the unit
##     of sampling. That is, points are summed across transects
##
## --------------------------- #

# NOTE in the final methods report the baseline period is stated as 2005-2021
source('CI_30_models_functions.R')
source('../R/functions.R')
source('CI_26_CC_baseline_models_functions.R') # includes the function CI_clean_inla_formula.R
CI_27_JU_baseline_models <- function() {
    
        load(paste0(DATA_PATH, 'processed/juv.df.RData')) # pre processing limits these data to post report_year 2006
         juv.df<-juv.df |>  filter(REPORT_YEAR<2022)  # AT added as states in the methods document that obs up to 2021 used in the baseline
         
        gbrmpa <- get(load(paste0(DATA_PATH, 'primary/gbrmpa.RData')))
        ## load(paste0(DATA_PATH, 'processed/points.analysis.data.transect.RData'))
        ## load(paste0(DATA_PATH, 'processed/spatial_lookup.RData'))
        ## load(paste0(DATA_PATH, 'primary/disturbances.RData'))

        bndry <- gbrmpa %>%
            st_transform(4326) %>%
            st_cast('POINT') %>%
            mutate(Longitude = st_coordinates(.)[,1],
                   Latitude = st_coordinates(.)[,2]) %>%
            select(Longitude, Latitude) %>%
            st_drop_geometry() %>%
            as.matrix() %>%
            inla.nonconvex.hull()


        fitModel_JU <- function(form, stack.est, spde, avail.area,family = 'poisson') {
          environment(form)<-environment()
           avail.area <- inla.stack.data(stack.est)$avail.area 
            inla(form,
                 ## offset = log(avail.area),
                 E = avail.area,
                 data = inla.stack.data(stack.est),
                 family= family,
                 control.predictor = list(compute = TRUE,
                                          link = 1,
                                          A = inla.stack.A(stack.est)
                                          ),
                 control.compute = list(config = TRUE),
                 verbose = FALSE)
        }


        juv.Baseline <- juv.df %>%
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

        ## Define the general formula

        juv.baseline <- juv.Baseline %>%
            mutate(
                ## STEP 1 - establish the full spatio-temporal grid #### 
                ## eventually replace with an expression that will provide unique
                ## long/lat/year combinations for all possible locations
                ## for now, I will use the locations defined in the data
                Full.grid = map(.x = data,
                                .f = ~ .x %>%
                                    dplyr::select(LONGITUDE, LATITUDE) %>%
                                    distinct() %>%
                                    tidyr::crossing(fYEAR = .x$fYEAR)
                                ),
                ## STEP 2 - establish a spatio-temporal grid of the actual data.class ####
                ## This is all observed locations for all monitored years
                Data.grid = map(.x = data,
                                .f = ~ .x %>%
                                    dplyr::select(LONGITUDE, LATITUDE) %>%
                                    distinct() %>%
                                    tidyr::crossing(fYEAR = .x$fYEAR)
                                ),
                ## STEP 3 - create a matrix of locations ####
                Coords = map(.x = Data.grid,
                             .f = ~ .x %>%
                                 dplyr::select(LONGITUDE, LATITUDE) %>%
                                 distinct() %>%
                                 as.matrix()
                             ),
                ## STEP 4 - create mesh #### 
                mesh = map(.x = Coords,
                           .f = ~ inla.mesh.2d(
                                   loc = .x,
                                   boundary = bndry,
                                   max.edge = c(0.5,3)*0.95,
                                   offset = c(0.95,3),
                                   cutoff = 0.95/20) ##40 secs with cutoff = 0.95/5
                           ),
                ## STEP 5 - generate the SPDE ####
                spde = map(.x = mesh, .f = ~ inla.spde2.matern(.x, alpha = 2)),
                ## define the spatial indicies for convenient access to model output 
                i.spatial = map(.x = spde,
                                .f = ~ inla.spde.make.index('spatial.field',
                                                            n.spde = .x$n.spde)),
                A.est = map2(.x = data, .y = mesh,
                             .f = ~ inla.spde.make.A(mesh = .y,
                                                     loc = as.matrix(.x[,c('LONGITUDE', 'LATITUDE')])
                                                     )
                             ),
                ## STEP 6 - generate the data stack ####
                stack.est = pmap(.l = list(data, A.est, i.spatial),
                                 .f = ~ inla.stack(data = list(y = ..1$value,
                                                               avail.area=..1$avail.area),
                                                   A=list(..2, 1,1,1,1),
                                                   effects = list(
                                                       ..3,
                                                       Intercept=rep(1, nrow(..1)),
                                                       list(fYEAR = ..1$fYEAR),
                                                       list(P_CODE = ..1$P_CODE),
                                                       list(Site = ..1$Site)
                                                   ),
                                                   tag = 'est'
                                                   )
                                 ),
                ## STEP 7 - define formula
                form = map(.x = data,
                           .f = ~ {
                               y ~ 0 +
                                   f(spatial.field, model = spde) +
                                   f(P_CODE, model = 'iid') +
                                   f(Site, model = 'iid') +
                                   f(fYEAR, model = 'iid') 
                           }
                           ),
                ## STEP 8 - remove formula terms for terms that have fewer than 2
                ## levels
                form = pmap(.l = list(form, data),
                            .f = ~ CI__clean_inla_formula(form = ..1, data = ..2) 
                            ),
                ## STEP 9 - define the family
                ## MIGHT WANT TO map this to a function that allows for a different
                ## family according to the response and/or depth
                family = 'poisson',
                ## STEP 10 - fit the model
                mod = pmap(.l = list(form, stack.est, spde, family),
                           .f = ~ fitModel_JU(..1, ..2, spde = ..3, family=..4) 
                           )
            )

  
           save(juv.baseline,
            # file=paste0(DATA_PATH, 'parameters/JU__baseline_mod.RData'))
           file=paste0(DATA_PATH, 'parameters/JU__baseline_mod_AT.RData')) # run with data through to 2021

################################################################################
####### PART 3 - get predictions for polygon samples ###########################
################################################################################
        load(file = paste0(DATA_PATH, "parameters/newdata_grid.RData"))  # created in CI_26_CC_baseline_models.R
        #load(file=paste0(DATA_PATH, 'parameters/JU__baseline_mod.RData'))
        load(file=paste0(DATA_PATH, 'parameters/JU__baseline_mod_AT.RData'))

        newdata <- newdata_grid %>%
            st_drop_geometry() %>%
            dplyr::select(BIOREGION, grid_ID, Longitude, Latitude, Weight) %>%
          #test influence of weights
          mutate(Weight=1) |>   
            mutate(grid_ID=as.factor(grid_ID),
                   BIOREGION=as.character(BIOREGION),
                   BIOREGION.agg=as.factor(case_when(BIOREGION %in% c("4", "3")~"4:3",
                                                     BIOREGION %in% c("35", "36")~"35:36",
                                                     !BIOREGION %in% c("4", "3", "35", "36")~BIOREGION))) %>%
            dplyr::rename(LONGITUDE = Longitude,
                          LATITUDE = Latitude)

        ## Predict to the level
        mods <- juv.baseline %>%
            mutate(Preds = pmap(.l = list(mod, mesh),
                                .f = ~ CI__make_predictions(newdata, mod = ..1, mesh = ..2)))

        ## Aggregate to BIOREGION level
        mods <- mods %>%
            mutate(bio.draws = map(.x = Preds,
                                   .f = ~ newdata %>%
                                       dplyr::select(grid_ID) %>%
                                       cbind(as.data.frame(t(.x))) %>%
                                       pivot_longer(cols=matches('[0-9]'), names_to='Rep') %>%
                                       right_join(newdata, by='grid_ID', relationship="many-to-many") 
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

        JUV_baseline <- mods %>%
            dplyr::select(bio.draws.mesh) %>%
            unnest(bio.draws.mesh)

      #  JUV_baseline |> group_by(BIOREGION.agg, .add=TRUE) |> tidybayes::summarise_draws(median, HDInterval::hdi)
        
        save(JUV_baseline, file = paste0(DATA_PATH, "parameters/JUV_baseline.RData"))

       
}
