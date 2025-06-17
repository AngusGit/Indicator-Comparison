
library(tidybayes)
CI_29_JU_consequence_models <- function() {
    
 
        source('functions.R')
        source('CI_29_JU_consequence_models_functions.R')
        load(file = paste0(DATA_PATH, "parameters/IPM_juv.RData"))

        ## Critical value from IPM
		
        Crit.value <- IPM_juv$mean[1] %>%  # just inshore at the moment
            round(1)                       # round to align with Manu's documentation
        
       CI__get_JUV_MA_site_data()
       load(file = paste0(DATA_PATH, 'processed/juv.ma.site.RData'))
       
       ## exclude 2022-24 data as not part of the original baseline
       
        ## Model to predict the density of Juvenile Acropora ~ MApLag
        ## The desire is to produce a a single value (distribution
        ## perhaps) of MApLag, that can be used to set a value of MAp
        ## that corresponds to the critical density of Acropora
        ## juveniles estimated by Manu for the Juvenile critical
        ## indicator. as Crit.value above with is 0.5 for inshore reefs
       

        ## Model Formula.
        form <- formula(Acropora ~ log(MApLag + 0.01) +
                            f(reef.depth.site, model = 'iid'))
        
        
        a <- juv.ma.site %>%
          filter(REPORT_YEAR<2022) %>%
            group_by(Shelf) %>%
            summarise(data = list(data.frame(reef.depth.site, 
                                             MApLag, Acropora, avail.area)),
                      .groups = "drop") %>%
					  # mutate(CV=case_when(Shelf=="Inshore" ~ Critical_Value[1])##MANU
            mutate(newdata = map(.x = data,
                                 .f = ~ .x %>%
                                     tidyr::expand(reef.depth.site = NA,
                                                   MApLag = seq(min(.$MApLag, na.rm = TRUE),
                                                                max(.$MApLag, na.rm = TRUE),
                                                                len = 100),
                                                   Acropora = NA, avail.area = 1) %>%
                                     distinct()
                                 )) %>%
            ## Modelling data
            mutate(fulldata = map2(.x = data, .y = newdata,
                                   .f = ~ .x %>% bind_rows(.y))) %>%
            ## indices
            mutate(Indices = map2(.x = data, .y = newdata,
                                  .f = ~ (1:nrow(.y)) + nrow(.x))) %>% 
            ## Fit model
            mutate(Mod = map(.x = fulldata,
                             .f = ~ {
                                 set.seed(123)
                               environment(form)<-environment() 
                                 inla(form,
                                      offset = log(.x$avail.area),
                                      data = .x, 
                                      family = 'nbinomial',
                                      verbose = F,
                                      control.family = list(link = 'log'),
                                      control.predictor = list(link = 1, compute = TRUE),
                                      control.compute = list(
                                          dic = TRUE, cpo = TRUE, waic = TRUE,
                                          config = TRUE)
                                      )
                             })) %>%
            ## Posterior Draws
            mutate(Draws = map(.x = Mod,
                               .f = ~ inla.posterior.sample(n = 1000, .x, seed = 123))) %>%
            ## Calculate MA value associated with Critical Acropora
            ## This involves inverting the regression equation so as
            ## to predict x from model
            ## The inverse is (y - beta0)/beta1 on the link scale
            mutate(MA = map(.x = Draws,
			# .y= CV,#MANU
                            .f = ~ {
                                wch <- str_which(.x[[1]][[2]] %>% dimnames() %>% `[[`(1),
                                                  "Intercept|MApLag")
                                betas <- sapply(.x, function(x) x[[2]][wch]) 
                                apply(betas, 2, function(x)
                                (exp((log(Crit.value) - x[1])/x[2])) + 0.01)## MANU
                            }),
                   MA.sum = map(.x = MA,
                                .f = ~ .x %>% median_hdci())
                   ) %>%
            ## The rest of this is not required, but it does provide a visual
            ## confirmation of the above calculation
            ## Cellmeans
            mutate(Cellmeans = map2(.x = Draws, .y = Indices,
                                   .f = ~ sapply(.x, function(x) x[[2]][.y]))) %>%
            ## MA dist
            mutate(MA.dist = map2(.x = newdata, .y = Cellmeans,
                                 .f = ~ .x %>%
                                     dplyr::select(MApLag) %>%
                                     cbind(exp(.y)) %>%
                                     pivot_longer(cols = matches('[0-9]'), names_to = 'draw'))) %>%
            mutate(model.pred = map(.x = MA.dist,
                                    .f = ~ .x %>%
                                        group_by(MApLag) %>%
                                        median_hdci(value) %>%
                                        rename(lower=`.lower`, upper=`.upper`, median=value) )) %>%
            ## Plots 
            mutate(gg = map2(.x = model.pred, .y = MA.sum,
                            .f = ~ .x %>%
                                ggplot(aes(y = median, x = MApLag)) +
                                geom_line()+
                                geom_ribbon(aes(ymin = lower, ymax = upper),
                                            show.legend = NA, color = 'blue',alpha = 0.4)+
                                scale_x_continuous('MA proportion')+
                                scale_y_continuous('Acropora density')+
                                geom_hline(yintercept = Crit.value,
                                           linetype = 2, color = 'red',
                                           linewidth = 0.5) +
                                {if (.y$ymax < 1)
                                    geom_vline(xintercept = .y$ymax,
                                               linetype = 2, color = 'orange',
                                               linewidth = 0.5)} +
                                         
                                        #scale_y_log10('Acropora density')+
                                theme_classic(base_size = 12)))
        threshold <- a[1,'MA.sum'][[1]][[1]]['ymax'] 
        map2(.x = a$Shelf,
             .y = a$gg,
             .f = ~ ggsave(filename = paste0(OUTPUT_PATH, "figures/JU_consequence_", .x, ".png"),
                           .y + geom_vline(xintercept = threshold[[1]],
                                           linetype = 2, colour = 'orange',
                                           linewidth = 0.5),
                           width = 6, height = 6)
             )
        map2(.x = a$Shelf,
             .y = a$gg,
             .f = ~ ggsave(filename = paste0(OUTPUT_PATH, "figures/JU_consequence_", .x, "_large.png"),
                           .y + geom_vline(xintercept = threshold[[1]],
                                           linetype = 2, colour = 'orange',
                                           linewidth = 0.5),
                           width = 6, height = 6, dpi = 300)
             )

        ## Although the above models were fit for Inshore and
        ## Offshore, we only want to use Inshore
        ma_from_juv_consequence <- a %>% filter(Shelf == 'Inshore') %>%
            dplyr::select(Shelf, MA.sum) %>%
            unnest(MA.sum)
        
        save(ma_from_juv_consequence,
              file = paste0(DATA_PATH, 'modelled/ma_from_juv_consequence.RData'))

       
}
