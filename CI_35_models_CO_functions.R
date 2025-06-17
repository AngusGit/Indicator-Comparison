CI_clear_models_CO_data <- function() {
    files <- list.files(path = paste0(DATA_PATH, "modelled"),
                        pattern = "CO.*|data_co.*",
                        full.names = TRUE)
    unlink(files)
}


CI_models_CO_prepare_data <- function() {
   
        load(file = paste0(DATA_PATH, 'processed/comp.transect.RData'))

        comp.transect <- comp.transect %>%
            mutate(P_CODE = ifelse(P_CODE == "RMRAP", "RM", P_CODE)) %>%
            mutate(DEPTH.f=factor(case_when(DEPTH>3~"deep slope",
                                            DEPTH<=3~"shallow slope")),
                   Site=factor(paste(REEF,DEPTH.f,SITE_NO)),
                   REEF.d=factor(paste(REEF, DEPTH.f)),
                   Transect=factor(paste(REEF,DEPTH.f,SITE_NO,TRANSECT_NO)),
                   fYEAR=factor(as.numeric(as.character(REPORT_YEAR)))) %>%
            suppressMessages() %>%
            suppressWarnings()
            
        ## Reef-level aggregation of cover data            
        coral.comms <- comp.transect %>%
            group_by(P_CODE, REEF.d, COMP_2021, VISIT_NO, REPORT_YEAR) %>%
            summarise(cover = mean(cover)) %>%
            ungroup() %>%
            group_by(REEF.d) %>%
            mutate(StartYR = min(REPORT_YEAR),
                   EndYR = max(REPORT_YEAR),
                   SurveyDuration = EndYR - StartYR) %>% 
            filter(EndYR >= 2016,
                   SurveyDuration >= 5) %>%
            droplevels() %>%
            pivot_wider(id_cols = everything(),
                        names_from = COMP_2021,
                        values_from = cover) %>%
            ## pivot longer according to all the taxa
            pivot_longer(-any_of(c("P_CODE", "REEF.d", "VISIT_NO", "REPORT_YEAR",
                                   "StartYR", "EndYR", "SurveyDuration")),
                         names_to = 'Taxa',
                         values_to = 'values') %>%
            ## replace missing values with 0's
            mutate(values = replace_na(values, 0)) %>%
            ## sqrt transform cover data
            mutate(values = sqrt(values)) %>% 
            mutate(REEF.d = factor(REEF.d)) %>%
            group_by(P_CODE, REEF.d, VISIT_NO, REPORT_YEAR) %>%
            ## relative abundance
            mutate(values = values/sum(values)) %>%
            ## pivot wider
            pivot_wider(id_cols = everything(),
                        names_from = Taxa,
                        values_from = values) %>%
            suppressMessages() %>%
            suppressWarnings()

        ## coral cover summary
        coral.cover <- comp.transect %>%
            group_by(REEF.d, GROUP_CODE, COMP_2021,
                     VISIT_NO, REPORT_YEAR) %>%
            summarise(cover = mean(cover)) %>%
            filter(GROUP_CODE == "HC") %>%     ## limit to just HC
            ungroup() %>%
            group_by(REEF.d, REPORT_YEAR) %>%
            summarise(cover = sum(cover)) %>%
            ungroup() %>%
            mutate(adj.cover = cover/100) %>%  ## adjusted cover for plotting later
            suppressMessages() %>%
            suppressWarnings()

        coral.comms <- coral.comms %>%
            left_join(coral.cover) %>%
            dplyr::relocate(cover, adj.cover, .after = SurveyDuration) %>% 
            suppressMessages() %>%
            suppressWarnings()
        
        save(coral.comms,
              file = paste0(DATA_PATH, 'modelled/coral.comms.RData'))
        
}

CI_models_CO_calc_mds <- function() {
    
          load(file = paste0(DATA_PATH, 'modelled/coral.comms.RData'))

        reef.nmds <- coral.comms %>%
            group_by(REEF.d) %>%
            summarise(data = list(cur_data_all()), .groups = "drop") %>%
            mutate(data = map(.x = data,
                              .f = ~ .x %>% arrange(REPORT_YEAR) %>%
                                  tibble::column_to_rownames(var = "REPORT_YEAR"))) %>%
            mutate(Composition_data = map(.x = data,
                                          .f = ~ .x %>% 
                                              dplyr::select(-any_of(c("P_CODE", "REEF.d",
                                                                      "VISIT_NO", "StartYR",
                                                                      "EndYR", "SurveyDuration",
                                                                      "cover", "adj.cover")))), 
                   NMDS = map(.x = Composition_data,
                              .f = ~ vegan::metaMDS(.x, k = 2, trace = FALSE) %>%
                                  suppressMessages() %>%
                                  suppressWarnings()
                              ),
                   Site_scores = map(.x = NMDS,
                                     .f = ~ .x %>% vegan::scores('sites') %>%
                                         as.data.frame() %>%
                                         tibble::rownames_to_column(var = 'Year')
                                     ),
                   Stress = map(.x = Composition_data,
                                .f = ~ {
                                    .x %>% vegan::metaMDS(k = 2, trace = FALSE, weakties = FALSE) %>%
                                        `[[`('stress') %>% 
                                        suppressMessages() %>%
                                        suppressWarnings()
                                })
                   ) %>%
            unnest(c(Site_scores, Stress)) %>%
            dplyr::select(-data, -Composition_data, -NMDS) %>% 
            group_by(REEF.d) %>%
            summarise(data = list(cur_data_all()), .groups = "drop") 

        save(reef.nmds,
              file = paste0(DATA_PATH, 'modelled/reef.nmds.RData'))
        
       
}



############# The following functions were supplied by the UQ team ##################################
## This version "timecutoff" is a two character vector that constrains
## the comparison matrix using year numbers.  testLOF now also reports
## a taxa table relative to the neighbors. This has to be run for each
## time point which complicates how the LOF scores are aggregated
## within aggLOF lapply k is now a *single* value only, rather than
## iterating over a bunch of k's.  ref.method to choose between
## year-based reference selection or row number based reference
## selection.  row based reference selection may be more generalisable
## given some reef surveys started later

LOF.Ref <- function(commData, k, timecutoff = NULL, method=method) {
  
  require(Rlof)
  
  aggLOF <- lapply(1:nrow(commData), function(n) {
    
      ## subset matrix. If timecutoff is NULL, get 1:n rows (no future
      ## rows) if is not NULL, contrain matrix to cutoff rows, add
      ## target row to end
      if (is.null(timecutoff)) {
          subMat <- commData[1:n,]
      } else if (!is.null(timecutoff) &
                 (!as.numeric(row.names(commData[n,])) %in%
                  c(timecutoff[1]:timecutoff[2]))) {
          ## start row is first row >= to start cutoff
          startRow <- which(as.numeric(row.names(commData)) >= timecutoff[1])[1]
          ## end row is the row that meets end cutoff
          endRow <- which(as.numeric(row.names(commData)) >= timecutoff[2])[1]
          ## but must be less than n
          endRow <- ifelse(endRow > n, n-1, endRow)
          
          subMat <- commData[unique(c(startRow:endRow, n)),]
      } else if (!is.null(timecutoff) &
                 (as.numeric(row.names(commData[n,])) %in%
                  c(timecutoff[1]:timecutoff[2]))) {
          subMat <- commData[n,]
      }
      
      ## run LOF function
      nullReturn <- list(data.frame(n=n, k=k,LOF=NA), NULL)
    
      if (nrow(subMat) < k+1) {
          LOFneigh <- nullReturn
          combLOF <- nullReturn
      } else if (nrow(subMat) >= k+1) {
          ## return the neighbor data for the supplementary table
          LOFneigh <- LOFveg(subMat, k=k, method=method, neighborData=TRUE)
          combLOF <- LOFveg(subMat, k=k, method=method, neighborData=FALSE)
      }

      if(sum(is.na(combLOF)) > (nrow(subMat)-2)) {
          LOFneigh <- nullReturn
          combLOF <- nullReturn
      }
    
      ## if we have no variation in LOF (all values = 0) we have a P of 1
      if(length(unique(combLOF))==1) {
          tempReturn <- nullReturn
          tempReturn[[1]]$LOF = 1
          return(tempReturn)}
    
      ## get the neighbor table ranking taxa by biggest change to
      ## smallest change k+1 row is the target
    
      if (nrow(subMat) >= k+1) {
          nData <- as.matrix(subMat[c(LOFneigh[1:k,ncol(LOFneigh)], nrow(subMat)),])
          
          ## sweep to subtract target row from neighbour rows
          nDiff <- - sweep(nData[-nrow(nData),], 2, nData[nrow(nData),], FUN = "-")
          nDiffDf <- t(apply(nDiff, 2, function(x) {
              cbind(mean = mean(x),
                    var = var(x)) }))
          colnames(nDiffDf) = c("meanDiff", "varDiff")
          nDiffDf <- as.data.frame(nDiffDf[order(abs(nDiffDf[,1]), decreasing=TRUE),])
                                        #nDiffDf <- nDiffDf[nDiffDf[,1] != 0,]
      } else if (nrow(subMat) < k+1) {
          nDiffDf <- data.frame(meanDiff=rep(NA, times=ncol(commData)),
                                varDiff=rep(NA, times=ncol(commData)))
          combLOF <- combLOF[[1]][,3L, drop=F]
      }
      
      return(list(data.frame(n=n, k=k, LOF=combLOF[nrow(combLOF),]),
                  nDiffDf))
  })
  
  ## return a two-entry list. 
  ## 1 is a dataframe of LOF scores
  ## 2 is a list of taxa neighbor comparison tables
  
    LOFscoreDf <- do.call('rbind', lapply(1:length(aggLOF), function(n) {
        aggLOF[[n]][[1]] }))
    LOFtables <- lapply(1:length(aggLOF), function(n) {
        aggLOF[[n]][[2]] })
    
    return(list(LOFscoreDf, LOFtables))
}

## LOF summary function adjust sigmoid.a value to change the
## sensitivity of LOF value threshold higher sigmoid.a value increases
## the sensitivity
minLOF <- function(lofList, sigmoid.a=3) {
  tData <- lofList[[1]]
  
  ## sigmoid scaling for LOF values
  inv.sig <- function(x, a, b) {
      (1/(1+exp(a*x-b))) / (1/(1+exp(a-b)))
  }
  
  tMin <- inv.sig(tData$LOF, sigmoid.a, 1)
  tMin[tMin>1] = 1

  return(tMin)
}


# LOF computation function
LOFveg <- function (data, k, method, neighborData=FALSE, ...) {
    if (is.null(k)) 
        stop("k is missing")
    if (!is.numeric(k)) 
        stop("k is not numeric")
    data <- as.matrix(data)
    if (!is.numeric(data)) 
        stop("the data contains non-numeric data type")
    v.k <- as.integer(k)
    if (max(v.k) >= dim(data)[1])
        stop("the maximum k value has to be less than the length of the data")
  

    ## TIM CUSTOMISATION ##

    ## substitute f.dist.to.knn() with a custom function that runs the
    ## same code but accepts vegdist instead of distmc (multicore
    ## version of dist) this allows for ecological dissimilarity
    ## indices
    distdata <- f.dist.to.knn.VEG(data, max(v.k), method=method, ...)

    if(neighborData) {
        return(distdata)
    }
  
    ## END CUSTOMISATION ##
  
    ## get dimensions and where distances start (at halfway)
    p <- dim(distdata)[2L]
    dist.start <- as.integer((dim(distdata)[1])/2)
    dist.end <- dim(distdata)[1]
  
    ## looped to allow for multiple k sizes
    m.lof <- sapply(v.k, function(ik) {
        lrddata <- f.reachability.VEG(distdata, ik)
        v.lof <- rep(0, p)
        for (i in 1:p) {
            nneigh <- sum(!is.na(distdata[c((dist.start + 1):dist.end), 
                                          i]) & (distdata[c((dist.start + 1):dist.end), 
                                                          i] <= distdata[(dist.start + ik), i]))
            v.lof[i] <- sum(lrddata[distdata[(1:nneigh), i]]/lrddata[i])/nneigh
        }
        ## any vlof of NaN is actually 0 (perfect reachability)
        v.lof[is.na(v.lof)] = 0
        v.lof
    })
  
    if (length(v.k) > 1) 
        colnames(m.lof) <- v.k
    return(m.lof)
}

f.dist.to.knn.VEG <- function(dataset, neighbors, method, ...) {
    require(vegan)
    
    m.dist <- as.matrix(vegdist(dataset, method=method, ...))
    
    num.col <- dim(m.dist)[2]
    
    ## get identities and distances to kth neighbours (allowing for
    ## ties)
    l.knndist <- lapply(c(1:num.col), function(i) {
        ## ordered vector of neighbor distances
        order.x <- sort(m.dist[-i,i])
        
        ## original code
    
        ## order points by distance to i
        order.x <- order(m.dist[, i])
    
        ## get distance of neighbours + 1 (accounting for 0 of self
        ## distance)
        kdist <- m.dist[, i][order.x[neighbors + 1]]
    
        ## how many distances are below the kth target-neighbour
        ## distance?  this accounts for other neighbor distances that
        ## are identical to the kth neighbour distance (so ties)
        numnei <- sum(m.dist[, i] <= kdist)
    
        ## return neighbor identities & distances
        data.frame(v.order = order.x[2:numnei], 
                   v.dist = m.dist[,i][order.x[2:numnei]])
    })
    rm(m.dist)
  
    ## this gets the max number of neighbors from l.knndist
    maxnum <- max(unlist(lapply(l.knndist, function(x) {
        dim(x)[1]
    })))
  
    ## this pads out the distance matrices with NAs to allow for some
    ## targets to have different numbers of kth comparisons.  returns
    ## a matrix with columns = results, rows = neighbor identities
    ## followed by distances. With no ties, rows should equal 2xk
    i <- numeric()
    knndist <- do.call("cbind", lapply(1:num.col, function(i) {
        len <- dim(l.knndist[[i]])[1]
        c(l.knndist[[i]]$v.order, rep(NA, (maxnum - len)), 
          l.knndist[[i]]$v.dist, rep(NA, (maxnum - len)))
    }))
    colnames(knndist) = paste0("result.", 1:ncol(knndist))
    
    return(knndist)
}

f.reachability.VEG <- function (distdata, k) {
    p <- dim(distdata)[2]
    lrd <- rep(0, p)
    dist.start <- as.integer((dim(distdata)[1])/2)
    dist.end <- dim(distdata)[1]
    
    for (i in 1:p) {
        ## crop out any NAs from both the identities and the distances
        numneigh <- sum(!is.na(distdata[c((dist.start + 1):dist.end), i]) & 
                        (distdata[c((dist.start + 1):dist.end), i] <= 
                         distdata[(dist.start + k), i]))
        ##numneigh = k if there aren't any NAs
    
        j <- c(1:numneigh)
        ## this puts the kth neighbour distance for each of target's
        ## neighbours beside the neigbhour distances for target
        temp <- rbind(distdata[dist.start + k, distdata[j, i]], 
                      distdata[dist.start + j, i])
        ## temp col 1 = distance to kth neighbour for each of target's
        ## closest k neighbours temp col 2 = target's distance to
        ## closest k neighbours
    
        ## reachability formula:  
        reach <- 1/(sum(apply(temp, 2, max))/numneigh)
    
        ## reachability fails when there are zero divisions, so
        ## replace Inf /0 failures with a zero (as a point completely
        ## overtop of neigbhours gets a "perfect" reachability of 0)
        reach <- ifelse(reach == Inf, 0, reach)
        
        lrd[i] <- reach
    }
    lrd
}

## LOF summary function
LOF.summary <- function(LOF.data, id.data, site.var, year.var) {
    out.df <- vector('list', length=length(LOF.data))
    
    for (i in 1:length(LOF.data)) {
        yr.placeholder <- id.data[[i]][,year.var]
        
        out.df[[i]] <- data.frame(Reef=unique(id.data[[i]][,site.var]),
                                  LOF=LOF.data[[i]],
                                  Year=yr.placeholder)
    }
    return(out.df)
}

## Murray's alternate version that works with nested data
LOF_summary <- function(data, LOF, k) {
    data %>% dplyr::select(REEF.d, REPORT_YEAR) %>%
        cbind(LOF = LOF, k = k)
}


## LOF table summary function
LOF.table.summary <- function(LOF.data, id.data, site.var, year.var) {
    tmp.tables <- vector('list', length=length(LOF.data))
    
    for (i in 1:length(LOF.data)) {
        tmp.tables[[i]] <- LOF.data[[i]][[2]]
        
        for (j in 1:length(tmp.tables[[i]])) {
            yr.placeholder <- as.numeric(levels(as.factor(id.data[[i]][,year.var]))[[j]])
            tmp.tables[[i]][[j]]$Reef <- unique(id.data[[i]][,site.var])
            tmp.tables[[i]][[j]]$Year <- yr.placeholder
            tmp.tables[[i]][[j]]$Taxon <- row.names(tmp.tables[[i]][[j]])
        }
        
        tmp.tables[[i]] <- as.data.frame(data.table::rbindlist(tmp.tables[[i]]))
        tmp.tables[[i]] <- tmp.tables[[i]][complete.cases(tmp.tables[[i]]),]
    }
    
    tmp.tables <- as.data.frame(data.table::rbindlist(tmp.tables))
    return(tmp.tables)
}

## Murray's alternate version that works with nested data
LOF_table_summary <- function(data, LOF, k) {
    dat <- LOF %>% `[[`(2)
    id <- data
    map_df(1:length(dat),
           .f = ~ dat[[.x]] %>%
               cbind(REPORT_YEAR = id$REPORT_YEAR[.x],
                     REEF.d = unique(id$REEF.d),
                     Taxon = rownames(dat[[.x]]),
                     k = k
                     ) %>%
               filter(!is.na(meanDiff),
                      !is.na(varDiff))
           )
}

############# End of the functions supplied by the UQ team ##########################################

CI_models_CO_calc_novelty <- function() {
    

        load(file = paste0(DATA_PATH, 'modelled/coral.comms.RData'))
        load(file = paste0(DATA_PATH, 'modelled/reef.nmds.RData'))
        ## Metrics computation ##
        ## Local outlier factor (LOF) metric ##
        ## ### ### ### ### ### ### ### ### ####  
        ## IMPORTANT: READ ANNOTATION IN THE NEXT LINE ###
        ## ### ### ### ### ### ### ### ### ### ### ### ### 

        ## The LOF function (in this helper script) will output NA
        ## deviation indicator values for the reference period
        ## (timecutoff argument). To see deviation indicator values
        ## for the reference period, mask this line and use
        ## gbrfLOF241122.R function in the next line instead.

        ## Novelty indicator (low LOF k value with no reference period
        ## to detect acute changes in community composition)
        GBR.LOF.novelty <- coral.comms %>%
            ungroup() %>%
            group_by(REEF.d) %>%
            summarise(data = list(cur_data_all()), .groups = "drop") %>%
            mutate(Composition_data = map(.x = data,
                                          .f = ~ .x %>% 
                                              dplyr::select(-any_of(c("P_CODE", "REEF.d",
                                                                      "VISIT_NO", "REPORT_YEAR",
                                                                      "StartYR", "EndYR",
                                                                      "SurveyDuration",
                                                                      "cover", "adj.cover")))),
                   ## Novelty indicator (low LOF k value with no
                   ## reference period to detect acute changes in
                   ## community composition)
                   Novelty = map(.x = Composition_data,
                                 .f = ~ LOF.Ref(commData=.x,
                                                k=3, 
                                                method='bray')),
                   ## adjust sigmoid.a value to change the sensitivity
                   ## of LOF value threshold.
                   Novelty.df = map(.x = Novelty,
                                    .f = ~ minLOF(.x, sigmoid.a = 3)),
                   ## replacement for LOF.summary
                   Novelty.df.sum = map2(.x = data, .y = Novelty.df,
                                        .f = ~ .x %>% LOF_summary(LOF = .y, k = 3)),
                   ## replacement for LOF.table.summary
                   Novelty.table = map2(.x = data, .y = Novelty,
                                         .f = ~ .x %>% LOF_table_summary(LOF = .y, k = 3))
                   )

        save(GBR.LOF.novelty,
              file = paste0(DATA_PATH, 'modelled/GBR.LOF.novelty.RData'))

       
}

CI_models_CO_calc_deviation <- function() {
    

        load(file = paste0(DATA_PATH, 'modelled/coral.comms.RData'))
        load(file = paste0(DATA_PATH, 'modelled/reef.nmds.RData'))

        ## Deviation indicator (high LOF k value with a reference
        ## period (here - 1995 to 2016) to detect changes in community
        ## composition from a previous understanding of reef state)
        GBR.LOF.deviation <- coral.comms %>%
            ungroup() %>%
            group_by(REEF.d) %>%
            summarise(data = list(cur_data_all()), .groups = "drop") %>%
            mutate(Composition_data = map(.x = data,
                                          .f = ~ .x %>% 
                                              dplyr::select(-any_of(c("P_CODE", "REEF.d",
                                                                      "VISIT_NO",
                                                                      "StartYR", "EndYR",
                                                                      "SurveyDuration",
                                                                      "cover", "adj.cover"))) %>%
                                          tibble::column_to_rownames(var = "REPORT_YEAR")),
                   Deviation = map(.x = Composition_data,
                                   .f = ~ LOF.Ref(commData=.x,
                                                  k=6, 
                                                  timecutoff = c(1995, 2016),
                                                  method='bray')),
                   ## adjust sigmoid.a value to change the sensitivity
                   ## of LOF value threshold.
                   Deviation.df = map(.x = Deviation,
                                    .f = ~ minLOF(.x, sigmoid.a = 3)),
                   ## replacement for LOF.summary
                   Deviation.df.sum = map2(.x = data, .y = Deviation.df,
                                        .f = ~ .x %>% LOF_summary(LOF = .y, k = 6)),
                   ## replacement for LOF.table.summary
                   Deviation.table = map2(.x = data, .y = Deviation,
                                         .f = ~ .x %>% LOF_table_summary(LOF = .y, k = 6))
                   )
        
                   
        save(GBR.LOF.deviation,
              file = paste0(DATA_PATH, 'modelled/GBR.LOF.deviation.RData'))

}



CI_models_CO_combine_metrics <- function() {
    

        load(file = paste0(DATA_PATH, 'modelled/GBR.LOF.novelty.RData'))
        load(file = paste0(DATA_PATH, 'modelled/GBR.LOF.deviation.RData'))
        load(file = paste0(DATA_PATH, 'modelled/reef.nmds.RData'))
        ## load(file = paste0(DATA_PATH, 'primary/sample.reef.RData'))
        load(file = paste0(DATA_PATH, 'processed/sample.reef.report.year.RData'))

        ## Get the spatio-temporal meta data for each reef/visit
        samp.reef <- sample.reef.report.year %>%
            group_by(REEF, DEPTH, VISIT_NO, REPORT_YEAR) %>%
            summarise(across(c(LATITUDE, LONGITUDE, Date), mean)) %>%
            ungroup() %>%
            mutate(DEPTH.f = factor(case_when(DEPTH >3 ~"deep slope",
                                              DEPTH <= 3 ~"shallow slope")),
                   REEF.d = factor(paste(REEF, DEPTH.f))) %>%
            dplyr::select(-DEPTH) %>%
            suppressMessages() %>%
            suppressWarnings()
            
           
        GBR.LOF <- GBR.LOF.novelty %>%
            dplyr::select(Novelty.df.sum) %>%
            unnest(cols = c(Novelty.df.sum)) %>%
            bind_rows(GBR.LOF.deviation %>%
                      dplyr::select(Deviation.df.sum) %>%
                      unnest(cols = c(Deviation.df.sum))
                      ) %>%
            left_join(reef.nmds %>%
                      dplyr::select(data) %>%
                      unnest(cols = c(data)) %>%
                      mutate(REPORT_YEAR = as.numeric(as.character(Year)))
                      ) %>%
            mutate(threshold_detection = ifelse(LOF > 0.5, FALSE, TRUE)) %>%
            group_by(REEF.d) %>%
            arrange(REPORT_YEAR) %>%
            ungroup() %>%
            left_join(samp.reef) %>%
            suppressMessages() %>%
            suppressWarnings()

        save(GBR.LOF,
              file = paste0(DATA_PATH, 'modelled/GBR.LOF.RData'))


        load(file = paste0(DATA_PATH, 'primary/video_codes.RData'))
        GBR.LOF.table <- GBR.LOF.novelty %>%
            dplyr::select(Novelty.table) %>%
            unnest(cols = c(Novelty.table)) %>%
            bind_rows(GBR.LOF.deviation %>%
                      dplyr::select(Deviation.table) %>%
                      unnest(cols = c(Deviation.table))
                      ) %>%
            left_join(samp.reef) %>%
            left_join(video_codes %>%
                      dplyr::select(COMP_2021, COMP_2021_DESCRIPTION, GROUP_DESC) %>%
                      distinct(),
                      by = c("Taxon" = "COMP_2021")) %>%
            mutate(GROUP_DESC = str_to_lower(GROUP_DESC)) %>% 
            suppressMessages() %>%
            suppressWarnings() 

        save(GBR.LOF.table,
              file = paste0(DATA_PATH, 'modelled/GBR.LOF.table.RData'))
        write_csv(GBR.LOF.table, file = paste0(TABS_PATH, "/Composition_change.csv"))
 
}

CI_models_CO_calc_distance <- function() {
   

        load(file = paste0(DATA_PATH, 'modelled/coral.comms.RData'))

        distance.df <- coral.comms %>%
            ungroup() %>%
            group_by(REEF.d) %>%
            nest() %>%
            mutate(Composition_data = map(.x = data,
                                          .f = ~ .x %>% 
                                              dplyr::select(-any_of(c("P_CODE", "REEF.d",
                                                                      "VISIT_NO",
                                                                      "StartYR", "EndYR",
                                                                      "SurveyDuration",
                                                                      "cover", "adj.cover"))) %>%
                                          tibble::column_to_rownames(var = "REPORT_YEAR")),
                   Dissimilarity = map(.x = Composition_data,
                                       .f = ~ as.matrix(vegdist(.x)))
                   ) 

        save(distance.df,
              file = paste0(DATA_PATH, 'modelled/distance.df.RData'))
       
}


CI_model_CI_standardise <- function() {
    

        load(file = paste0(DATA_PATH, 'modelled/GBR.LOF.RData'))
        mods <- GBR.LOF %>%
            group_by(REEF.d) %>%
            summarise(data = list(cur_data_all()), .groups = "drop") %>%
            mutate(Scores = map(.x = data,
                                 .f = ~ .x %>%
                                     mutate(fYEAR = factor(REPORT_YEAR),
                                            .draw = 1,
                                            Metric = ifelse(k == 3, 'Critical', 'Reference')
                                            ) %>%
                                     dplyr::select(fYEAR, REEF.d, .draw,
                                                   .value = LOF,
                                                   REEF, DEPTH.f,
                                                   Metric)
                                ),
                   Summary = map(.x = Scores,
                                 .f = ~ .x %>%
                                     dplyr::select(-any_of(c(
                                                        "P_CODE",
                                                        "Model",
                                                        "value",
                                                        "baseline",
                                                        "DEPTH.f"))) %>%
                                     group_by(fYEAR, REEF, REEF.d, Metric) %>%
                                     ## summarise_draws(median = ~ median(.x, na.rm = TRUE),
                                     ##                 mean = ~ mean(.x, na.rm = TRUE),
                                     ##                 sd = ~ sd(.x, na.rm = TRUE),
                                     ##                 HDInterval::hdi,
                                     ##                 `p<0.5` = ~ mean(.x < 0.5, na.rm = TRUE)
                                     ##                 )
                                     summarise(variable = '.value',
                                               tidybayes::median_hdci(.value, na.rm = TRUE),
                                               mean = mean(.value, na.rm = TRUE),
                                               sd = sd(.value, na.rm = TRUE),
                                               `p<0.5` = mean(.value < 0.5, na.rm = TRUE)
                                               ) %>%
                                     dplyr::select(fYEAR, REEF, REEF.d, Metric, variable,
                                                   median = y, mean, sd,
                                                   lower = ymin, upper = ymax, `p<0.5`)
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
              file = paste0(DATA_PATH, 'modelled/CO__scores_reef_year.RData'))
 
}

CI_models_CO_aggregation <- function(level = 'NRM') {
  
  
  mods <- get(load(file = paste0(DATA_PATH, 'modelled/CO__scores_reef_year.RData')))
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
    summarise(.value = mean(.value, na.rm=TRUE)) %>%
    ungroup() |> 
    group_by(fYEAR, !!sym(level), Metric) |> 
    summarise(Score=median(.value, na.rm=TRUE),
              q95=quantile(.value, 0.95, na.rm=TRUE)) |> 
    mutate(BelowPar=ifelse(q95<0.5, 0,1),
           Indicator="Composition") |> 
    ungroup()
  
  save(mods,
       file = paste0(DATA_PATH, 'modelled/CO__scores_', level,'_year.RData'))
  
  
}



############# Plotting functions supplied by the UQ team ############################################
theme_Publication <- function(base_size=14, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.8, "cm"),
            legend.spacing = unit(0, "cm"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#1153aa","#ff9927","#4ac14a","#ee2b1a",
                                                            "#662506","#bd25e2","#fc807f","#9638a4","#ffe319", "#b1b1b1")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#1153aa","#ff9927","#4ac14a","#ee2b1a",
                                                              "#662506","#bd25e2","#fc807f","#9638a4","#ffe319", "#b1b1b1")), ...)
  
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

titleGrob <- function(label, margin=unit(c(b=5, l=0, t=2, r=0), "line"), ..., debug=FALSE){
  library(gtable)
  tg <- textGrob(label, ...)
  w <- grobWidth(tg)
  h <- grobHeight(tg)
  
  g <- gtable("title",
              widths = unit.c(margin[2], w, margin[4]),
              heights = unit.c(margin[3], h, margin[1]), respect = FALSE)
  if(debug){
    rg <- rectGrob()
    pos <- expand.grid(x=1:3, y=1:3)[-5,]
    g <- gtable_add_grob(g, list(rg, rg, rg, rg, rg, rg, rg, rg), t=pos$y, l=pos$x)
  }
  gtable_add_grob(g, tg, t=2, l = 2)
}

gg_circle <- function(r, xc, yc, color="black", fill=NA, ...) {
  x <- xc + r*cos(seq(0, pi, length.out=100))
  ymax <- yc + r*sin(seq(0, pi, length.out=100))
  ymin <- yc + r*sin(seq(0, -pi, length.out=100))
  annotate("ribbon", x=x, ymin=ymin, ymax=ymax, color=color, fill=fill, ...)
}

Comb.plot <- function(df, reef.name, pool.k) {
  Dissimilarity.outplot <- ggplot() +
    geom_line(data=coral.cover[coral.cover$REEF==reef.name,],
              aes(x=REPORT_YEAR, y=adj.cover), color='black', size=1) +
    geom_line(data=GBR.LOF[GBR.LOF$Reef==reef.name &
                             GBR.LOF$k==3,],
              aes(x=Year, y=LOF), color='coral', size=1) +
    geom_line(data=GBR.LOF[GBR.LOF$Reef==reef.name &
                             GBR.LOF$k==6,],
              aes(x=Year, y=LOF), color='dodgerblue3', size=1) +
    geom_hline(yintercept=0.5, linetype=2) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(name='Coral community indicator index',
                       sec.axis=sec_axis(~.*100, name='Coral cover')) +
    theme_Publication() +
    theme(panel.border=element_rect(color='black',fill=NA, size=1),
          # axis.title.x=element_blank(),
          # axis.title.y=element_blank(),
          plot.title=element_text(hjust=0)) +
    xlab('Year') +
    ylab('Coral community indicator index') +
    ggtitle('(a)')
  
  NMDS.outplot <- ggplot() +
    geom_segment(aes(x=unique(df[df$Reef==reef.name,
                                 c('Year','NMDS1','NMDS2')])[1L,'NMDS1'], 
                     xend=((unique(df[df$Reef==reef.name,
                                      c('Year','NMDS1','NMDS2')])[1L,'NMDS1']) +
                             (unique(df[df$Reef==reef.name,
                                        c('Year','NMDS1','NMDS2')])[2L,'NMDS1']))/2,
                     y=unique(df[df$Reef==reef.name,
                                 c('Year','NMDS1','NMDS2')])[1L,'NMDS2'],
                     yend=((unique(df[df$Reef==reef.name,
                                      c('Year','NMDS1','NMDS2')])[1L,'NMDS2']) +
                             (unique(df[df$Reef==reef.name,
                                        c('Year','NMDS1','NMDS2')])[2L,'NMDS2']))/2),
                 arrow=arrow(type='closed', length=unit(5,'mm'))) +
    geom_path(data=unique(df[df$Reef==reef.name,
                             c('Year','NMDS1','NMDS2')]),
              aes(x=NMDS1, y=NMDS2)) +
    geom_point(data=df[df$Reef==reef.name &
                         df$Year %in% c(1995:2016),],
               aes(x=NMDS1, y=NMDS2), color='forestgreen', size=7, alpha=0.5) +
    geom_point(data=unique(df[df$Reef==reef.name,
                              c('Year','NMDS1','NMDS2')]),
               aes(x=NMDS1, y=NMDS2)) +
    geom_point(data=df[df$Reef==reef.name &
                         df$k==3 &
                         df$threshold_detection=='Yes',],
               aes(x=NMDS1, y=NMDS2), color='coral', size=5) +
    geom_point(data=df[df$Reef==reef.name &
                         df$k==6 &
                         df$threshold_detection=='Yes',],
               aes(x=NMDS1, y=NMDS2), color='dodgerblue3', size=3) +
    geom_text_repel(data=unique(df[df$Reef==reef.name,
                                   c('Year','NMDS1','NMDS2')]),
                    aes(x=NMDS1, y=NMDS2, label=Year)) +
    theme_Publication() +
    theme(panel.border=element_rect(color='black',fill=NA, size=1),
          # axis.title.x=element_blank(),
          # axis.title.y=element_blank(),
          plot.title=element_text(hjust=0)) +
    xlab('NMDS1') +
    ylab('NMDS2') +
    ggtitle('(b)')
  
  outplot <- (Dissimilarity.outplot | NMDS.outplot)
  
  outplot <- outplot +
    plot_annotation(title=reef.name) &
    theme(text=element_text(size=15, face='bold'))
  
  plot(outplot)
  
  ggsave(filename=paste0('Figures/',
                         reef.name, '.png'),
         outplot, height=150, width=300,
         units='mm', device='png', dpi=300)
  
  return(outplot)
}

############# End of the Plotting functions supplied by the UQ team #################################

CI_models_CO_trajectory <- function() {
    

        load(file = paste0(DATA_PATH, 'modelled/coral.comms.RData'))

        save(distance.df,
              file = paste0(DATA_PATH, 'modelled/distance.df.RData'))
 
}
