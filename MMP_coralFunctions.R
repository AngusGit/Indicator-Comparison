loadPackages = function() {
  require(scales)  
  require(brms)
  require(zoo)
  require(mgcv)
  require(gmodels)
  require(lubridate)
  require(broom)
  require(coda)
  require(knitr)
  require(kableExtra)
  require(gtable)  #temperature plots
  require(TeachingDemos) #temperature plots
  require(maptools) #temperature plots
  require(cowplot) # plot_grid for subregion indicator panel
  require(gridExtra)# regional disease boxplots
  require(grid) # Discharge plots
  require(MuMIn) # discharge to index scores
  require(lme4) #discharge to index scores
  require(lmerTest) #discharge to index scores
  require(rgdal) #maps
  require(rgeos) #maps
  require(gamlss) # scores to env analysis
  require(betareg) #scores to env analysis
  require(gdata)   #scores to env plotting
  require(raster)
  #require(mapping)
  require(ggmap)
  
  require(tidyverse)
  
  
 

  
  
}

rollMean = function(x, width=4) {
    #print(x)
    #x = x %>% arrange(-yr)
    #print(x %>% dplyr::select(REEF,yr,change.score.linear))
    dat <- NULL
    for (i in unique(x$yr)) {
      xx = x %>% filter(yr<(i+1) & yr>(i-width))
      curYear <- max(xx$yr, na.rm=TRUE)
      firstYear <- curYear-3
      DateRange <- paste(firstYear, '-',curYear)
        cs.linear.mean = mean(xx$change.score.linear, na.rm=TRUE)
        dat <- rbind(dat,data.frame(yr=max(xx$yr),DateRange=DateRange,NumYearsInRange=nrow(xx),
                                    Change.score.linear.mean=cs.linear.mean                
        ))
    }
     dat
}

MMP_generateGrades <- function(x) {
  ifelse(is.na(x),'NA',ifelse(x>=0.805, 'A', ifelse(x>=0.605, 'B', ifelse(x>=0.405, 'C',  ifelse(x>=0.205, 'D', 'E'))))) #updated 2018 to reflect grading based on rounding to whole percentage points eg: score of 0.204 rounds to 0.20.
}

rs = function(dat) {
    s = vector(length=nrow(dat))
    for (i in 1:nrow(dat)){
        if (is.na(dat[i,1])) {
            s[i]=NA
        } else {
            s[i]=with(dat, scales::rescale(x[i], from=c(f1[i],f2[i]), to=c(t1[i],t2[i])))
        }
    }
    s
}


boot_signal <-function(Var) {
  Confidence=scales::rescale(sqrt(Var), from=c(0,sqrt(var(c(rep(0,5000),rep(1,5000))))), to=c(1,0))
  Signal = ifelse(Confidence>=0.8,5,ifelse(Confidence>=0.6,4,ifelse(Confidence>=0.4,3,ifelse(Confidence>=0.2,2,1))))
  Signal
}

# a function to estimate confidence intervals from the following boot strapped distributions
conf <- function(x) {
  ## Tabuate the values
  tt<-table(x)
  ## Sum up the tabulated values divided by the minimum
  ## This should give the number of values that wetn into making up
  ## the non-bootstrap sample
  l <- sum(tt/min(tt))
  #d <- unique(x)
  #l <- length(d)
  s <- NULL
  for (i in 1:1000) {
    s <- c(s, mean(sample(x,size=l, replace=TRUE), na.rm=TRUE))
  }
  #mean(s, na.rm=TRUE)
  q<-quantile(s, na.rm=TRUE, p=c(0.025,0.975))
  c(q[1],q[2])
}

boot_sum <- function(x) {
  require(scales)
  Mean=mean(x$Boot, na.rm=TRUE)
  Var=var(x$Boot, na.rm=TRUE)
  Confidence=scales::rescale(sqrt(Var), from=c(0,sqrt(var(c(rep(0,5000),rep(1,5000))))), to=c(1,0))
  Signal = ifelse(Confidence>0.805,5,ifelse(Confidence>0.605,4,ifelse(Confidence>0.405,3,ifelse(Confidence>0.205,2,1)))) #updated 2018 to reflect change in grades at scores rounded to two decimal places 0.205 = 0.20 etc
  Grade=MMP_generateGrades(Mean)
  data.frame(Mean=Mean, Var=Var, Signal=Signal,Grade=Grade)
}

trafficLightPalette <- (c('#FF0000','#FFC000','#FFFF00','#92D050','#00B050'))




MMP_coralWorm <- function(dat, region=NULL, Subregion=NULL, title=FALSE, ytitle=FALSE, guide=FALSE, subLabels=FALSE, indicators=FALSE, nrow=1, strip=FALSE, point.size=4, line.size=0.5) {
    require(lubridate)
    nFac=1
    if(!is.null(region)) {
        if(class(region)!='logical') {
            dat <- dat %>% filter(Region==region)
            if (indicators!=FALSE) indicators <- indicators %>% filter(Region==region)
            ti=region
        } 
    }
    if(!is.null(Subregion)) {
        if(class(Subregion)!='logical') {
            dat <- dat %>% filter(subregion==Subregion) %>% droplevels
            if (indicators!=FALSE) indicators <- indicators %>% filter(subregion==Subregion)
            ti=Subregion
        }
    }
    load('data/primary/CoralDateRanges.RData')
    s=data.frame(Date=CoralDateRanges[[1]], subLabels=subLabels)
  
    lims <- rev(LETTERS[1:5]) 
    
    dtt <- filter(dat, Year>=year(round_date(CoralDateRanges[[1]], unit='year'))) %>% droplevels
    
    p <- ggplot(dtt, aes(y=Score, x=as.Date(paste0(Year,'-01-01')))) +
        ylim(0,1)+
        xlim(CoralDateRanges$min, CoralDateRanges$max)+
        geom_hline(yintercept=seq(0.205,0.805,by=0.2), color='grey80', linetype='dashed')+ #changed 2018 to reflect score changes at 0.205, 0.405 etc. rather than 0.2, 0.4 
        scale_x_date('',limits=c(CoralDateRanges[[1]],CoralDateRanges[[2]]))+
        scale_y_continuous('', breaks=seq(0,1,by=0.2), limits=c(0,1),expand=c(0,0))+
        scale_fill_manual('', values=trafficLightPalette, limits=lims, guide=FALSE)+
        theme_classic(10)+
        theme(strip.background=element_blank(),
              strip.text=element_blank(),
              axis.title.x=element_blank(),
              plot.margin=unit(c(0.2,0.5,0,1),"lines"),
              panel.spacing=unit(c(1),"lines"),
              axis.title.y=element_text(vjust=1, size=rel(1.25)),
              axis.line.x=element_line(),axis.line.y=element_line()
                                        #axis.title.y=element_text(vjust=2)
              )
    if (indicators!=FALSE) {
        dt <- filter(indicators, Year>=year(round_date(CoralDateRanges[[1]], unit='year')))
        p <- p+
            geom_line(data=dt, aes(y=Score, x=as.Date(paste0(Year,'-01-01')),color=Indicator, size=Insicator),size=line.size)+
            scale_color_manual('Indicator',breaks=c('Cover change', 'Coral composition', 'Coral cover', 'Juvenile density', 'Macroalgal cover','Overall'),
                               values=c('cyan','darkblue','purple','pink','brown','black'))+
            scale_size_manual('Indicator',breaks=c('Cover change', 'Coral composition', 'Coral cover', 'Juvenile density', 'Macroalgal cover','Overall'),
                              values=c(rep(0.5,5),0.75))+
            scale_shape_manual('Indicator',
                               breaks=c('Cover change', 'Coral composition', 'Coral cover', 'Juvenile density', 'Macroalgal cover','Overall'),
                               values=c(22,23,24,25,11,19))
        p=p+geom_blank(data=dt, aes(y=Score, x=as.Date(paste0(Year,'-01-01')), shape='Overall', linetype='Overall', color='Overall', size='Overall'))+
            theme(legend.key.width=unit(1,'cm'))
        
    } 
    p <- p+geom_linerange(aes(ymin=Lower, ymax=Upper))+
        geom_line(size=0.75) + geom_point(data=dtt,aes(fill=Grade), shape=21,size=point.size)+
        scale_x_date('',limits=c(CoralDateRanges[[1]],CoralDateRanges[[2]]))
    
    if(!is.null(region)) {
        if(region==TRUE) {
            p = p + facet_wrap(~Region, scales='free', nrow=nrow)
            nFac <- length(unique(dat$Region))
            if(subLabels==TRUE) {
                subLabels=paste0(letters[1:nFac],')')
                s <- data.frame(Region=unique(dat$Region), Date=CoralDateRanges[[1]], subLabels=subLabels)
            }
        }
    }
    if(!is.null(Subregion)) {
        if(Subregion==TRUE) {
            p = p + facet_wrap(~subregion, scales='free', nrow=nrow)
            nFac <- length(unique(dat$subregion))
            if(subLabels==TRUE) {
                subLabels=paste0(letters[1:nFac],')')
                s <- data.frame(subregion=unique(dat$subregion), Date=CoralDateRanges[[1]], subLabels=subLabels)
            }
        }
    }
    if(subLabels!=FALSE) {
        p = p + geom_text(data=s, aes(y=Inf, x=Date, label=subLabels),hjust=-0.2, vjust=1)
    }
    ## It only makes sense to have a title on a single figure, so
    ## region and subregion must not be both TRUE or NULL
    print(title)
    if (title==TRUE || is.character(title)) {
        if (is.character(title)) ti=title
        p <- p+ggtitle(ti)
    }
    if (class(ytitle)=='expression') {
        p<-p+scale_y_continuous(ytitle, breaks=seq(0,1,by=0.2),limits=c(0,1))
    } else if(ytitle!=FALSE) {
        if(ytitle==TRUE) ytitle='Score'
        p<-p+scale_y_continuous(ytitle, breaks=seq(0,1,by=0.2),limits=c(0,1))
    }
    if(strip==TRUE) p <- p+theme(strip.text=element_text())
    p
}
###

MMP_coralWorm.comp <- function(dat, region=NULL, Subregion=NULL, title=FALSE, ytitle=FALSE, guide=FALSE, subLabels=FALSE, nrow=1, strip=FALSE, point.size=4, line.size=0.5) {
  require(lubridate)
  nFac=1
  if(!is.null(region)) {
    if(class(region)!='logical') {
      dat <- dat %>% filter(Region==region)
      ti=region
    } 
  }
  if(!is.null(Subregion)) {
    if(class(Subregion)!='logical') {
      dat <- dat %>% filter(subregion==Subregion) %>% droplevels
      ti=Subregion
    }
  }
  load('data/primary/CoralDateRanges.RData')
  s=data.frame(Date=CoralDateRanges[[1]], subLabels=subLabels)
  
  lims <- rev(LETTERS[1:5]) 
  
  dtt <- filter(dat, Year>=year(round_date(CoralDateRanges[[1]], unit='year'))) %>% droplevels
  
  p <- ggplot(dtt, aes(y=Score, x=as.Date(paste0(Year,'-01-01')))) +
    ylim(0,1)+
    xlim(CoralDateRanges$min, CoralDateRanges$max)+
    geom_hline(yintercept=seq(0.205,0.805,by=0.2), color='grey80', linetype='dashed')+ 
    scale_x_date('',limits=c(CoralDateRanges[[1]],CoralDateRanges[[2]]))+
    scale_y_continuous('', breaks=seq(0,1,by=0.2), limits=c(0,1),expand=c(0,0))+
    scale_fill_manual('', values=trafficLightPalette, limits=lims, guide=FALSE)+
    theme_classic(10)+
    theme(strip.background=element_blank(),
          strip.text=element_blank(),
          axis.title.x=element_blank(),
          plot.margin=unit(c(0.2,0.5,0,1),"lines"),
          panel.spacing=unit(c(1),"lines"),
          axis.title.y=element_text(vjust=1, size=rel(1.25)),
          axis.line.x=element_line(),axis.line.y=element_line()
          #axis.title.y=element_text(vjust=2)
    )
  
 
  p <- p+
    geom_linerange(data=dtt %>% filter(indicators=='no.comp'), aes(ymin=Lower, ymax=Upper, x=as.Date(paste0(Year,'-03-01'))), alpha=0.5, color='blue')+
    geom_linerange(data=dtt %>% filter(indicators=='all'), aes(ymin=Lower, ymax=Upper), alpha=0.5)+
    geom_line(data=dtt %>% filter(indicators=='no.comp'),aes(y=Score, x=as.Date(paste0(Year,'-03-01'))), size=line.size, color='blue') +
    geom_line(data=dtt %>% filter(indicators=='all'), size=line.size) + 
    geom_point(data=dtt %>% filter(indicators=='no.comp'), aes(y=Score, x=as.Date(paste0(Year,'-03-01')),fill=Grade), shape=21,size=point.size) +
    geom_point(data=dtt %>% filter(indicators=='all'),aes(fill=Grade), shape=21,size=point.size)+
   
    scale_x_date('',limits=c(CoralDateRanges[[1]],CoralDateRanges[[2]]))
  
     
 if(!is.null(region)) {
    if(region==TRUE) {
      p = p + facet_wrap(~Region, scales='free', nrow=nrow)
      nFac <- length(unique(dat$Region))
      if(subLabels==TRUE) {
        subLabels=paste0(letters[1:nFac],')')
        s <- data.frame(Region=unique(dat$Region), Date=CoralDateRanges[[1]], subLabels=subLabels)
      }
    }
  }
  if(!is.null(Subregion)) {
    if(Subregion==TRUE) {
      p = p + facet_wrap(~subregion, scales='free', nrow=nrow)
      nFac <- length(unique(dat$subregion))
      if(subLabels==TRUE) {
        subLabels=paste0(letters[1:nFac],')')
        s <- data.frame(subregion=unique(dat$subregion), Date=CoralDateRanges[[1]], subLabels=subLabels)
      }
    }
  }
  if(subLabels!=FALSE) {
    p = p + geom_text(data=s, aes(y=Inf, x=Date, label=subLabels),hjust=-0.2, vjust=1)
  }
  ## It only makes sense to have a title on a single figure, so
  ## region and subregion must not be both TRUE or NULL
  print(title)
  if (title==TRUE || is.character(title)) {
    if (is.character(title)) ti=title
    p <- p+ggtitle(ti)
  }
  if (class(ytitle)=='expression') {
    p<-p+scale_y_continuous(ytitle, breaks=seq(0,1,by=0.2),limits=c(0,1))
  } else if(ytitle!=FALSE) {
    if(ytitle==TRUE) ytitle='Score'
    p<-p+scale_y_continuous(ytitle, breaks=seq(0,1,by=0.2),limits=c(0,1))
  }
  if(strip==TRUE) p <- p+theme(strip.text=element_text())
  p
}
###

MMP_coralWormBLANK <- function(dat, region=NULL, Subregion=NULL, title=FALSE, ytitle=FALSE, guide=FALSE, subLabels=FALSE, indicators=FALSE, nrow=1, strip=FALSE, point.size=4, line.size=0.5) {
  require(lubridate)
  nFac=1
  if(!is.null(region)) {
    if(class(region)!='logical') {
      dat <- dat %>% filter(Region==region)
      if (indicators!=FALSE) indicators <- indicators %>% filter(Region==region)
      ti=region
    } 
  }
  if(!is.null(Subregion)) {
    if(class(Subregion)!='logical') {
      dat <- dat %>% filter(subregion==Subregion) %>% droplevels
      if (indicators!=FALSE) indicators <- indicators %>% filter(subregion==Subregion)
      ti=Subregion
    }
  }
  load('data/primary/CoralDateRanges.RData')
  s=data.frame(Date=CoralDateRanges[[1]], subLabels=subLabels)
  
  lims <- rev(LETTERS[1:5]) 
  
  dtt <- filter(dat, Year>=year(round_date(CoralDateRanges[[1]], unit='year'))) %>% droplevels
  
  p <- ggplot(dtt, aes(y=Score, x=as.Date(paste0(Year,'-01-01')))) +
    ylim(0,1)+
    xlim(CoralDateRanges$min, CoralDateRanges$max)+
    #geom_hline(yintercept=seq(0.205,0.805,by=0.2), color='grey80', linetype='dashed')+ #changed 2018 to reflect score changes at 0.205, 0.405 etc. rather than 0.2, 0.4 
    scale_x_date('',limits=c(CoralDateRanges[[1]],CoralDateRanges[[2]]))+
    scale_y_continuous('', breaks=seq(0,1,by=0.2), limits=c(0,1),expand=c(0,0))+
    scale_fill_manual('', values=trafficLightPalette, limits=lims, guide=FALSE)+
    theme_classic(10)+
    theme(strip.background=element_blank(),
          strip.text=element_blank(),
          axis.title.x=element_blank(),
          plot.margin=unit(c(0.2,0.5,0,1),"lines"),
          panel.spacing=unit(c(1),"lines"),
          axis.title.y=element_text(vjust=1, size=rel(1.25)),
          axis.line.x=element_line(),axis.line.y=element_line()
          #axis.title.y=element_text(vjust=2)
    )
  if (indicators!=FALSE) {
    dt <- filter(indicators, Year>=year(round_date(CoralDateRanges[[1]], unit='year')))
    p <- p+
      geom_line(data=dt, aes(y=Score, x=as.Date(paste0(Year,'-01-01')),color=Indicator, size=Insicator),size=line.size)+
      scale_color_manual('Indicator',breaks=c('Cover change', 'Coral composition', 'Coral cover', 'Juvenile density', 'Macroalgal cover','Overall'),
                         values=c('cyan','darkblue','purple','pink','brown','black'))+
      scale_size_manual('Indicator',breaks=c('Cover change', 'Coral composition', 'Coral cover', 'Juvenile density', 'Macroalgal cover','Overall'),
                        values=c(rep(0.5,5),0.75))+
      scale_shape_manual('Indicator',
                         breaks=c('Cover change', 'Coral composition', 'Coral cover', 'Juvenile density', 'Macroalgal cover','Overall'),
                         values=c(22,23,24,25,11,19))
    p=p+geom_blank(data=dt, aes(y=Score, x=as.Date(paste0(Year,'-01-01')), shape='Overall', linetype='Overall', color='Overall', size='Overall'))+
      theme(legend.key.width=unit(1,'cm'))
    
  } 
  p <- p+geom_linerange(aes(ymin=Lower, ymax=Upper))+
    geom_line(size=0.75) + geom_point(data=dtt,aes(fill=Grade), shape=21,size=point.size)+
    scale_x_date('',limits=c(CoralDateRanges[[1]],CoralDateRanges[[2]]))
  
  if(!is.null(region)) {
    if(region==TRUE) {
      p = p + facet_wrap(~Region, scales='free', nrow=nrow)
      nFac <- length(unique(dat$Region))
      if(subLabels==TRUE) {
        subLabels=paste0(letters[1:nFac],')')
        s <- data.frame(Region=unique(dat$Region), Date=CoralDateRanges[[1]], subLabels=subLabels)
      }
    }
  }
  if(!is.null(Subregion)) {
    if(Subregion==TRUE) {
      p = p + facet_wrap(~subregion, scales='free', nrow=nrow)
      nFac <- length(unique(dat$subregion))
      if(subLabels==TRUE) {
        subLabels=paste0(letters[1:nFac],')')
        s <- data.frame(subregion=unique(dat$subregion), Date=CoralDateRanges[[1]], subLabels=subLabels)
      }
    }
  }
  if(subLabels!=FALSE) {
    p = p + geom_text(data=s, aes(y=Inf, x=Date, label=subLabels),hjust=-0.2, vjust=1)
  }
  ## It only makes sense to have a title on a single figure, so
  ## region and subregion must not be both TRUE or NULL
  print(title)
  if (title==TRUE || is.character(title)) {
    if (is.character(title)) ti=title
    p <- p+ggtitle(ti)
  }
  if (class(ytitle)=='expression') {
    p<-p+scale_y_continuous(ytitle, breaks=seq(0,1,by=0.2),limits=c(0,1))
  } else if(ytitle!=FALSE) {
    if(ytitle==TRUE) ytitle='Score'
    p<-p+scale_y_continuous(ytitle, breaks=seq(0,1,by=0.2),limits=c(0,1))
  }
  if(strip==TRUE) p <- p+theme(strip.text=element_text())
  p
}


#########################################################################
## The following function calculates the expected hard coral and other ## ## coral cover based on the Bayesian model and the estK, Acr, Ohc and  ##
## SC.  It calculates the Mean, Median and lower/upper confidence      ##
## bounds of what would be expected...                                 ##
## Parameters:                                                         ##
##    - mod.acr   the Bayesian Acropora growth model                   ##
##    - mod.ohc   the Bayesian Other coral growth model                ##
##    - dat       the specific set of reef/years that are not          ##
##                 disturbed etc                                       ##
## Returns:                                                            ##
##    - Mean      mean expected growth rate                            ##
##    - Median    median expected growth rate                          ##
##    - lower     lower interval of expected growth rate               ##
##    - upper     upper interval of expected growth rate               ##
##    - Posterior all samples                                          ##
#########################################################################

addHCOHC <- function(mod.acr, mod.ohc, dat,include.region=TRUE, include.reef=FALSE){
  require(plyr)
  require(coda)
  rAc <- mod.acr$BUGSoutput$sims.list[['rAc']]
  alpha <- mod.acr$BUGSoutput$sims.list[['alpha']]
 
  #wAcMa <- mod.acr$BUGSoutput$sims.list[['wAcMa']]
  wAcMa <- rep(0, length(alpha))
  res <- mod.acr$BUGSoutput$sims.list[['res']]
  vRegion <- mod.acr$BUGSoutput$sims.list[['vary_region']]#[,3] #Changed
  regMatrix <- model.matrix(~-1+NRM_REGION, dat)#[,1] #Changed
  #reefMatrix <- model.matrix(~-1+REEF, dat)
  
  pq.acr <- matrix(NA,length(alpha), length(dat$Acr2))
  for (i in 1:nrow(alpha)) {
    rAC <- rAc[i]
    if (include.region==TRUE) rAC <- alpha[i] + t(cbind(vRegion)[i,] %*% t(cbind(regMatrix))) #+ t(vReef[i,] %*% t(reefMatrix))
    if (include.reef==TRUE) rAC <- alpha[i] + t(cbind(vRegion)[i,] %*% t(cbind(regMatrix))) + t(vReef[i,] %*% t(reefMatrix))
    pq.acr[i,] <- with(dat, exp(rAC +log(Acr) -((rAC/log(estK))*log(Acr+Ohc+SC+wAcMa[i]*MA)) +res[1]))
    
  }
  # other coral model
  alpha <- mod.ohc$BUGSoutput$sims.list[['alpha']]
  rAc <- mod.ohc$BUGSoutput$sims.list[['rAc']] #this is ok the rate is for ohc
  
  wAcMa <- rep(0, length(alpha))
  res <- mod.ohc$BUGSoutput$sims.list[['res']]
  vRegion <- mod.ohc$BUGSoutput$sims.list[['vary_region']][,3] 
  regMatrix <- model.matrix(~-1+NRM_REGION, dat)[,1] 
  pq.ohc <- matrix(NA,length(alpha), length(dat$Acr2)) # Acr2 simply used for convenience
  for (i in 1:nrow(alpha)) {
    rAC <- rAc[i]
    if (include.region==TRUE) rAC <- alpha[i] + t(cbind(vRegion)[i,] %*% t(cbind(regMatrix))) 
    if (include.reef==TRUE) rAC <- alpha[i] + t(cbind(vRegion)[i,] %*% t(cbind(regMatrix))) + t(vReef[i,] %*% t(reefMatrix)) 
    pq.ohc[i,] <- with(dat, exp(rAC +log(Ohc) -((rAC/log(estK))*log(Acr+Ohc+SC+wAcMa[i]*MA)) +res[1]))
  }
  
  pq<-pq.acr + pq.ohc
  
  list(
    Summary=adply(pq, 2, function(x) {
          data.frame(Mean=mean(x,na.rm=TRUE), Median=median(x,na.rm=TRUE), HPDinterval(as.mcmc(x)))
    }),
    Posterior=pq
  )
}

#########

scaleChangeScores <- function(data,a, fromBoundary=FALSE) {

  data1 <- NULL
  for (i in 1:length(data$HC2_adj)) {
      data1 <- rbind(data1,scaleChangeScore(data[i,],a[,i], fromBoundary=fromBoundary))
  }
  data1
}


## ###################################################
## Function scales observed change in Hard coral cover (data$HC_adj)  relative to 
## Predicted change in cover (data$Median)
## scores 0 if there has been no positive change
## scores 1 if HC_adj is greater than twice the upper CI of the predicted change (mid point of vgood categorisation)
## scaled between 0.1 and 0.4, as equating to a v-poor to poor categorisation, for positive change less than lower 95%CI of Median
## scaled between 0.4 and 0.6, as equating to a moderate categorisation, for positive change between upper and lower 95%CI of Median
## scaled between 0.6 and 0.9, as equating to a good to v good categorisation, for positive change between upper and 2*upper 95%CI of Median
########################################################

scaleChangeScore <- function(data,a, fromBoundary=FALSE) {
  require(scales)
  ##Linear difference
  if (fromBoundary==TRUE) {
    if (data$HC > data$HC2_adj) {  ## changed from >= to > between visit 12 and visit 13 reports
        data$change.score.linear<-0
    } else {
      if (data$HC2_adj<=data$lower) {
        
        #deal with values lower than boundary
          data$change.score.linear <- rescale(data$lower-data$HC2_adj, from=c(0, data$lower-data$HC), to=c(0.4,0.1))
      
    }else{
        if (data$lower == 0 & data$upper==0) {
            data$change.score.linear=0.5    
        
    } else {
        if (data$HC2_adj>data$lower & data$HC2_adj<=data$upper) {
          
          #values between lower and upper
          data$change.score.linear <- rescale(data$upper-data$HC2_adj, from=c(0,data$upper-data$lower), to=c(0.6,0.4))      
          
    } else {
        if (data$HC2_adj>data$upper & data$HC2_adj<=(data$upper+(data$upper-data$Median))) { # changed < to <= between visit 12 and visit 13 reports
          
          #values between upper and cap
          data$change.score.linear <- rescale(((data$upper*2)-data$Median)-data$HC2_adj, from=c(0,data$upper-data$Median), to=c(0.9,0.6))
          
      } else {   
        if (data$HC2_adj>=(2*(data$upper)-data$Median)) {  
          
          #cap at 2 times the predicted change
          data$change.score.linear<-1
             }    
            }
            }
            }
            }
    }
  }
  data
}


#########################################################################
## The following function takes the coral data and (using the relevant ##
## depth growth model), it calculates the expected amount of coral for ##
## non-disturbance years.  Based on a comparison of observed coral     ##
## cover to expected coral cover, I then generate a number of possible ##
## scores.                                                             ##
## Parameters:                                                         ##
## Returns:                                                            ##
##   - data.frame    a data frame contraining the coral data with the  ##
##                   expected coral cover and scores                   ##
#########################################################################

MMP_combineHCOHC <- function(fromBoundary=FALSE, include.region=TRUE) {
  #load(file='data/processed/Stage2/acr.jags.RData')
  #load(file='data/processed/Stage2/oth.jags.RData')
  load(file='data/perpetual/acr.jags.2.RData')
  load(file='data/perpetual/oth.jags.2.RData')
  load(file='data/perpetual/acr.jags.5.RData')
  load(file='data/perpetual/oth.jags.5.RData')
  load(file='data/processed/change.reef.RData')
    
      
  change.scores<-NULL
  for (pcode in levels(change.reef$P_CODE)) {
    data.i <- subset(change.reef, P_CODE==pcode & !is.na(Date))
    for (ireef in levels(droplevels(data.i$REEF))) {
      data.ii <- subset(data.i, REEF==ireef)
      for (depth in unique(data.ii$DEPTH)) {
        data.iii <- subset(data.ii, DEPTH==depth)
        data.iii <- arrange(data.iii, VISIT_NO)
        data.iv <- data.iii
        data.iv$Dt.num <- decimal_date(data.iv$Date) 
        dt.num <- max(data.iv$Dt.num,na.rm=TRUE)  
        waterYear <- max(data.iv$yr, na.rm=TRUE) 
        DateRange <- paste(min(data.iv$yr, na.rm=TRUE), '-',max(data.iv$yr, na.rm=TRUE))
        vp<-subset(data.iv, DISTURBANCE2 %in% c("n","d")) 
 print(pcode)
 print(ireef)
 print(depth)
        if (dim(vp)[1]!=0) {
          if(depth==2) {
            a <- addHCOHC(acr.jags.2, oth.jags.2, vp,include.region=include.region)
            vp <- cbind(vp,subset(a[['Summary']], select=-X1))
          }
          if(depth==5) {
                          
            a<-addHCOHC(acr.jags.5, oth.jags.5, vp, include.region=include.region)
            vp <- cbind(vp,subset(a[['Summary']], select=-X1))
          }
          vp$DEPTH <- depth
          vp <- scaleChangeScores(vp,a[[2]], fromBoundary=fromBoundary)
          change.scores <- rbind(change.scores, vp)
        }
      
      }
    }
  }
  change.scores
}


MMP_analysisSedCoralLMER <- function(Variable, VariableName,region=NULL, Subregion=NULL, title=FALSE, ytitle=FALSE, guide=FALSE, subLabels=FALSE) {
  load(file='data/processed/allCoralData.reef.RData')
  require(lme4)
  dat = allCoralData.reef
  nFac=1
  if(!is.null(region)) {
    dat <- dat %>% filter(Region==region)
    ti=region
  }
  if(!is.null(Subregion)) {
    dat <- dat %>% filter(subregion==Subregion)
    ti=Subregion
  }
  data.i<-droplevels(dat)
  data.i$REEF <- with(data.i, interaction(P_CODE,REEF))
  data.i <- arrange(data.i, REEF,Dt.num)
  iVar=Variable
  ##fit the linear mixed effect model
  newdata <- NULL
  err <- 0
  #rm(i.gamm)
  eval(parse(text=paste0('data.ii <- droplevels(filter(data.i, !is.na(',iVar,')))')))
  eval(parse(text=paste0('data.means = data.ii %>% group_by(Region,subregion,Year) %>% summarize(',Variable,'=mean(',iVar,',na.rm=TRUE))')))
  if (Variable=='gen.composition') {
    i.data = data.ii %>% filter_(.dots=paste('!is.na(',iVar,')')) %>% mutate(Year=factor(Year)) %>%
      mutate_(.dots=setNames(paste(iVar),'Resp'))
    tryCatch(
      i.lmer <- lmer(Resp ~ Year + (1|REEF), data=i.data),
      error=function(x) err<<-1
    )
  } else {
    i.data = data.ii %>% filter_(.dots=paste('!is.na(',iVar,')')) %>% mutate(Year=factor(Year)) %>%
      mutate_(.dots=setNames(paste0('asin(sqrt(',iVar,'/100))'),'Resp'))
    tryCatch(
      i.lmer <-lmer(Resp ~ Year + (1|REEF), data=i.data),
      error=function(x) err<<-1
    )
  }
  newdata=data.frame(Year=levels(i.lmer@frame$Year))
  
  predFun <- function(fit) {
    predict(fit,newdata, re.form=~0)
  }
  bb <- bootMer(i.lmer,nsim=1000,FUN=predFun,seed=123)$t
  if (Variable!='gen.composition') {
    bb=100*sin(bb)^2
  }
  
  Xmat <- model.matrix(~Year, data=newdata)
  fit=bb %*% Xmat 
  dissim.matrix = data.frame(Mean=apply(fit,2,mean),lower=apply(fit,2,quantile, p=0.025),upper=apply(fit,2,quantile, p=0.975) )
  pred = data.frame(fit=apply(bb,2,mean),lower=apply(bb,2,quantile, p=0.025),upper=apply(bb,2,quantile, p=0.975))
  
  newdata = cbind(newdata, pred)
  newdata = newdata %>% mutate(Date=as.Date(paste0(Year,'-06-01')))
  data.i = data.i %>% mutate(Date=as.Date(lubridate:::date_decimal(Dt.num)))
  
  remove.missing=missing_span(data.i, iVar)
  eval(parse(text=paste0('data.i = data.i %>% filter(!is.na(',iVar,'1))')))
  load('data/primary/CoralDateRanges.RData')
  s=data.frame(Date=CoralDateRanges[[1]], subLabels=subLabels)
 
  p <- ggplot(newdata, aes(y=fit, x=Date)) +
    xlim(CoralDateRanges$min, CoralDateRanges$max)+
    geom_ribbon(aes(ymin=lower, ymax=upper, x=Date), fill='blue',alpha=0.3)+
    geom_line(aes(y=fit, x=Date), color='blue')+
    eval(parse(text=paste0('geom_line(data=data.i, aes(y=',iVar,'1, x=Date,linetype=factor(P_CODE),group=interaction(REEF)), alpha=0.3)')))+  
    scale_linetype_manual(values=c('solid','dashed'), guide=FALSE)+
    scale_x_date('', limits=c(CoralDateRanges[[1]],CoralDateRanges[[2]]))+
    scale_y_continuous(VariableName)+
    theme_classic(10)+
    theme(strip.background=element_blank(),
          strip.text=element_blank(),
          axis.title.x=element_blank(),
          plot.margin=unit(c(0.2,0.5,0,1),"lines"),
          panel.spacing=unit(c(0),"lines"),
          axis.title.y=element_text(vjust=1, size=rel(1.25)),
          axis.line.x=element_line(),axis.line.y=element_line()
    )
  if(nrow(remove.missing)>0) eval(parse(text=paste0('p <- p + geom_line(data=remove.missing, aes(y=',iVar,'1, x=Date,linetype=factor(P_CODE),group=interaction(REEF)), alpha=0.3, linetype="dashed")')))
  if(!is.null(region)) {
    if(region==TRUE) {
      p = p + facet_wrap(~Region, scales='free', nrow=1)
      nFac <- length(unique(dat$Region))
      if(subLabels==TRUE) {
        subLabels=paste0(letters[1:nFac],')')
        s <- data.frame(Region=unique(dat$Region), Date=CoralDateRanges[[1]], subLabels=subLabels)
      }
    }
  }
  if(!is.null(Subregion)) {
    if(Subregion==TRUE) {
      p = p + facet_wrap(~subregion, scales='free', nrow=1)
      nFac <- length(unique(dat$subregion))
      if(subLabels==TRUE) {
        subLabels=paste0(letters[1:nFac],')')
        s <- data.frame(s=unique(dat$subregion), Date=CoralDateRanges[[1]], subLabels=subLabels)
      }
    }
  }
  if(subLabels!=FALSE) {
    p = p + geom_text(data=s, aes(y=Inf, x=Date, label=subLabels),hjust=-0.2, vjust=1)
  }
  ## It only makes sense to have a title on a single figure, so
  ## region and subregion must not be both TRUE or NULL
  if (title==TRUE) p <- p+ggtitle(ti)
  if (class(ytitle)=='expression') {
    p<-p+scale_y_continuous(ytitle)
  } else if(ytitle!=FALSE) {
    if(ytitle==TRUE) ytitle='Mean'
    p<-p+scale_y_continuous(ytitle)
  }
  newdata=newdata %>% cbind(Var=as.character(VariableName),dat %>% select(Region, subregion) %>% distinct())
  list(g=p, m=data.means, p=newdata, dissim.matrix=dissim.matrix)
}


##################### lmer model including random reef and depth

MMP_analysisSedCoralLMER_depth <- function(Variable, VariableName,region=NULL, Subregion=NULL, title=FALSE, ytitle=FALSE, guide=FALSE, subLabels=FALSE) {
  load(file='data/processed/allCoralData.RData')
  load(file='data/primary/CoralDateRanges.RData')
  require(lme4)
  dat = allCoralData %>%filter(Year>=year(CoralDateRanges$min))
  nFac=1
  if(!is.null(region)) {
    dat <- dat %>% filter(Region==region)
    ti=region
  }
  if(!is.null(Subregion)) {
    dat <- dat %>% filter(subregion==Subregion)
    ti=Subregion
  }
  data.i<-droplevels(dat)
  data.i <- arrange(data.i, ReefDepth,Dt.num)
  iVar=Variable
  ##fit the linear mixed effect model
  newdata <- NULL
  err <- 0
  
  eval(parse(text=paste0('data.ii <- droplevels(filter(data.i, !is.na(',iVar,')))')))
  eval(parse(text=paste0('data.means = data.ii %>% group_by(Region,subregion,Year) %>% summarize(',Variable,'=mean(',iVar,',na.rm=TRUE))')))
  if (Variable=='gen.composition') {
    i.data = data.ii %>% filter_(.dots=paste('!is.na(',iVar,')')) %>% mutate(Year=factor(Year)) %>%
      mutate_(.dots=setNames(paste(iVar),'Resp'))
    tryCatch(
      
      i.lmer <- lmer(Resp ~ Year + (1|ReefDepth), data=i.data),
      error=function(x) err<<-1
    )
  } else {
    i.data = data.ii %>% filter_(.dots=paste('!is.na(',iVar,')')) %>% mutate(Year=relevel(factor(Year),ref="2010")) %>%
      mutate_(.dots=setNames(paste0('asin(sqrt(',iVar,'/100))'),'Resp'))
    tryCatch(
     i.lmer <-lmer(Resp ~ Year + (1|ReefDepth), data=i.data),
      #i.lmer <-glmer(Resp ~ Year + (1|ReefDepth), data=i.data, family=binomial),
      error=function(x) err<<-1
    )
  }
  newdata=data.frame(Year=levels(i.lmer@frame$Year))
 
  
  predFun <- function(fit) {
    predict(fit,newdata, re.form=~0)
  }
  bb <- bootMer(i.lmer,nsim=1000,FUN=predFun,seed=123)$t
  if (Variable!='gen.composition') {
    bb=100*sin(bb)^2
  }
  
  Xmat <- model.matrix(~Year, data=newdata)
  fit=bb %*% Xmat 
  dissim.matrix = data.frame(Mean=apply(fit,2,mean),lower=apply(fit,2,quantile, p=0.025),upper=apply(fit,2,quantile, p=0.975) )
  pred = data.frame(fit=apply(bb,2,mean),lower=apply(bb,2,quantile, p=0.025),upper=apply(bb,2,quantile, p=0.975))
  
  newdata = cbind(newdata, pred)
  newdata = newdata %>% mutate(Date=as.Date(paste0(Year,'-06-01')))
  if (Variable=="Change.score.linear.mean") {
    newdata=newdata %>% filter(year(Date)>2006)
  } else {newdata=newdata
  }
  data.i = data.i %>% mutate(Date=as.Date(lubridate:::date_decimal(Dt.num)))
  
   remove.missing=missing_span_depth(data.i, iVar)
   eval(parse(text=paste0('data.i = data.i %>% filter(!is.na(',iVar,'))')))
 
  s=data.frame(Date=CoralDateRanges[[1]], subLabels=subLabels)
  
  p <- ggplot(newdata, aes(y=fit, x=Date)) +
    xlim(CoralDateRanges$min, CoralDateRanges$max)+
    geom_ribbon(aes(ymin=lower, ymax=upper, x=Date), fill='blue',alpha=0.3)+
    geom_line(aes(y=fit, x=Date), color='blue')+
    eval(parse(text=paste0('geom_line(data=data.i, aes(y=',iVar,', x=Date,linetype=factor(DEPTH),group=interaction(ReefDepth)), alpha=0.3)')))+  
    scale_linetype_manual(values=c('solid','dashed'), guide=FALSE)+
    scale_x_date('', limits=c(CoralDateRanges[[1]],CoralDateRanges[[2]]))+
    scale_y_continuous(VariableName)+
    theme_classic(10)+
    theme(strip.background=element_blank(),
          strip.text=element_blank(),
          axis.title.x=element_blank(),
          plot.margin=unit(c(0.2,0.5,0,1),"lines"),
          panel.spacing=unit(c(0),"lines"),
          axis.title.y=element_text(vjust=1, size=rel(1.25)),
          axis.line.x=element_line(),axis.line.y=element_line()
    )
  #if(nrow(remove.missing)>0) eval(parse(text=paste0('p <- p + geom_line(data=remove.missing, aes(y=',iVar,'1, x=Date,linetype=factor(P_CODE),group=interaction(ReefDepth)), alpha=0.3, linetype="dashed")')))
  if(!is.null(region)) {
    if(region==TRUE) {
      p = p + facet_wrap(~Region, scales='free', nrow=1)
      nFac <- length(unique(dat$Region))
      if(subLabels==TRUE) {
        subLabels=paste0(letters[1:nFac],')')
        s <- data.frame(Region=unique(dat$Region), Date=CoralDateRanges[[1]], subLabels=subLabels)
      }
    }
  }
  if(!is.null(Subregion)) {
    if(Subregion==TRUE) {
      p = p + facet_wrap(~subregion, scales='free', nrow=1)
      nFac <- length(unique(dat$subregion))
      if(subLabels==TRUE) {
        subLabels=paste0(letters[1:nFac],')')
        s <- data.frame(s=unique(dat$subregion), Date=CoralDateRanges[[1]], subLabels=subLabels)
      }
    }
  }
  if(subLabels!=FALSE) {
    p = p + geom_text(data=s, aes(y=Inf, x=Date, label=subLabels),hjust=-0.2, vjust=1)
  }
  ## It only makes sense to have a title on a single figure, so
  ## region and subregion must not be both TRUE or NULL
  if (title==TRUE) p <- p+ggtitle(ti)
  if (class(ytitle)=='expression') {
    p<-p+scale_y_continuous(ytitle)
  } else if(ytitle!=FALSE) {
    if(ytitle==TRUE) ytitle='Mean'
    p<-p+scale_y_continuous(ytitle)
  }
  newdata=newdata %>% cbind(Var=as.character(VariableName),dat %>% dplyr::select(Region, subregion) %>% distinct())
  #list(g=p, m=data.means, p=newdata, dissim.matrix=dissim.matrix)
  list(g=p, p=newdata)
}
#####################
#lmer model including random reef and depth for point
#####################
MMP_analysisPointCoverLMER_depth <- function(Variable, VariableName,PointPool,region=NULL, Subregion=NULL, title=FALSE, ytitle=FALSE, guide=FALSE, subLabels=FALSE) {
  load(file='data/processed/allCoralData.RData')
  load(file='data/primary/CoralDateRanges.RData')
  require(lme4)
  require(emmeans)
  dat = allCoralData %>%filter(Year>=year(CoralDateRanges$min))
  nFac=1
  if(!is.null(region)) {
    dat <- dat %>% filter(Region==region)
    ti=region
  }
  if(!is.null(Subregion)) {
    dat <- dat %>% filter(subregion==Subregion)
    ti=Subregion
  }
  
  iVar=Variable
  iPool=PointPool
  
  data.i<-droplevels(dat) %>%
    arrange(ReefDepth,Dt.num) %>%
    mutate(ivar=eval(parse(text=iVar)),
           total = eval(parse(text=iPool)))
 
  ##fit the linear mixed effect model
  #newdata <- NULL

  data.ii <- data.i %>%
    filter(!is.na(ivar)) %>%
    droplevels

  #eval(parse(text=paste0('data.means = data.ii %>% group_by(Region,subregion,Year) %>% summarize(',Variable,'=mean(',iVar,',na.rm=TRUE))')))

  i.data = data.ii %>% mutate(Year=relevel(factor(Year),ref="2010"))
  
    i.lmer <-glmer(cbind(ivar,total-ivar) ~ Year + (1|ReefDepth),
                   data=i.data, family=binomial,nAGQ=9)

  newdata<-   emmeans(i.lmer, ~Year, type='response') %>%
     as.data.frame %>%
    mutate(Date=as.Date(paste0(Year,'-06-01')),
                        prob=prob*100,
                        asymp.LCL=asymp.LCL*100,
                        asymp.UCL=asymp.UCL*100)
  
  data.i = data.i %>% mutate(Date=as.Date(lubridate:::date_decimal(Dt.num)))
  
  remove.missing=missing_span_depth(data.i, iVar)
  eval(parse(text=paste0('data.i = data.i %>% filter(!is.na(',iVar,'))')))
  
  s=data.frame(Date=CoralDateRanges[[1]], subLabels=subLabels)    
  
 p<- ggplot(newdata, aes(y=prob, x=Date)) +
    xlim(CoralDateRanges$min, CoralDateRanges$max)+
    geom_ribbon(aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='blue', alpha=0.3) +
    geom_line(color='blue')+
    geom_line(data=data.i, aes(y=(ivar/total)*100, x=Date,linetype=factor(DEPTH),group=interaction(ReefDepth)), alpha=0.3)+
    scale_linetype_manual(values=c('solid','dashed'), guide=FALSE)+
    scale_x_date('', limits=c(CoralDateRanges[[1]],CoralDateRanges[[2]]))+
    scale_y_continuous(VariableName)+
    theme_classic(10)+
    theme(strip.background=element_blank(),
          strip.text=element_blank(),
          axis.title.x=element_blank(),
          plot.margin=unit(c(0.2,0.5,0,1),"lines"),
          panel.spacing=unit(c(0),"lines"),
          axis.title.y=element_text(vjust=1, size=rel(1.25)),
          axis.line.x=element_line(),axis.line.y=element_line()
    )
 # if(nrow(remove.missing)>0) eval(parse(text=paste0('p <- p + geom_line(data=remove.missing, aes(y=',iVar,', x=Date,linetype=factor(P_CODE),group=interaction(ReefDepth)), alpha=0.3, linetype="dashed")')))
  if(!is.null(region)) {
    if(region==TRUE) {
      p = p + facet_wrap(~Region, scales='free', nrow=1)
      nFac <- length(unique(dat$Region))
      if(subLabels==TRUE) {
        subLabels=paste0(letters[1:nFac],')')
        s <- data.frame(Region=unique(dat$Region), Date=CoralDateRanges[[1]], subLabels=subLabels)
      }
    }
  }
  if(!is.null(Subregion)) {
    if(Subregion==TRUE) {
      p = p + facet_wrap(~subregion, scales='free', nrow=1)
      nFac <- length(unique(dat$subregion))
      if(subLabels==TRUE) {
        subLabels=paste0(letters[1:nFac],')')
        s <- data.frame(s=unique(dat$subregion), Date=CoralDateRanges[[1]], subLabels=subLabels)
      }
    }
  }
  if(subLabels!=FALSE) {
    p = p + geom_text(data=s, aes(y=Inf, x=Date, label=subLabels),hjust=-0.2, vjust=1)
  }
  ## It only makes sense to have a title on a single figure, so
  ## region and subregion must not be both TRUE or NULL
  if (title==TRUE) p <- p+ggtitle(ti)
  if (class(ytitle)=='expression') {
    p<-p+scale_y_continuous(ytitle)
  } else if(ytitle!=FALSE) {
    if(ytitle==TRUE) ytitle='Mean'
    p<-p+scale_y_continuous(ytitle)
  }
  newdata=newdata %>% cbind(Var=as.character(VariableName),dat %>% dplyr::select(Region, subregion) %>% distinct())
  #list(g=p, m=data.means, p=newdata)
  list(g=p, p=newdata)
}
####################
#lmer model including random reef and depth for point with no reef lines
#####################
MMP_analysisPointCoverLMER_noReef <- function(Variable, VariableName,PointPool,region=NULL, Subregion=NULL, title=FALSE, ytitle=FALSE, guide=FALSE, subLabels=FALSE) {
  load(file='data/processed/allCoralData.RData')
  load(file='data/primary/CoralDateRanges.RData')
  require(lme4)
  require(emmeans)
  dat = allCoralData %>%filter(Year>=year(CoralDateRanges$min))
  nFac=1
  if(!is.null(region)) {
    dat <- dat %>% filter(Region==region)
    ti=region
  }
  if(!is.null(Subregion)) {
    dat <- dat %>% filter(subregion==Subregion)
    ti=Subregion
  }
  
  iVar=Variable
  iPool=PointPool
  
  data.i<-droplevels(dat) %>%
    arrange(ReefDepth,Dt.num) %>%
    mutate(ivar=eval(parse(text=iVar)),
           total = eval(parse(text=iPool)))
  
  ##fit the linear mixed effect model
  #newdata <- NULL
  
  data.ii <- data.i %>%
    filter(!is.na(ivar)) %>%
    droplevels
  
  #eval(parse(text=paste0('data.means = data.ii %>% group_by(Region,subregion,Year) %>% summarize(',Variable,'=mean(',iVar,',na.rm=TRUE))')))
  
  i.data = data.ii %>% mutate(Year=relevel(factor(Year),ref="2010"))
  
  i.lmer <-glmer(cbind(ivar,total-ivar) ~ Year + (1|ReefDepth),
                 data=i.data, family=binomial,nAGQ=9)
  
  newdata<-   emmeans(i.lmer, ~Year, type='response') %>%
    as.data.frame %>%
    mutate(Date=as.Date(paste0(Year,'-06-01')),
           prob=prob*100,
           asymp.LCL=asymp.LCL*100,
           asymp.UCL=asymp.UCL*100)
  
  data.i = data.i %>% mutate(Date=as.Date(lubridate:::date_decimal(Dt.num)))
  
  remove.missing=missing_span_depth(data.i, iVar)
  eval(parse(text=paste0('data.i = data.i %>% filter(!is.na(',iVar,'))')))
  
  s=data.frame(Date=CoralDateRanges[[1]], subLabels=subLabels)    
  
  p<- ggplot(newdata, aes(y=prob, x=Date)) +
    xlim(CoralDateRanges$min, CoralDateRanges$max)+
    geom_ribbon(aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='blue', alpha=0.3) +
    geom_line(color='blue')+
   # geom_line(data=data.i, aes(y=ivar/total, x=Date,linetype=factor(DEPTH),group=interaction(ReefDepth)), alpha=0.3)+
    scale_linetype_manual(values=c('solid','dashed'), guide=FALSE)+
    scale_x_date('', limits=c(CoralDateRanges[[1]],CoralDateRanges[[2]]))+
    scale_y_continuous(VariableName)+
    theme_classic(10)+
    theme(strip.background=element_blank(),
          strip.text=element_blank(),
          axis.title.x=element_blank(),
          plot.margin=unit(c(0.2,0.5,0,1),"lines"),
          panel.spacing=unit(c(0),"lines"),
          axis.title.y=element_text(vjust=1, size=rel(1.25)),
          axis.line.x=element_line(),axis.line.y=element_line()
    )
  #if(nrow(remove.missing)>0) eval(parse(text=paste0('p <- p + geom_line(data=remove.missing, aes(y=',iVar,', x=Date,linetype=factor(P_CODE),group=interaction(ReefDepth)), alpha=0.3, linetype="dashed")')))
  if(!is.null(region)) {
    if(region==TRUE) {
      p = p + facet_wrap(~Region, scales='free', nrow=1)
      nFac <- length(unique(dat$Region))
      if(subLabels==TRUE) {
        subLabels=paste0(letters[1:nFac],')')
        s <- data.frame(Region=unique(dat$Region), Date=CoralDateRanges[[1]], subLabels=subLabels)
      }
    }
  }
  if(!is.null(Subregion)) {
    if(Subregion==TRUE) {
      p = p + facet_wrap(~subregion, scales='free', nrow=1)
      nFac <- length(unique(dat$subregion))
      if(subLabels==TRUE) {
        subLabels=paste0(letters[1:nFac],')')
        s <- data.frame(s=unique(dat$subregion), Date=CoralDateRanges[[1]], subLabels=subLabels)
      }
    }
  }
  if(subLabels!=FALSE) {
    p = p + geom_text(data=s, aes(y=Inf, x=Date, label=subLabels),hjust=-0.2, vjust=1)
  }
  ## It only makes sense to have a title on a single figure, so
  ## region and subregion must not be both TRUE or NULL
  if (title==TRUE) p <- p+ggtitle(ti)
  if (class(ytitle)=='expression') {
    p<-p+scale_y_continuous(ytitle)
  } else if(ytitle!=FALSE) {
    if(ytitle==TRUE) ytitle='Mean'
    p<-p+scale_y_continuous(ytitle)
  }
  # newdata=newdata %>% cbind(Var=as.character(VariableName),dat %>% dplyr::select(Region, subregion) %>% distinct())
  # list(g=p, m=data.means, p=newdata)
  # list(g=p, p=newdata)
  list(g=p)
}

#############################################
#lmer model including random reef and depth for juvenile counts
#####################
MMP_analysisPoissonGLMER_depth <- function(Variable, VariableName,Offset,region=NULL, Subregion=NULL, title=FALSE, ytitle=FALSE, guide=FALSE, subLabels=FALSE) {
  load(file='data/processed/juv.counts.RData')
  load(file='data/primary/CoralDateRanges.RData')
  require(lme4)
  require(emmeans)
  dat = juv.counts %>%filter(Year>=year(CoralDateRanges$min))
  nFac=1
  if(!is.null(region)) {
    dat <- dat %>% filter(Region==region)
    ti=region
  }
  if(!is.null(Subregion)) {
    dat <- dat %>% filter(subregion==Subregion)
    ti=Subregion
  }
  
  iVar=Variable
  iOff=Offset
  
  data.i<-droplevels(dat) %>%
    arrange(ReefDepth,Dt.num) %>%
    mutate(ivar=eval(parse(text=iVar)),
           ioff = eval(parse(text=iOff)))
  
  ##fit the linear mixed effect model
 # newdata <- NULL
  
  data.ii <- data.i %>%
    filter(!is.na(ivar)) %>%
    droplevels
  
  #eval(parse(text=paste0('data.means = data.ii %>% group_by(Region,subregion,Year) %>% summarize(',Variable,'=mean(',iVar,',na.rm=TRUE))')))
  
  i.data = data.ii %>% mutate(Year=relevel(factor(Year),ref="2010"))
  
  i.lmer <-glmer(ivar ~ Year + offset(log(ioff)) +(1|ReefDepth/SITE_NO),
                 data=i.data,
                 family=poisson,nAGQ=1)
  
  newdata<- emmeans(i.lmer, ~Year, type='response', offset=0) %>%
    as.data.frame %>%
    mutate(Date=as.Date(paste0(Year,'-06-01')))
 
  #density data to plot 
  data.reef = data.i %>% mutate(Date=as.Date(lubridate:::date_decimal(Dt.num))) %>%
    mutate(juv.density=ivar/ioff) %>%
    group_by(Year,REEF,DEPTH,ReefDepth) %>%
    summarise(juv.d=mean(juv.density),
              Date=min(Date))
    
   s=data.frame(Date=CoralDateRanges[[1]], subLabels=subLabels)    
  
  p<- ggplot(newdata, aes(y=rate, x=Date)) +
    xlim(CoralDateRanges$min, CoralDateRanges$max)+
    geom_ribbon(aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='blue', alpha=0.3) +
    geom_line(color='blue')+
    geom_line(data=data.reef, aes(y=juv.d, x=Date,linetype=factor(DEPTH),group=interaction(ReefDepth)), alpha=0.3)+
    scale_linetype_manual(values=c('solid','dashed'), guide=FALSE)+
    scale_x_date('', limits=c(CoralDateRanges[[1]],CoralDateRanges[[2]]))+
    scale_y_continuous(VariableName)+
    theme_classic(10)+
    theme(strip.background=element_blank(),
          strip.text=element_blank(),
          axis.title.x=element_blank(),
          plot.margin=unit(c(0.2,0.5,0,1),"lines"),
          panel.spacing=unit(c(0),"lines"),
          axis.title.y=element_text(vjust=1, size=rel(1.25)),
          axis.line.x=element_line(),axis.line.y=element_line()
    )
  #if(nrow(remove.missing)>0) eval(parse(text=paste0('p <- p + geom_line(data=remove.missing, aes(y=',iVar,', x=Date,linetype=factor(P_CODE),group=interaction(ReefDepth)), alpha=0.3, linetype="dashed")')))
  if(!is.null(region)) {
    if(region==TRUE) {
      p = p + facet_wrap(~Region, scales='free', nrow=1)
      nFac <- length(unique(dat$Region))
      if(subLabels==TRUE) {
        subLabels=paste0(letters[1:nFac],')')
        s <- data.frame(Region=unique(dat$Region), Date=CoralDateRanges[[1]], subLabels=subLabels)
      }
    }
  }
  if(!is.null(Subregion)) {
    if(Subregion==TRUE) {
      p = p + facet_wrap(~subregion, scales='free', nrow=1)
      nFac <- length(unique(dat$subregion))
      if(subLabels==TRUE) {
        subLabels=paste0(letters[1:nFac],')')
        s <- data.frame(s=unique(dat$subregion), Date=CoralDateRanges[[1]], subLabels=subLabels)
      }
    }
  }
  if(subLabels!=FALSE) {
    p = p + geom_text(data=s, aes(y=Inf, x=Date, label=subLabels),hjust=-0.2, vjust=1)
  }
  ## It only makes sense to have a title on a single figure, so
  ## region and subregion must not be both TRUE or NULL
  if (title==TRUE) p <- p+ggtitle(ti)
  if (class(ytitle)=='expression') {
    p<-p+scale_y_continuous(ytitle)
  } else if(ytitle!=FALSE) {
    if(ytitle==TRUE) ytitle='Mean'
    p<-p+scale_y_continuous(ytitle)
  }
   newdata=newdata %>% cbind(Var=as.character(VariableName),dat %>% dplyr::select(Region, subregion) %>% distinct())
   #list(g=p, m=data.means, p=newdata)
   list(g=p,  p=newdata)
}

  ##############################

# no reef lines
MMP_analysisSedCoralLMER_noReef <- function(Variable, VariableName,region=NULL, Subregion=NULL, title=FALSE, ytitle=FALSE, guide=FALSE, subLabels=FALSE) {
  load(file='data/processed/allCoralData.RData')
  load(file='data/primary/CoralDateRanges.RData')
  require(lme4)
  dat = allCoralData %>%filter(Year>=year(CoralDateRanges$min))
  nFac=1
  if(!is.null(region)) {
    dat <- dat %>% filter(Region==region)
    ti=region
  }
  if(!is.null(Subregion)) {
    dat <- dat %>% filter(subregion==Subregion)
    ti=Subregion
  }
  data.i<-droplevels(dat)
  data.i <- arrange(data.i, ReefDepth,Dt.num)
  iVar=Variable
  ##fit the linear mixed effect model
  newdata <- NULL
  err <- 0
  
  eval(parse(text=paste0('data.ii <- droplevels(filter(data.i, !is.na(',iVar,')))')))
  eval(parse(text=paste0('data.means = data.ii %>% group_by(Region,subregion,Year) %>% summarize(',Variable,'=mean(',iVar,',na.rm=TRUE))')))
  if (Variable=='gen.composition') {
    i.data = data.ii %>% filter_(.dots=paste('!is.na(',iVar,')')) %>% mutate(Year=factor(Year)) %>%
      mutate_(.dots=setNames(paste(iVar),'Resp'))
    tryCatch(
      #i.lmer <- lmer(Resp ~ Year + (Year|ReefDepth), data=i.data),
      i.lmer <- lmer(Resp ~ Year + (1|ReefDepth), data=i.data),
      error=function(x) err<<-1
    )
  } else {
    i.data = data.ii %>% filter_(.dots=paste('!is.na(',iVar,')')) %>% mutate(Year=factor(Year)) %>%
      mutate_(.dots=setNames(paste0('asin(sqrt(',iVar,'/100))'),'Resp'))
    tryCatch(
      #i.lmer <-lmer(Resp ~ Year + (Year|ReefDepth), data=i.data),
      i.lmer <-lmer(Resp ~ Year + (1|ReefDepth), data=i.data),
      error=function(x) err<<-1
    )
  }
  newdata=data.frame(Year=levels(i.lmer@frame$Year))
  
  predFun <- function(fit) {
    predict(fit,newdata, re.form=~0)
  }
  bb <- bootMer(i.lmer,nsim=1000,FUN=predFun,seed=123)$t
  if (Variable!='gen.composition') {
    bb=100*sin(bb)^2
  }
  
  Xmat <- model.matrix(~Year, data=newdata)
  fit=bb %*% Xmat 
  dissim.matrix = data.frame(Mean=apply(fit,2,mean),lower=apply(fit,2,quantile, p=0.025),upper=apply(fit,2,quantile, p=0.975) )
  pred = data.frame(fit=apply(bb,2,mean),lower=apply(bb,2,quantile, p=0.025),upper=apply(bb,2,quantile, p=0.975))
  
  newdata = cbind(newdata, pred)
  newdata = newdata %>% mutate(Date=as.Date(paste0(Year,'-06-01')))
  data.i = data.i %>% mutate(Date=as.Date(lubridate:::date_decimal(Dt.num)))
  
  
  remove.missing=missing_span_depth(data.i, iVar)
  eval(parse(text=paste0('data.i = data.i %>% filter(!is.na(',iVar,'))')))
  load('data/primary/CoralDateRanges.RData')
  s=data.frame(Date=CoralDateRanges[[1]], subLabels=subLabels)
  
  p <- ggplot(newdata, aes(y=fit, x=Date)) +
    xlim(CoralDateRanges$min, CoralDateRanges$max)+
    geom_ribbon(aes(ymin=lower, ymax=upper, x=Date), fill='blue',alpha=0.3)+
    geom_line(aes(y=fit, x=Date), color='blue')+
    scale_x_date('', limits=c(CoralDateRanges[[1]],CoralDateRanges[[2]]))+
    scale_y_continuous(VariableName, limits=c(0,60))+
    theme_classic(8)+
    theme(strip.background=element_blank(),
          strip.text=element_blank(),
          axis.title.x=element_blank(),
          plot.margin=unit(c(0.2,0.5,0,1),"lines"),
          panel.spacing=unit(c(0),"lines"),
          axis.title.y=element_text(vjust=1, size=rel(1)),
          axis.line.x=element_line(),axis.line.y=element_line()
    )
  
  if(subLabels!=FALSE) {
    p = p + geom_text(data=s, aes(y=Inf, x=Date, label=subLabels),hjust=-0.1, vjust=1)
  }
  ## It only makes sense to have a title on a single figure, so
  ## region and subregion must not be both TRUE or NULL
  if (title==TRUE) p <- p+ggtitle(ti)
  if (class(ytitle)=='expression') {
    p<-p+scale_y_continuous(ytitle)
  } else if(ytitle!=FALSE) {
    if(ytitle==TRUE) ytitle='Mean'
    p<-p+scale_y_continuous(ytitle)
  }
  
  newdata=newdata %>% cbind(Var=as.character(Variable),dat %>% dplyr::select(Region) %>% distinct())
  list(g=p, m=data.means, p=newdata, dissim.matrix=dissim.matrix)
}


#This is a clever function for isolating a missing span..
missing_span_depth <- function(data, iVar) {
  data %>% group_by(P_CODE, Region, subregion, ReefDepth) %>% do({
    wch<-which(is.na(.[,iVar]))
    wch<-wch[wch>1 & wch<nrow(.)]
    if (length(wch)==0) w<-data.frame(NULL)
    else {
      w<-.[c(wch-1,wch+1),]
      eval(parse(text=paste0('w<-filter(w,!is.na(',iVar,'))')))
      w<-data.frame(w)
    }
    print(w)
    w        
  })%>% as.data.frame()
}

missing_span <- function(data, iVar) {
  data %>% group_by(P_CODE, Region, subregion, REEF) %>% do({
    wch<-which(is.na(.[,iVar]))
    wch<-wch[wch>1 & wch<nrow(.)]
    if (length(wch)==0) w<-data.frame(NULL)
    else {
      w<-.[c(wch-1,wch+1),]
      eval(parse(text=paste0('w<-filter(w,!is.na(',iVar,'))')))
      w<-data.frame(w)
    }
    #print(w)
    w        
  })%>% as.data.frame()
}

RC_boot_aggregate = function(boot=NULL, df, size=10, seed=set.seed(123),
                             grouping_cols=colnames(df)[!colnames(df) %in% c('Score','Grade')],
                             over='') {
  if(is.null(boot)) boot=df
  ## Bootstrap
  set.seed(seed)
  x2=boot %>%
    group_by(.dots=grouping_cols) %>%
    RC_bootstrap(., size=size, group=over) %>% ungroup
  
  ## Confidence intervals
  x3=x2 %>%
    group_by(.dots=grouping_cols) %>%
    summarize(Lower=quantile(Score, p=0.025, na.rm=TRUE),
              Upper=quantile(Score, p=0.975, na.rm=TRUE),
              Boot.Mean=mean(Score, na.rm=TRUE)
    ) %>%
    full_join(df) %>%
    ungroup
  list(dist=x2, sum=x3)
}

RC_bootstrap = function(df, size=1000, group='Measure') {
  if (group=='' | is.null(group)) {
    grp=NULL
  } else if (is.list(group)) {
    grp=group
  } else {
    grp=group
  }
  
  df %>% do({
    x=.
    xx=replicate(n=size, tapply(x$Score, x[[grp]], function(x) sample(x, size=1)))
    w = tapply(x$Weight, x[[grp]], mean)
    xxx = apply(xx, 2, function(x) weighted.mean(x, w = w, na.rm=TRUE))
    data.frame(Score=xxx)
  })
}




RC_aggregate = function(df, grouping_cols=colnames(df)[!colnames(df) %in% c('Score','Grade')]) {
  df %>%
    group_by(.dots=grouping_cols) %>%
    summarize(Score=weighted.mean(Score, w=Weight, na.rm=TRUE)) %>% mutate(Grade=MMP_generateGrades(Score))
}

RC_boot_accumulate = function(boot=NULL, df, size=10, seed=set.seed(123),
                              grouping_cols=colnames(df)[!colnames(df) %in% c('Score','Grade')]) {
  if(is.null(boot)) boot=df
  ## Bootstrap
  set.seed(seed)
  x2 = boot %>% group_by(.dots=grouping_cols) %>%
    do({
      x = .
      xx=replicate(n=size, sample(x$Score, replace=TRUE))
      data.frame(Score=colMeans(xx))
    })
  x3 = x2 %>% 
    group_by(.dots=grouping_cols) %>%
    summarize(Lower=quantile(Score, p=0.025),
              Upper=quantile(Score, p=0.975),
              Boot.mean=mean(Score)) %>%
    full_join(df) %>% ungroup
  list(dist=x2, sum=x3)
}
