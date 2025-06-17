source('MMP_coralFunctions.R')
#loadPackages()
library(tidyverse)
size =100
seed=123
bootstrap=TRUE

#############################################################################
# This scipt uses reported data from 2020 to provide index scores with and without the 
# Composition indicator
#
#
################################################################################

load(file="../../../../Visit 20 2024/code/data/processed/final/coral.RData")

coral <- coral  %>% rename(Region=NRM_REGION, Year=yr) %>%
    mutate(subregion=recode_factor(subregion,
                                   'Daintree'='Barron Daintree',
                                   'Johnstone'='Johnstone Russell Mulgrave',
                                   'Tully'='Tully Herbert',
                                   'Burdekin'='Burdekin',
                                   'Proserpine'='Mackay Whitsunday',
                                   'Fitzroy'='Fitzroy'))

## ---- Bootstrapp Aggregations---------------------------------------------------------------------------------------------------------------------------------------------------

## Prepare for bootstrapping aggregation
coral.melt <- coral %>% filter(Year>2005) %>%
    dplyr::select(P_CODE,Region,subregion,REEF,DEPTH,Year,Date,Change.score.linear.mean,ma.score,CoralCover.score,juv.score,composition.score) %>%
    dplyr::rename('Macroalgal cover'=ma.score,
                  'Coral cover'=CoralCover.score,
                  'Juvenile density'=juv.score,
                  'Cover change'=Change.score.linear.mean,
                  'Coral composition'=composition.score) %>%
  gather(key=Indicator,value=Score,-P_CODE,-Region,-subregion,-REEF,-DEPTH,-Year,-Date) %>%
  filter(!is.na(Score)) %>%
  mutate(ReefDepth=paste0(REEF,DEPTH)) %>%
  ungroup

coral.melt.no.comp <- coral %>% filter(Year>2005) %>%
  dplyr::select(P_CODE,Region,subregion,REEF,DEPTH,Year,Date,Change.score.linear.mean,ma.score,CoralCover.score,juv.score) %>%
  dplyr::rename('Macroalgal cover'=ma.score,
                'Coral cover'=CoralCover.score,
                'Juvenile density'=juv.score,
                'Cover change'=Change.score.linear.mean) %>%
  gather(key=Indicator,value=Score,-P_CODE,-Region,-subregion,-REEF,-DEPTH,-Year,-Date) %>%
  filter(!is.na(Score)) %>%
  mutate(ReefDepth=paste0(REEF,DEPTH)) %>%
  ungroup
## ---- Output data to enable sub regional pairwise comparison of changes in index and indicator scores between focal years
indicator.reef.depth<-coral.melt
#save(indicator.reef.depth, file='data/processed/indicator.reef.depth.RData')

indicator.reef.depth.no.comp<-coral.melt.no.comp
#save(indicator.reef.depth.no.comp, file='data/processed/indicator.reef.depth.no.comp.RData')

###Index/ReefDepth================================================================

#arithmetic mean.
index.reef.depth = indicator.reef.depth %>% mutate(Weight=1) %>% # alternatively indicator level weights could be applied by merging in a dataframe containing weights. 
  RC_aggregate(grouping_cols=c('Year','Region','subregion', 'REEF','DEPTH','ReefDepth')) %>%
  ungroup
# bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.index.reef.depth = indicator.reef.depth %>% mutate(Weight=1) %>% 
    RC_boot_accumulate(df=index.reef.depth, size=size, seed=seed,
                       grouping_cols=c('Year','Region','subregion', 'REEF','DEPTH', 'ReefDepth')) 
  save(boot.index.reef.depth, file='../data/processed/boot.index.reef.depth.RData')
}
#output index at reef depth for pairwise comparisons
save(index.reef.depth, file='../data/processed/index.reef.depth.RData')

index.reef.depth.no.comp = indicator.reef.depth.no.comp %>% mutate(Weight=1) %>% # alternatively indicator level weights could be applied by merging in a dataframe containing weights. 
  RC_aggregate(grouping_cols=c('Year','Region','subregion', 'REEF','DEPTH','ReefDepth')) %>%
  ungroup
# bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.index.reef.depth.no.comp = indicator.reef.depth.no.comp %>% mutate(Weight=1) %>% 
    RC_boot_accumulate(df=index.reef.depth.no.comp, size=size, seed=seed,
                       grouping_cols=c('Year','Region','subregion', 'REEF','DEPTH', 'ReefDepth')) 
  save(boot.index.reef.depth.no.comp, file='../data/processed/boot.index.reef.depth.no.comp.RData')
}
#output index at reef depth for pairwise comparisons with no composition
save(index.reef.depth.no.comp, file='../data/processed/index.reef.depth.no.comp.RData')

### Indicator/Subregion/Depth ==============================================================

# arithmetic mean
indicator.subregion.depth = indicator.reef.depth %>% mutate(Weight=1) %>%
  RC_aggregate(grouping_cols=c('Year','Region','subregion', 'DEPTH','Indicator')) %>%
  ungroup

#bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.indicator.subregion.depth = indicator.reef.depth %>% mutate(Weight=1) %>% 
    RC_boot_accumulate(df=indicator.subregion.depth, size=size, seed=seed,
                       grouping_cols=c('Year','Region','subregion', 'DEPTH', 'Indicator')) 
  save(boot.indicator.subregion.depth, file='../data/processed/boot.indicator.subregion.depth.RData')
}

indicator.subregion.depth.no.comp = indicator.reef.depth.no.comp %>% mutate(Weight=1) %>%
  RC_aggregate(grouping_cols=c('Year','Region','subregion', 'DEPTH','Indicator')) %>%
  ungroup

#bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.indicator.subregion.depth.no.comp = indicator.reef.depth.no.comp %>% mutate(Weight=1) %>% 
    RC_boot_accumulate(df=indicator.subregion.depth.no.comp, size=size, seed=seed,
                       grouping_cols=c('Year','Region','subregion', 'DEPTH', 'Indicator')) 
  save(boot.indicator.subregion.depth.no.comp, file='../data/processed/boot.indicator.subregion.depth.no.comp.RData')
}


## Index/Subregion/Depth ==============================================================
# arithmetic mean
index.subregion.depth = indicator.subregion.depth %>% mutate(Weight=1) %>% 
  RC_aggregate(grouping_cols=c('Year','Region','subregion', 'DEPTH')) %>%
  ungroup

#bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.index.subregion.depth = boot.indicator.subregion.depth[['dist']] %>% mutate(Weight=1) %>% 
    RC_boot_aggregate(df=index.subregion.depth, size=size, seed=seed,
                      grouping_cols=c('Year','Region','subregion', 'DEPTH'),
                      over='Indicator') 
  save(boot.index.subregion.depth, file='../data/processed/boot.index.subregion.depth.RData')
}

index.subregion.depth.no.comp = indicator.subregion.depth.no.comp %>% mutate(Weight=1) %>% # alternatively indicator level weights could be applied by merging in a dataframe containing weights. 
  RC_aggregate(grouping_cols=c('Year','Region','subregion', 'DEPTH')) %>%
  ungroup

#bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.index.subregion.depth.no.comp = boot.indicator.subregion.depth.no.comp[['dist']] %>% mutate(Weight=1) %>% 
    RC_boot_aggregate(df=index.subregion.depth.no.comp, size=size, seed=seed,
                      grouping_cols=c('Year','Region','subregion', 'DEPTH'),
                      over='Indicator') 
  save(boot.index.subregion.depth.no.comp, file='../data/processed/boot.index.subregion.depth.no.comp.RData')
}


### Indicator/Subregion ==============================================================

# arithmetic mean
indicator.subregion = indicator.reef.depth %>% mutate(Weight=1) %>%
  RC_aggregate(grouping_cols=c('Year','Region','subregion','Indicator')) %>%
  ungroup

#bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.indicator.subregion = indicator.reef.depth %>% mutate(Weight=1) %>% 
    RC_boot_accumulate(df=indicator.subregion, size=size, seed=seed,
                       grouping_cols=c('Year','Region','subregion', 'Indicator')) 
  save(boot.indicator.subregion, file='../data/processed/boot.indicator.subregion.RData')
}
write.csv(boot.indicator.subregion[['sum']] %>% filter(Year==max(Year)), file='../output/tables/subregion.indicator.scores.current.year.csv')

indicator.subregion.no.comp = indicator.reef.depth.no.comp %>% mutate(Weight=1) %>%
  RC_aggregate(grouping_cols=c('Year','Region','subregion','Indicator')) %>%
  ungroup

#bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.indicator.subregion.no.comp = indicator.reef.depth.no.comp %>% mutate(Weight=1) %>% 
    RC_boot_accumulate(df=indicator.subregion.no.comp, size=size, seed=seed,
                       grouping_cols=c('Year','Region','subregion', 'Indicator')) 
  save(boot.indicator.subregion.no.comp, file='../data/processed/boot.indicator.subregion.no.comp.RData')
}
write.csv(boot.indicator.subregion.no.comp[['sum']] %>% filter(Year==max(Year)), file='../output/tables/subregion.indicator.scores.current.year.no.comp.csv')


## Index/Subregion ==============================================================
# arithmetic mean
index.subregion = indicator.subregion %>% mutate(Weight=1) %>% # alternatively indicator level weights could be applied by merging in a dataframe containing weights. 
  RC_aggregate(grouping_cols=c('Year','Region','subregion')) %>%
  ungroup

#bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.index.subregion = boot.indicator.subregion[['dist']] %>% mutate(Weight=1) %>% 
    RC_boot_aggregate(df=index.subregion, size=size, seed=seed,
                      grouping_cols=c('Year','Region','subregion'),
                      over='Indicator') 
  save(boot.index.subregion, file='../data/processed/boot.index.subregion.RData')
}

# Output current year subregional index scores for appendix table
write.csv(boot.index.subregion[['sum']] %>% filter(Year==max(Year)), file='../output/tables/subregion.index.scores.current.year.csv')

# arithmetic mean
index.subregion.no.comp = indicator.subregion.no.comp %>% mutate(Weight=1) %>% # alternatively indicator level weights could be applied by merging in a dataframe containing weights. 
  RC_aggregate(grouping_cols=c('Year','Region','subregion')) %>%
  ungroup

#bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.index.subregion.no.comp = boot.indicator.subregion.no.comp[['dist']] %>% mutate(Weight=1) %>% 
    RC_boot_aggregate(df=index.subregion.no.comp, size=size, seed=seed,
                      grouping_cols=c('Year','Region','subregion'),
                      over='Indicator') 
  save(boot.index.subregion.no.comp, file='../data/processed/boot.index.subregion.no.comp.RData')
}

# Output current year subregional index scores for appendix table
write.csv(boot.index.subregion.no.comp[['sum']] %>% filter(Year==max(Year)), file='../output/tables/subregion.index.scores.current.year.no.compcsv')



### Indicator/Region ==============================================================

# arithmetic mean
indicator.region = indicator.reef.depth %>% mutate(Weight=1) %>%
  RC_aggregate(grouping_cols=c('Year','Region','Indicator')) %>%
  ungroup

#bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.indicator.region = indicator.reef.depth %>% mutate(Weight=1) %>% 
    RC_boot_accumulate(df=indicator.region, size=size, seed=seed,
                       grouping_cols=c('Year','Region','Indicator')) 
  save(boot.indicator.region, file='../data/processed/boot.indicator.region.RData')
}

indicator.region.no.comp = indicator.reef.depth.no.comp %>% mutate(Weight=1) %>%
  RC_aggregate(grouping_cols=c('Year','Region','Indicator')) %>%
  ungroup

#bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.indicator.region.no.comp = indicator.reef.depth.no.comp %>% mutate(Weight=1) %>% 
    RC_boot_accumulate(df=indicator.region.no.comp, size=size, seed=seed,
                       grouping_cols=c('Year','Region','Indicator')) 
  save(boot.indicator.region.no.comp, file='../data/processed/boot.indicator.region.no.comp.RData')
}

## Index/Region ==============================================================
# arithmetic mean
index.region = indicator.region %>% mutate(Weight=1) %>% 
  RC_aggregate(grouping_cols=c('Year','Region')) %>%
  ungroup

#bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.index.region = boot.indicator.region[['dist']] %>% mutate(Weight=1) %>% 
    RC_boot_aggregate(df=index.region, size=size, seed=seed,
                      grouping_cols=c('Year','Region'),
                      over='Indicator') 
  save(boot.index.region, file='../data/processed/boot.index.region.RData')
}


index.region.no.comp = indicator.region.no.comp %>% mutate(Weight=1) %>% 
  RC_aggregate(grouping_cols=c('Year','Region')) %>%
  ungroup

#bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.index.region.no.comp = boot.indicator.region.no.comp[['dist']] %>% mutate(Weight=1) %>% 
    RC_boot_aggregate(df=index.region.no.comp, size=size, seed=seed,
                      grouping_cols=c('Year','Region'),
                      over='Indicator') 
  save(boot.index.region.no.comp, file='../data/processed/boot.index.region.no.comp.RData')
}


## Indicator/GBR ========================================================================
region.weights=rbind(data.frame(Region="Wet Tropics", Weight=0.209),
                     data.frame(Region="Burdekin", Weight=0.092),
                     data.frame(Region="Mackay Whitsunday", Weight=0.381),
                     data.frame(Region="Fitzroy", Weight=0.318))

# arithmetic mean
indicator.gbr = indicator.region %>% left_join(region.weights) %>%
  RC_aggregate(grouping_cols=c('Year','Indicator')) %>%
  ungroup

#bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.indicator.gbr = boot.indicator.region[['dist']] %>% left_join(region.weights) %>% 
    RC_boot_aggregate(df=indicator.gbr, size=size, seed=seed,
                      grouping_cols=c('Year','Indicator'),
                      over='Region') 
  save(boot.indicator.gbr, file='../data/processed/boot.indicator.gbr.RData')
}
load(file='../data/processed/boot.indicator.gbr.RData')


# arithmetic mean
indicator.gbr.no.comp = indicator.region.no.comp %>% left_join(region.weights) %>%
  RC_aggregate(grouping_cols=c('Year','Indicator')) %>%
  ungroup

#bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.indicator.gbr.no.comp = boot.indicator.region.no.comp[['dist']] %>% left_join(region.weights) %>% 
    RC_boot_aggregate(df=indicator.gbr.no.comp, size=size, seed=seed,
                      grouping_cols=c('Year','Indicator'),
                      over='Region') 
  save(boot.indicator.gbr.no.comp, file='../data/processed/boot.indicator.gbr.no.comp.RData')
}


## Index/GBR ==============================================================
# arithmetic mean
index.gbr = indicator.gbr %>% mutate(Weight=1) %>% ## NOTE weights applied already at indicator level
  RC_aggregate(grouping_cols=c('Year')) %>%
  ungroup

#bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.index.gbr = boot.indicator.gbr[['dist']] %>% mutate(Weight=1,
                                                           indicators='all') %>% ## NOTE weights applied already at indicator level
    RC_boot_aggregate(df=index.gbr, size=size, seed=seed,
                      grouping_cols=c('Year','indicators'),
                      over='Indicator') 
  save(boot.index.gbr, file='../data/processed/boot.index.gbr.RData')
}


index.gbr.no.comp = indicator.gbr.no.comp %>% mutate(Weight=1) %>% ## NOTE weights applied already at indicator level
  RC_aggregate(grouping_cols=c('Year')) %>%
  ungroup

#bootstrap distribution and summary
if(bootstrap==TRUE) {
  boot.index.gbr.no.comp = boot.indicator.gbr.no.comp[['dist']] %>% mutate(Weight=1,
                                                                           indicators='no.comp') %>% ## NOTE weights applied already at indicator level
    RC_boot_aggregate(df=index.gbr.no.comp, size=size, seed=seed,
                      grouping_cols=c('Year', 'indicators'),
                      over='Indicator') 
  save(boot.index.gbr.no.comp, file='../data/processed/boot.index.gbr.no.comp.RData')
}

subregion.index.comparison<-boot.index.subregion[['sum']] %>%
  mutate(indicators="all") %>%
  rbind(boot.index.subregion.no.comp[['sum']] %>% mutate(indicators='no.comp'))

region.index.comparison<-boot.index.region[['sum']] %>%
  mutate(indicators="all") %>%
  rbind(boot.index.region.no.comp[['sum']] %>% mutate(indicators='no.comp'))


gbr.index.comparison<-boot.index.gbr[['sum']] %>% rbind(boot.index.gbr.no.comp[['sum']])

ggplot(gbr.index.comparison %>% filter(Year>2005), aes(y=Score, x=Year)) +
  geom_linerange(data=gbr.index.comparison %>% filter(indicators=='all'), aes(ymin=Lower, ymax=Upper, x=Year), size=1.1)+
  geom_linerange(data=gbr.index.comparison %>% filter(indicators=='no.comp'), aes(ymin=Lower, ymax=Upper, x=Year+0.1), size=1.1, color='blue')+
  geom_line(data=gbr.index.comparison %>% filter(indicators=='all'), size=1.5) +
  geom_point(data=gbr.index.comparison %>% filter(indicators=='all'), aes(y=Score, x=Year,fill=Grade),shape=21,size=5,color='grey') +
  geom_line(data=gbr.index.comparison %>% filter(indicators=='no.comp'),aes(y=Score, x=Year+0.1), color='blue', size=1.5) +
  geom_point(data=gbr.index.comparison %>% filter(indicators=='no.comp'), aes(y=Score, x=Year+0.1,fill=Grade),shape=21,size=5,color='grey') +
  scale_fill_manual('Grade', breaks=LETTERS[1:5], values=rev(trafficLightPalette), limits=rev(LETTERS[1:5])) +
  scale_y_continuous('Score', limits = c(0.2, 0.6)) +
  scale_x_continuous('', breaks=function(limits) pretty(limits, 5)) +
  theme_bw() +
  theme(strip.background=element_blank())

theme_regions<-theme_classic(12)+
  theme(strip.background=element_blank(),
        strip.text=element_blank(),
        axis.title.x=element_blank(),
        plot.margin=unit(c(0.2,0.5,0,1),"lines"),
        panel.spacing=unit(c(1),"lines"),
        axis.title.y=element_text(vjust=1, size=rel(1.25)),
        axis.line.x=element_line(),axis.line.y=element_line())


gbr.index.comp<-MMP_coralWorm.comp(gbr.index.comparison, title=FALSE, strip=TRUE,
                                      ytitle=expression(Coral~index),subLabels=FALSE,nrow=1, line.size=0.75, point.size=5)+
  theme_regions

ggsave(filename='output/figures/index.with.without.composition.gbr.png', gbr.index.comp, width=7, height=7, dpi=300)

region.index.comp<-MMP_coralWorm.comp(region.index.comparison, region=TRUE, title=FALSE, strip=TRUE,
                                      ytitle=expression(Coral~index),subLabels=FALSE,line.size=0.75, point.size=3, nrow=2)

ggsave(filename='output/figures/index.with.without.composition.region.png', region.index.comp, width=7, height=7, dpi=300)+
  theme_regions


subregion.index.comp <- MMP_coralWorm.comp(subregion.index.comparison, Subregion=TRUE,title=FALSE,strip=TRUE,
                                             ytitle=expression(Coral~index),subLabels=FALSE,
                                             nrow=3, line.size=0.75, point.size=5)

ggsave(filename='output/figures/index.with.without.composition.subregion.png', subregion.index.comp, width=15, height=22.5, dpi=300)
