
library(tidyverse)
spatial_lookup <- get(load(file = "../data/processed/spatial_lookup.RData"))

depth<-spatial_lookup |> 
  dplyr::select(REEF.d, DEPTH.f) |> 
  unique()

j<-get(load(file="../data/modelled/JU__scores_reef_year.RData"))

ju<-j |> 
  dplyr::select(Summary)  |> 
  unnest(Summary)  |> 
  dplyr::select(-c(variable, BIOREGION.agg)) |> 
  mutate(Indicator="Juvenile",
         Metric=case_match(Metric, "Total"~"distance.metric",
                           "Acropora"~"consequence.metric"))
  
m<-get(load(file="../data/modelled/MA__scores_reef_year.RData")) |> 
  ungroup()

ma<-m |> dplyr::select(Summary)  |> 
  unnest(Summary) |> 
  filter(variable==".value") |> 
  dplyr::select(-c(variable, BIOREGION.agg)) |> 
  mutate(Indicator="Macroalgae")|> 
  left_join(depth)
  
c<-get(load(file="../data/modelled/CC__scores_reef_year.RData")) |> 
  ungroup()

cc<-c |> dplyr::select(Summary)  |> 
  unnest(Summary) |> 
  filter(variable==".value") |> 
  dplyr::select(-c(variable, BIOREGION.agg)) |> 
  mutate(Indicator="CoralCover",
         Metric=case_match(Metric, "rescale.dist.metric"~"distance.metric",
                           "pcb.rescale.dist.metric"~"consequence.metric")) |> 
  left_join(depth)

.o<-get(load(file="../data/modelled/CO__scores_reef_year.RData")) |> 
  ungroup()

co<-.o |> dplyr::select(Summary)  |> 
  unnest(Summary) |> 
  left_join(spatial_lookup |> 
              dplyr::select(REEF.d, DEPTH.f, Shelf) |> 
              unique()) |> 
  filter(Shelf=="Inshore") |> 
  dplyr::select(-variable, -Shelf) |> 
  mutate(Indicator="Composition",
         Metric=case_match(Metric, "Reference"~"distance.metric",
                           "Critical"~"consequence.metric")) 

.r<-get(load(file="../data/modelled/RPI__scores_reef_year.RData")) |> 
  ungroup()

rp<-.r |> dplyr::select(Summary)  |> 
  unnest(Summary) |> 
  left_join(depth) |> 
  dplyr::select(-c(variable, BIOREGION.agg)) |> 
  mutate(Indicator="Performance",
         Metric=case_match(Metric, "reference"~"distance.metric",
                           "critical"~"consequence.metric")) 
### infill performance of disturbance years

rp.fill<-rp |> 
  group_by(REEF, REEF.d, Metric, DEPTH.f,Indicator) |> 
      arrange(as.numeric(as.character(fYEAR)))  |> 
     tidyr::fill(median, mean, sd, lower, upper, 'p<0.5')  


Framework_indicator_scores_reef<-rbind(ju, ma, cc, co, rp.fill) 



save(Framework_indicator_scores_reef, file="../data/final/Framework_indicator_scores_reef.RData")


# NRM Framework_indicator_scores

j.nrm<-get(load(file="../data/modelled/JU__scores_NRM_year.RData")) |> 
  mutate(Metric=case_match(Metric, "Total"~"baseline", "Acropora"~"consequence"))

ma.nrm<-get(load(file="../data/modelled/MA__scores_NRM_year.RData")) |> 
  mutate(Metric=case_match(Metric, "distance.metric"~"baseline", "consequence.metric"~"consequence"))

cc.nrm<-get(load(file="../data/modelled/CC__scores_NRM_year.RData")) |> 
  mutate(Metric=case_match(Metric, "rescale.dist.metric"~"baseline", "pcb.rescale.dist.metric"~"consequence"))
  
co.nrm<-get(load(file="../data/modelled/CO__scores_NRM_year.RData")) |> 
  mutate(Metric=case_match(Metric, "Reference"~"baseline", "Critical"~"consequence"))

p.nrm<-get(load(file="../data/modelled/RPI__scores_NRM_year.RData")) |> 
  mutate(Metric=case_match(Metric, "reference"~"baseline", "critical"~"consequence"))


Framework_indicator_scores_NRM<-rbind(j.nrm, ma.nrm, cc.nrm, co.nrm, p.nrm) 

save(Framework_indicator_scores_NRM, file="../data/final/Framework_indicator_scores_NRM.RData")

Framework_index_scores_NRM_all<-Framework_indicator_scores_NRM |> 
  filter(Metric=='baseline') |> 
  group_by(fYEAR,NRM) |> 
  summarise(Score=mean(Score, na.rm=TRUE)) |> 
  ungroup()

Framework_index_scores_NRM_exc_composition<-Framework_indicator_scores_NRM |> 
  filter(Metric=='baseline' & ! Indicator=="Composition") |> 
  group_by(fYEAR,NRM) |> 
  summarise(Score_exc_comp=mean(Score, na.rm=TRUE)) |> 
  ungroup()

Framework_index_scores_NRM_exc_composition_combined<-Framework_indicator_scores_NRM |> 
  filter( ! Indicator=="Composition") |> 
  group_by(fYEAR,NRM,Indicator) |> 
  summarise(Score=mean(Score, na.rm=TRUE)) |> 
  ungroup() |> 
  group_by(fYEAR,NRM) |> 
  summarise(Score_exc_comp=mean(Score, na.rm=TRUE)) |> 
  ungroup()

MMP_generateGrades <- function(x) {
  ifelse(is.na(x),'NA',ifelse(x>=0.805, 'A', ifelse(x>=0.605, 'B', ifelse(x>=0.405, 'C',  ifelse(x>=0.205, 'D', 'E'))))) #updated 2018 to reflect grading based on rounding to whole percentage points eg: score of 0.204 rounds to 0.20.
}


# framework categories NRM
cats<-Framework_indicator_scores_NRM |> 
  filter(Metric=="baseline") |> 
  pivot_wider(id_cols=c(fYEAR, NRM),names_from = Indicator, values_from = BelowPar) |>   # not here that BelowPar=0 
  mutate(
         level2=Juvenile+Macroalgae+Composition,
         Juv_MA=Juvenile+Macroalgae,
         cat_grade=ifelse(CoralCover==1 & Performance==1 & level2==3, 'A',
                      ifelse(CoralCover==1 & Performance==1 & level2>0, 'B',
                             ifelse(CoralCover==1 & Performance==1 & level2==0, 'C',
                                    ifelse(CoralCover==1 & Performance==0 & level2==3, 'B',
                                           ifelse(CoralCover==1 & Performance==0 & level2>0, 'C',
                                                  ifelse(CoralCover==1 & Performance==1 & level2==0, 'D',
                                                         ifelse(CoralCover==0 & Performance==1 & level2==3, 'B',
                                                                ifelse(CoralCover==0 & Performance==1 & level2>0, 'C',
                                                                       ifelse(CoralCover==0 & Performance==1 & level2==0, 'D',
                                                                              ifelse(CoralCover==0 & Performance==0 & level2==3, 'D', 'E')))))))))),
        cat_grade_exc_comp=ifelse(CoralCover==1 & Performance==1 & Juv_MA==2, 'A',
                               ifelse(CoralCover==1 & Performance==1 & Juv_MA>0, 'B',
                                      ifelse(CoralCover==1 & Performance==1 & Juv_MA==0, 'C',
                                             ifelse(CoralCover==1 & Performance==0 & Juv_MA==2, 'B',
                                                    ifelse(CoralCover==1 & Performance==0 & Juv_MA>0, 'C',
                                                           ifelse(CoralCover==1 & Performance==1 & Juv_MA==0, 'D',
                                                                  ifelse(CoralCover==0 & Performance==1 & Juv_MA==2, 'B',
                                                                         ifelse(CoralCover==0 & Performance==1 & Juv_MA>0, 'C',
                                                                                ifelse(CoralCover==0 & Performance==1 & Juv_MA==0, 'D',
                                                                                       ifelse(CoralCover==0 & Performance==0 & Juv_MA==2, 'D', 'E')))))))))))
                                          

Framework_index_scores_NRM<- Framework_index_scores_NRM_all |> 
  left_join(Framework_index_scores_NRM_exc_composition) |> 
  mutate(f.grade=MMP_generateGrades(Score),
         f.grade.exc.comp=MMP_generateGrades(Score_exc_comp)) |> 
  left_join(cats)

save(Framework_index_scores_NRM, file="../data/final/Framework_index_scores_NRM.RData")

save(Framework_index_scores_NRM_exc_composition_combined, file="../data/final/Framework_index_scores_NRM_exc_composition_combined.RData")


