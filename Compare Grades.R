
library(tidyverse)
library(patchwork)




MMP_generateGrades <- function(x) {
  ifelse(is.na(x),'NA',ifelse(x>=0.805, 'A', ifelse(x>=0.605, 'B', ifelse(x>=0.405, 'C',  ifelse(x>=0.205, 'D', 'E'))))) #updated 2018 to reflect grading based on rounding to whole percentage points eg: score of 0.204 rounds to 0.20.
}

fw<-get(load(file="../data/final/Framework_indicator_scores_reef.RData")) |> 
  mutate(DEPTH=case_match(DEPTH.f, 
                          "shallow slope" ~"2",
                          "deep slope" ~ "5"),
         REEF=ifelse(REEF %in% c("Border Island",
                     "Fitzroy Island" ,
                     "Havannah Island" ,
                     "Hayman Island",
                     "Langford and Bird Isles",
                     "Pandora Reef"),
                     case_match(REEF, 
                         "Border Island" ~ "Border",
                         "Fitzroy Island" ~ "Fitzroy West LTMP",
                         "Havannah Island" ~ "Havannah North",
                         "Hayman Island" ~  "Hayman",
                         "Langford and Bird Isles" ~ "Langford",
                         "Pandora Reef" ~ "Pandora North"), REEF)
          ) 
fw.grades<- fw |> 
  dplyr::select(-c(DEPTH.f, REEF.d, mean, sd, lower, upper, 'p<0.5')) %>%
  group_by(fYEAR, REEF, DEPTH) |> 
  pivot_wider(names_from = Metric, 
              values_from = 'median') |> 
  mutate(combined.metric=(consequence.metric+distance.metric)/2)|> 
  mutate(fw.ref.grade=MMP_generateGrades(distance.metric),
         fw.con.grade=MMP_generateGrades(consequence.metric),
         fw.comb.grade=MMP_generateGrades(combined.metric),
         fYEAR=as.numeric(as.character(fYEAR)))

  
###
#summary of baselines 

cc.baselines <- get(load(file = paste0(DATA_PATH,
                                    'modelled/CC__baseline_posteriors.RData')))
j.baselines <- get(load(file = paste0(DATA_PATH,
                                       'modelled/JU__baseline_posteriors.RData')))
cc.baselines.median<-cc.baselines |>  
  group_by(BIOREGION.agg, DEPTH.f) |>  
  summarise(median=median(value))

j.baselines.median<-j.baselines |> 
  filter(Taxa=='Total')|>  
  group_by(BIOREGION.agg, DEPTH.f) |>  
  summarise(median=median(value))
#MMP indicator scores

mmp<-get(load(file="../../../../Visit 20 2024/code/data/processed/final/coral.RData"))

mmp.grades<-mmp |> 
  filter(!is.na(HC)) |>  # removes years not sampled by the mmp for which scores and grades would  carried forward 
  #remove Green as not classified as Inshore in the Framework
  filter(!REEF=='Green') |> 
  dplyr::select(REEF,DEPTH,yr,Change.score.linear.mean, CoralCover.score, juv.score, ma.score, composition.score) |> 
  rename (Juvenile =  juv.score,
          Macroalgae = ma.score,
          CoralCover = CoralCover.score,
          Composition = composition.score, 
          Performance = Change.score.linear.mean,
          fYEAR=yr) |> 
pivot_longer(cols=c("Performance", "CoralCover", "Juvenile", "Macroalgae", "Composition"), 
             names_to="Indicator",
             values_to="mmp.score") |> 
  mutate(mmp.grade=MMP_generateGrades(mmp.score),
         fYEAR=as.numeric(fYEAR))


observed.values<-mmp |> 
  filter(!is.na(HC)) |>  # removes years not sampled by the mmp for which scores and grades would  carried forward 
  #remove Green as not classified as Inshore in the Framework
  filter(!REEF=='Green') |> 
  dplyr::select(REEF,DEPTH,yr,HC, SC, CoralCover,  MAprop, MA, juv5.d.nf) |> 
  mutate(fYEAR=as.numeric(yr))

mmp.j5<-read.csv(file="../../../../Visit 20 2024/code/data/primary/juv5mmp.csv") |> 
  rename(total.juv=SUM.ABUNDANCE.)
ltmp.j5<-read.csv(file="../../../../Visit 20 2024/code/data/primary/juv5.ltmp.csv") |> 
  rename(total.juv=SUM.ABUNDANCE.)
j5nf<-read.csv(file="../../../../Visit 20 2024/code/data/primary/juv5nofungia.csv")|> 
  rename(nf=SUM.ABUNDANCE.)
j5nf.ltmp<-read.csv(file="../../../../Visit 20 2024/code/data/primary/juv5nofungia.ltmp.csv")|> 
  rename(nf=SUM.ABUNDANCE.)

obs.juv<-mmp.j5 |> 
  left_join(j5nf) |> 
  rbind(ltmp.j5 |> left_join(j5nf.ltmp)) |> 
  mutate(p.fun=(total.juv-nf)/total.juv,
         DEPTH=as.factor(DEPTH)) |> 
  left_join(mmp |> dplyr::select(REEF,DEPTH,VISIT_NO,yr) |> distinct()) |> 
  group_by(REEF,DEPTH,yr) |> 
  summarise(total.juv=sum(total.juv),
            nf=sum(nf),
            p.fun=mean(p.fun)) |> 
  ungroup() |> 
    mutate(fYEAR=as.numeric(yr)) |> 
  left_join(observed.values |> 
              dplyr::select(REEF,DEPTH,fYEAR,juv5.d.nf))
             

fw.mmp<-mmp.grades |> 
  left_join(fw.grades) 

########################
#Coral cover
########################

#grades
CC<-fw.mmp |> 
  filter(Indicator=='CoralCover') |> 
  group_by(mmp.grade, fw.ref.grade) |> 
  summarise(obs=n()) |> 
  ungroup() |> 
  mutate(fw.ref.grade=factor(fw.ref.grade, levels=c('E','D','C','B','A')),
    mmp.grade=factor(mmp.grade, levels=c('E','D','C','B','A')))

CC.con<-fw.mmp |> 
  filter(Indicator=='CoralCover') |> 
  group_by(mmp.grade, fw.con.grade) |> 
  summarise(obs=n()) |> 
  ungroup() |> 
  mutate(fw.con.grade=factor(fw.con.grade, levels=c('E','D','C','B','A')),
         mmp.grade=factor(mmp.grade, levels=c('E','D','C','B','A')))

CC.comb<-fw.mmp |> 
  filter(Indicator=='CoralCover') |> 
  group_by(mmp.grade, fw.comb.grade) |> 
  summarise(obs=n()) |> 
  ungroup() |> 
  mutate(fw.comb.grade=factor(fw.comb.grade, levels=c('E','D','C','B','A')),
         mmp.grade=factor(mmp.grade, levels=c('E','D','C','B','A')))

 
cc.plot<- ggplot(CC, aes(x=mmp.grade, y=fw.ref.grade)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = obs), colour="black") +
  scale_fill_gradient(low = "white", high = "green") +
  geom_segment(x = "A", y = "A", xend = "E", yend = "E", colour = "red")+
  labs(x="Index grade", y="Framework grade (baseline)") +
  geom_text(aes(label = obs), color = "black", size = 3)+
  theme_bw(base_size=10)+
  theme(legend.position="none") 

cc.plot.con<- ggplot(CC.con, aes(x=mmp.grade, y=fw.con.grade)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = obs), colour="black") +
  scale_fill_gradient(low = "white", high = "green") +
  geom_segment(x = "A", y = "A", xend = "E", yend = "E", colour = "red")+
  labs(x="Index grade", y="Framework grade (consequence)") +
  geom_text(aes(label = obs), color = "black", size = 3)+
  theme_bw(base_size=10)+
  theme(legend.position="none") 

cc.plot.comb<- ggplot(CC.comb, aes(x=mmp.grade, y=fw.comb.grade)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = obs), colour="black") +
  scale_fill_gradient(low = "white", high = "green") +
  geom_segment(x = "A", y = "A", xend = "E", yend = "E", colour = "red")+
  labs(x="Index grade", y="Framework grade (combined)") +
  geom_text(aes(label = obs), color = "black", size = 3)+
  theme_bw(base_size=10) +
  theme(legend.position="none") 

cc.plots<-(cc.plot|cc.plot.con|cc.plot.comb)+
  plot_annotation(tag_levels = 'A')
                      
ggsave(cc.plots, file="../output/figures/CoralCover.grades.png", width=16, height=6, units="cm", scale=1.2)

### scores
Scores.cc.wide<-fw.mmp |> filter(Indicator=="CoralCover") |> 
  left_join(mmp |> 
              mutate(fYEAR=as.numeric(yr),
                     SC_proportion=SC/CoralCover)) 

Scores.cc <-Scores.cc.wide |> 
  dplyr::select(REEF,DEPTH,fYEAR,SC_proportion,HC,CoralCover,mmp.score, distance.metric) |> 
  group_by(REEF,DEPTH,fYEAR,SC_proportion,HC) |> 
  pivot_longer(cols=c("mmp.score","distance.metric"), names_to="metric", values_to="Score") |> 
  mutate(Metric=ifelse(metric=="distance.metric", "Framework", "Index")) |> 
  ungroup()

hc.scores_v_cover<- ggplot(Scores.cc, aes(x=HC, y=Score, colour=Metric)) +                         
  geom_point(aes(size=SC_proportion), shape=21) +
  geom_vline(xintercept=c(75), linetype=c("dotted"))+
  labs(x="Hard coral cover (%)", y="Metric score")+
  theme_bw(base_size=10)+
  theme(legend.position="none") 

cc.scores_v_cover<- ggplot(Scores.cc, aes(x=CoralCover, y=Score, colour=Metric)) +                         
  geom_point(aes(size=SC_proportion), shape=21) +
  geom_vline(xintercept=c(75), linetype=c("dotted"))+
  labs(x="Total coral cover (%)", y="Metric score")+
  theme_bw(base_size=10)

cc.scores_v_cover<-(hc.scores_v_cover|cc.scores_v_cover)+
  plot_annotation(tag_levels = 'A')+ 
  plot_layout(guides = "collect")


ggsave(cc.scores_v_cover, file="../output/figures/CoralCover.scores.png", width=16, height=6, units="cm", scale=1.2)

#######################
#Macroalgae
#######################
ma<-fw.mmp |> 
  filter(Indicator=='Macroalgae') |> 
  group_by(mmp.grade, fw.ref.grade) |> 
  summarise(obs=n()) |> 
  ungroup() |> 
  mutate(fw.ref.grade=factor(fw.ref.grade, levels=c('E','D','C','B','A')),
         mmp.grade=factor(mmp.grade, levels=c('E','D','C','B','A')))

ma.con<-fw.mmp |> 
  filter(Indicator=='Macroalgae') |> 
  group_by(mmp.grade, fw.con.grade) |> 
  summarise(obs=n()) |> 
  ungroup() |> 
  mutate(fw.con.grade=factor(fw.con.grade, levels=c('E','D','C','B','A')),
         mmp.grade=factor(mmp.grade, levels=c('E','D','C','B','A')))

ma.comb<-fw.mmp |> 
  filter(Indicator=='Macroalgae') |> 
  group_by(mmp.grade, fw.comb.grade) |> 
  summarise(obs=n()) |> 
  ungroup() |> 
  mutate(fw.comb.grade=factor(fw.comb.grade, levels=c('E','D','C','B','A')),
         mmp.grade=factor(mmp.grade, levels=c('E','D','C','B','A')))

ma.plot<- ggplot(ma, aes(x=mmp.grade, y=fw.ref.grade)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = obs), colour="black") +
  scale_fill_gradient(low = "white", high = "green") +
  geom_segment(x = "A", y = "A", xend = "E", yend = "E", colour = "red")+
  labs(x="Index grade", y="Framework grade (baseline)") +
  geom_text(aes(label = obs), color = "black", size = 3)+
  theme_bw(base_size=10)+
  theme(legend.position="none") 

ma.plot.con<- ggplot(ma.con, aes(x=mmp.grade, y=fw.con.grade)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = obs), colour="black") +
  scale_fill_gradient(low = "white", high = "green") +
  geom_segment(x = "A", y = "A", xend = "E", yend = "E", colour = "red")+
  labs(x="Index grade", y="Framework grade (consequence)") +
  geom_text(aes(label = obs), color = "black", size = 3)+
  theme_bw(base_size=10)+
  theme(legend.position="none") 

ma.plot.comb<- ggplot(ma.comb, aes(x=mmp.grade, y=fw.comb.grade)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = obs), colour="black") +
  scale_fill_gradient(low = "white", high = "green") +
  geom_segment(x = "A", y = "A", xend = "E", yend = "E", colour = "red")+
  labs(x="Index grade", y="Framework grade (combined)") +
  geom_text(aes(label = obs), color = "black", size = 3)+
  theme_bw(base_size=10) +
  theme(legend.position="none") 

ma.plots<-(ma.plot|ma.plot.con|ma.plot.comb)+
  plot_annotation(tag_levels = 'A')

ggsave(ma.plots, file="../output/figures/Macroalgae.grades.png", width=16, height=6, units="cm", scale=1.2)

### scores
Scores.ma.wide<-fw.mmp |> filter(Indicator=="Macroalgae") |> 
  left_join(observed.values)

Scores.ma <-Scores.ma.wide |> 
  dplyr::select(REEF,DEPTH,fYEAR,MA,MAprop,mmp.score, distance.metric) |> 
  group_by(REEF,DEPTH,fYEAR,MA,MAprop) |> 
  pivot_longer(cols=c("mmp.score","distance.metric"), names_to="metric", values_to="Score") |> 
  mutate(Metric=ifelse(metric=="distance.metric", "Framework", "Index")) |> 
  ungroup()

ma.scores_v_maprop<- ggplot(Scores.ma, aes(x=MAprop, y=Score, colour=Metric)) +                         
  geom_jitter(height=0.002, shape=21) +
 labs(x="Macroalgae proportion", y="Metric score",title='A: Framework (baseline)')+
  geom_vline(xintercept=c(23,25), linetype='dotted')+
  theme_bw(base_size=10)+
  theme(legend.position="none",
        plot.title = element_text(size = 10))  

Scores.ma.con <-Scores.ma.wide |> 
  dplyr::select(REEF,DEPTH,fYEAR,MA,MAprop,mmp.score, consequence.metric) |> 
  group_by(REEF,DEPTH,fYEAR,MA,MAprop) |> 
  pivot_longer(cols=c("mmp.score","consequence.metric"), names_to="metric", values_to="Score") |> 
  mutate(Metric=ifelse(metric=="consequence.metric", "Framework", "Index")) |> 
  ungroup()

ma.scores.con_v_maprop<- ggplot(Scores.ma.con, aes(x=MAprop, y=Score, colour=Metric)) +                         
  geom_jitter(height=0.002, shape=21) +
  labs(x="Macroalgae proportion", y="Metric score", title='B: Framework (consequence)')+
  geom_vline(xintercept=c(23,25), linetype='dotted')+
  theme_bw(base_size=10)+
  theme(legend.position="none",
        plot.title = element_text(size = 10))  

Scores.ma.com <-Scores.ma.wide |> 
  dplyr::select(REEF,DEPTH,fYEAR,MA,MAprop,mmp.score, combined.metric) |> 
  group_by(REEF,DEPTH,fYEAR,MA,MAprop) |> 
  pivot_longer(cols=c("mmp.score","combined.metric"), names_to="metric", values_to="Score") |> 
  mutate(Metric=ifelse(metric=="combined.metric", "Framework", "Index")) |> 
  ungroup()

ma.scores.com_v_maprop<- ggplot(Scores.ma.com, aes(x=MAprop, y=Score, colour=Metric)) +                         
  geom_jitter(height=0.002, shape=21) +
  labs(x="Macroalgae proportion", y="Metric score", title='C: Framework (combined)')+
  labs(x="Macroalgae proportion", y="Metric score")+
  geom_vline(xintercept=c(23,25), linetype='dotted')+
  theme_bw(base_size=10)+
  theme(legend.position="right",
        plot.title = element_text(size = 10)) 

ma.scores<-(ma.scores_v_maprop|ma.scores.con_v_maprop|ma.scores.com_v_maprop)+
  plot_layout(guides = "collect")


ggsave(ma.scores, file="../output/figures/MA.scores.png", width=16, height=6, units="cm", scale=1.2)

############
# Juveniles
############

j<-fw.mmp |> 
  filter(Indicator=='Juvenile' & fYEAR>2006) |>  # the framework estimates start in 2007
  filter(!(REEF=='Low Isles' & fYEAR=='2015')) |>  # juvs not counted in this visit but mmp carried through the score
  group_by(mmp.grade, fw.ref.grade) |> 
  summarise(obs=n()) |> 
  ungroup() |> 
  mutate(fw.ref.grade=factor(fw.ref.grade, levels=c('E','D','C','B','A')),
         mmp.grade=factor(mmp.grade, levels=c('E','D','C','B','A')))

j.con<-fw.mmp |> 
  filter(Indicator=='Juvenile' & fYEAR>2006) |> 
  filter(!(REEF=='Low Isles' & fYEAR=='2015')) |>
  group_by(mmp.grade, fw.con.grade) |> 
  summarise(obs=n()) |> 
  ungroup() |> 
  mutate(fw.con.grade=factor(fw.con.grade, levels=c('E','D','C','B','A')),
         mmp.grade=factor(mmp.grade, levels=c('E','D','C','B','A')))

j.comb<-fw.mmp |> 
  filter(Indicator=='Juvenile' & fYEAR>2006) |> 
  filter(!(REEF=='Low Isles' & fYEAR=='2015')) |>
  group_by(mmp.grade, fw.comb.grade) |> 
  summarise(obs=n()) |> 
  ungroup() |> 
  mutate(fw.comb.grade=factor(fw.comb.grade, levels=c('E','D','C','B','A')),
         mmp.grade=factor(mmp.grade, levels=c('E','D','C','B','A')))

j.plot<- ggplot(j, aes(x=mmp.grade, y=fw.ref.grade)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = obs), colour="black") +
  scale_fill_gradient(low = "white", high = "green") +
  geom_segment(x = "A", y = "A", xend = "E", yend = "E", colour = "red")+
  labs(x="Index grade", y="Framework grade (baseline)") +
  geom_text(aes(label = obs), color = "black", size = 3)+
  theme_bw(base_size=10)+
  theme(legend.position="none") 

j.plot.con<- ggplot(j.con, aes(x=mmp.grade, y=fw.con.grade)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = obs), colour="black") +
  scale_fill_gradient(low = "white", high = "green") +
  geom_segment(x = "A", y = "A", xend = "E", yend = "E", colour = "red")+
  labs(x="Index grade", y="Framework grade (consequence)") +
  geom_text(aes(label = obs), color = "black", size = 3)+
  theme_bw(base_size=10)+
  theme(legend.position="none") 

j.plot.comb<- ggplot(j.comb, aes(x=mmp.grade, y=fw.comb.grade)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = obs), colour="black") +
  scale_fill_gradient(low = "white", high = "green") +
  geom_segment(x = "A", y = "A", xend = "E", yend = "E", colour = "red")+
  labs(x="Index grade", y="Framework grade (combined)") +
  geom_text(aes(label = obs), color = "black", size = 3)+
  theme_bw(base_size=10) +
  theme(legend.position="none") 

j.plots<-(j.plot|j.plot.con|j.plot.comb)+
  plot_annotation(tag_levels = 'A')

ggsave(j.plots, file="../output/figures/juvenile.grades.png", width=16, height=6, units="cm", scale=1.2)

### scores

load("../data/processed/juv.df.RData")
# aggregate juv.data to reef
j.data<-juv.df |> 
  mutate(DEPTH=case_match(DEPTH.f, 
                          "shallow slope" ~"2",
                          "deep slope" ~ "5"),
         REEF=ifelse(REEF %in% c("Border Island",
                                 "Fitzroy Island" ,
                                 "Havannah Island" ,
                                 "Hayman Island",
                                 "Langford and Bird Isles",
                                 "Pandora Reef"),
                     case_match(REEF, 
                                "Border Island" ~ "Border",
                                "Fitzroy Island" ~ "Fitzroy West LTMP",
                                "Havannah Island" ~ "Havannah North",
                                "Hayman Island" ~  "Hayman",
                                "Langford and Bird Isles" ~ "Langford",
                                "Pandora Reef" ~ "Pandora North"), REEF)
  ) |> 
  group_by(REEF,DEPTH,REEF.d,DEPTH.f,fYEAR) |> 
  summarise(total.juv=sum(total.juv, na.rm=TRUE),
            Acropora=sum(Acropora, na.rm=TRUE),
            avail.area=sum(avail.area)) |> 
  mutate(tot.d=total.juv/avail.area,
         p.acr=Acropora/total.juv,
         fYEAR=as.numeric(as.character(fYEAR))) |> 
left_join(obs.juv) |>  
   ungroup()

Scores.j.wide<-fw.mmp |> filter(Indicator=="Juvenile") |> 
  left_join(j.data) |> 
  filter(fYEAR>2006) |> 
  filter(!(REEF=='Low Isles' & fYEAR=='2015'))   # juvs not counted in this visit but mmp carried through the score

Scores.j.base <-Scores.j.wide |> 
  dplyr::select(REEF,DEPTH,fYEAR,p.fun, tot.d, mmp.score, distance.metric) |> 
  group_by(REEF,DEPTH,fYEAR,p.fun,tot.d) |> 
  pivot_longer(cols=c("mmp.score","distance.metric"), names_to="metric", values_to="Score") |> 
  mutate(Metric=ifelse(metric=="distance.metric", "Framework", "Index")) |> 
  ungroup()

j.scores.ref_v_den<- ggplot(Scores.j.base, aes(x=tot.d, y=Score, fill=Metric)) +                         
  geom_point(shape=21) +
  labs(x="Juvenile density", 
       y="Metric score", 
       title='A: Framework (baseline)')+
  theme_bw(base_size=10)+
  theme(plot.title = element_text(size = 10)) 
 
Scores.j.con <-Scores.j.wide |> 
  dplyr::select(REEF,DEPTH,fYEAR,p.acr, tot.d, mmp.score, consequence.metric) |> 
  group_by(REEF,DEPTH,fYEAR,p.acr,tot.d) |> 
  pivot_longer(cols=c("mmp.score","consequence.metric"), names_to="metric", values_to="Score") |> 
  mutate(Metric=ifelse(metric=="consequence.metric", "Framework", "Index")) |> 
  ungroup()

j.scores.con_v_cover<- ggplot(Scores.j.con, aes(x=tot.d, y=Score, fill=Metric)) +                         
  geom_point( shape=21) +
  labs(x="Juvenile density", 
       y="Metric score", 
       title='B: Framework (consequence)')+
  theme_bw(base_size=10)+
  theme(plot.title = element_text(size = 10))
  

j.score.plots<-(j.scores.ref_v_den|j.scores.con_v_cover)+
  #plot_annotation(tag_levels = 'A')+ 
  plot_layout(guides = "collect")


ggsave(j.score.plots, file="../output/figures/juv.scores.png", width=16, height=6, units="cm", scale=1.2)

                                
j.score.dif<-Scores.j.wide |>
  mutate(score.dif=mmp.score-distance.metric,
         score.dif.con=mmp.score-consequence.metric)
  
j.scores.dif.ref<- ggplot(j.score.dif, aes(x=score.dif, y=p.fun)) +                         
  geom_point() +
 
  labs(y="Proportion of juveniles Fungia", 
       x="Metric score difference", 
       title='Framework (baseline)')+
  theme_bw(base_size=10)+
  theme(plot.title = element_text(size = 10)) 

j.scores.dif.con<- ggplot(j.score.dif, aes(x=score.dif, y=p.acr)) +                         
  geom_point() +
  labs(y="Proportion of juvenile Acropora", 
       x="Metric score difference", 
       title='B: Framework (consequence)')+
  theme_bw(base_size=10)+
  theme(plot.title = element_text(size = 10)) 

j.score.diferences<-(j.scores.dif.ref|j.scores.dif.con) #+
#plot_annotation(tag_levels = 'A') 
#plot_layout(guides = "collect")


ggsave(j.scores.dif.ref, file="../output/figures/juv.scores.diff.png", width=8, height=6, units="cm", scale=1.2)


########################
#Performance
########################

change<-fw.mmp |> 
  filter(!is.na(mmp.score) & Indicator=='Performance' )  
 #dim(change) 1089

change.ref<-change |> filter(!is.na(distance.metric))
#dim(change.ref) 745

change.con<-change.ref |> filter(!is.na(consequence.metric))
#dim(change.con) 568

P<-change.ref |> 
  group_by(mmp.grade, fw.ref.grade) |> 
  summarise(obs=n()) |> 
  ungroup() |> 
  mutate(fw.ref.grade=factor(fw.ref.grade, levels=c('E','D','C','B','A')),
         mmp.grade=factor(mmp.grade, levels=c('E','D','C','B','A')))

P.con<-change.con |> 
  group_by(mmp.grade, fw.con.grade) |> 
  summarise(obs=n()) |> 
  ungroup() |> 
  mutate(fw.con.grade=factor(fw.con.grade, levels=c('E','D','C','B','A')),
         mmp.grade=factor(mmp.grade, levels=c('E','D','C','B','A')))

P.comb<-change.con |> 
  group_by(mmp.grade, fw.comb.grade) |> 
  summarise(obs=n()) |> 
  ungroup() |> 
  mutate(fw.comb.grade=factor(fw.comb.grade, levels=c('E','D','C','B','A')),
         mmp.grade=factor(mmp.grade, levels=c('E','D','C','B','A')))

p.plot<- ggplot(P, aes(x=mmp.grade, y=fw.ref.grade)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = obs), colour="black") +
  scale_fill_gradient(low = "white", high = "green") +
  geom_segment(x = "A", y = "A", xend = "E", yend = "E", colour = "red")+
  labs(x="Index grade", y="Framework grade (baseline)") +
  geom_text(aes(label = obs), color = "black", size = 3)+
  theme_bw(base_size=10)+
  theme(legend.position="none") 

p.plot.con<- ggplot(P.con, aes(x=mmp.grade, y=fw.con.grade)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = obs), colour="black") +
  scale_fill_gradient(low = "white", high = "green") +
  geom_segment(x = "A", y = "A", xend = "E", yend = "E", colour = "red")+
  labs(x="Index grade", y="Framework grade (consequence)") +
  geom_text(aes(label = obs), color = "black", size = 3)+
  theme_bw(base_size=10)+
  theme(legend.position="none") 

p.plot.comb<- ggplot(P.comb, aes(x=mmp.grade, y=fw.comb.grade)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = obs), colour="black") +
  scale_fill_gradient(low = "white", high = "green") +
  geom_segment(x = "A", y = "A", xend = "E", yend = "E", colour = "red")+
  labs(x="Index grade", y="Framework grade (combined)") +
  geom_text(aes(label = obs), color = "black", size = 3)+
  theme_bw(base_size=10) +
  theme(legend.position="none") 

p.plots<-(p.plot|p.plot.con|p.plot.comb)+
  plot_annotation(tag_levels = 'A')

ggsave(p.plots, file="../output/figures/Performance.grades.png", width=16, height=6, units="cm", scale=1.2)

########################
#Composition
########################

comp<-fw.mmp |> 
  filter(!is.na(mmp.score) & Indicator=='Composition' )  
#dim(change) 1089

comp.ref<-comp |> filter(!is.na(distance.metric))
#dim(change.ref) 745


CO<-comp.ref |> 
  group_by(mmp.grade, fw.ref.grade) |> 
  summarise(obs=n()) |> 
  ungroup() |> 
  mutate(fw.ref.grade=factor(fw.ref.grade, levels=c('E','D','C','B','A')),
         mmp.grade=factor(mmp.grade, levels=c('E','D','C','B','A')))

CO.con<-comp.ref |> 
  group_by(mmp.grade, fw.con.grade) |> 
  summarise(obs=n()) |> 
  ungroup() |> 
  mutate(fw.con.grade=factor(fw.con.grade, levels=c('E','D','C','B','A')),
         mmp.grade=factor(mmp.grade, levels=c('E','D','C','B','A')))

CO.plot<- ggplot(CO, aes(x=mmp.grade, y=fw.ref.grade)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = obs), colour="black") +
  scale_fill_gradient(low = "white", high = "green") +
  geom_segment(x = "A", y = "A", xend = "E", yend = "E", colour = "red")+
  labs(x="Index grade", y="Framework grade (baseline)") +
  geom_text(aes(label = obs), color = "black", size = 3)+
  theme_bw(base_size=10)+
  theme(legend.position="none") 

CO.plot.con<- ggplot(CO.con, aes(x=mmp.grade, y=fw.con.grade)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = obs), colour="black") +
  scale_fill_gradient(low = "white", high = "green") +
  geom_segment(x = "A", y = "A", xend = "E", yend = "E", colour = "red")+
  labs(x="Index grade", y="Framework grade (consequence)") +
  geom_text(aes(label = obs), color = "black", size = 3)+
  theme_bw(base_size=10)+
  theme(legend.position="none") 


CO.plots<-(CO.plot|CO.plot.con)+
  plot_annotation(tag_levels = 'A')

ggsave(CO.plots, file="../output/figures/composition.grades.png", width=12, height=6, units="cm", scale=1.2)

