library(tidyverse)

load(file="../data/final/Framework_index_scores_NRM.RData")
load(file="../data/final/Framework_index_scores_NRM_exc_composition_combined.RData")

# distribution of mean scores within categorical grades

trafficLightPalette <- (c('#FF0000','#FFC000','#FFFF00','#92D050','#00B050'))
lims <- LETTERS[1:5]

grade.plot.framework<-ggplot(Framework_index_scores_NRM |>  filter(!is.na(cat_grade)), aes(x = cat_grade, y = Score, fill=cat_grade))+
  geom_violin()+
  scale_fill_manual(values=rev(trafficLightPalette), limits=lims)+
  labs(x="Framework categorised grade", y="Mean indicator score") +
  scale_y_continuous(breaks=c(0.2,0.4,0.6,0.8), limits=c(0.1,0.9))+
  geom_text(data = count(Framework_index_scores_NRM |> filter(!is.na(cat_grade)),cat_grade), aes(y = 0.1, label = n))+
  geom_hline(yintercept=c(0.2,0.4,0.6,0.8), linetype="dashed", colour='black', linewidth=0.5)+
  theme_bw(base_size=10)+
  theme(legend.position="none") 

# distribution of grades when the composition score is excluded
# also drop NA values due to pre huvenile years and years where performance could not be estimated
grade.plot.framework_noComp<-ggplot(Framework_index_scores_NRM |> filter(!is.na(cat_grade_exc_comp)), aes(x = cat_grade_exc_comp, y = Score_exc_comp, fill=cat_grade_exc_comp))+
  geom_violin()+
  scale_fill_manual(values=rev(trafficLightPalette), limits=lims)+
  labs(x="Framework categorised grade", y="Mean indicator score") +
  scale_y_continuous(breaks=c(0.2,0.4,0.6,0.8), limits=c(0.1,0.9))+
  geom_text(data = count(Framework_index_scores_NRM |> filter(!is.na(cat_grade_exc_comp)),cat_grade_exc_comp), aes(y = 0.1, label = n))+
  geom_hline(yintercept=c(0.2,0.4,0.6,0.8), linetype="dashed", colour='black', linewidth=0.5)+
  theme_bw(base_size=10)+
  theme(legend.position="none") 

framework.index.scores<-(grade.plot.framework|grade.plot.framework_noComp)+
  plot_annotation(tag_levels = 'A')

ggsave(framework.index.scores, file="../output/figures/Framework.grades.png", width=16, height=6, units="cm", scale=1.2)





# next bring in the coral index scores
# mmp.index<-get(load(file="../../../../Visit 20 2024/code/data/processed/boot.index.region.RData"))
mmp.index<-get(load(file="../data/processed/boot.index.region.RData"))


mmp<-mmp.index$sum

framework<- Framework_index_scores_NRM |> 
  rename('Region'='NRM', f.score=Score) |> 
  mutate(Year=as.numeric(as.character(fYEAR))) |> 
  dplyr::select(Year, Region, f.score,cat_grade,Score_exc_comp,cat_grade_exc_comp)

mmp.frame<-mmp |> 
  left_join(framework)

mmp.framework.regional.scores<-ggplot(mmp.frame |> filter(!is.na(cat_grade)), aes(x = f.score, y = Score, fill=cat_grade))+
  geom_point(shape=21, size=3)+
  labs(x="", y="Coral Index score", fill="Framework Category") +
  scale_fill_manual(values=rev(trafficLightPalette), limits=lims)+
  scale_y_continuous(breaks=c(0.2,0.4,0.6,0.8), limits=c(0.19,0.81))+
  scale_x_continuous(breaks=c(0.2,0.4,0.6,0.8), limits=c(0.19,0.81))+
  geom_abline(intercept=0, slope=1)+
  geom_hline(yintercept=c(0.2,0.4,0.6,0.8), linetype="dashed", colour='black', linewidth=0.5)+
  geom_vline(xintercept=c(0.2,0.4,0.6,0.8), linetype="dashed", colour='black', linewidth=0.5)+
  theme_bw(base_size=10)+
  theme(legend.position="bottom") 

mmp.framework.regional.scores.trend<-ggplot(mmp.frame |> filter(Year>2016), aes(x = Year))+
  scale_y_continuous(breaks=c(0.2,0.4,0.6), limits=c(0.19,0.70))+
  geom_line(aes(y=Score, color="Coral Index")) +
  geom_line(aes(y=f.score, color="Framework")) +
  geom_point(aes(y=Score, fill=Grade), shape=21,size=3) +
  geom_point(aes(y=f.score, fill=cat_grade), shape=21,size=3) +
  scale_fill_manual(values=rev(trafficLightPalette), limits=lims)+
  labs(x="", y="Mean indicator score", colour='Scoring method') +
  geom_hline(yintercept=c(0.205,0.405,0.605), linetype="dashed", colour='black', linewidth=0.5)+
  geom_vline(xintercept=c(), linetype="dashed", colour='black', linewidth=0.5)+
  theme_bw(base_size=10)+
  theme(legend.position="right") +
  facet_wrap(vars(Region))
ggsave(mmp.framework.regional.scores.trend, file="../output/figures/regional.grade.comparison.png", width=16, height=16, units="cm", scale=1.2)




#### MMP Index scores excluding composition
mmp.index.no.comp<-get(load(file="../data/processed/boot.index.region.no.comp.RData"))
mmp.no.comp<-mmp.index.no.comp$sum

framework.comb.ex.comp<-Framework_index_scores_NRM_exc_composition_combined |> 
  rename('Region'='NRM', f.score_exc_comp=Score_exc_comp) |> 
  mutate(Year=as.numeric(as.character(fYEAR))) |> 
  dplyr::select(Year, Region, f.score_exc_comp)

mmp.frame.no.comp<-mmp.no.comp |> 
  left_join(framework.comb.ex.comp)


mmp.framework.regional.scores.trend.no.comp<-ggplot(mmp.frame.no.comp |> filter(Year>2016), aes(x = Year))+
  scale_y_continuous(breaks=c(0.2,0.4,0.6), limits=c(0.19,0.70))+
  geom_line(aes(y=Score, color="Coral Index")) +
  geom_line(aes(y=f.score_exc_comp, color="Framework \n combined \n metrics")) +
  geom_point(aes(y=Score, fill=Grade), shape=21,size=3) +
  geom_point(aes(y=f.score_exc_comp), shape=21,size=3) +
  scale_fill_manual(values=rev(trafficLightPalette), limits=lims)+
  labs(x="", y="Mean indicator score excluding composition", colour='Scoring method') +
  geom_hline(yintercept=c(0.205,0.405,0.605), linetype="dashed", colour='black', linewidth=0.5)+
  geom_vline(xintercept=c(), linetype="dashed", colour='black', linewidth=0.5)+
  theme_bw(base_size=10)+
  theme(legend.position="right") +
  facet_wrap(vars(Region))

ggsave(mmp.framework.regional.scores.trend.no.comp, file="../output/figures/regional.grade.comparison.comb.no.comp.png", width=16, height=16, units="cm", scale=1.2)
