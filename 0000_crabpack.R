#devtools::install_github("AFSC-Shellfish-Assessment-Program/crabpack")
library(crabpack)
library(dplyr)
library(ggplot2)
library(reshape2)



## Connect to Oracle
#channel <- crabpack::get_connected()

## Pull specimen data
specimen_data <- crabpack::get_specimen_data(species = "SNOW",
                                             region = "EBS",
                                             years = c(1982:2024),
                                             channel = 'API')


test_101<-crabpack::calc_bioabund(crab_data = specimen_data,
                                  species = "SNOW",
                                  region = "EBS",
                                  sex='male',
                                  crab_category = c("preferred_male"))

test_101$mod_abnd<-test_101$TOTAL_AREA/max(test_101$TOTAL_AREA)

ggplot(test_101)+
  geom_line(aes(x=YEAR,y=ABUNDANCE))+
  geom_line(aes(x=YEAR,y=ABUNDANCE/mod_abnd),col=2)

#==compare to existing data pulls
#==crabpack
AB_Sizegroup <- crabpack::calc_bioabund(crab_data = specimen_data,
                                        species = "SNOW",
                                        region = "EBS",
                                        sex='male',
                                        crab_category = c("legal_male","all_categories","large_male","preferred_male"))

EBSdata	<-read.csv("data/survey/EBSCrab_AB_Sizegroup.csv")
names(EBSdata)[1]<-"YEAR"

SurveyTotals<- filter(EBSdata,SIZE_CLASS_MM=='TOTAL')%>%
  group_by(YEAR,SEX) %>%
  summarise(Totals = sum(ABUNDANCE,na.rm=T))

#==make upper and lower CIs
EBSdata$upper_bio_ci<-exp(log(EBSdata$BIOMASS_MT)+1.96*(sqrt(log(1+EBSdata$BIOMASS_MT_CV^2))))
EBSdata$lower_bio_ci<-exp(log(EBSdata$BIOMASS_MT)-1.96*(sqrt(log(1+EBSdata$BIOMASS_MT_CV^2))))
EBSdata$upper_abnd_ci<-exp(log(EBSdata$ABUNDANCE)+1.96*(sqrt(log(1+EBSdata$ABUNDANCE_CV^2))))
EBSdata$lower_abnd_ci<-exp(log(EBSdata$ABUNDANCE)-1.96*(sqrt(log(1+EBSdata$ABUNDANCE_CV^2)))) 

#==make df for comparison
plot_males<-c("MALE_GE78","MALE_GE95","MALE_GE102","MALE_TOTAL")
replacement_values <- c("legal_male","large_male","preferred_male","all_categories")
tmp_old<-filter(EBSdata,SIZE_GROUP%in%plot_males)
for(x in 1:length(plot_males))
  tmp_old$SIZE_GROUP[tmp_old$SIZE_GROUP==plot_males[x]]<-replacement_values[x]

tmp_old$data<-"AKFIN"
colnames(tmp_old)[5]<-"CATEGORY"
AB_Sizegroup$data<-"crabpack"
use_em<-intersect(colnames(tmp_old),colnames(AB_Sizegroup))

comp_ugh<-rbind(tmp_old%>%select(all_of(use_em)),
AB_Sizegroup%>%select(all_of(use_em)))

#==no 'all_categories' here
ggplot(comp_ugh)+
  geom_line(aes(x=YEAR,y=BIOMASS_MT,col=data))+
  facet_wrap(~CATEGORY)

#==========================================================
# compare female size comps and index of abundance
#=============================================================
#==existing code
# this file is from EBC Crab, Large data Download, Abundance/Biomass, with MATURITY_SIZE_1MM_SHELLCON
EBSdata_in<-read.csv("data/survey/EBSCrab_Abundance_Biomass_female.csv",skip=7)
EBSdata<-EBSdata_in[EBSdata_in$SIZE_CLASS_MM>=25 & EBSdata_in$SEX=="FEMALE",]

female_total<- filter(EBSdata)%>%
  group_by(SURVEY_YEAR) %>%
  summarise(Totals = sum(ABUNDANCE))

Years				<-sort(unique(EBSdata$SURVEY_YEAR))
LengthBins			<-seq(25,135,5)
mid_pts<-seq(27.5,132.5,5)
FemaleMature		<-matrix(nrow=length(Years),ncol=length(LengthBins))
FemaleImmature		<-matrix(nrow=length(Years),ncol=length(LengthBins))

for(x in 1:length(Years))
{
  temp<-EBSdata[EBSdata$SURVEY_YEAR==Years[x],]
  for(y in 1:(length(LengthBins)-1))
  {
    FemaleMature[x,y]	<-sum(temp$ABUNDANCE[ temp$MATURITY=="MATURE" & temp$SIZE_CLASS_MM>=LengthBins[y] & temp$SIZE_CLASS_MM<LengthBins[y+1] ])
    FemaleImmature[x,y]	<-sum(temp$ABUNDANCE[ temp$MATURITY=="IMMATURE" & temp$SIZE_CLASS_MM>=LengthBins[y] & temp$SIZE_CLASS_MM<LengthBins[y+1] ])
  }
}

FemaleMature<-FemaleMature
FemaleImmature<-FemaleImmature

nal_fem_mat<-round(FemaleMature,4)[,1:(length(LengthBins)-1)]
sc_fem_mat<-sweep(nal_fem_mat,1,apply(nal_fem_mat,1,sum,na.rm=T),FUN="/")

nal_fem_imm<-round(FemaleImmature,4)[,1:(length(LengthBins)-1)]
sc_fem_imm<-sweep(nal_fem_imm,1,apply(nal_fem_imm,1,sum,na.rm=T),FUN="/")

#==========================
#==calc female biomasses===
#==========================
fem_bios<-EBSdata%>%
  group_by(SURVEY_YEAR,MATURITY) %>%
  summarise(biomass = round(sum(BIOMASS_MT,na.rm=T),2))


#===========================
# crabpack
#====================

#==indices of abundance
mat_fem_snow <- crabpack::calc_bioabund(crab_data = specimen_data,
                                        species = "SNOW",
                                        region = "EBS",
                                        crab_category = "mature_female",
                                        female_maturity = "morphometric",
                                        size_min = 25,
                                        bin_1mm = TRUE)

nal_fem_mat_cp<-mat_fem_snow %>%
  group_by(YEAR,group = cut(SIZE_1MM , breaks = seq(0, max(SIZE_1MM), 5))) %>%
  summarise(n = sum(ABUNDANCE),
            data='crabpack')


midpts_cp<-seq(22.5,102.5,5)
tmp<-strsplit(as.character(nal_fem_mat_cp$group),',')
nal_fem_mat_cp$mid_pts<-NA
for(x in 1:length(tmp))
{
  in1<-as.numeric(substr(unlist(tmp[x])[1],2, nchar(unlist(tmp[x])[1])))
  in2<-as.numeric(substr(unlist(tmp[x])[2],1, nchar(unlist(tmp[x])[2])-1))
  nal_fem_mat_cp$mid_pts[x]<-sum(in1,in2)/2
}


imm_fem_snow <- crabpack::calc_bioabund(crab_data = specimen_data,
                                        species = "SNOW",
                                        region = "EBS",
                                        crab_category = "immature_female",
                                        female_maturity = "morphometric",
                                        size_min = 25,
                                        bin_1mm = TRUE)

nal_fem_im_cp<-imm_fem_snow %>%
  group_by(YEAR,group = cut(SIZE_1MM , breaks = seq(0, max(SIZE_1MM), 5))) %>%
  summarise(n = sum(ABUNDANCE),
            data='crabpack')

midpts_cp<-seq(22.5,102.5,5)
tmp<-strsplit(as.character(nal_fem_im_cp$group),',')
nal_fem_im_cp$mid_pts<-NA
for(x in 1:length(tmp))
{
 in1<-as.numeric(substr(unlist(tmp[x])[1],2, nchar(unlist(tmp[x])[1])))
 in2<-as.numeric(substr(unlist(tmp[x])[2],1, nchar(unlist(tmp[x])[2])-1))
 nal_fem_im_cp$mid_pts[x]<-sum(in1,in2)/2
}
 
 
rownames(nal_fem_mat)<-Years
colnames(nal_fem_mat)<-mid_pts
nal_fem_mat_og<-melt(nal_fem_mat)
nal_fem_mat_og$data<-"OG"
colnames(nal_fem_mat_og)<-c("YEAR","mid_pts","n",'data')
rownames(nal_fem_imm)<-Years
colnames(nal_fem_imm)<-mid_pts
nal_fem_imm_og<-melt(nal_fem_imm)
nal_fem_imm_og$data<-"OG"
colnames(nal_fem_imm_og)<-c("YEAR","mid_pts","n",'data')

comp_fem_mat_plt<-rbind(nal_fem_mat_cp[,c(1,5,3,4)],nal_fem_mat_og)
comp_fem_imm_plt<-rbind(nal_fem_im_cp[,c(1,5,3,4)],nal_fem_imm_og)

fem_mat_comp_sc<-ggplot(comp_fem_mat_plt)+
  geom_line(aes(x=mid_pts,y=n,col=data))+
  facet_wrap(~YEAR)+theme_bw()+xlim(25,75)

fem_imm_comp_sc<-ggplot(comp_fem_imm_plt)+
  geom_line(aes(x=mid_pts,y=n,col=data))+
  facet_wrap(~YEAR)+theme_bw()+xlim(25,75)

comp_imm_tots<-filter(comp_fem_imm_plt,mid_pts>25)%>%
  group_by(YEAR,data)%>%
  summarize(tot=sum(n))

comp_mat_tots<-filter(comp_fem_mat_plt,mid_pts>25)%>%
  group_by(YEAR,data)%>%
  summarize(tot=sum(n))

comp_mat_tots_pl<-ggplot(comp_mat_tots)+
  geom_line(aes(x=YEAR,y=tot,col=data))+
  theme_bw()

comp_imm_tots_pl<-ggplot(comp_imm_tots)+
  geom_line(aes(x=YEAR,y=tot,col=data))+
  theme_bw()

png("plots/data_comp/fem_mat_comp_sc.png",height=8,width=8,res=400,units='in')
print(fem_mat_comp_sc)
dev.off()

png("plots/data_comp/fem_imm_comp_sc.png",height=8,width=8,res=400,units='in')
print(fem_imm_comp_sc)
dev.off()

png("plots/data_comp/comp_imm_tots_pl.png",height=8,width=8,res=400,units='in')
print(comp_imm_tots_pl)
dev.off()

png("plots/data_comp/comp_mat_tots_pl.png",height=8,width=8,res=400,units='in')
print(comp_mat_tots_pl)
dev.off()

#====female biomass===
EBSdata_in<-read.csv("data/survey/EBSCrab_Abundance_Biomass_female.csv",skip=7)
EBSdata<-EBSdata_in[EBSdata_in$SIZE_CLASS_MM>=25 & EBSdata_in$SEX=="FEMALE",]

fem_bios<-EBSdata%>%
  group_by(SURVEY_YEAR,MATURITY) %>%
  summarise(biomass = round(sum(BIOMASS_MT,na.rm=T),2))
fem_bios$data<-"OG"
colnames(fem_bios)<-c("Year","Maturity","Biomass","data")
mat_fem_snow_ind <- crabpack::calc_bioabund(crab_data = specimen_data,
                                        species = "SNOW",
                                        region = "EBS",
                                        crab_category = c("immature_female","mature_female"),
                                        female_maturity = "morphometric")

comp_bio<-mat_fem_snow_ind[,c(2,6,10)]
comp_bio$CATEGORY[comp_bio$CATEGORY=="immature_female"]<-"IMMATURE"
comp_bio$CATEGORY[comp_bio$CATEGORY=="mature_female"]<-"MATURE"
comp_bio$data<-"crabpack"
colnames(comp_bio)<-c("Year","Maturity","Biomass","data")

plot_fem_bio<-rbind(fem_bios,comp_bio)
fem_bio_comp<-ggplot(plot_fem_bio)+
  geom_line(aes(x=Year,y=Biomass,col=data),lwd=1.2,alpha=.8)+
  facet_wrap(~Maturity)+
  theme_bw()

png("plots/data_comp/comp_fem_bio.png",height=5,width=8,res=400,units='in')
print(fem_bio_comp)
dev.off()


#==========================================================
# compare male size comps and index of abundance
#=============================================================
#==crabpack
#==indices of abundance

male_maturity_data <- crabpack::get_male_maturity(species = "SNOW",
                                                  region = "EBS",
                                                  channel = 'API')

male_snow <- crabpack::calc_bioabund(crab_data = specimen_data,
                                        species = "SNOW",
                                        region = "EBS",
                                        sex='male',
                                        shell_condition="all_categories",
                                        crab_category = c("large_male","small_male"),
                                        size_min = 25,
                                        bin_1mm = TRUE)

new_male_snow<-filter(male_snow,SHELL_TEXT=='new_hardshell')%>%
  group_by(YEAR,group = cut(SIZE_1MM , breaks = seq(0, max(SIZE_1MM), 5))) %>%
  summarise(n = sum(ABUNDANCE),
            data='crabpack')

midpts_cp<-seq(22.5,152.5,5)
tmp<-strsplit(as.character(new_male_snow$group),',')
new_male_snow$mid_pts<-NA
for(x in 1:length(tmp))
{
  in1<-as.numeric(substr(unlist(tmp[x])[1],2, nchar(unlist(tmp[x])[1])))
  in2<-as.numeric(substr(unlist(tmp[x])[2],1, nchar(unlist(tmp[x])[2])-1))
  new_male_snow$mid_pts[x]<-sum(in1,in2)/2
}

old_male_snow<-filter(male_snow,SHELL_TEXT%in%c('oldshell','very_oldshell'))%>%
  group_by(YEAR,group = cut(SIZE_1MM , breaks = seq(0, max(SIZE_1MM), 5))) %>%
  summarise(n = sum(ABUNDANCE),
            data='crabpack')

midpts_cp<-seq(22.5,152.5,5)
tmp<-strsplit(as.character(old_male_snow$group),',')
old_male_snow$mid_pts<-NA
for(x in 1:length(tmp))
{
  in1<-as.numeric(substr(unlist(tmp[x])[1],2, nchar(unlist(tmp[x])[1])))
  in2<-as.numeric(substr(unlist(tmp[x])[2],1, nchar(unlist(tmp[x])[2])-1))
  old_male_snow$mid_pts[x]<-sum(in1,in2)/2
}

#==translate to mature/immature based on 
male_maturity_data <- crabpack::get_male_maturity(species = "SNOW",
                                                  region = "EBS",
                                                  channel = 'API')


#==compare newshell/oldshell males

EBSdata_in<-read.csv("data/survey/EBSCrab_Abundance_Biomass_male.csv",skip=7)
EBSdata<-EBSdata_in[EBSdata_in$SIZE_CLASS_MM>=25 & EBSdata_in$SEX=="MALE" &EBSdata_in$SURVEY_YEAR>1981,]

male_total<- filter(EBSdata)%>%
  group_by(SURVEY_YEAR) %>%
  summarise(Totals = sum(ABUNDANCE))

old_shell<- filter(EBSdata,SHELL_CONDITION>=3)%>%
  group_by(SURVEY_YEAR) %>%
  summarise(Abundance = sum(ABUNDANCE))
old_shell$shell<-"Old shell"
new_shell<- filter(EBSdata,SHELL_CONDITION<3)%>%
  group_by(SURVEY_YEAR) %>%
  summarise(Abundance = sum(ABUNDANCE))
new_shell$shell<-"New shell"


Years				<-sort(unique(EBSdata$SURVEY_YEAR))
LengthBins	<-seq(25,135,5)
mid_pts<-seq(27.5,132.5,5)
MaleNew			<-matrix(nrow=length(Years),ncol=length(mid_pts))
MaleOld			<-matrix(nrow=length(Years),ncol=length(mid_pts))
MaleNewMature		<-matrix(nrow=length(Years),ncol=length(mid_pts))
MaleOldMature		<-matrix(nrow=length(Years),ncol=length(mid_pts))
MaleNewImmature		<-matrix(nrow=length(Years),ncol=length(mid_pts))
MaleOldImmature		<-matrix(0,nrow=length(Years),ncol=length(mid_pts))
NewShellIndex		<-3

for(x in 1:length(Years))
{
  temp<-EBSdata[EBSdata$SURVEY_YEAR==Years[x],]
  for(y in 1:(length(LengthBins)-1))
  {
    MaleNew[x,y]	<-sum(temp$ABUNDANCE[temp$SHELL_CONDITION<NewShellIndex & temp$SEX=="MALE"  & temp$SIZE_CLASS_MM>=LengthBins[y] & temp$SIZE_CLASS_MM<LengthBins[y+1] ])
    MaleOld[x,y]	<-sum(temp$ABUNDANCE[temp$SHELL_CONDITION>=NewShellIndex & temp$SEX=="MALE"  & temp$SIZE_CLASS_MM>=LengthBins[y] & temp$SIZE_CLASS_MM<LengthBins[y+1] ])
  }
}

rownames(MaleNew)<-Years
colnames(MaleNew)<-mid_pts
nal_new_mal_og<-melt(MaleNew)
nal_new_mal_og$data<-"OG"
colnames(nal_new_mal_og)<-c("YEAR","mid_pts","n",'data')

rownames(MaleOld)<-Years
colnames(MaleOld)<-mid_pts
nal_old_mal_og<-melt(MaleOld)
nal_old_mal_og$data<-"OG"
colnames(nal_old_mal_og)<-c("YEAR","mid_pts","n",'data')

old_male_snow
new_male_snow

comp_new_mal_plt<-rbind(old_male_snow[,c(1,5,3,4)],nal_old_mal_og)
comp_old_mal_plt<-rbind(new_male_snow[,c(1,5,3,4)],nal_new_mal_og)

new_male_plot_sc<-ggplot(comp_new_mal_plt)+
  geom_line(aes(x=mid_pts,y=n,col=data))+
  facet_wrap(~YEAR)+theme_bw()

old_male_plot_sc<-ggplot(comp_old_mal_plt)+
  geom_line(aes(x=mid_pts,y=n,col=data))+
  facet_wrap(~YEAR)+theme_bw()

comp_new_tots_pl_m<-filter(comp_new_mal_plt,mid_pts>25)%>%
  group_by(YEAR,data)%>%
  summarize(tot=sum(n))

comp_old_tots<-filter(comp_old_mal_plt,mid_pts>25)%>%
  group_by(YEAR,data)%>%
  summarize(tot=sum(n))

new_male_pl_tot<-ggplot(comp_new_tots_pl_m)+
  geom_line(aes(x=YEAR,y=tot,col=data))+
  theme_bw()
old_male_pl_tot<-ggplot(comp_old_tots)+
  geom_line(aes(x=YEAR,y=tot,col=data))+
  theme_bw()

png("plots/data_comp/new_male_plot_sc.png",height=8,width=8,res=400,units='in')
print(new_male_plot_sc)
dev.off()

png("plots/data_comp/old_male_plot_sc.png",height=8,width=8,res=400,units='in')
print(old_male_plot_sc)
dev.off()

png("plots/data_comp/new_male_pl_tot.png",height=8,width=8,res=400,units='in')
print(new_male_pl_tot)
dev.off()

png("plots/data_comp/old_male_pl_tot.png",height=8,width=8,res=400,units='in')
print(old_male_pl_tot)
dev.off()

#================================================
# compare mature biomasses that go into 