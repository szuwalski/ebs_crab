
library(maps)
library("rnaturalearth")
library(interp)
library(RColorBrewer)
library(reshape2) # for melt
library(mgcv)  
library(PBSmapping)
library(mapdata)    #some additional hires data
library(maptools)   #useful tools such as reading shapefiles
library(mapproj)
library(ggplot2)
library(patchwork)
library(dplyr)
library(reshape2)
library(ggridges)

#==shared map data
world <- ne_countries(scale = "medium", returnclass = "sf")
lon_1<- -179
lon_2<- -155
lat_1<- 53
lat_2<- 62.2
nmiSurv<-140350 # surveyed area

#==============================
# EBS data
#==============================
opilio<-read.csv("C:/Users/cody.szuwalski/Work/snow_index/data/EBSCrab_Haul.csv",header=T,skip=5)
opilio<-filter(opilio,AKFIN_SURVEY_YEAR>1981)
years<-unique(opilio$AKFIN_SURVEY_YEAR)
bins<-seq(25,135,5)
mid_pts<-seq(27.5,132.5,5)
loc_dat<-filter(opilio,HAUL_TYPE==3&nchar(GIS_STATION)<5)%>%
  group_by(GIS_STATION)%>%
  summarize(mid_lat=mean(MID_LATITUDE,na.rm=T),
            mid_lon=mean(MID_LONGITUDE,na.rm=T))
stations<-unique(loc_dat$GIS_STATION)

png("plots/stations_location.png",height=8,width=10,res=400,units='in')
ggplot(loc_dat)+
  geom_text(aes(x=mid_lon,y=mid_lat,label=GIS_STATION),size=3)
dev.off()

bbrkc<-read.csv("C:/Users/cody.szuwalski/Work/bbrkc_down/data/EBSCrab_Haul.csv",header=T,skip=5)
tanner<-read.csv("C:/Users/cody.szuwalski/Work/tanner_simple/data/EBSCrab_Haul.csv",header=T,skip=5)
bkc<-read.csv("C:/Users/cody.szuwalski/Work/smbkc/data/EBSCrab_Haul.csv",header=T,skip=5)

#=======================================
# plot maps of distribution by species
#=======================================
snow_samp<-opilio%>%
  group_by(GIS_STATION)%>%
  summarize(avg_dens=sum(SAMPLING_FACTOR,na.rm=T))
snow_dat<-merge(snow_samp,loc_dat,by='GIS_STATION')
snow_dat<-filter(snow_dat,avg_dens!="NaN")

snow_dat[snow_dat$avg_dens>15000,2]<-15000
snow_dat$in_alpha<-snow_dat$avg_dens/max(snow_dat$avg_dens)
use_snow<-filter(snow_dat,in_alpha>0.1)
use_snow$species<-"Snow"
use_snow$color<-"#3da550"
use_snow$plot<-1

#===BBRKC   
  bbrkc_samp<-bbrkc%>%
    group_by(GIS_STATION)%>%
    summarize(avg_dens=sum(SAMPLING_FACTOR,na.rm=T))
  bbrkc_dat<-merge(bbrkc_samp,loc_dat,by='GIS_STATION')
  bbrkc_dat<-filter(bbrkc_dat,avg_dens>0)
  bbrkc_dat[bbrkc_dat$avg_dens>1000,2]<-1000
  bbrkc_dat$in_alpha<-bbrkc_dat$avg_dens/max(bbrkc_dat$avg_dens)
  use_bbrkc<-filter(bbrkc_dat,in_alpha>0.01)
  use_bbrkc$species<-"Red king"
  use_bbrkc$color<-"#ff5050"
  use_bbrkc$plot<-1
    
  #===tanner   
  tan_samp<-tanner%>%
    group_by(GIS_STATION)%>%
    summarize(avg_dens=sum(SAMPLING_FACTOR,na.rm=T))
  tan_dat<-merge(tan_samp,loc_dat,by='GIS_STATION')
  tan_dat<-filter(tan_dat,avg_dens!="NaN")
  
  tan_dat[tan_dat$avg_dens>6000,2]<-6000
  tan_dat$in_alpha<-tan_dat$avg_dens/max(tan_dat$avg_dens)
  use_tan<-filter(tan_dat,in_alpha>0.1)
  use_tan$species<-"Tanner"
  use_tan$color<-"#ff8f38"
  use_tan$plot<-2   
  
  #===bkc 
  bkc_samp<-bkc%>%
    group_by(GIS_STATION)%>%
    summarize(avg_dens=sum(SAMPLING_FACTOR,na.rm=T))
  bkc_dat<-merge(bkc_samp,loc_dat,by='GIS_STATION')
  bkc_dat<-filter(bkc_dat,avg_dens!="NaN")
  
  bkc_dat[bkc_dat$avg_dens>100,2]<-100
  bkc_dat$in_alpha<-bkc_dat$avg_dens/max(bkc_dat$avg_dens)
  use_bkc<-filter(bkc_dat,in_alpha>0.1)
  use_bkc$species<-"Blue king"
  use_bkc$color<-"#0034c3"
  use_bkc$plot<-2

in_dat<-rbind(use_snow,use_bbrkc,use_tan,use_bkc)


p<-ggplot() + 
  geom_tile(data=in_dat, aes(x = mid_lon, y = mid_lat, fill = species,alpha=in_alpha),width=.5,height=.25) +
  scale_fill_manual(values=c("#0034c3","#ff5050","#3da550","#ff8f38"))+
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  geom_point(data=loc_dat,aes(x=mid_lon,y=mid_lat),size=.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA),
        legend.position='none')+
  labs(fill="Species")+
  guides(alpha='none')+
  facet_wrap(~plot,ncol=1)+ylab("Latitude")+xlab("Longitude")

png("plots/species_dist.png",height=8,width=6,res=400,units='in')
print(p)
dev.off()

#==add stations with crab and total crab observed
pres<-opilio%>%
  group_by(AKFIN_SURVEY_YEAR,GIS_STATION)%>%
  summarize(presence=sum(!is.na(SAMPLING_FACTOR)))
snow_st<-pres%>%
  group_by(AKFIN_SURVEY_YEAR)%>%
  summarize(stations=sum(presence>0))




#==add the catch and value from each species as stacked barplots
ref<-read.csv('data/Assessment_Summary_Data.csv')
ref<-ref[-c(3,5,10),]
catches<-read.csv("data/catches.csv")
catches<-catches[,-c(4,6,11)]
colnames(catches)[2:ncol(catches)]<-ref$Stock.Name
melted<-melt(catches,id.var="Year")
melted$speces<-NULL
colnames(melted)[2]<-"Stock"
tmmp<-strsplit(as.character(melted$Stock),split="-")
for(x in 1:nrow(melted))
  melted$species[x]<-tmmp[[x]][1]

melted$species[melted$species=="Snow crab " | melted$species=="Southern Tanner crab "]<-"Snow/Tanner crab"
stocks<-unique(melted$Stock)
catch_plot<-filter(melted,Stock%in%stocks[c(1,2,5,7,8,9)])

png("plots/catches.png",height=6,width=8,res=350,units='in') 
catch_out<-ggplot(catch_plot)+
  geom_area(aes(x=Year,y=value,fill=Stock),lwd=2)+
  theme_bw()+
  ylab("Catch (metric tons)")+
  theme(legend.position=c(.75,.72),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  xlim(1975,2024)+
  scale_y_continuous(position = "right")+
  scale_fill_manual(values=c("#0034c377","#0034c3","#ff5050","#ff505077","#3da550","#ff8f38"))
dev.off()

cat_val<-read.csv("data/FOSS_landings.csv")

val_out<-ggplot(cat_val)+
  geom_area(aes(x=Year,y=Dollars/10000000,fill=NMFS.Name),lwd=2)+
  theme_bw()+
  ylab("Value ($US 10,000,000)")+
  theme(legend.position='none')+
  xlim(1975,2024)+
  scale_y_continuous(position = "right")+
  scale_fill_manual(values=c("#ff5050","#3da550","#ff8f38"))




png("plots/fig_1.png",height=8,width=12,res=400,units='in')
print(p)+(catch_out/val_out)
dev.off()
