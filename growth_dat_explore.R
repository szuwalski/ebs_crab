library(ggplot2)
library(reshape2)
install.packages("fpc") 
library(fpc)

grow_dat<-read_excel("data/HistoricTagging_RKC_Bering.xlsx", sheet = "Sheet1")
grow_dat<-filter(grow_dat,!is.na(Growth_CL))

ggplot()+
  geom_point(data=grow_dat,aes(x=Rel_Length,y=Growth_CL))+
  theme_bw()

grow_dat$time_lib<-as.Date(grow_dat$Rec_Date,format="%m/%d/%Y")-as.Date(grow_dat$Rel_Date,format="%m/%d/%Y")
plot_lib<-as.numeric(grow_dat$time_lib)
plot_lib[plot_lib<0]<-0
plot_lib[plot_lib>3000]<-3000
hist(plot_lib)
grow_dat$plot_lib<-plot_lib
grow_dat$plot_lib_bin<-cut(plot_lib, breaks=seq(0,max(plot_lib,na.rm=T),365))

png("plots/bbrkc_growth.png",height=8,width=8,res=400,units='in')

ggplot()+
  geom_point(data=grow_dat,aes(x=Rel_Length,y=Growth_CL),alpha=0.1)+
  theme_bw()+facet_wrap(~plot_lib_bin)+
  geom_hline(yintercept=0,lty=2)

dev.off()

bins<-unique(grow_dat$plot_lib_bin)
keeper<-NULL
for(x in 1:(length(bins)-2))
{
  tmp<-filter(grow_dat,plot_lib_bin==bins[x])
  in_eps<-3.5
  if(x==4)
    in_eps<-3
  tmp$clust<-dbscan((tmp[,c(26,20)]),MinPts=20,eps=in_eps)$cluster
  # ggplot(tmp)+
  #   geom_point(aes(x=Rel_Length,y=Growth_CL,col=as.factor(clust)),alpha=.2)+
  #   theme_bw()
  tmp$time<-x
  keeper<-rbind(tmp,keeper)
}

png("plots/bbrkc_growth.png",height=8,width=8,res=400,units='in')
ggplot()+
  geom_point(data=keeper,aes(x=Rel_Length,y=Growth_CL,col=as.factor(clust)),alpha=0.1)+
  theme_bw()+facet_wrap(~time)+
  geom_hline(yintercept=0,lty=2)+
  theme(legend.position='none')
dev.off()

library(dplyr)
avg_molt_inc<-keeper%>%
  dplyr::group_by(clust,time)%>%
  dplyr::summarize(mean_inc=mean(Growth_CL),
            sd_ind=sd(Growth_CL))

avg_molt_inc$molt_n<-c(rep(NA,6),
                      0,2,1,0,1,2,
                      1,3,2,1,2,3)


