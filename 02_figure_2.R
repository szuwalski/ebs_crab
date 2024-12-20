#==pull in output from models
#==pull in pictures
#==plot time series with uncertainty
#==call out important periods in time serieslibrary(dplyr)
library(reshape2)
library(ggplot2)
library(mgcv)  
library(dplyr)
library(ggridges)
library(png)
library(PBSmodelling)
library(patchwork)

annotation_custom2 <-   function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf,ymax = Inf, data){ layer(data = data, 
                                                                                                      stat = StatIdentity, 
                                                                                                      position = PositionIdentity,
                                                                                                      geom = ggplot2:::GeomCustomAnn,
                                                                                                      inherit.aes = TRUE, 
                                                                                                      params = list(grob = grob,xmin = xmin, xmax = xmax,
                                                                                                                    ymin = ymin, ymax = ymax))}


#==add PIBKC, PIRKC
rep_files<-c("/models/bbrkc/rkc.rep",
             "/models/pirkc/rkc.rep",
             "/models/smbkc/bkc.rep",
             "/models/pibkc/bkc.rep")

species<-c("BBRKC","PIRKC","SMBKC","PIBKC")
maturity<-c(120,120,105,120)
outs_in<-list(list())

for(x in 1:length(rep_files))
  outs_in[[x]]<-readList(paste(getwd(),rep_files[x],sep=""))

#==need a flag for the type of mortality to split for snow + tanner?
#==plot snow and tanner together and then the kind crabs together?
#==snow and tanner with recruitment and immature mortality; fishing and mature mortality?
all_dat<-NULL
for(x in 1:length(outs_in))
{
years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)

#names(outs_in[[x]])
outs_in[[x]]$pred_pop_num
tmp_ssn<-apply(outs_in[[x]]$pred_pop_num[,which(outs_in[[x]]$sizes>maturity[x])],1,sum)
tmp_ssn<-tmp_ssn/max(tmp_ssn)
plot_dat<-data.frame(values=c(outs_in[[x]]$recruits/max(outs_in[[x]]$recruits),
                              outs_in[[x]]$`natural mortality`[,1],
                              outs_in[[x]]$est_fishing_mort,
                              tmp_ssn),
                     Year=c(years-5,rep(years,3)),
                     process=c(rep("Recruitment",length(years)),
                               rep("Other mortality",length(years)),
                               rep("Fishing mortality",length(years)),
                               rep("Spawner abundance",length(years))))
plot_dat$species<-species[x]
    all_dat<-rbind(all_dat,plot_dat)                
}
in_col<-c("#ff5050","#0034c377","#ff505077","#0034c3")
all_dat$values[all_dat$process=='Fishing mortality' & all_dat$values==0]<-NA
king_proc<-ggplot()+
  geom_line(data=filter(all_dat,process!="Spawner abundance"),aes(x=Year,y=values,col=species),lwd=1.2)+
  facet_grid(rows=vars(process),cols=vars(species),scales='free_y')+
  theme_bw()+ylab("")+
  scale_color_manual(values=in_col)+
  theme(legend.position='none',
        axis.line = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

king_proc_agg<-ggplot()+
  geom_line(data=all_dat,aes(x=Year,y=values,col=species),lwd=1.2)+
  facet_grid(rows=vars(process),scales='free_y')+
  theme_bw()+ylab("")+
  scale_color_manual(values=in_col)+
  theme(legend.position='none',
        axis.line = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

png("plots/fig_agg_king.png",height=8,width=5,res=400,units='in')
print(king_proc_agg)
dev.off()


##==cross correlation within a stock of processes
library(forecast)
in_max_lag<-8
tmp<-filter(all_dat,species=="BBRKC")
ggCcf(filter(tmp,process=="Recruitment")$values,
      filter(tmp,process=="Spawner abundance")$values,lwd=2,lag.max=in_max_lag)+theme_bw()
ggCcf(filter(tmp,process=="Other mortality")$values,
      filter(tmp,process=="Recruitment")$values,lwd=2,lag.max=in_max_lag)+theme_bw()


#==need to put the CVs in the .REP files and pull here
div_n<-c(max(outs_in[[1]]$n_obs),max(outs_in[[2]]$n_obs),max(outs_in[[3]]$n_obs),max(outs_in[[4]]$n_obs))
ind_dat<-NULL

for(x in 1:length(outs_in))
{
years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
df_1<-data.frame(pred=outs_in[[x]]$numbers_pred/div_n[x],
                 obs=outs_in[[x]]$n_obs/div_n[x],
                 year=years,
                 ci_dn=(outs_in[[x]]$n_obs/div_n[x]) /  exp(1.96*sqrt(log(1+outs_in[[x]]$sigma_numbers^2))),
                 ci_up=(outs_in[[x]]$n_obs/div_n[x]) *  exp(1.96*sqrt(log(1+outs_in[[x]]$sigma_numbers^2))))

df_1$species<-species[x]
df_1$color<-in_col[x]
ind_dat<-rbind(ind_dat,df_1)
}

ind_dat$ci_up[ind_dat$ci_up>1.5]<-1.5

king_abnd<-ggplot(data=ind_dat)+
  geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=pred,col=species),lwd=1.5)+
  theme_bw()+
  scale_color_manual(values=in_col)+
  ylab("Relative abundance")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  facet_grid(~species,scales='free_y')+xlim(1970,2023)+
  xlab("")+ guides(color="none")

png("plots/fig_2.png",height=6,width=8,res=400,units='in')
king_abnd/king_proc + plot_layout(heights = c(2, 3))
dev.off()


#=========================
# plot spawner abundance vs recruits
#===============================
srr_dat<-filter(all_dat,process%in%c("Recruitment","Spawner abundance"))
in_dat<-data.frame(dcast(data=srr_dat,formula=Year+species~process , value.var = "values"))

king_srr<-ggplot()+
  geom_point(data=in_dat,aes(x=Spawner.abundance,y=Recruitment,col=species),size=2)+
  facet_wrap(~species)+
  theme_bw()+ylab("")+
  scale_color_manual(values=in_col)+
  theme(legend.position='none',
        axis.line = element_line(colour = "black"),
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ylab("Recruits")+xlab("Relative spawner abundance")+
  expand_limits(x=0)

png("plots/king_srr.png",height=6,width=6,res=400,units='in')
print(king_srr)
dev.off()

#============================================
# need to pull spawning biomass in there too
#============================================
all_dat<-NULL
for(x in 1:length(outs_in))
{
  
  years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  plot_dat<-data.frame(values=c(outs_in[[x]]$recruits/max(outs_in[[x]]$recruits),
                                outs_in[[x]]$`natural mortality`[,1],
                                outs_in[[x]]$est_fishing_mort,
                                pred=outs_in[[x]]$numbers_pred/div_n[x]),
                       Year=c(years-5,rep(years,3)),
                       process=c(rep("Recruitment",length(years)),
                                 rep("Other mortality",length(years)),
                                 rep("Fishing mortality",length(years)),
                                 rep("Abundance",length(years))))
  plot_dat$species<-species[x]
  all_dat<-rbind(all_dat,plot_dat)                
}

library(GGally)
all_dat$species_process<-paste(all_dat$species,"_",substring(all_dat$process,1,1),sep="")
casted<-dcast(all_dat,Year~species_process,value.var="values")[,-1]
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p+theme_bw()
}

p1 = ggpairs(casted, lower = list(continuous = my_fn))


# Correlation matrix plot
p2 <- ggcorr(casted, label = TRUE, label_round = 2)
g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
p<-ncol(casted)
for (k1 in 1:(p-1)) {
  for (k2 in (k1+1):p) {
    plt <- getPlot(p1,k1,k2) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    p1 <- putPlot(p1,plt,k1,k2)
    idx <- idx+1
  }
}

png("plots/rkc_cors.png",height=13,width=13,res=400,units='in')
print(p1)
dev.off()





casted<-dcast(filter(all_dat,process%in%c("Abundance","Recruitment",),Year~species_process,value.var="values"))
p1<-ggpairs(casted[,-1])+ theme_grey(base_size = 8)
p2 <- ggcorr(casted[,-1], label = TRUE)

abund<-filter(all_dat,process%in%c("Abundance"))[,c(1,2,4)]
colnames(abund)<-c("Abundance","Year","Species")
rec<-filter(all_dat,process%in%c("Recruitment"))[,c(1,2,4)]
colnames(rec)<-c("Recruitment","Year","Species")
plot_srr<-merge(rec,abund,by=c("Year"))
plot_srr <- rec %>% right_join(abund, by=c("Year","Species"))
ggplot(plot_srr)+
  geom_point(aes(x=Abundance,y=Recruitment))+
  geom_smooth(aes(x=Abundance,y=Recruitment),span=1)+
  theme_bw()+facet_wrap(~Species)

