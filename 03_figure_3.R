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


#==color coordinate the estimates with figure 1
rep_files<-c("/models/snow/snow_down.rep",
             "/models/tanner/tanner.rep")

species<-c("Snow", "Tanner")
outs_in<-list(list())

for(x in 1:length(rep_files))
  outs_in[[x]]<-readList(paste(getwd(),rep_files[x],sep=""))

#==need a flag for the type of mortality to split for snow + tanner?
#==plot snow and tanner together and then the kind crabs together?
#==snow and tanner with recruitment and immature mortality; fishing and mature mortality?
all_dat_imm<-NULL
for(x in 1:length(outs_in))
{
  years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  plot_dat<-data.frame(values=c(outs_in[[x]]$recruits/max(outs_in[[x]]$recruits),
                                outs_in[[x]]$`natural mortality`[,1]),
                       Year=c(years-5,rep(years,1)),
                       process=c(rep("Recruitment",length(years)),
                                 rep("Other mortality (imm)",length(years))))
  plot_dat$species<-species[x]
  all_dat_imm<-rbind(all_dat_imm,plot_dat)                
}


all_dat_imm$process<-as.character(all_dat_imm$process)
all_dat_imm$process<-factor(all_dat_imm$process, levels=c("Recruitment", "Other mortality (imm)"))

all_dat_mat<-NULL
for(x in 1:length(outs_in))
{
  years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  plot_dat<-data.frame(values=c(outs_in[[x]]$`mature natural mortality`[,1],
                                outs_in[[x]]$est_fishing_mort),
                       Year=c(rep(years,2)),
                       process=c(rep("Other mortality (mat)",length(years)),
                                 rep("Fishing mortality",length(years))))
  plot_dat$species<-species[x]
  all_dat_mat<-rbind(all_dat_mat,plot_dat)                
}
in_col<-c("#3da550","#ff8f38","#3da550","#ff8f38")
all_dat_mat$process<-as.character(all_dat_mat$process)
all_dat_mat$process<-factor(all_dat_mat$process, levels=c("Other mortality (mat)","Fishing mortality"))
all_dat_mat$values[all_dat_mat$process=='Fishing mortality' & all_dat_mat$values==0]<-NA
imm_proc<-ggplot()+
  geom_line(data=all_dat_imm,aes(x=Year,y=values,col=species),lwd=1.2)+
  facet_grid(rows=vars(process),cols=vars(species),scales='free_y')+
  theme_bw()+ylab("")+  scale_color_manual(values=in_col)+
  theme(legend.position='none',
        axis.line = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

mat_proc<-ggplot()+
  geom_line(data=all_dat_mat,aes(x=Year,y=values,col=species),lwd=1.2)+
  facet_grid(rows=vars(process),cols=vars(species),scales='free_y')+
  theme_bw()+ylab("")+  scale_color_manual(values=in_col)+
  theme(legend.position='none',
        axis.line = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

chion_proc_agg<-ggplot()+
  geom_line(data=rbind(all_dat_mat,all_dat_imm),aes(x=Year,y=values,col=species),lwd=1.2)+
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

png("plots/fig_agg_chion.png",height=8,width=5,res=400,units='in')
print(chion_proc_agg)
dev.off()


#==need to put the CVs in the .REP files and pull here
#==immature indices
div_n<-c(max(outs_in[[1]]$imm_n_obs),max(outs_in[[2]]$imm_n_obs))
ind_dat_imm<-NULL

for(x in 1:length(outs_in))
{
  years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  df_1<-data.frame(pred=outs_in[[x]]$imm_numbers_pred/div_n[x],
                   obs=outs_in[[x]]$imm_n_obs/div_n[x],
                   year=years,
                   ci_dn=(outs_in[[x]]$imm_n_obs/div_n[x]) /  exp(1.96*sqrt(log(1+0.15^2))),
                   ci_up=(outs_in[[x]]$imm_n_obs/div_n[x]) *  exp(1.96*sqrt(log(1+0.15^2))))
  
  df_1$species<-species[x]
  df_1$color<-in_col[x]
  ind_dat_imm<-rbind(ind_dat_imm,df_1)
}


imm_abnd<-ggplot(data=ind_dat_imm)+
  geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=pred,col=color),lwd=1.5,alpha=.8)+
  theme_bw()+
  scale_color_manual(values=in_col)+
  ylab("Relative abundance")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  facet_grid(~species,scales='free_y')+xlim(1970,2023)+
  xlab("")+ guides(color="none")+ggtitle("IMMATURE")


#==mature indices
div_n<-c(max(outs_in[[1]]$mat_n_obs),max(outs_in[[2]]$mat_n_obs))
ind_dat_mat<-NULL
for(x in 1:length(outs_in))
{
  years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  df_1<-data.frame(pred=outs_in[[x]]$mat_numbers_pred/div_n[x],
                   obs=outs_in[[x]]$mat_n_obs/div_n[x],
                   year=years,
                   ci_dn=(outs_in[[x]]$mat_n_obs/div_n[x]) /  exp(1.96*sqrt(log(1+0.15^2))),
                   ci_up=(outs_in[[x]]$mat_n_obs/div_n[x]) *  exp(1.96*sqrt(log(1+0.15^2))))
  
  df_1$species<-species[x]
  df_1$color<-in_col[x]
  ind_dat_mat<-rbind(ind_dat_mat,df_1)
}

mat_abnd<-ggplot(data=ind_dat_mat)+
  geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=pred,col=color),lwd=1.5,alpha=.8)+
  theme_bw()+
  scale_color_manual(values=in_col)+
  ylab("Relative abundance")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),axis.line.y=element_blank())+
  facet_grid(~species,scales='free_y')+xlim(1970,2023)+
  xlab("")+ guides(color="none")+ggtitle("MATURE")

png("plots/fig_3.png",height=6,width=9,res=400,units='in')
(imm_abnd/imm_proc) | (mat_abnd/mat_proc)
dev.off()


#============================================
# correlations between time series
# need to pull spawning biomass in there too
#============================================
div_n_imm<-c(max(apply(outs_in[[1]]$pred_imm_pop_num,1,sum)),max(apply(outs_in[[2]]$pred_imm_pop_num,1,sum)))
div_n_mat<-c(max(apply(outs_in[[1]]$pred_mat_pop_num,1,sum)),max(apply(outs_in[[2]]$pred_mat_pop_num,1,sum)))

all_dat<-NULL
for(x in 1:length(outs_in))
{
  years<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  plot_dat<-data.frame(values=c(outs_in[[x]]$recruits/max(outs_in[[x]]$recruits),
                                outs_in[[x]]$`natural mortality`[,1],
                                outs_in[[x]]$`mature natural mortality`[,1],
                                outs_in[[x]]$est_fishing_mort,
                                apply(outs_in[[x]]$pred_imm_pop_num,1,sum)/div_n_imm[x],
                                apply(outs_in[[x]]$pred_mat_pop_num,1,sum)/div_n_mat[x]),
                       Year=c(years-5,rep(years,5)),
                       process=c(rep("Rec",length(years)),
                                 rep("M_imm",length(years)),
                                 rep("M_mat",length(years)),
                                 rep("F",length(years)),
                                 rep("N_imm",length(years)),
                                 rep("N_mat",length(years))))
  plot_dat$species<-species[x]
  all_dat<-rbind(all_dat,plot_dat)                
}

library(GGally)
casted<-dcast(all_dat,Year~species+process,value.var="values")[,-1]
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p+theme_bw()
}

p1 = ggpairs(casted, lower = list(continuous = my_fn))


#============================================
# SR relationsihp
#============================================

chion_srr_dat<-dcast(filter(all_dat,process%in%c("Rec","N_mat")),species+Year~process,value.var="values")
chion_srr<-ggplot()+
  geom_point(data=chion_srr_dat,aes(x=N_mat,y=Rec,col=species),size=2)+
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

png("plots/chion_srr.png",height=4,width=6,res=400,units='in')
print(chion_srr)
dev.off()

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

png("plots/chion_cors.png",height=13,width=13,res=400,units='in')
print(p1)
dev.off()