#==read in par and input
library(reshape2)
library(ggplot2)
library(mgcv)  
library(dplyr)
library(ggridges)
library(png)
library(PBSmodelling)
library(patchwork)


#==QUESTION
#==what would the historical trajectory look like with no fishing
#==constant M, varying M, density dependent M

rep_files<-c("/models/bbrkc/rkc.rep")
maturity<-c(120)
outs_in<-list(list())

for(x in 1:length(rep_files))
  outs_in[[x]]<-readList(paste(getwd(),rep_files[x],sep=""))

names(outs_in[[1]])
in_outs<-outs_in[[1]]
pop_dy<-function(in_outs,
                 in_f,
                 proj_yr,
                 mort='constant',
                 disc_surv=0.80,
                 steepness,
                 srr_r0,
                 srr_s0,
                 n0,
                 min_m=0.1,
                 max_m=0.5,
                 rec_type="obs",
                 rec_min_yr=NA,
                 rec_max_yr=NA)
{

  yr_sq<-seq(in_outs$styr,in_outs$endyr)
  years<-length(yr_sq)
  n_at_size<-matrix(ncol=length(in_outs$sizes),nrow=years)
  n_at_size[1,]<-in_outs$"numbers_pred"[1]*in_outs$"pred numbers at size"[1,]
  MMB<-rep(0,years)
  vuln_n<-rep(0,years)
  if(rec_type=="obs")
   recruits<-in_outs$"recruits"
  if(rec_type=="mean")
   recruits<-rep(mean(in_outs$"recruits"[match(rec_min_yr,yr_sq):match(rec_max_yr,yr_sq)]),length(in_outs$"recruits"))
  
  ret_catch<-matrix(ncol=length(in_outs$sizes),nrow=years)
  
  #==fill nat M
  if(mort=='constant')
    nat_m<-matrix(median(in_outs$`natural mortality`),ncol=length(in_outs$sizes),nrow=years) 
 
  if(mort=='historical')
    nat_m<-in_outs$`natural mortality`

  if(mort=='sample')
  {
    nat_m<-matrix(ncol=length(in_outs$sizes),nrow=years)
    for(x in 1:years)
      nat_m[x,]<-sample(in_outs$`natural mortality`,size=1)
  }
  
  if(mort=='density')
    nat_m<-matrix(median(in_outs$`natural mortality`),ncol=length(in_outs$sizes),nrow=years)   
  
    f_mort<-rep(in_f,years)
  if(in_f=="est")
    f_mort<-in_outs$est_fishing_mort
  if(in_f=='no_fish')
    f_mort<-rep(0,years)
  
  for(x in 1:(years-1))
  {
    tmp_n <- n_at_size[x,]*exp(-0.17*nat_m[x,]) 
    
    vuln_n[x]<-sum(tmp_n*in_outs$ret_fish_sel)
    # fishery
    tmp_catch <- tmp_n * (1-exp(-f_mort[x]*in_outs$total_fish_sel))
    ret_catch[x,] <- tmp_catch*in_outs$ret_fish_sel
    tmp_n <- tmp_n *exp(-f_mort[x]*in_outs$total_fish_sel)
    tmp_n <- tmp_n + (tmp_catch*(1-in_outs$ret_fish_sel))*disc_surv
    
    # growth
    trans_n <- c(in_outs$size_trans %*% (tmp_n*in_outs$prob_molt))
    tmp_n <- trans_n + tmp_n*(1-in_outs$prob_molt)
    
    # recruits
    MMB[x]<-sum(tmp_n[in_outs$sizes>120])
    # num   <- 4*srr_r0*steepness*MMB[x]
    # denom <- (1-steepness)*srr_r0*(srr_s0/srr_r0)+(5*steepness-1)*MMB[x]
    # recruits[x]<- num/denom 
    
    tmp_n[1]<-0.3*recruits[x]
    tmp_n[2]<-0.5*recruits[x]
    tmp_n[3]<-0.2*recruits[x]
    
    n_at_size[x+1,] <- tmp_n*exp(-0.83*nat_m[x,]) 
  }
  list(MMB=MMB,recruits=recruits,yield=apply(ret_catch,1,sum),
       nat_m=nat_m,abund=apply(n_at_size,1,sum),
       years=seq(in_outs$styr,in_outs$endyr),
       vuln_n=vuln_n)
}


no_fish<-pop_dy(in_outs=outs_in[[1]],
                in_f=0,
                mort="historical",
                disc_surv=0.80,
                steepness=1,
                srr_r0=100000,
                srr_s0=100)

obs_fish<-pop_dy(in_outs=outs_in[[1]],
                in_f='est',
                mort="historical",
                disc_surv=0.80,
                steepness=1,
                srr_r0=100000,
                srr_s0=100)

no_fish_m<-pop_dy(in_outs=outs_in[[1]],
                in_f=0,
                mort="constant",
                disc_surv=0.80,
                steepness=1,
                srr_r0=100000,
                srr_s0=100)

obs_fish_m<-pop_dy(in_outs=outs_in[[1]],
                 in_f='est',
                 mort="constant",
                 disc_surv=0.80,
                 steepness=1,
                 srr_r0=100000,
                 srr_s0=100)

plotter<-data.frame(abund=c(no_fish$abund,
                   obs_fish$abund,
                   no_fish_m$abund,
                   obs_fish_m$abund,
                   no_fish$vuln_n,
                   obs_fish$vuln_n,
                   no_fish_m$vuln_n,
                   obs_fish_m$vuln_n),
           year=rep(obs_fish$years,8),
           fishing=rep(c(rep("No fishing",length(obs_fish$years)),
                      rep("Historical fishing mortality",length(obs_fish$years))),4),
           quantity=rep(c(rep("Abundance",4*length(obs_fish$years)),
                        rep("Vulnerable N",4*length(obs_fish$years)))),
           mortality=rep(c(rep("Variable mortality",2*length(obs_fish$years)),
                         rep("Constant mortality",2*length(obs_fish$years))),2))

png("plots/alt_histories.png",height=5,width=7,res=350,units='in') 
ggplot(plotter)+
  geom_line(aes(x=year,y=abund,col=fishing),lwd=1.2)+
  facet_grid(quantity~mortality,scales='free_y')+theme_bw()+
  theme(legend.position=c(.8,.8))
dev.off()


