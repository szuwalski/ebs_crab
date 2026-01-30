ret_dat<-read.csv("C:\\Users\\cody.szuwalski\\Work\\ebs_crab\\data\\tanner catch\\RetainedCatch.ABs.DirectedAndIncidental.csv")
tot_dat<-read.csv("C:\\Users\\cody.szuwalski\\Work\\ebs_crab\\data\\tanner catch\\TotalCatchABs.CrabFisheries.csv")

ret_sc<-read.csv("C:\\Users\\cody.szuwalski\\Work\\ebs_crab\\data\\tanner catch\\RetainedCAtchZCs.CrabFisheries.Scaled.Binned.csv")

library(dplyr)
library(ggplot2)
library(reshape2)

ggplot(filter(tot_dat,area=='all EBS'&sex=='male'))+
  geom_line(aes(x=year,y=abundance,col=fishery,group=fishery))

tot_dat_noshell<-tot_dat%>%
  group_by(year,sex,fishery,area)%>%
  summarize(tot_n=sum(abundance))

ggplot(filter(tot_dat_noshell,area=='all EBS'&sex=='male'))+
  geom_line(aes(x=year,y=tot_n,col=fishery,group=fishery))


all_ret<-filter(ret_dat,area=='all EBS')%>%
  group_by(year)%>%
  summarize(tot_n=sum(abundance))

all_tot<-filter(tot_dat,sex=="male")%>%
  group_by(year)%>%
  summarize(tot_n=sum(abundance))

#==are the sums of total TCF same as sums of ret?



