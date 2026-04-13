h_out<-c(1,1,0,1,0,0,0,1,0,0,1,1,0,0,0,0,0,1,0,0)
t_out<-c(1,1,0,1,0,1,1,1,0,0,1,1,0,0,1,0,0,0,1,0)

indat<-data.frame(outcome=c(h_out,t_out),start=c(rep('heads',length(h_out)),rep('tails',length(t_out))))

library(ggplot2)
ggplot(indat)+
  geom_histogram(aes(x=outcome,col=start,fill=start),position = "dodge",binwidth = 0.5)+
  theme_bw()

