library(ggplot2)
library(Rmisc)
library(reshape2)
library(dplyr)
library(pracma)
library(ggbeeswarm)
library(dabestr)
library(tidyr)
library(data.table)
library(Cairo)
library(ggthemes)
library(DescTools)

setwd('C:/Users/optim/Desktop/ArrestManuscript/OngoingFigureData/')

# Figure 7G - Repeat testing
setwd('D:/Imaging/03132019_xa266_gcamp6s')

df <- read.csv2('Caimaging_nobase_normd.csv',header=T,sep='\t')


L = df$Cell == 7

df6 <- t(df[L,])

time <- seq(0,60,1)
groups <- c(rep("Flow1",61),rep('Flow2',61),rep('NoStim',61))

dat <- data.frame(rep(time,3),groups,as.numeric(df6[-1:-2,]))
colnames(dat) <- c("time","groups","Fl")

ggplot(dat,aes(x=time,colour=groups,y=Fl))+
  geom_line(stat="identity",show.legend=F,size=1)+
  xlab("Time (s)") +
  ylab("GCaMP6s fluorescence") +
  scale_x_continuous(breaks=seq(0,60,10),limits=c(0,61))+
  #scale_fill_manual(values=c("gray20","red","cyan"))+
  ylim(0,3) + 
  theme(text=element_text(size=16,family="Serif")) + 
  theme_classic()

ggsave("xa266_gcamp6s_neuron10.eps", plot = last_plot(), device = "eps",
       scale = 1,width=10,height=10,units="cm",
       dpi = 300, limitsize = TRUE)


#################################################################################################################
######################### xa266 masked pERK/tERK labeling #######################################################

df <- read.table('maskvals.csv',header=T)
my.data <- as_tibble(df)

xa334.shared <- 
  df %>%
  dabest(Group,PT,
         idx = c("C","T"),
         paired=F)

xa334.shared.effsize <- cohens_d(xa334.shared)

plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="pERK/tERK ratio",
     rawplot.ylim = c(0,1.5),
     effsize.ylabel="Effect size")

svg(filename="Fig7C_ptERK_estimplot.svg",
    width=4,
    height=5,
    pointsize=14)
plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="pERK/tERK ratio",
     rawplot.ylim = c(0,1.5),
     effsize.ylabel="Effect size")
dev.off()


## Fig 6E - xa334 calcium imaging heatmap 

responsives <- read.table('Fig6E_calciumresponders.csv',sep=",",header=T)

time <- responsives$time
cellnum <- rep(seq(1,ncol(responsives)-1),length(time))
cellnum <- sort(cellnum,decreasing=F)

dffmtx <- melt(as.matrix(responsives[,-1]))
dff <- dffmtx$value

timevar <- rep(time,ncol(responsives)-1)

caresp <- data.frame(cellnum,timevar,dff)

ggplot(caresp,aes(x=timevar,y=cellnum))+
  geom_raster(aes(fill=dff),interpolate=F)+
  scale_fill_viridis_c(option = "magma",limits = c(0,max(caresp$dff)),
                       direction = 1)+
  theme_minimal()

ggsave("Fig6E_responders.svg",width=30,height=5)



nonresponsives <- read.table('Fig6E_calciumnonresponders.csv',sep="\t",header=T)

time <- nonresponsives$time
cellnum <- rep(seq(1,ncol(nonresponsives)-1),length(time))
cellnum <- sort(cellnum,decreasing=F)

dffmtx_non <- melt(as.matrix(nonresponsives[,-1]))
dff_non <- dffmtx_non$value

timevar <- rep(time,ncol(nonresponsives)-1)

caresp <- data.frame(cellnum,timevar,dff_non)

ggplot(caresp,aes(x=timevar,y=cellnum))+
  geom_raster(aes(fill=dff_non),interpolate=F)+
  scale_fill_viridis_c(option = "magma",limits = c(0,max(caresp$dff)),
                       direction = 1)

ggsave("Fig6E_nonresponders.svg",width=30,height=20)



data <- read.table('Fig7I_y334cochr.csv',header=T,sep=',')

dspeed <- 0.1*(data$Blk24 - data$baseline)
dspeed <- 1 - ((data$baseline - data$Blk24)/data$baseline)

df <- data.frame(data$Group,data$baseline,dspeed)
colnames(df) <- c('group','baseline','dspeed')
df <- as_tibble(df)


ggplot(df,aes(x=group,y=dspeed))+
  stat_summary(fun.data = quantiles_95, geom="boxplot",fill=c("grey","red"),weight=4)+
  geom_beeswarm(color="grey20",size=3,cex=2)+
  #scale_x_discrete(labels = c("a0p5" = "14","a1p0" = "20","a1p5" = "23","a2p0" = "26"))+
  ylab("Speed re. baseline")+
  xlab("")+
  ylim(0,2)+
  theme_classic()+
  theme(axis.text=element_text(size=18,colour='black',family="Arial",face='plain'))+
  theme(axis.title = element_text(size = 20,colour='black',family="Arial",face='plain'))

ggsave('xa334_cochr_boxplot.svg',width=4,height=7)

t <- aov(df$dspeed ~ df$group)
summary(t)

xa334.shared <- 
  df %>%
  dabest(group,dspeed,
         idx = c("c","r"),
         paired=F)

xa334.shared.effsize <- cohens_d(xa334.shared)

plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Speed (re. baseline)",
     rawplot.ylim = c(0,2),
     effsize.ylabel="Effect size")


svg(filename="Fig7I_xa334cochr_estimplot.svg",
    width=4,
    height=5,
    pointsize=14)
plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Speed (re. baseline)",
     rawplot.ylim = c(0,2),
     effsize.ylabel="Effect size")
dev.off()

# Figure 7I - xa334cochr time series

data <- read.table('Fig7I_xa334cochr_timeseries.csv',header=F,sep=',')

dat <- as.matrix(data[,2:49])
dat <- t(dat)

time <- c(seq(-115,0,by=5),seq(15,130,by=5))
groups <- c(rep("Control",48),rep("cochr",48))
newdat <- cbind(rep(time,2),groups,melt(as.numeric(dat[,1:2])),+
                  melt(as.numeric(dat[,3:4])))

colnames(newdat) <- c("time","group","mean","SEM")


ggplot(newdat,aes(x=time,y=mean,colour=group,fill=group))+
  geom_ribbon(aes(ymin=mean-SEM,ymax=mean+SEM),show.legend=F)+
  #geom_point(size=1, shape=21)+
  geom_line(show.legend=F)+
  scale_color_manual(values=c("red","black"))+
  scale_fill_manual(values=c("red","gray"))+
  xlim(-120,120)+
  ylim(0,2.5)+
  theme(text=element_text(size=16,family="Serif")) + 
  xlab("Time (s)") +
  ylab("Speed (mm/s)") +
  theme_classic()

ggsave("Fig7I_xa334cochr_timeseries.svg", plot = last_plot(),
       scale = 1,width=15,height=10,units="cm",
       dpi = 300, limitsize = TRUE)