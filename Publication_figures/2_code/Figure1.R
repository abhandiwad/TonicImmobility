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

quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

##########################################################################################################

#### Initial behavior

## Figure 1A - Tracks

df <- read.table('Fig1A_movementtraces.csv',header=T,sep=',')

pre <- df[1:190,]
post <- df[192:490,]

pre$mmspeed = movavg(pre$mmspeed,5,type='m')

ggplot(pre,aes(x=x,y=y,colour=mmspeed))+
  geom_path(lineend='round',size=3)+
  scale_color_viridis_c(option = "magma",limits = c(0,max(pre$mmspeed)),direction = 1)+
  scale_y_reverse()+
  xlab('')+
  ylab('')+
  theme_classic()

ggsave('Fig1A_v2_pre.svg', plot = last_plot(),
       scale = 1,width=13,height=10,units="cm",
       dpi = 300, limitsize = TRUE)

post$mmspeed = movavg(post$mmspeed,5,type='m')

ggplot(post,aes(x=x,y=y,colour=mmspeed))+
  geom_path(lineend='round',size=3)+
  scale_color_viridis_c(option = "magma",limits = c(0,max(pre$mmspeed)),direction = 1)+
  scale_y_reverse()+
  xlab('')+
  ylab('')+
  theme_classic()

ggsave('Fig1A_v2_post.svg', plot = last_plot(),
       scale = 1,width=13,height=10,units="cm",
       dpi = 300, limitsize = TRUE)

# Individual tracks

data <- read.table('Fig1B_rasterdata_1min.csv',header = TRUE, sep = ",")
dato <- data.matrix(data[,-1])
mm <- movavg(dato,5)

time <- c(seq(-55,0,by=5),seq(15,120,by=5))
gg <- melt(t(data[,-1]))
temp <- length(gg$value)/34
g<- rep(seq(1,temp,1),34)
g<- t(matrix(g,nrow=temp))
g <- melt(g)

df2 <- data.frame(rep(time,temp),g$value,gg$value)
colnames(df2) <- c("time","cellnum","dff")

ggplot(df2,aes(time,cellnum))+
  geom_raster(aes(fill=dff),interpolate=F)+
  #scale_fill_distiller(type="seq",palette = "Magma",direction = 1)+
  scale_fill_viridis_c(option = "magma")+
  scale_x_continuous(breaks=seq(-60,120,by=30))+
  theme_classic()+
  ylab(" ") +
  xlab("Time (s)") +
  theme_classic()+
  theme(text=element_text(size=18,colour='black',family="Arial",face='plain')) + 
  theme(axis.text=element_text(size=12,colour='black',family="Arial",face='plain'))


ggsave('Fig1B_v2_TI_initial_raster_single.svg', plot = last_plot(),
       scale = 1,width=15,height=10,units="cm",
       dpi = 300, limitsize = TRUE)


time <- c(seq(-120,0,by=5),seq(15,215,by=5))
groups <- rep("t1",66)
for (j in 2:106){
  groups <- append(groups,rep(paste("t",as.character(j),sep=""),66))
}

newdat <- cbind(time,groups,melt(as.numeric(data[-1,])))

colnames(newdat) <- c("time","groups","dist")

ggplot(newdat,aes(x=time,y=dist,color=groups))+
  geom_line(show.legend=FALSE)+
  scale_x_continuous(breaks=seq(-120,220,30)) +
  xlab("Time (s)") +
  ylab("Normalized displacement") +
  theme_classic()

ggsave("InitialBehavior.eps", plot = last_plot(), device = "eps",
       scale = 1,width=15,height=7,units="cm",
       dpi = 300, limitsize = TRUE)

# Individual tracks histogram

data <- read.table('TImovementdata.csv')

mtx <- data.matrix(data)

diffmtx <- diff(mtx,lag=1)

thetapre <- abs(atan2d(diffmtx[,2],diffmtx[,1]))
thetapost <- abs(atan2d(diffmtx[,4],diffmtx[,3]))

# Fig 1E - Kinematic analysis - Scoots and turns

data <- read.table('Fig1E_swimkinematics.csv',sep = ',', header=T)
data <- as_tibble(data)
data$Displacement <- 0.1*data$Displacement


# Responsiveness
# Box plot

byscoot <- data[which(data$Bout == "Fwd"),]


dat <- byscoot %>%
  group_by(Fish) %>%
  summarise(dspeed = mean(dspeed),baseline = mean(baseline))

ggplot(byscoot,aes(x=Cond,y=Response))+
  stat_summary(fun.data = quantiles_95, geom="boxplot",fill=c("#e3dbe0","#ff894a"),weight=4)+
  geom_beeswarm(color="grey20",size=3)+
  #scale_x_discrete(labels = c("a0p5" = "14","a1p0" = "20","a1p5" = "23","a2p0" = "26"))+
  xlab("")+
  theme_classic()+
  theme(axis.text=element_text(size=18,colour='black',family="Arial",face='plain'))+
  theme(axis.title = element_text(size = 20,colour='black',family="Arial",face='plain'))

#6e4aff
# Estimation plot
multi.paired <- 
  data %>%
  dabest(Cond,Response,
         idx = list(c("Base","Post"),c("Basel","Postf")),
         paired=T, id.col = Fish)

multi.paired.meandiff <- cohens_d(multi.paired)



CairoSVG(file="Fig1E_responsiveness.svg",
         width=4,
         height=8,
         pointsize=12)
plot(multi.paired.meandiff, effsize.ylabel = "Effect Size", 
     color.column = Cond, palette = c("#D936A0","#D936A0","#4C48D9","#4C48D9"), 
     slopegraph = TRUE)
dev.off()


# Initial angle

multi.paired <- 
  data %>%
  dabest(Cond,Angle,
         idx = list(c("Base","Post"),c("Basel","Postf")),
         paired=T, id.col = Fish)

multi.paired.meandiff <- cohens_d(multi.paired)

plot(multi.paired.meandiff, effsize.ylabel = "Effect Size", 
     color.column = Cond, palette = c("#D936A0","#D936A0","#4C48D9","#4C48D9"), 
     slopegraph = TRUE)

CairoSVG(file="Fig1E_angle.svg",
         width=4,
         height=8,
         pointsize=12)
plot(multi.paired.meandiff, effsize.ylabel = "Effect Size", 
     color.column = Cond, palette = c("#D936A0","#D936A0","#4C48D9","#4C48D9"), 
     slopegraph = TRUE)
dev.off()


# Displacement

multi.paired <- 
  data %>%
  dabest(Cond,Displacement,
         idx = list(c("Base","Post"),c("Basel","Postf")),
         paired=T, id.col = Fish)

multi.paired.meandiff <- cohens_d(multi.paired)

plot(multi.paired.meandiff, effsize.ylabel = "Effect Size", 
     color.column = Cond, palette = c("#D936A0","#D936A0","#4C48D9","#4C48D9"), 
     slopegraph = TRUE)

CairoSVG(file="Fig1E_displacement.svg",
         width=4,
         height=8,
         pointsize=12)
plot(multi.paired.meandiff, effsize.ylabel = "Effect Size", 
     color.column = Cond, palette = c("#D936A0","#D936A0","#4C48D9","#4C48D9"), 
     slopegraph = TRUE)
dev.off()



#####################################################################################################################
################## Multiple stimulus amplitudes (Fig 1) ############################################################

data <- read.table('Fig1C_increasingstimuli.csv',header = TRUE, sep = ",")

data <- t(data)
time <- c(seq(-115,0,by=5),seq(15,220,by=5))
groups <- c(rep("0.5V",66),rep("1V",66),rep("1.5V",66),rep("2V",66))
newdat <- cbind(rep(time,4),groups,melt(as.numeric(data[-1,1:4])),+
                  melt(as.numeric(data[-1,6:9])))

colnames(newdat) <- c("time","geno","mean","SEM")

# Plot time series

ggplot(newdat,aes(x=time,y=mean,colour=geno,fill=geno))+
  geom_ribbon(aes(ymin=mean-SEM,ymax=mean+SEM),show.legend=F)+
  #geom_point(size=1, shape=21)+
  geom_line(show.legend=F)+
  scale_color_manual(values=c("#fee5d9","#fcae91","#fb6a4a","#cb181d"))+
  scale_fill_manual(values=c("#fee5d9","#fcae91","#fb6a4a","#cb181d"))+
  scale_x_continuous(breaks=seq(-120,220,30))+
  coord_cartesian(xlim=c(-60,120),ylim = c(0,1.2))+
  #ylim(0,1.2)+
  xlab("Time (s)") +
  ylab("Speed (re. baseline)") +
  theme_classic()+
  theme(text=element_text(size=16,family="Arial",color="Black"))+
  theme(axis.text=element_text(size=14,family="Arial",color="Black"))

ggsave("Fig1C_multiplestims.svg", plot = last_plot(),
       scale = 1,width=27,height=18,units="cm",
       dpi = 300, limitsize = TRUE)


# Figure 1C - Multiple stims change in speed

data <- read.table('Fig1C_MultipleStimAmplitudes.csv',header=T,sep=',')

data$displ <- paste0(substr(data$displ,12,18),substr(data$displ,27,29))
matx <- data.matrix(data)
rec <- matx[,26:ncol(matx)]
matx <- matx[,2:27]

base <- 0.1*rowMeans(matx[,1:24])
ti <- 0.1*data$Blk24

dspeed <- ti - base

df <- data.frame(data$displ,base,dspeed)
colnames(df) <- c('fish','baseline','dspeed')
df <- as_tibble(df)

df <- df[which(df$baseline >= 0.5),]

dat <- df %>%
  group_by(fish) %>%
  summarise(dspeed = mean(dspeed),baseline = mean(baseline))

dat$group <- substr(dat$fish,1,4)


ggplot(dat,aes(x=group,y=dspeed))+
  stat_summary(fun.data = quantiles_95, geom="boxplot",fill=c("gray80","#fee5d9","#fcae91","#fb6a4a","#cb181d"),weight=2)+
  geom_beeswarm(color="grey60",size=3)+
  scale_x_discrete(labels = c("a000" = "0", "a0p5" = "14","a1p0" = "20","a1p5" = "23","a2p0" = "26"))+
  ylab("")+
  xlab("")+
  theme_classic()+
  theme(axis.text=element_text(size=18,colour='black',family="Arial",face='plain'))

ggsave("Fig1D_MultipleStims_changeinspeed.svg", plot = last_plot(),
       scale = 1,width=10,height=20,units="cm",
       dpi = 300, limitsize = TRUE)

g = aov(dat$dspeed~dat$group)
summary(g)
dat$group <- as.factor(dat$group)
DunnettTest(dspeed ~ group, data=dat ,control="a000")

# Multiple stims time to recovery (tau)

setwd('D:/TransferInhibition/RealTimeTracking/InitialBehavior/09162016_9well_multiplestimamplitudes/anadata8e_00-0')
x <- {}
y <- {}

data <- read.table('a2p0.csv',header=T)

timeseq <- seq(5,210,5)
TI <- data[,26:ncol(data)]
#TI <- matrix(melt(TI)$value,nrow=nrow(TI))

recov <- {}

for (i in 1:nrow(TI)){
  
  halfway <- TI[i,1] + (1-TI[i,1])/2
  idx <- which(TI[i,]>halfway)
  if (isEmpty(idx)){
    idx = length(timeseq)
  }
  
  if (is.nan(idx) == F){
    recov <- append(recov,timeseq[idx[1]])
  }
  
}

x <- append(x,rep(4,length(recov)))
y <- append(y,recov)

df <- data.frame(factor(x),y)
colnames(df) <- cbind('Stim','dt')

res.aov <- aov(dt ~ Stim, data=df)
summary(res.aov)

# Multiple stims recovery figure

df <- read.table('multiplestims_timetorecovery.csv',header=T,sep=',')

s <- group_by(df,Stim)

avs <- summarize(s,avs = mean(dt),ste=std_err(dt))
avs$Stim <- c(0.5,1.0,1.5,2.0)

ggplot(avs,aes(Stim,avs))+
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymin=avs-ste,ymax=avs+ste),width=0.1)+
  ylim(0,27)+
  xlab("Voltage output") +
  ylab("T") +
  theme(text=element_text(size=16,family="Serif")) + 
  theme_classic()

setwd('D:/Manuscripts/SEBA paper/OngoingFigureData')

ggsave("MultipleStims_timetorecovery.eps", plot = last_plot(), device = "eps",
       scale = 1,width=8,height=14,units="cm",
       dpi = 300, limitsize = TRUE)

####################### Multiple consecutive stimuli #################################################################
######################################################################################################################

data <- read.table('initial_consecutive_stims.csv',header = TRUE, sep = "\t")

data <- t(data)
time <- c(seq(-120,-5,by=5),seq(15,220,by=5),seq(225,340,by=5),seq(360,565,by=5),seq(570,685,by=5),seq(700,905,by=5))
groups <- c(rep("cnt",198),rep("abl",198))
newdat <- cbind(rep(time,2),groups,melt(as.numeric(data[-1,1:2])),+
                  melt(as.numeric(data[-1,4:5])))

colnames(newdat) <- c("time","geno","mean","SEM")

# Plot time series

ggplot(newdat,aes(x=time,y=mean,colour=geno))+
  geom_ribbon(aes(ymin=mean-SEM,ymax=mean+SEM),fill='grey60',show.legend=F)+
  #geom_point(size=1, shape=21)+
  geom_line(show.legend=F)+
  scale_color_manual(values=c("red","grey50"))+
  #scale_fill_manual(values=c("red","grey50"))+
  scale_x_continuous(breaks=seq(-120,920,60)) +
  ylim(0,1.2)+
  xlab("Time (s)") +
  ylab("Normalized displacement") +
  theme(text=element_text(size=16,family='Arial'))+
  theme_classic()

ggsave("ConsecutiveStims.eps", plot = last_plot(), device = "eps",
       scale = 1,width=30,height=10,units="cm",
       dpi = 300, limitsize = TRUE)


## Figure 1F - Balance data
# Replicate 2 - 10/19 
# Changes: [1] Done in triplicate and avgd, [2] Recorded for 15 sec

da <- read.table('Fig1F_tibalance_v3.csv',header=T,sep=',')
code <- read.table('Fig1F_TIbalance_coding.csv',header=F,sep=",")
da <- as_tibble(da)
da$Time <- (8*(da$Frame-1800)/1000)

temp <- unique(da$Plate)
group <- {}
for (i in 1:length(temp)){
  group[which(da$Plate == temp[i])] <- code$V1[i]
}

da$group <- group
da$Cond <- substr(da$group,1,1)

dat <- da %>% 
  group_by(Time,Cond,Plate) %>%
  summarise(punb = mean(punb))


dsumms <- dat %>% 
  group_by(Time,Cond) %>%
  summarise(mean_punbal = mean(punb),se_punbal = std_err(punb))

ggplot(dat,aes(x=Time,color=Cond))+
  geom_line(aes(y=punb,group=Plate),alpha=0.3,show.legend=F)+
  geom_line(data=dsumms,aes(y=mean_punbal),size=1.5,show.legend=F)+
  geom_ribbon(data=dsumms,aes(ymin=mean_punbal-se_punbal,
                              ymax=mean_punbal+se_punbal,
                              fill=Cond),size=0.2,alpha=0.5,show.legend=F)+
  scale_color_manual(values=c("grey30","red"))+
  scale_fill_manual(values=c("grey30","red"))+
  scale_x_continuous(breaks=seq(0,14,2)) +
  ylab("Proportion of unbalanced fish")+
  xlab("Time (s)")+
  theme_classic()+
  theme(text=element_text(size=20,family="Arial",color='black',))+
  theme(axis.text=element_text(size=18,family="Arial",color='black'))

ggsave("Fig1F_balance_reflex.svg", plot = last_plot(), 
       scale = 1,width=25,height=20,units="cm",
       dpi = 300, limitsize = TRUE)

temp <- unique(dat$Plate)
time <- unique(da$Time)
plate <- {}
Cond <- {}
tau <- {}

for (i in 1:length(temp)){
  plate[i] <- temp[i]
  
  platenum <- dat[which(dat$Plate == temp[i]),]
  
  Cond[i] <- platenum$Cond[1]
  
  if (platenum$punb[1] > 0){
    
    timeidx <- which(platenum$punb <= 0.05*platenum$punb[1])
    if (isempty(timeidx)){
      tau[i] <- time[length(time)]
    } else {
      tau[i] <- time[timeidx[1]]
    }
    
  } else {
    
    tau[i] <- 0
    
  }
}

punbtau <- data.frame(Cond,plate,tau)

punb.unpaired <- 
  punbtau %>%
  dabest(Cond,tau,
         idx = c("c","t"),
         paired=F, id.col = Cond)

punb.unpaired.meandiff <- cohens_d(punb.unpaired)

CairoSVG(file="Fig1F_TimetoBalanced.svg",
         width=4,
         height=8,
         pointsize=12)
plot(punb.unpaired.meandiff, rawplot.ylabel = "Time to balanced (s)", effsize.ylabel = "Effect Size", 
     color.column = Cond, palette = c("grey","red"), show.legend = F,
     slopegraph = F)
dev.off()

summ <- aov(tau ~ Cond)
summary(summ)

#Fig 1G - SLC during arrest

da <- read.table('Fig1G_TI_SLC.csv',header=T,sep=',')
da <- as_tibble(da)

dat <- da %>%
  group_by(Fishnum,Cond) %>%
  summarise(mSLC = mean(SLC))


dsumms <- dat %>% 
  group_by(Cond) %>%
  summarise(mean_SLC = mean(mSLC,na.rm=T),se_SLC = std_err(mSLC))

ggplot(dat,aes(x=Cond))+
  #geom_line(aes(y=mSLC,group=Fishnum),alpha=0.3,show.legend=F)+
  geom_line(data=dsumms,aes(y=mean_SLC),size=1,show.legend=F)+
  geom_point(data=dsumms,aes(y=mean_SLC),size=3,show.legend=F)+
  geom_errorbar(data=dsumms,aes(ymin=mean_SLC-se_SLC,
                                ymax=mean_SLC+se_SLC),
                size=1,width=0.2,show.legend=F)+
  scale_x_continuous(breaks = c(1,2,3),labels = c("Base","TI","Recovery")) +
  ylab("Startle response (%)")+
  ylim(0,80)+
  xlab("")+
  theme_classic()+
  theme(text=element_text(size=16,family="sans",color='black'))


ggsave("Fig1G_SLC.svg", plot = last_plot(), device = "svg",
       scale = 1,width=9,height=15,units="cm",
       dpi = 300, limitsize = TRUE)

# Figure 1G- O-bend during TI

a <- read.table('Fig1G_TI_OBend.csv',header=T,sep=',')
da <- as_tibble(a)

dat <- da %>%
  group_by(Fishnum,Cond) %>%
  summarise(mObend = mean(Obend))

dsumms <- dat %>% 
  group_by(Cond) %>%
  summarise(mean_obend = mean(mObend,na.rm=T),se_obend = std_err(mObend))

ggplot(dat,aes(x=Cond))+
  #geom_line(aes(y=mSLC,group=Fishnum),alpha=0.3,show.legend=F)+
  geom_line(data=dsumms,aes(y=mean_obend),size=1,show.legend=F)+
  geom_point(data=dsumms,aes(y=mean_obend),size=3,show.legend=F)+
  geom_errorbar(data=dsumms,aes(ymin=mean_obend-se_obend,
                                ymax=mean_obend+se_obend),
                size=1,width=0.2,show.legend=F)+
  scale_x_continuous(breaks = c(1,2,3),labels = c("Base","TI","Recovery")) +
  ylab("O-bend response (%)")+
  ylim(0,80)+
  xlab("")+
  theme_classic()+
  theme(text=element_text(size=16,family="sans",color='black'))


ggsave("Fig1G_obend.svg", plot = last_plot(), device = "svg",
       scale = 1,width=9,height=15,units="cm",
       dpi = 300, limitsize = TRUE)


SLC <- data.frame(Genotype=c("Baseline","Post-stimulus","Recovery"),NormDisp = c(64.5,5.76,46.93),Err=c(12.17,5.49,13.17))
LLC <- data.frame(Genotype=c("Baseline","Post-stimulus","Recovery"),NormDisp = c(65.44,54.54,49.45),Err=c(18.17,10.69,14.08))
obend <- data.frame(Genotype=c("Baseline","Post-stimulus","Recovery"),NormDisp = c(45.28,14.28,35.45),Err=c(9.6,5.9,9.29))

escapes <- cbind(rep('SLC',3),rep('LLC',3),rep('Obd',3))
escapes <- melt(escapes)

df <- rbind(SLC,LLC,obend)
df <- cbind(escapes$value,df)

df$Genotype <- factor(df$Genotype, as.character(df$Genotype))

colnames(df) <- c('Esc','Cond','Val','Sterr')

ggplot(df, aes(x=Esc,y=Val,fill=Cond)) + 
  geom_bar(position="dodge",stat='identity',show.legend=FALSE) +
  geom_errorbar(aes(ymin=Val-Sterr, ymax=Val+Sterr),position="dodge", width = 0.2) +
  xlab("Condition") +
  ylab("Percent response") +
  scale_fill_manual(values=c("gray20","red","cyan"))+
  ylim(0,80) + 
  theme(text=element_text(size=16,family="Serif")) + 
  theme_classic()

ggsave("escapes.eps", plot = last_plot(), device = "eps",
       scale = 1,width=15,height=10,units="cm",
       dpi = 300, limitsize = TRUE)

## SLC stats

df <- read.table('D:/TransferInhibition/HighSpeedTracking/TransferInhibitionStartleInduction/SLCduringSEBA.csv',sep=',',header=T)

rmanova <- aov(SLC~factor(Cond)+Error(factor(Fishnum)), data = df)

summary(rmanova)

## Obend stats

df <- read.table('D:/TransferInhibition/HighSpeedTracking/TransferInhibitionStartleInduction/ObendduringSEBA.csv',sep=',',header=T)

rmanova <- aov(Obend~factor(Cond)+Error(factor(Fishnum)), data = df)

summary(rmanova)

# SLC: Acoustic vs. eshock

A <- data.frame(Genotype=c("Baseline","Post-stimulus"),NormDisp = c(63.5,7.26),Err=c(10.17,6.49))
E <- data.frame(Genotype=c("Baseline","Post-stimulus"),NormDisp = c(40.74,47.91),Err=c(3.28,10.19))

escapes <- cbind(rep('Acoustic',2),rep('Electric',2))
escapes <- melt(escapes)

df <- rbind(A,E)
df <- cbind(escapes$value,df)

df$Genotype <- factor(df$Genotype, as.character(df$Genotype))

df$Genotype <- factor(df$Genotype, as.character(df$Genotype))

colnames(df) <- c('Esc','Cond','Val','Sterr')

ggplot(df, aes(x=Esc,y=Val,fill=Cond)) + 
  geom_bar(position="dodge",stat='identity',show.legend=FALSE) +
  geom_errorbar(aes(ymin=Val-Sterr, ymax=Val+Sterr),position="dodge", width = 0.2) +
  xlab("Condition") +
  ylab("Percent response") +
  scale_fill_manual(values=c("gray20","red"))+
  ylim(0,80) + 
  theme(text=element_text(size=16,family="Serif")) + 
  theme_classic()

ggsave("SLC_acoustic_electric.eps", plot = last_plot(), device = "eps",
       scale = 1,width=10,height=15,units="cm",
       dpi = 300, limitsize = TRUE)

# Kinematic analysis for SLCs [09072019 data]

data <- read.table('SEBA_SLC_angle.csv',header=T)

group <- c(rep(1,times=nrow(data)),rep(2,times=nrow(data)))
group <- as.character(group)
resp <- c(data$Pre,data$Blk06)

df <- data.frame(group,resp)
colnames(df) <- c('group','resp')

ggplot(df, aes(x=group,y=resp)) + 
  geom_boxplot(fill=c('red','blue')) +
  geom_beeswarm(size=2,fill='gray80')+
  xlab("") +
  ylab("Bend angle") +
  ylim(0,200) + 
  theme(text=element_text(size=16,family="Serif")) + 
  theme_classic()

ggsave("slc_angle.eps", plot = last_plot(), device = "eps",
       scale = 1,width=7,height=15,units="cm",
       dpi = 300, limitsize = TRUE)

data <- read.table('SEBA_SLC_c2angle.csv',header=T)

group <- c(rep(1,times=nrow(data)),rep(2,times=nrow(data)))
group <- as.character(group)
resp <- c(data$Pre,data$Blk06)

df <- data.frame(group,resp)
colnames(df) <- c('group','resp')

ggplot(df, aes(x=group,y=resp)) + 
  geom_boxplot(fill=c('red','blue')) +
  geom_beeswarm(size=2,fill='gray80')+
  xlab("") +
  ylab("Counterbend angle") +
  ylim(0,150) + 
  theme(text=element_text(size=16,family="Serif")) + 
  theme_classic()

ggsave("slc_c2angle.eps", plot = last_plot(), device = "eps",
       scale = 1,width=7,height=15,units="cm",
       dpi = 300, limitsize = TRUE)

data <- read.table('SEBA_SLC_disp.csv',header=T)

group <- c(rep(1,times=nrow(data)),rep(2,times=nrow(data)))
group <- as.character(group)
resp <- c(data$Pre,data$Blk06)

df <- data.frame(group,resp)
colnames(df) <- c('group','resp')

ggplot(df, aes(x=group,y=resp)) + 
  geom_boxplot(fill=c('red','blue')) +
  geom_beeswarm(size=2,fill='gray80')+
  xlab("") +
  ylab("Displacement") +
  ylim(0,125) + 
  theme(text=element_text(size=16,family="Serif")) + 
  theme_classic()

ggsave("slc_displ.eps", plot = last_plot(), device = "eps",
       scale = 1,width=7,height=15,units="cm",
       dpi = 300, limitsize = TRUE)

data <- read.table('SEBA_SLC_latency.csv',header=T)

group <- c(rep(1,times=nrow(data)),rep(2,times=nrow(data)))
group <- as.character(group)
resp <- c(data$Pre,data$Blk06)

df <- data.frame(group,resp)
colnames(df) <- c('group','resp')

ggplot(df, aes(x=group,y=resp)) + 
  geom_boxplot(fill=c('red','blue')) +
  geom_beeswarm(size=2,fill='gray80')+
  xlab("") +
  ylab("Latency") +
  ylim(0,15) + 
  theme(text=element_text(size=16,family="Serif")) + 
  theme_classic()

ggsave("slc_latency.eps", plot = last_plot(), device = "eps",
       scale = 1,width=7,height=15,units="cm",
       dpi = 300, limitsize = TRUE)


# Kinematics - O-bend [0516 data]

data <- read.table('SEBA_obend_angle.csv',header=T)

group <- c(rep(1,times=nrow(data)),rep(2,times=nrow(data)))
group <- as.character(group)
resp <- c(data$pre,data$post)

df <- data.frame(group,resp)
colnames(df) <- c('group','resp')

ggplot(df, aes(x=group,y=resp)) + 
  geom_boxplot(fill=c('red','blue')) +
  geom_beeswarm(size=2,fill='gray80')+
  xlab("") +
  ylab("Bend angle") +
  ylim(0,200) + 
  theme(text=element_text(size=16,family="Serif")) + 
  theme_classic()

ggsave("obend_angle.eps", plot = last_plot(), device = "eps",
       scale = 1,width=7,height=15,units="cm",
       dpi = 300, limitsize = TRUE)

data <- read.table('SEBA_obend_c2angle.csv',header=T)

group <- c(rep(1,times=nrow(data)),rep(2,times=nrow(data)))
group <- as.character(group)
resp <- c(data$pre,data$post)

df <- data.frame(group,resp)
colnames(df) <- c('group','resp')

ggplot(df, aes(x=group,y=resp)) + 
  geom_boxplot(fill=c('red','blue')) +
  geom_beeswarm(size=2,fill='gray80')+
  xlab("") +
  ylab("Counterbend angle") +
  ylim(0,50) + 
  theme(text=element_text(size=16,family="Serif")) + 
  theme_classic()

ggsave("obend_c2angle.eps", plot = last_plot(), device = "eps",
       scale = 1,width=7,height=15,units="cm",
       dpi = 300, limitsize = TRUE)

data <- read.table('SEBA_obend_disp.csv',header=T)

group <- c(rep(1,times=nrow(data)),rep(2,times=nrow(data)))
group <- as.character(group)
resp <- c(data$pre,data$post)

df <- data.frame(group,resp)
colnames(df) <- c('group','resp')

ggplot(df, aes(x=group,y=resp)) + 
  geom_boxplot(fill=c('red','blue')) +
  geom_beeswarm(size=2,fill='gray80')+
  xlab("") +
  ylab("Displacement") +
  ylim(0,50) + 
  theme(text=element_text(size=16,family="Serif")) + 
  theme_classic()

ggsave("obend_displ.eps", plot = last_plot(), device = "eps",
       scale = 1,width=7,height=15,units="cm",
       dpi = 300, limitsize = TRUE)

data <- read.table('SEBA_obend_latency.csv',header=T)

group <- c(rep(1,times=nrow(data)),rep(2,times=nrow(data)))
group <- as.character(group)
resp <- c(data$pre,data$post)

df <- data.frame(group,resp)
colnames(df) <- c('group','resp')

ggplot(df, aes(x=group,y=resp)) + 
  geom_boxplot(fill=c('red','blue')) +
  geom_beeswarm(size=2,fill='gray80')+
  xlab("") +
  ylab("Latency") +
  ylim(0,300) + 
  theme(text=element_text(size=16,family="Serif")) + 
  theme_classic()

ggsave("obend_latency.eps", plot = last_plot(), device = "eps",
       scale = 1,width=7,height=15,units="cm",
       dpi = 300, limitsize = TRUE)