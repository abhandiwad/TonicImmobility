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


df <- read.table('Fig3c_xa239xc43.csv',header=T,sep=',')

df$dSpeed <- 0.1*df$dSpeed
df$base <- 0.1*df$Base

xa334.shared <- 
  df %>%
  dabest(Group,dSpeed,
         idx = c("control","r"),
         paired=F)

xa334.shared.effsize <- cohens_d(xa334.shared)

plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Change in speed (mm/s)",
     rawplot.ylim= c(-3,1),
     effsize.ylabel="Effect size")

svg(filename="Fig3C_xa239xc43_estimplot.svg",
    width=4,
    height=5,
    pointsize=14)
plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Change in speed (mm/s)",
     rawplot.ylim= c(-3,1),
     effsize.ylabel="Effect size")
dev.off()

# ggplot(df, aes(x=Group,y=dSpeed)) + 
#   geom_boxplot(fill = c('gray80','magenta'),color='black',size=1)+
#   geom_beeswarm(size=4,colour='gray20',alpha=0.96,cex=4)+
#   ylab("Change in speed (mm/s)") +
#   xlab(" ") +
#   scale_x_discrete(labels=c("Control","Ablated"))+
#   ylim(-2.5,0.5) + 
#   theme_classic()+
#   theme(text=element_text(size=24,colour='black',family="Arial",face='bold')) + 
#   theme(axis.text=element_text(size=18,colour='black',family="Arial",face='plain'))
# 
# ggsave("Fig3c_xa239xc43_t0_boxplot.svg", plot = last_plot(), 
#        scale = 1,width=10,height=20,units="cm",
#        dpi = 300, limitsize = TRUE)

t.test(df$dSpeed ~ df$Group)

### xa239;xc121 [Negative intersect]

df <- read.table('xa239xc121_normd.csv',header=T)

ggplot(df, aes(x=Group,y=1-Blk24)) + 
  geom_boxplot(fill = c('gray80','magenta'))+
  geom_beeswarm(dotsize=2,fill='gray80')+
  xlab("Genotype") +
  ylab("Normalized displacement") +
  ylim(-1,1) + 
  theme(text=element_text(size=16,family="Serif")) + 
  theme_classic()

ggsave("xa239xc121_t0_boxplot.svg", plot = last_plot(),
       scale = 1,width=7,height=15,units="cm",
       dpi = 300, limitsize = TRUE)



###############################################################################################
####################### Figure 3G - Primary sensory  ##########################################

data <- read.table('Fig3G_dark.csv',header=T,sep=',')
data$displ <- paste(substr(data$displ,13,18),substr(data$displ,26,28))
matx <- data.matrix(data)
rec <- matx[,26:ncol(matx)]
matx <- matx[,2:27]

base <- rowMeans(matx[,1:24])
ti <- data$Blk24

dspeed <- 0.1*(ti - base)

df <- data.frame(data$displ,base,dspeed)
colnames(df) <- c('fish','baseline','dspeed')
df <- as_tibble(df)

dat <- df %>%
  group_by(fish) %>%
  summarise(dspeed = mean(dspeed),baseline = mean(baseline))

dat$group <- substr(dat$fish,1,2)

t <- aov(dat$dspeed ~ dat$group)
summary(t)

ggplot(dat,aes(x=group,y=dspeed))+
  geom_boxplot(fill = c('gray80','red'),color='black',size=1)+
  geom_beeswarm(size=4,colour='gray20',alpha=0.96,cex=4)+
  ylab("Change in speed (mm/s)") +
  xlab(" ") +
  scale_x_discrete(labels=c("Control","Ablated"))+
  ylim(-5,1.5) + 
  theme_classic()+
  theme(text=element_text(size=24,colour='black',family="Arial",face='bold')) + 
  theme(axis.text=element_text(size=18,colour='black',family="Arial",face='plain'))

ggsave("Fig3G_dark.svg", plot = last_plot(),
       scale = 1,width=8,height=15,units="cm",
       dpi = 300, limitsize = TRUE)



data <- read.table('Fig3G_SAG.csv',header=T,sep=',')
data <- data[seq(1,nrow(data),3),]
data$displ <- paste(substr(data$displ,13,18),substr(data$displ,26,28))
matx <- data.matrix(data)
rec <- matx[,26:ncol(matx)]
matx <- matx[,2:27]

base <- rowMeans(matx[,1:24])
ti <- data$Blk24

dspeed <- 0.1*(ti - base)

df <- data.frame(data$displ,base,dspeed)
colnames(df) <- c('fish','baseline','dspeed')
df <- as_tibble(df)

dat <- df %>%
  group_by(fish) %>%
  summarise(dspeed = mean(dspeed),baseline = mean(baseline))

dat$group <- substr(dat$fish,1,2)

t <- aov(dat$dspeed ~ dat$group)
summary(t)

ggplot(dat,aes(x=group,y=dspeed))+
  geom_boxplot(fill = c('gray80','red'),color='black',size=1)+
  geom_beeswarm(size=4,colour='gray20',alpha=0.96,cex=4)+
  ylab("Change in speed (mm/s)") +
  xlab(" ") +
  scale_x_discrete(labels=c("Control","Ablated"))+
  ylim(-5,1.5) + 
  theme_classic()+
  theme(text=element_text(size=24,colour='black',family="Arial",face='bold')) + 
  theme(axis.text=element_text(size=18,colour='black',family="Arial",face='plain'))

ggsave("Fig3G_SAG.svg", plot = last_plot(),
       scale = 1,width=8,height=15,units="cm",
       dpi = 300, limitsize = TRUE)

data <- read.table('Fig3G_Neo2.trk',header=T,sep=',')
data$displ <- paste(substr(data$displ,13,18),substr(data$displ,26,28))
matx <- data.matrix(data)
rec <- matx[,26:ncol(matx)]
matx <- matx[,2:27]

base <- rowMeans(matx[,1:24])
ti <- data$Blk24

dspeed <- 0.1*(ti - base)

df <- data.frame(data$displ,base,dspeed)
colnames(df) <- c('fish','baseline','dspeed')
df <- as_tibble(df)

dat <- df %>%
  group_by(fish) %>%
  summarise(dspeed = mean(dspeed),baseline = mean(baseline))

dat$group <- substr(dat$fish,1,2)

t <- aov(dat$dspeed ~ dat$group)
summary(t)

ggplot(dat,aes(x=group,y=dspeed))+
  geom_boxplot(fill = c('gray80','red'),color='black',size=1)+
  geom_beeswarm(size=4,colour='gray20',alpha=0.96,cex=4)+
  ylab("Change in speed (mm/s)") +
  xlab(" ") +
  scale_x_discrete(labels=c("Control","Ablated"))+
  ylim(-5,1.5) + 
  theme_classic()+
  theme(text=element_text(size=24,colour='black',family="Arial",face='bold')) + 
  theme(axis.text=element_text(size=18,colour='black',family="Arial",face='plain'))

ggsave("Fig3G_Neo.svg", plot = last_plot(),
       scale = 1,width=8,height=15,units="cm",
       dpi = 300, limitsize = TRUE)

data <- read.table('Fig3G_SAG_Neo2.csv',header=T,sep=',')
data$displ <- paste(substr(data$displ,13,18),substr(data$displ,26,28))
matx <- data.matrix(data)
rec <- matx[,26:ncol(matx)]
matx <- matx[,2:27]

base <- rowMeans(matx[,1:24])
ti <- data$Blk24

dspeed <- 0.1*(ti - base)

df <- data.frame(data$displ,base,dspeed)
colnames(df) <- c('fish','baseline','dspeed')
df <- as_tibble(df)

dat <- df %>%
  group_by(fish) %>%
  summarise(dspeed = mean(dspeed),baseline = mean(baseline))

dat$group <- substr(dat$fish,1,2)

t <- aov(dat$dspeed ~ dat$group)
summary(t)

ggplot(dat,aes(x=group,y=dspeed))+
  geom_boxplot(fill = c('gray80','red'),color='black',size=1)+
  geom_beeswarm(size=4,colour='gray20',alpha=0.96,cex=4)+
  ylab("Change in speed (mm/s)") +
  xlab(" ") +
  scale_x_discrete(labels=c("Control","Ablated"))+
  ylim(-5,1.5) + 
  theme_classic()+
  theme(text=element_text(size=24,colour='black',family="Arial",face='bold')) + 
  theme(axis.text=element_text(size=18,colour='black',family="Arial",face='plain'))

ggsave("Fig3G_SAG_neo.svg", plot = last_plot(),
       scale = 1,width=8,height=15,units="cm",
       dpi = 300, limitsize = TRUE)

data <- read.table('Fig3G_Trig.csv',header=T,sep=',')
data$displ <- paste(substr(data$displ,13,18),substr(data$displ,26,28))
matx <- data.matrix(data)
rec <- matx[,26:ncol(matx)]
matx <- matx[,2:27]

base <- rowMeans(matx[,1:24])
ti <- data$Blk24

dspeed <- 0.1*(ti - base)

df <- data.frame(data$displ,base,dspeed)
colnames(df) <- c('fish','baseline','dspeed')
df <- as_tibble(df)

dat <- df %>%
  group_by(fish) %>%
  summarise(dspeed = mean(dspeed),baseline = mean(baseline))

dat$group <- substr(dat$fish,1,2)

t <- aov(dat$dspeed ~ dat$group)
summary(t)

ggplot(dat,aes(x=group,y=dspeed))+
  geom_boxplot(fill = c('gray80','red'),color='black',size=1)+
  geom_beeswarm(size=4,colour='gray20',alpha=0.96,cex=4)+
  ylab("Change in speed (mm/s)") +
  xlab(" ") +
  scale_x_discrete(labels=c("Control","Ablated"))+
  ylim(-5,1.5) + 
  theme_classic()+
  theme(text=element_text(size=24,colour='black',family="Arial",face='bold')) + 
  theme(axis.text=element_text(size=18,colour='black',family="Arial",face='plain'))

ggsave("Fig3G_Trig.svg", plot = last_plot(),
       scale = 1,width=8,height=15,units="cm",
       dpi = 300, limitsize = TRUE)