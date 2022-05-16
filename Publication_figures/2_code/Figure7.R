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

# Figure 7E - Displacement traces

data <- read.table('Fig5E_RoL1_ablated_displ.csv',header=T,sep=',')

data$Group <- paste(substr(data$Group,13,18),substr(data$Group,26,28))
matx <- data.matrix(data)
rec <- matx[,26:ncol(matx)]
matx <- matx[,2:27]

base <- rowMeans(matx[,1:24])
ti <- data$Blk25

dspeed <- 0.1*(ti - base)

df <- data.frame(data$Group,base,dspeed)
colnames(df) <- c('fish','baseline','dspeed')
df <- as_tibble(df)

dat <- df %>%
  group_by(fish) %>%
  summarise(dspeed = mean(dspeed),baseline = mean(baseline))

dat$group <- substr(dat$fish,1,2)

t <- aov(dat$dspeed ~ dat$group)
summary(t)

t <- aov(dat$baseline ~ dat$group)
summary(t)


ggplot(dat,aes(x=group,y=dspeed))+
  stat_summary(fun.data = quantiles_95, geom="boxplot",fill=c("grey","red"),weight=4)+
  geom_beeswarm(color="grey20",size=3,cex=2)+
  ylab("Speed re. baseline")+
  xlab("")+
  ylim(-5,2)+
  theme_classic()+
  theme(axis.text=element_text(size=18,colour='black',family="Arial",face='plain'))+
  theme(axis.title = element_text(size = 20,colour='black',family="Arial",face='plain'))

ggsave('Fig6E_rol1_dspeed.svg',width=3,height=5.5)


dat$baseline <- 0.1*dat$baseline
ggplot(dat,aes(x=group,y=baseline))+
  stat_summary(fun.data = quantiles_95, geom="boxplot",fill=c("grey","red"),weight=4)+
  geom_beeswarm(color="grey20",size=3,cex=2)+
  ylab("Baseline speed (mm/s)")+
  xlab("")+
  ylim(0,6)+
  theme_classic()+
  theme(axis.text=element_text(size=18,colour='black',family="Arial",face='plain'))+
  theme(axis.title = element_text(size = 20,colour='black',family="Arial",face='plain'))

ggsave('Fig6F_rol1_baseline.svg',width=3,height=5.5)

xa334.shared <- 
  dat %>%
  dabest(group,dspeed,
         idx = c("r0","r1"),
         paired=F)

xa334.shared.effsize <- cohens_d(xa334.shared)

plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Change in speed (mm/s)",
     rawplot.ylim = c(-5,2),
     effsize.ylabel="Effect size")


svg(filename="Fig5E_Rol1_dspeed_estimplot.svg",
    width=4,
    height=5,
    pointsize=14)
plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Speed (re. baseline)",
     rawplot.ylim = c(-5,2),
     effsize.ylabel="Effect size")
dev.off()


t <- aov(dat$baseline ~ dat$group)
summary(t)

dat$baseline <- 0.1*dat$baseline

xa334.shared <- 
  dat %>%
  dabest(group,baseline,
         idx = c("r0","r1"),
         paired=F)

xa334.shared.effsize <- cohens_d(xa334.shared)

plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Baseline speed (mm/s)",
     rawplot.ylim = c(0,6),
     effsize.ylabel="Effect size")


svg(filename="Fig5E_Rol1_baseline_estimplot.svg",
    width=4,
    height=5,
    pointsize=14)
plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Baseline speed (mm/s)",
     rawplot.ylim = c(0,6),
     effsize.ylabel="Effect size")
dev.off()



# Figure 7F - Rol1 kinematics - turns + scoots

data <- read.table('Fig5E_rol1kinematics.csv',header=T,sep='\t')

ggplot(data,aes(x=File,y=AllTurns))+
  stat_summary(fun.data = quantiles_95, geom="boxplot",fill=c("grey70","red"))+
  geom_beeswarm(color="grey20",size=3, linewidth=2,cex=3)+
  #scale_x_discrete(labels = c("a0p5" = "14","a1p0" = "20","a1p5" = "23","a2p0" = "26"))+
  ylab("Turn initiation (%)")+
  xlab("")+
  theme_classic()+
  theme(axis.text=element_text(size=18,colour='black',family="Arial",face='plain'))+
  theme(axis.title = element_text(size = 20,colour='black',family="Arial",face='plain'))

# Turns

xa334.shared <- 
  data %>%
  dabest(File,AllTurns,
         idx = c("con","rol"),
         paired=F)

xa334.shared.effsize <- cohens_d(xa334.shared)

plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Turn initiation (%)",
     rawplot.ylim = c(0,50),
     effsize.ylabel="Effect size")


svg(filename="Fig5E_Rol1_turns_estimplot.svg",
    width=4,
    height=5,
    pointsize=14)
plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Turn initiation (%)",
     rawplot.ylim = c(0,50),
     effsize.ylabel="Effect size")
dev.off()

# Forward swims

xa334.shared <- 
  data %>%
  dabest(File,Scoot,
         idx = c("con","rol"),
         paired=F)

xa334.shared.effsize <- cohens_d(xa334.shared)

plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Forward swim initiation (%)",
     rawplot.ylim = c(0,80),
     effsize.ylabel="Effect size")


svg(filename="Fig5E_Rol1_forwardswims_estimplot.svg",
    width=4,
    height=5,
    pointsize=14)
plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Forward swim initiation (%)",
     rawplot.ylim = c(0,80),
     effsize.ylabel="Effect size")
dev.off()