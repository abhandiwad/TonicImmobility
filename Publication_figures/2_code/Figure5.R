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


df <- read.table('Fig5c_ethanoltreated.csv',header=T,sep=',')
df$dSpeed <- 0.1*df$dSpeed

xa334.shared <- 
  df %>%
  dabest(Group,dSpeed,
         idx = c("control","treated"),
         paired=F)

xa334.shared.effsize <- cohens_d(xa334.shared)

plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Change in speed (mm/s)",
     rawplot.ylim = c(-5,3),
     effsize.ylabel="Effect size")

svg(filename="Fig5C_ethanol_estimplot.svg",
    width=4,
    height=5,
    pointsize=14)
plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Change in speed (mm/s)",
     rawplot.ylim = c(-5,3),
     effsize.ylabel="Effect size")
dev.off()
# ggplot(df, aes(x=Group,y=dSpeed)) + 
#   geom_boxplot(fill = c('gray80','magenta'),color='black',size=1)+
#   geom_beeswarm(size=4,colour='gray20',alpha=0.96,cex=4)+
#   ylab("Change in speed (mm/s)") +
#   xlab(" ") +
#   scale_x_discrete(labels=c("Control","EtOH treated"))+
#   ylim(-5,5) + 
#   theme_classic()+
#   theme(text=element_text(size=24,colour='black',family="Arial",face='bold')) + 
#   theme(axis.text=element_text(size=18,colour='black',family="Arial",face='plain'))
# 
# ggsave("Fig5c_ethanoltreat_beeswarm.svg", plot = last_plot(), 
#        scale = 1,width=10,height=20,units="cm",
#        dpi = 300, limitsize = TRUE)

########### Washout

df <- read.table('Fig5c_ethanolrecov.csv',header=T,sep=',')
df$dSpeed <- 0.1*df$dSpeed

xa334.shared <- 
  df %>%
  dabest(Group,dSpeed,
         idx = c("control","treated"),
         paired=F)

xa334.shared.effsize <- cohens_d(xa334.shared)

plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Change in speed (mm/s)",
     rawplot.ylim = c(-5,3),
     effsize.ylabel="Effect size")

svg(filename="Fig5C_ethanolrecov_estimplot.svg",
    width=4,
    height=5,
    pointsize=14)
plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Change in speed (mm/s)",
     rawplot.ylim = c(-4,3),
     effsize.ylabel="Effect size")
dev.off()


## Figure 5D - GABA

setwd('C:/Users/optim/Downloads/2020_1115_GABA/anadata8e_00-0')

dat <- read.table('2020_1115_GABA.trk',header=T,sep=',')

dat$id <- paste0(substr(dat$Group,13,15),substr(dat$Group,nchar(dat$Group)-3,nchar(dat$Group)))
dat <- as_tibble(dat)
dat$ti_raw <- (dat$Blk24+dat$Blk25)/2 - dat$baseline
dat$ti_raw <- 0.1*dat$ti_raw
dat$ti_norm <- dat$ti_raw/dat$baseline

timeseq <- seq(5,210,5)
TI <- data.matrix(dat[,26:ncol(dat)])

recov <- {}

for (i in 1:nrow(TI)){
  
  halfway <- TI[i,1] + (dat$baseline[i]-TI[i,1])/2
  idx <- which(TI[i,] > halfway)
  if (isempty(idx)){
    idx = length(timeseq)
  }
  
  if (is.nan(idx) == F){
    recov <- append(recov,timeseq[idx[1]])
  }
  
}

dat$recov <- recov
dat <- dat[which(dat$baseline >= 2),]

df<- dat %>% 
  group_by(id) %>%
  summarise(base = median(baseline),disp24 = median(Blk24), ti_raw = median(ti_raw),ti_norm = median(ti_norm),ti_recov = median(recov,na.rm=T))

df$cond <- substr(df$id,1,1)

gaba.shared <- 
  df %>%
  dabest(cond,ti_raw,
         idx = c("c","d","e","f"),
         paired=F)

gaba.shared.meandiff <- gaba.shared %>% mean_diff()

setwd('C:/Users/optim/Desktop/ArrestManuscript/OngoingFigureData/')

svg(filename="Fig5D_GABA.svg",
    width=5,
    height=4,
    pointsize=12)
plot(gaba.shared.meandiff, rawplot.type = "swarmplot",
     rawplot.ylabel="Change in speed (mm/s)",
     effsize.ylabel="Mean difference")
dev.off()


ggplot(df,aes(x=cond,y=ti_raw))+
  geom_boxplot()+
  geom_beeswarm()+
  ylab("Change in speed (mm/s)") +
  xlab(" ") +
  scale_x_discrete(labels=c("Control","25 uM","50 uM","100 uM"))+
  ylim(-5,2.5) + 
  theme_classic()+
  theme(text=element_text(size=24,colour='black',family="Arial",face='bold')) + 
  theme(axis.text=element_text(size=18,colour='black',family="Arial",face='plain'))

g <- aov(df$ti_recov ~ df$cond)
summary(g)

g <- aov(df$ti_norm ~ df$cond)
summary(g)

g <- aov(df$ti_raw ~ df$cond)
summary(g)

setwd('C:/Users/optim/Desktop/ArrestManuscript/OngoingFigureData/')

ggsave("Fig5D_GABAseries.svg", plot = last_plot(),
       scale = 1,width=13,height=15,units="cm",
       dpi = 300, limitsize = TRUE)

# Figure 6F - ICL

data <- read.table('Fig5F_xa334ICL.csv',header=T,sep=',')

data$displ <- paste0(substr(data$displ,13,21),substr(data$displ,29,32))
matx <- data.matrix(data)
rec <- matx[,26:ncol(matx)]
matx <- matx[,2:27]

base <- rowMeans(matx[,1:24])
ti <- data$Blk24

dspeed <- 0.1*(ti - base)
base <- 0.1*base

df <- data.frame(data$displ,base,dspeed)
colnames(df) <- c('fish','baseline','dspeed')
df <- as_tibble(df)

dat <- df %>%
  group_by(fish) %>%
  summarise(dspeed = mean(dspeed),baseline = mean(baseline))

dat$group <- substr(dat$fish,7,8)

ICL <- rbind(dat[which(dat$group == 'c0'),],dat[which(dat$group == 'd0'),])

t <- aov(ICL$baseline ~ ICL$group)
summary(t)

ICL <- ICL[which(ICL$baseline >= 3 & ICL$baseline <= 4),]

t <- aov(ICL$dspeed ~ ICL$group)
summary(t)


ICLmut <- rbind(dat[which(dat$group == 'e0'),],dat[which(dat$group == 'f0'),])

t <- aov(ICLmut$baseline ~ ICLmut$group)
summary(t)

ICLmut <- ICLmut[which(ICLmut$baseline >= 3 & ICLmut$baseline <= 4),]

t <- aov(ICLmut$dspeed ~ ICLmut$group)
summary(t)

xa334.shared <- 
  dat %>%
  dabest(group,dspeed,
         idx = c("r0","r1"),
         paired=F)

xa334.shared.effsize <- cohens_d(xa334.shared)

plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Change in speed (mm/s)",
     rawplot.ylim = c(-3,2),
     effsize.ylabel="Effect size")

svg(filename="Fig2C_xa266_estimplot.svg",
    width=4,
    height=5,
    pointsize=14)
plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Change in speed (mm/s)",
     rawplot.ylim = c(-3,2),
     effsize.ylabel="Effect size")
dev.off()