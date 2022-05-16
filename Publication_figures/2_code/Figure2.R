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


setwd('D:/TransferInhibition/')
data <- read.table('CircuitBreakingScreen.csv',header = TRUE, sep = "\t")

ggplot(data,aes(x=EffectSize))+
  geom_histogram(position="identity", alpha=1,binwidth=0.06)+
  scale_x_continuous(breaks=seq(-2,1,0.5))+
  scale_y_continuous(breaks=seq(0,5,1))+
  xlab('Normalized inhibition effect size')+
  theme(text=element_text(size=16,family="Serif"))+ 
  theme_classic()

setwd('C:/Users/Bhandiwad/Desktop/OngoingFigureData')

ggsave("CircuitBreakingScreen_Histogram.eps", plot = last_plot(), device = "eps",
       scale = 1,width=15,height=10,units="cm",
       dpi = 300, limitsize = TRUE)


############## Single line t0 information ###################################

# Figure 2c - xa239 (y318)

data <- read.table('Fig2c_xa239_dSpeed.csv',header=T,sep=',')

data$displ <- paste0(substr(data$displ,13,18),substr(data$displ,26,28))
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

svg(filename="Fig2C_xa239_estimplot.svg",
    width=4,
    height=5,
    pointsize=14)
plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Change in speed (mm/s)",
     rawplot.ylim = c(-3,2),
     effsize.ylim = c(-2,4),
     effsize.ylabel="Effect size")
dev.off()



# Figure 2d - xa239 time to recovery

timeseq <- seq(5,210,5)
recov <- {}

for (i in 1:nrow(rec)){
  
  halfway <- rec[i,1] + (dat$baseline-rec[i,1])/2 
  idx <- which(rec[i,]>halfway)
  if (isempty(idx)){
    idx = length(timeseq)
  }
  
  if (is.nan(idx) == F){
    recov <- append(recov,timeseq[idx[1]])
  }
  
}

df <- as_tibble(data.frame(data$displ,recov))
colnames(df) <- c('fish','recov')

d <- ddply(df,'fish',summarise,tau=mean(recov))
d$group <- substr(d$fish,1,2)

xa334.shared <- 
  d %>%
  dabest(group,tau,
         idx = c("r0","r1"),
         paired=F)

xa334.shared.effsize <- cohens_d(xa334.shared)

plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Half recovery time (s)",
     rawplot.ylim= c(0,200),
     effsize.ylim = c(-1.5,3),
     effsize.ylabel="Effect size")

svg(filename="Fig2C_xa239_tau_estimplot.svg",
    width=4,
    height=5,
    pointsize=14)
plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Half recovery time (s)",
     rawplot.ylim= c(0,200),
     effsize.ylim = c(-1.5,3),
     effsize.ylabel="Effect size")
dev.off()

t <- aov(d$tau ~ d$group)
summary(t)

# Figure 2c - xa334 (y334)

data <- read.table('Fig2c_xa334_dSpeed.csv',header=T,sep=',')

data$Group <- paste0(substr(data$Group,13,18),substr(data$Group,26,28))
matx <- data.matrix(data)
rec <- matx[,26:ncol(matx)]
matx <- matx[,2:27]

base <- rowMeans(matx[,1:24])
ti <- data$Blk25

dspeed <- 0.1*(ti - base)
base <- 0.1*base

df <- data.frame(data$Group,base,dspeed)
colnames(df) <- c('fish','baseline','dspeed')
df <- as_tibble(df)

dat <- df %>%
  group_by(fish) %>%
  summarise(dspeed = mean(dspeed),baseline = mean(baseline))

dat$group <- substr(dat$fish,1,2)

xa334.shared <- 
  dat %>%
  dabest(group,dspeed,
         idx = c("r0","r1"),
         paired=F)

xa334.shared.effsize <- cohens_d(xa334.shared)



svg(filename="Fig2C_xa334_estimplot.svg",
    width=4,
    height=5,
    pointsize=14)
plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Change in speed (mm/s)",
     rawplot.ylim = c(-3,2),
     effsize.ylabel="Effect size")
dev.off()

t <- aov(dat$dspeed ~ dat$group)
summary(t)

# ggplot(dat,aes(x=group,y=dspeed))+
#   geom_boxplot(fill = c('gray80','magenta'),color='black',size=1)+
#   geom_beeswarm(size=4,colour='gray20',alpha=0.96,cex=4)+
#   ylab("Change in speed (mm/s)") +
#   xlab(" ") +
#   scale_x_discrete(labels=c("Control","Ablated"))+
#   ylim(-2.5,1.5) + 
#   theme_classic()+
#   theme(text=element_text(size=24,colour='black',family="Arial",face='bold')) + 
#   theme(axis.text=element_text(size=18,colour='black',family="Arial",face='plain'))
# 
# ggsave("Fig2c_xa334_beeswarm.svg", plot = last_plot(), 
#        scale = 1,width=10,height=20,units="cm",
#        dpi = 300, limitsize = TRUE)

# Figure 2d - xa334 time to recovery

timeseq <- seq(5,210,5)
recov <- {}

for (i in 1:nrow(rec)){
  
  halfway <- rec[i,1] + (dat$baseline-rec[i,1])/2 
  #halfway <- TI[i,1] - TI[i,1]/2 # for xa266 ablated
  idx <- which(rec[i,]>halfway)
  if (isempty(idx)){
    idx = length(timeseq)
  }
  
  if (is.nan(idx) == F){
    recov <- append(recov,timeseq[idx[1]])
  }
  
}

df <- as_tibble(data.frame(data$Group,recov))
colnames(df) <- c('fish','recov')

d <- ddply(df,'fish',summarise,tau=mean(recov))
d$group <- substr(d$fish,1,2)

d <- d[which(is.na(d$tau)==F),]

xa334.shared <- 
  d %>%
  dabest(group,tau,
         idx = c("r0","r1"),
         paired=F)

xa334.shared.effsize <- cohens_d(xa334.shared)

plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Half recovery time (s)",
     rawplot.ylim= c(0,200),
     effsize.ylim = c(-1.5,1),
     effsize.ylabel="Effect size")

svg(filename="Fig2C_xa334_tau_estimplot.svg",
    width=4,
    height=5,
    pointsize=14)
plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Change in speed (mm/s)",
     rawplot.ylim= c(0,200),
     effsize.ylabel="Effect size")
dev.off()

t <- aov(d$tau ~ d$group)
summary(t)

# Figure 2c - xa266 (y405)

data <- read.table('Fig2c_xa266_dSpeed.csv',header=T,sep=',')

data$displ <- paste(substr(data$displ,13,18),substr(data$displ,26,28))
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

dat$group <- substr(dat$fish,1,2)

t <- aov(dat$dspeed ~ dat$group)
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



## Fig 2d - xa266 time to recovery

timeseq <- seq(5,210,5)
recov <- {}

for (i in 1:nrow(rec)){
  
  halfway <- rec[i,1] + (dat$baseline-rec[i,1])/2 
  
  if (halfway<= dat$baseline) {
    idx <- which(rec[i,]>halfway)
    if (isempty(idx)){
      idx = length(timeseq)
    }
    
    if (is.nan(idx) == F){
      recov <- append(recov,timeseq[idx[1]])
    }
    
  } else {recov <- append(recov,timeseq[idx[1]])}
}

df <- as_tibble(data.frame(data$displ,recov))
colnames(df) <- c('fish','recov')

d <- ddply(df,'fish',summarise,tau=mean(recov))
d$group <- substr(d$fish,1,2)



xa334.shared <- 
  d %>%
  dabest(group,tau,
         idx = c("r0","r1"),
         paired=F)

xa334.shared.effsize <- cohens_d(xa334.shared)

plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Half recovery time (s)",
     rawplot.ylim= c(0,200),
     effsize.ylabel="Effect size")

svg(filename="Fig2C_xa266_tau_estimplot.svg",
    width=4,
    height=5,
    pointsize=14)
plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Half recovery time (s)",
     rawplot.ylim= c(0,200),
     effsize.ylabel="Effect size")
dev.off()

t <- aov(d$tau ~ d$group)
summary(t)
