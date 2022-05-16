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


df <- read.table('Fig4c_xa266laser_dspeed.csv',header=T,sep=',')
df$dSpeed <- 0.1*df$dSpeed

xa334.shared <- 
  df %>%
  dabest(Group,dSpeed,
         idx = c("c","r"),
         paired=F)

xa334.shared.effsize <- cohens_d(xa334.shared)

plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Change in speed (mm/s)",
     rawplot.ylim = c(-3,2),
     effsize.ylabel="Effect size")

svg(filename="Fig4C_xa266_laser_estimplot.svg",
    width=4,
    height=5,
    pointsize=14)
plot(xa334.shared.effsize, rawplot.type = "swarmplot",
     rawplot.ylabel="Change in speed (mm/s)",
     rawplot.ylim = c(-3,2),
     effsize.ylabel="Effect size")
dev.off()


# Statistical analysis

t.test(dSpeed ~ Group, data=df)

# Summary stats

summstats <- df %>%
  group_by(Group) %>%
  summarise(avg = mean(Blk24),std = std_err(Blk24))


## xa266 laser baseline locomotion

df <- read.table('xa266laser_baseline.csv',header=T)

ggplot(df, aes(x=Group,y=Avg)) + 
  #geom_bar(stat="identity") +
  #geom_errorbar(aes(ymin=NormDisp - Err, ymax=NormDisp + Err), width = 0.3) +
  geom_boxplot(fill = c('gray80','magenta'))+
  geom_dotplot(dotsize=0.5,binwidth=1,binaxis='y',stackdir='center',fill='gray80')+
  ylab("Baseline speed") +
  ylim(0,40) + 
  theme(text=element_text(size=16,family="Serif")) + 
  theme_classic()

ggsave("xa266laser_baseline.eps", plot = last_plot(), device = "eps",
       scale = 1,width=10,height=15,units="cm",
       dpi = 300, limitsize = TRUE)