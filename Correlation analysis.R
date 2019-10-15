require(tidyr)
require(stringr)
require(ggplot2)
require(ncf)
require(grDevices)
require(WaveletComp)
require(biwavelet)
require(survival)

#mSynch is the function to estimate the mean (cross-)correlation in a spatiotemporal dataset as discussed in Bjornstad et al. (1999). The function requires multiple observations at each location.
mSynchs <- data.frame("1"=rep(NA,14),"2"=rep(NA,14),"3"=rep(NA,14),"4"=rep(NA,14),"5"=rep(NA,14),
                      "6"=rep(NA,14),"7"=rep(NA,14),"8"=rep(NA,14),"9"=rep(NA,14),"10"=rep(NA,14))
corrs <- data.frame("1"=rep(NA,22500),"2"=rep(NA,22500),"3"=rep(NA,22500),"4"=rep(NA,22500),"5"=rep(NA,22500),"6"=rep(NA,22500),
                    "7"=rep(NA,22500),"8"=rep(NA,22500),"9"=rep(NA,22500),"10"=rep(NA,22500),"11"=rep(NA,22500),"12"=rep(NA,22500),
                    "13"=rep(NA,22500),"14"=rep(NA,22500))


for (i in 1:14){
  for (j in 1:50){ #10
    time_series <- read.csv(paste0(directory,'time_series',i,'_',j-1,'.csv'),header=TRUE)
    time_series <- as.data.frame(time_series[2:501,])
    time_series_corr <- sapply(time_series, as.numeric)
    corr <- cor(time_series_corr,method="pearson")
    time_series_corr_transpose <- t(time_series_corr)
    mean_corr <- mSynch(time_series_corr_transpose,na.rm=TRUE)
    mSynchs[i,j] <- mean_corr$real
    cat(i,sum(time_series_corr),"\n")
  }
}

mSynchs <- as.data.frame(cbind(seq(1,14,1),mSynchs))
names(mSynchs) <- c("Incubation_Period","Correlation1","Correlation2","Correlation3","Correlation4","Correlation5")
mSynchs$mean <- rowMeans(mSynchs[,2:11],na.rm=TRUE)
write.csv(mSynchs,paste0(directory,'mSynchs.csv'))
mSynchs <- read.csv(paste0(directory,'mSynchs.csv'))
ggplot(mSynchs,aes(x=Incubation_Period,y=mean)) + geom_point(aes(x=Incubation_Period,y=mean,color=factor(Incubation_Period))) + 
  theme(panel.background=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=14)) + 
  geom_smooth(method='lm',formula=y~x) + xlab("Incubation Period") + ylab("Correlation") + theme(legend.position = "none")
summary(lm(mSynchs$mean~mSynchs$Incubation_Period))
# 
# x_coord <- rep(seq(-12.92,-10.4,0.28),15)
# y_coord<- c(rep(9.9,10),rep(9.7,10),rep(9.5,10),rep(9.3,10),rep(9.1,10),rep(8.9,10),rep(8.7,10),rep(8.5,10),rep(8.3,10),rep(8.1,10),rep(7.9,10),rep(7.7,10),rep(7.5,10),rep(7.3,10),rep(7.1,10))

x_coord <- rep(seq(1,10,1),15)
y_coord<- rep(seq(1,15,1),10)
Correlogs <- data.frame("1"=rep(NA,10),"2"=rep(NA,10),"3"=rep(NA,10),"4"=rep(NA,10),"5"=rep(NA,10),"6"=rep(NA,10),"7"=rep(NA,10),"8"=rep(NA,10),
                        "9"=rep(NA,10),"10"=rep(NA,10),"11"=rep(NA,10),"12"=rep(NA,10),"13"=rep(NA,10),"14"=rep(NA,10))
for (i in 1:14){
  Correlog <- data.frame("1"=rep(NA,10),"2"=rep(NA,10),"3"=rep(NA,10),"4"=rep(NA,10),"5"=rep(NA,10),"6"=rep(NA,10),"7"=rep(NA,10),"10"=rep(NA,10),"10"=rep(NA,10),"10"=rep(NA,10))
  for (j in 1:50){ #10
    time_series <- read.csv(paste0(directory,'time_series',i,'_',j-1,'.csv'),header=TRUE)
    time_series_corr <- as.data.frame(time_series[2:501,])
    time_series_corr <- sapply(time_series_corr, as.numeric)
    time_series_corr_transpose <- t(time_series_corr)
    mean_corr <- correlog.nc(x_coord,y_coord,time_series_corr_transpose,increment=2,na.rm=TRUE)
    Correlog[,j] <- mean_corr$correlation
    cat(i,j,"\n")
  }
  Correlogs[,i] <- apply(Correlog,1,median,na.rm=TRUE)
}

Correlogs_df <- data.frame(cbind(rep(as.vector(mean_corr$mean.of.class),14),
                                 c(Correlogs$X1,Correlogs$X2,Correlogs$X3,Correlogs$X4,Correlogs$X5,Correlogs$X6,Correlogs$X7,Correlogs$X8,
                                   Correlogs$X9,Correlogs$X10,Correlogs$X11,Correlogs$X12,Correlogs$X13,Correlogs$X14),
                                 rep(seq(1,14,1),each=10)))
names(Correlogs_df) <- c("Distance","Correlation","Incubation")

ggplot(Correlogs_df2,aes(x=Distance,y=Correlation))+ 
  geom_line(aes(x=Distance,y=Correlation,group=Incubation,color=factor(Incubation)))+xlab("Distance (steps on lattice)")+
  labs(color="Inc. Period") +theme(panel.background=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=14)) 
