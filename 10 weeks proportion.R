
# time_series top row is row, second row is column

chiefdoms <- c(1,1,2,2,3,3,3,6,6,6,
               1,1,2,2,3,3,3,6,6,6,
               4,4,4,4,5,5,5,6,6,6,
               4,4,4,4,5,5,5,6,6,6,
               7,9,9,9,10,10,10,11,11,11,
               8,9,9,9,10,10,10,11,11,11,
               8,9,9,9,10,10,10,12,12,12,
               8,9,9,9,10,10,10,12,12,12,
               13,13,13,13,15,15,15,16,16,16,
               13,13,13,13,15,15,15,16,16,16,
               13,13,14,14,15,15,15,16,16,16,
               13,13,14,14,15,15,15,16,16,16,
               17,17,17,17,18,18,18,19,19,19,
               17,17,17,17,18,18,18,19,19,19,
               17,17,17,17,18,18,18,19,19,19)

districts <- c(1,1,1,1,2,2,2,2,2,2,
              1,1,1,1,2,2,2,2,2,2,
              1,1,1,1,2,2,2,2,2,2,
              1,1,1,1,2,2,2,2,2,2,
              3,1,1,1,5,5,5,5,5,5,
              4,1,1,1,5,5,5,5,5,5,
              4,1,1,1,5,5,5,5,5,5,
              4,1,1,1,5,5,5,5,5,5,
              6,6,6,6,6,6,6,8,8,8,
              6,6,6,6,6,6,6,8,8,8,
              6,6,6,6,6,6,6,8,8,8,
              6,6,6,6,6,6,6,8,8,8,
              7,7,7,7,7,7,7,8,8,8,
              7,7,7,7,7,7,7,8,8,8,
              7,7,7,7,7,7,7,8,8,8)

lattice_summary_perc <- data.frame("1"=rep(NA,150),"2"=rep(NA,150),"3"=rep(NA,150),"4"=rep(NA,150),"5"=rep(NA,150),"6"=rep(NA,150),
                                    "7"=rep(NA,150),"8"=rep(NA,150),"9"=rep(NA,150),"10"=rep(NA,150),"11"=rep(NA,150),"12"=rep(NA,150),
                                    "13"=rep(NA,150),"14"=rep(NA,150))
chiefdom_summary_perc <- data.frame("1"=rep(NA,19),"2"=rep(NA,19),"3"=rep(NA,19),"4"=rep(NA,19),"5"=rep(NA,19),"6"=rep(NA,19),
                                 "7"=rep(NA,19),"8"=rep(NA,19),"9"=rep(NA,19),"10"=rep(NA,19),"11"=rep(NA,19),"12"=rep(NA,19),
                                 "13"=rep(NA,19),"14"=rep(NA,19))
district_summary_perc <- data.frame("1"=rep(NA,8),"2"=rep(NA,8),"3"=rep(NA,8),"4"=rep(NA,8),"5"=rep(NA,8),"6"=rep(NA,8),
                                 "7"=rep(NA,8),"8"=rep(NA,8),"9"=rep(NA,8),"10"=rep(NA,8),"11"=rep(NA,8),"12"=rep(NA,8),
                                 "13"=rep(NA,8),"14"=rep(NA,8))
for (i in 1:14){
  lattice_summary <- data.frame("1"=rep(NA,150),"2"=rep(NA,150),"3"=rep(NA,150),"4"=rep(NA,150),"5"=rep(NA,150),"6"=rep(NA,150),
                                 "7"=rep(NA,150),"8"=rep(NA,150),"9"=rep(NA,150),"10"=rep(NA,150))
  chiefdom_summary <- data.frame("1"=rep(NA,19),"2"=rep(NA,19),"3"=rep(NA,19),"4"=rep(NA,19),"5"=rep(NA,19),"6"=rep(NA,19),
                            "7"=rep(NA,19),"8"=rep(NA,19),"9"=rep(NA,19),"10"=rep(NA,19))
  district_summary <- data.frame("1"=rep(NA,8),"2"=rep(NA,8),"3"=rep(NA,8),"4"=rep(NA,8),"5"=rep(NA,8),"6"=rep(NA,8),
                                        "7"=rep(NA,8),"8"=rep(NA,8),"9"=rep(NA,8),"10"=rep(NA,8))
  for (j in 1:10){
    time_series <- read.csv(paste0(directory,'time_series',i,'_',j-1,'.csv'),header=TRUE)
    time_series <- time_series[2:nrow(time_series),]
    total <- apply(time_series,2,sum)
    ten_weeks_lattice_perc <- rep(NA,150)
    ten_weeks_chiefdom_perc <- rep(NA,19)
    ten_weeks_district_perc <- rep(NA,8)
    for (l in 1:150){
      time_series_lattice_aggregate <- time_series[,l]
      start <- as.vector(which(time_series_lattice_aggregate!=0))[1]
      if (!is.na(start)){
        end <- start+70
        ten_weeks_lattice_perc[l] <- sum(time_series_lattice_aggregate[start:end])/sum(time_series_lattice_aggregate)
      } else{
        ten_weeks_lattice_perc[l] <- 0
      }
    }
    for (k in 1:19){
      time_series_chiefdom <- time_series[,which(chiefdoms==k)]
      if (k!=7){
        time_series_chiefdom_aggregate <- apply(time_series_chiefdom,1,sum)
      } else{
        time_series_chiefdom_aggregate <- time_series_chiefdom
      }
      start <- as.vector(which(time_series_chiefdom_aggregate!=0))[1]
      end <- start+70
      ten_weeks_chiefdom_perc[k] <- sum(time_series_chiefdom_aggregate[start:end])/sum(time_series_chiefdom_aggregate)
    }
    for (d in 1:8){
      time_series_district <- time_series[,which(districts==d)]
      if (d!=3){
        time_series_district_aggregate <- apply(time_series_district,1,sum)
      } else{
        time_series_district_aggregate <- time_series_district
      }
      start <- as.vector(which(time_series_district_aggregate!=0))[1]
      end <- start+70
      ten_weeks_district_perc[d] <- sum(time_series_district_aggregate[start:end])/sum(time_series_district_aggregate)
    }
    lattice_summary[,j] <- ten_weeks_lattice_perc
    chiefdom_summary[,j] <- ten_weeks_chiefdom_perc
    district_summary[,j] <- ten_weeks_district_perc
  } 
  lattice_summary_perc[,i] <- apply(lattice_summary,1,mean)
  chiefdom_summary_perc[,i] <- apply(chiefdom_summary,1,mean)
  district_summary_perc[,i] <- apply(district_summary,1,mean)
}

chiefdom_av <- as.vector(apply(chiefdom_summary_perc,2,mean))
district_av <- as.vector(apply(district_summary_perc,2,mean))
lattice_av <- as.vector(apply(lattice_summary_perc,2,mean))

summary_perc <- as.data.frame(cbind(rep(seq(1,14,1),3),c(lattice_av,chiefdom_av,district_av),rep(c("Lattice","Chiefdom","District"),each=14)))
names(summary_perc) <- c("Inc. Period","Percent","Scale")

summary_perc <- summary_perc[summary_perc$Scale!="Chiefdom2",]

prop_cases <- ggplot(summary_perc,aes(x=Inc..Period,y=Percent,group=Scale)) + geom_line(aes(linetype=Scale))+
  scale_linetype_manual(values=c("dashed", "solid"))+
  theme(panel.background=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=14)) + 
  xlab("Inc. Period") + ylab("Proportion of cases reported by 10 weeks") + labs(linetype="Level of Analysis") + ylim(c(0,1))+
  theme(panel.background=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=16), legend.position = "bottom")
#ggsave("propcases.tiff",prop_cases,height=3,width=2.6,dpi=300)

