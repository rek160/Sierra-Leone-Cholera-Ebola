# Overlap function  
FUN <- function(x){
  sqrt(prod(x))
}

# create overlap dataframe to store overlap for each time point across 14 incubation periods
overlap <- array(NA,dim=c(500,14))

total_pop <- 45000
#for (j in 1:14){ #14
cat("j",j,"\n")
for (t in 1:100){ #500
  #cat("t",t,"\n")
  pis_summary <- c()
  for (i in 1:50){ #nsim
    for (k in 1:50){ #nsim
      pis <- array(NA,dim=c(2,150))
      if (i<k){
        cat(i,k,"\n")
        active_agents_1 <- read.csv(paste0(directory,'active_agents',j,'_',i-1,'.csv'),header=TRUE)
        active_agents_1 <- as.data.frame(active_agents_1[2:501,])
        active_agents_2 <- read.csv(paste0(directory,'active_agents',j,'_',k-1,'.csv'),header=TRUE)
        active_agents_2 <- as.data.frame(active_agents_2[2:501,])
        for (x in 1:150){
          pis[1,x] <- active_agents_1[t,x]/sum(active_agents_1[t,])
          pis[2,x] <- active_agents_2[t,x]/sum(active_agents_2[t,])
        }
        sim_active <- apply(pis,2,FUN)
        sim_active_sum <- ifelse(all(is.na(sim_active)),NA,sum(sim_active,na.rm=TRUE))
        pis_summary <- c(pis_summary,sim_active_sum)
      }
    }
  }
  overlap[t,j] <- mean(pis_summary,na.rm=TRUE)
  
}

write.csv(overlap,paste0('overlap2_',j,'.csv'))