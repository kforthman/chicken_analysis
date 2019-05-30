library(ggplot2)
simFits <- read.csv("/Volumes/T1000/Analysis/kforthman/Chicken Task/Chicken_code/Data/Simulated_Analysis/adaptivityModelFits.csv")

simTypes <- unique(substr(simFits$subID,1,13))
simAvg <- list()

for(i in 1:length(simTypes)){
  this_simType <- simFits[which(substr(simFits$subID,1,13) == simTypes[i]),]
  H_subj_EST_avg <- mean(this_simType$H_subj_EST)
  H_subj_EST_stder<- sd(this_simType$H_subj_EST)/sqrt(10)
  simAvg <- rbind(simAvg, cbind(this_simType[1,1:5], H_subj_EST_avg, H_subj_EST_stder))
}

simAvg_df <- as.data.frame(simAvg)

for(this_sigma in c(30,110,160)){
  for(this_Htrue in c(0.05, 0.95)){
    ggp <- ggplot(subset(simAvg_df, N == 150 & H_true == this_Htrue & sigma == this_sigma), aes(x=H_subj, y=H_subj_EST_avg)) + 
      geom_errorbar(aes(ymin=H_subj_EST_avg-H_subj_EST_stder, ymax=H_subj_EST_avg+H_subj_EST_stder), width=.1) +
      geom_point()+
      ylim(0,1)+
      geom_abline(intercept = 0)+
      ggtitle(paste0("H_true = ", this_Htrue, " & sigma = ", this_sigma))
    print(ggp)
  }
}
