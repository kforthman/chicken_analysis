library(ggplot2)
library(stringr)

if(!file.exists('/Volumes/T1000/Analysis/kforthman/Chicken_Task/Chicken_code/Output/Simulated_nl_Analysis/simAvg2.csv')){
  simFits <- read.csv("/Volumes/T1000/Analysis/kforthman/Chicken_Task/Chicken_code/Output/Simulated_nl_Analysis/adaptivityModelFits.csv")
  #'../Output/Simulated_nl_Analysis/subj_info.mat'
  
  simFits$pct <- simFits$pct/100
  
  getid <- function(i){
    return(paste(
      simFits$N[i],
      simFits$sigma[i],
      simFits$H_true[i],
      simFits$H_subj[i],
      simFits$noise[i],
      simFits$lapse[i],
      sep = '_'
    ))
  }
  
  simFits$subID <- getid(1:dim(simFits)[1])
  
  simTypes <- unique(simFits$subID)
  simAvg <- list()
  
  for(i in 1:length(simTypes)){
    this_simType <- simFits[which(simFits$subID == simTypes[i]),]
    
    H_subj_EST_avg <- mean(this_simType$H_subj_EST)
    H_subj_EST_stder<- sd(this_simType$H_subj_EST)/sqrt(10)
    
    noise_EST_avg <- mean(this_simType$noise_EST)
    noise_EST_stder<- sd(this_simType$noise_EST)/sqrt(10)
    
    lapse_EST_avg <- mean(this_simType$lapse_EST)
    lapse_EST_stder<- sd(this_simType$lapse_EST)/sqrt(10)
    
    pct_avg <- mean(this_simType$pct)
    pct_stder<- sd(this_simType$pct)/sqrt(10)
    
    simAvg <- rbind(simAvg, cbind(this_simType[1,c(1:5,7:8)], 
                                  H_subj_EST_avg, H_subj_EST_stder, 
                                  noise_EST_avg, noise_EST_stder,
                                  lapse_EST_avg, lapse_EST_stder,
                                  pct_avg, pct_stder))
  }
  
  simAvg_df <- as.data.frame(simAvg)
  
  write.csv(simAvg_df, file = '/Volumes/T1000/Analysis/kforthman/Chicken_Task/Chicken_code/Output/Simulated_nl_Analysis/simAvg2.csv')
}else{
  simAvg_df <- read.csv("/Volumes/T1000/Analysis/kforthman/Chicken_Task/Chicken_code/Output/Simulated_nl_Analysis/simAvg2.csv")
}


## Get MAE

Hsubj_mae <- matrix(nrow=12*13*6, ncol = 6)
index = 1
for(this_sigma in c(110,160)){
  for(this_Htrue in c(0.05, 0.95)){
    for(this_N in c(50,100,150)){
      for(this_noise in c(0.01, seq(0.5, 6, by = 0.5))){
        for(this_lapse in seq(0, 0.25, by = 0.05)){
          this_subset <- subset(simAvg_df, N == this_N & H_true == this_Htrue & sigma == this_sigma & abs(lapse-this_lapse) < 0.001 & noise == this_noise)
          this_mae <- mean(abs(this_subset$H_subj - this_subset$H_subj_EST_avg))
          Hsubj_mae[index,] <- c(this_N, this_Htrue, this_sigma, this_noise, this_lapse, this_mae)
          index = index+1
        }
      }
    }
  }
}

Hsubj_mae <- as.data.frame(Hsubj_mae)
colnames(Hsubj_mae) <- c("N", "H_true", "sigma", "noise", "lapse", "mae")
Hsubj_mae$N <- as.character(Hsubj_mae$N)

for(this_sigma in c(110,160)){
  for(this_Htrue in c(0.05, 0.95)){
    
    this_subset <- subset(Hsubj_mae, H_true == this_Htrue & sigma == this_sigma)
    
    filename = paste0("/Volumes/T1000/Analysis/kforthman/Chicken_Task/Chicken_code/Figures/", "MAEplot_Htrue-", this_Htrue, "_sigma-", this_sigma, ".png")
    png(filename = filename, width = 3000, height = 1000)
    
    ggp <- ggplot(this_subset, aes(x = lapse, y = mae, group = N, color = N)) +
      geom_point() + 
      geom_line(size = 3) + scale_color_manual(values=c("red", "darkred", "maroon1")) +
      ylim(0,.3) +
      facet_grid(.~noise, scales = "free") +
      labs(title = paste0("H_true = ", this_Htrue, " & sigma = ", this_sigma), subtitle = "noise", y = "MAE of the estimate of H_subj")+
      theme(plot.subtitle = element_text(hjust = 0.5))
    
    print(ggp)
    
    dev.off()
  }
}


# ## Create all the plots
# 
# for(this_sigma in c(110,160)){
#   for(this_Htrue in c(0.05, 0.95)){
#     for(this_N in c(50,100,150)){
#       for(this_noise in c(0.01, seq(0.5, 6, by = 0.5))){
#         for(this_lapse in seq(0, 0.25, by = 0.05)){
#           filename = paste0("/Volumes/T1000/Analysis/kforthman/Chicken_Task/Chicken_code/Figures/Simulation_results/H_subj/", "N-", str_pad(this_N, 3, pad = "0"), "_Htrue-", this_Htrue, "_sigma-", this_sigma, "_noise-", format(this_noise, nsmall = 2), "_lapse-", format(this_lapse, nsmall = 2), ".png")
#           png(filename = filename)
#           ggp <- ggplot(subset(simAvg_df, N == this_N & H_true == this_Htrue & sigma == this_sigma & noise == this_noise & abs(lapse-this_lapse) < 0.001), aes(x=H_subj, y=H_subj_EST_avg)) +
#             geom_errorbar(aes(ymin=H_subj_EST_avg-H_subj_EST_stder, ymax=H_subj_EST_avg+H_subj_EST_stder), width=.1) +
#             geom_point()+
#             ylim(0,1)+
#             geom_abline(intercept = 0)+
#             ggtitle(paste0("H_true = ", this_Htrue, " & sigma = ", this_sigma, " & N = ", this_N, " & noise = ", this_noise, " & lapse = ", this_lapse))
#           print(ggp)
#           dev.off()
#         }
#       }
#     }
#   }
# }
# 
# 
# for(this_sigma in c(110,160)){
#   for(this_Htrue in c(0.05, 0.95)){
#     for(this_N in c(50,100,150)){
#       for(this_H_subj in seq(0.05, .95, by = 0.05)){
#         for(this_lapse in seq(0, 0.25, by = 0.05)){
#           filename = paste0("/Volumes/T1000/Analysis/kforthman/Chicken_Task/Chicken_code/Figures/Simulation_results/noise/", "N-", str_pad(this_N, 3, pad = "0"), "_Htrue-", this_Htrue, "_sigma-", this_sigma, "_Hsubj-", format(this_H_subj, nsmall = 2), "_lapse-", format(this_lapse, nsmall = 2), ".png")
#           png(filename = filename)
#           ggp <- ggplot(subset(simAvg_df, N == this_N & H_true == this_Htrue & sigma == this_sigma & abs(H_subj-this_H_subj) < 0.001 & abs(lapse-this_lapse) < 0.001), aes(x=noise, y=noise_EST_avg)) +
#             geom_errorbar(aes(ymin=noise_EST_avg-noise_EST_stder, ymax=noise_EST_avg+noise_EST_stder), width=.1) +
#             geom_point()+
#             ylim(0,6)+
#             geom_abline(intercept = 0)+
#             ggtitle(paste0("H_true = ", this_Htrue, " & sigma = ", this_sigma, " & N = ", this_N, " & H_subj = ", this_H_subj, " & lapse = ", this_lapse))
#           print(ggp)
#           dev.off()
#         }
#       }
#     }
#   }
# }
# 
# for(this_sigma in c(110,160)){
#   for(this_Htrue in c(0.05, 0.95)){
#     for(this_N in c(50,100,150)){
#       for(this_H_subj in seq(0.05, .95, by = 0.05)){
#         for(this_noise in c(0.01, seq(0.5, 6, by = 0.5))){
#           filename = paste0("/Volumes/T1000/Analysis/kforthman/Chicken_Task/Chicken_code/Figures/Simulation_results/lapse/", "N-", str_pad(this_N, 3, pad = "0"), "_Htrue-", this_Htrue, "_sigma-", this_sigma, "_Hsubj-", format(this_H_subj, nsmall = 2), "_noise-", format(this_noise, nsmall = 2), ".png")
#           png(filename = filename)
#           subset(simAvg_df, N == this_N & H_true == this_Htrue & sigma == this_sigma &  abs(H_subj-this_H_subj) < 0.001 & noise == this_noise)
#           ggp <- ggplot(subset(simAvg_df, N == this_N & H_true == this_Htrue & sigma == this_sigma &  abs(H_subj-this_H_subj) < 0.001 & noise == this_noise), aes(x=lapse, y=lapse_EST_avg)) + 
#             geom_errorbar(aes(ymin=lapse_EST_avg-lapse_EST_stder, ymax=lapse_EST_avg+lapse_EST_stder), width=.1) +
#             geom_point()+
#             ylim(0,.45)+
#             geom_abline(intercept = 0)+
#             ggtitle(paste0("H_true = ", this_Htrue, " & sigma = ", this_sigma, " & N = ", this_N, " & H_subj = ", this_H_subj, " & noise = ", this_noise))
#           print(ggp)
#           dev.off()
#         }
#       }
#     }
#   }
# }