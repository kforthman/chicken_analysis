library(ggplot2)
library(randomcoloR)
library(psych)

# Set working directory
setwd("/Volumes/T1000/Analysis/kforthman/Chicken Task/Chicken_code")

# Read in the data.
data <- read.csv("Data/Analysis/adaptivityModelFits.csv")
# Change percent to decimal
data[,'pct'] = data[,'pct']/100

# Now load in the pattern information
pat_data <- read.csv("Data/blockPatterns.csv")
pat_T1 <- pat_data[which(pat_data[,'t']=='T1'),]

# These are the h and sig for each block, in block order.
sig <- c(33,140,33,140,33,140)
h <- c(0.05,0.05,0.2,0.2,0.95,0.95)
# Below are some clearer definitions of what each h and sig entails
sigdef <- c("Small Egg Distribution", "Large Egg Distribution","Small Egg Distribution", "Large Egg Distribution","Small Egg Distribution", "Large Egg Distribution")
hdef <- c("Very Infrequent Switching","Very Infrequent Switching", "Infrequent Switching", "Infrequent Switching", "Very Frequent Switching", "Very Frequent Switching")

# how the pattern will be represented in the plot
pattern <- c("turquoise3", "seagreen4", "slateblue", "maroon")


# Block-wise subjective (fit) H vs objective H. Solid line is a least-squares fit. 
png(paste0("Figures/H_Default_Comparison_to_H_true.png"))
plot(data$H_true, data$H_subjective, xlim = c(0,1), ylim = c(0,1))
tru2def_fit <- lm(data$H_subjective~data$H_true)
abline(tru2def_fit)
dev.off()

t_cat <- levels(unique(data$trial))
sig_cat <- unique(data$sigma)
for(i in 1:length(t_cat)){
  for(j in 1:length(sig_cat)){
    index <- which(data$trial == t_cat[i] & data$sigma == sig_cat[j])
    png(paste0("Figures/H_Default_Comparison_to_H_true-", t_cat[i], "-", sig_cat[j],".png"))
    plot(jitter(data$H_true[index]), data$H_subjective[index], xlim = c(0,1), ylim = c(0,1), 
         main = paste0("Trial: ", t_cat[i], ", Sigma: ", sig_cat[j]),
         xlab = "True H",
         ylab = "Subjective H")
    tru2def_fit <- lm(data$H_subjective[index]~data$H_true[index])
    abline(tru2def_fit)
    #text(.9,0,paste0("Adj R2 = ", round(signif(summary(tru2def_fit)$adj.r.squared, 5), 3)))
    dev.off()
  }
}

# Block-wise subjective (fit) H vs objective H.
id <- levels(unique(data$ID))
png(paste0("Figures/H_Default_Comparison_to_H_true-spaghetti.png"))
plot(data$H_true, data$H_subjective, xlim = c(0,1), ylim = c(0,1))
for(i in 1:length(id)){
  index <- which(data$ID == id[i])
  tru2def_fit <- lm(data$H_subjective[index]~data$H_true[index])
  abline(tru2def_fit, col = randomColor(1))
}
dev.off()

# The following plots show the default H compared at T1 and T2. The wheels marking each point indicate the percent correct. The length of the + spokes indicate percent correct at T1 and the length of the x spokes indicate percent correct at T2. If the spokes reach the edge of the wheel, that means the participant had ~100% accuracy. The pink line is the least-squares fit. The dashed line is just a 45 degree line.
for(i in 1:6){
  
  t1 <- which(data[,2] == 'T1' & data[,3]==i)
  t2 <- which(data[,2] == 'T2' & data[,3]==i)
  
  fit <- lm(data[t2,6]~data[t1,6])
  
  icc <- ICC(cbind(data[t1,6],data[t2,6]))
  
  gg <- ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    stat_smooth(method = "lm", se = F, colour = "pink") + # Linear regression
    geom_abline(intercept = 0, colour = "slateblue4", linetype = "dashed") + # 45 degree line
    geom_point(colour = pattern[pat_T1[,'pattern']], size = data[t1,'pct']*4, pch = 3) + # wheels
    geom_point(colour = pattern[pat_T1[,'pattern']], size = data[t2,'pct']*4, pch = 4) +
    geom_point(colour = pattern[pat_T1[,'pattern']], size = 6, pch = 1)+ scale_fill_continuous(guide = guide_legend()) +
    # geom_point(colour = rgb(data[t1,11], 0, data[t2,11]), pch = 16, size = 5) + #color
    labs(title = paste0("Block ", i, ": H = ", h[i], ", sigma = ", sig[i]),
         caption = paste0("Adj R2 = ", round(signif(summary(fit)$adj.r.squared, 5), 3),
                          ", Intercept =", round(signif(fit$coef[[1]],5 ), 3),
                          ", Slope =", round(signif(fit$coef[[2]], 5), 3),
                          ", P =", round(signif(summary(fit)$coef[2,4], 5), 3),
                          ",\nICC (single absolute) = ", round(icc$results$ICC[1],3),
                          ", ICC (single fixed) = ", round(icc$results$ICC[3], 3)),
         subtitle = paste0(sigdef[i], ", ", hdef[i]))+
    xlim(0,1)+
    ylim(0,1)+
    xlab("T1 Subjective H")+
    ylab("T2 Subjective H")  + theme(legend.position = "bottom")
  
  png(paste0("Figures/H_Default_Comparison-Block_", i, ".png"),
      bg = "transparent")
  plot(gg)
  dev.off()
  
  
  #cool color palette graph
  # x <- rep(seq(0,1,length.out = 100),100)
  # y <- seq(0,1,length.out = 10000)
  # plot(x,y)
  # dummy <- data.frame(x,y)
  # 
  # gg <- ggplot(dummy, aes_string(x = x, y = y)) +
  #   geom_point(colour = rgb(x, 0, y), pch = 16, size = 5) + #color
  #   xlim(0,1)+
  #   ylim(0,1)+
  #   xlab("T1 Accuracy")+
  #   ylab("T2 Accuracy")
  # 
  # plot(gg)
  
  #plot just pct
  # fit <- lm(data[t2,11]~data[t1,11])
  # 
  # gg <- ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  #   geom_point(colour = "turquoise4") +
  #   stat_smooth(method = "lm", se = F, colour = "turquoise") +
  #   labs(title = paste0("Block ", i, ": H = ", h[i], ", sigma = ", sig[i]),
  #        caption = paste0("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
  #                         ", Intercept =",signif(fit$coef[[1]],5 ),
  #                         ", Slope =",signif(fit$coef[[2]], 5),
  #                         ", P =",signif(summary(fit)$coef[2,4], 5)),
  #        subtitle = paste0(sigdef[i], ", ", hdef[i]))+
  #   xlim(0,100)+
  #   ylim(0,100)+
  #   xlab("T1 Percent Correct")+
  #   ylab("T2 Percent Correct")
  # 
  # plot(gg)
}

# Distribution of parameter estimates across subjects
data <- read.csv("Data/Analysis/adaptivityModelFits.csv")
hist(as.matrix(data["noise_in_DV"]), breaks=30)
hist(as.matrix(data["lapse_rate"]), breaks=30)
hist(as.matrix(data["H_subjective"]))


