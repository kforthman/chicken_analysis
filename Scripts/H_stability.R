data <- read.csv("/Volumes/T1000/Analysis/kforthman/Chicken_Task/Chicken_code/Output/Analysis/adaptivityModelFits_stabilityTest.csv")

H_stability <- matrix(nrow = 0,ncol = 1)
H_range <- matrix(nrow = 0,ncol = 1)
for(i in 1:dim(data)[1]){
  H <- as.numeric(data[i,6:16])
  H_stability[i] <- var(H)
  H_range[i] <- max(H)-min(H)
}
hist(H_stability, breaks = 30)
hist(H_range, breaks = 30)

noise <- as.numeric(data[,"noise_in_DV_5"])
lapse <- as.numeric(data[,"lapse_rate_5"])

plot(noise, lapse)#, xlim = c(0,.5), ylim = c(0,.1))
nl_rela <- lm(lapse~noise)
abline(nl_rela)
summary(nl_rela)
cor(noise, lapse)

plot(H_range, lapse)#, xlim = c(0,.5), ylim = c(0,.1))
rl_rela <- lm(lapse~H_range)
abline(rl_rela)
summary(rl_rela)
cor(H_range, lapse)

sig <- c(33,140,33,140,33,140)
h <- c(0.05,0.05,0.2,0.2,0.95,0.95)

stab_block <- matrix(nrow = length(H_stability), ncol = 6)
for(b in 1:6){
  ind <- which(data["block"] == b)
  stab_block[,b] <- H_stability[ind]
  #hist(H_stability[ind], breaks = 30, main = paste0("Block: ", b, ", Sigma: ", sig[b], ", H: ", h[b]))
  #hist(H_range[ind], breaks = 30, main = b)
}

colnames(stab_block) <- c(1,2,3,4,5,6)
boxplot(stab_block)

for(b in 1:6){
  ind <- which(data["block"] == b)
  H <- as.numeric(data[ind,"H_subjective_5"])
  plot(lapse[ind], H, main = paste0("Block: ", b, ", Sigma: ", sig[b], ", H: ", h[b]))#, xlim = c(0,.5), ylim = c(0,.1))
  abline(h[b], 0, col = "red")
  abline(0.5, 0, col = "blue")
}

