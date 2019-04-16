# First, set the parameters.
# h = .05, .2, .95
h <- seq(.05,.95,0.05)
s <- seq(30,200,10)

pc_mat <- matrix(nrow = length(h), ncol = length(s))
mag_mat <- matrix(nrow = length(h), ncol = length(s))
know_mat <- matrix(nrow = length(h), ncol = length(s))
llr_mat <- matrix(nrow = length(h), ncol = length(s))
psi_mat <- matrix(nrow = length(h), ncol = length(s))

pc_mat_avg <- matrix(0, nrow = length(h), ncol = length(s))
mag_mat_avg <- matrix(0, nrow = length(h), ncol = length(s))
know_mat_avg <- matrix(0, nrow = length(h), ncol = length(s))
llr_mat_avg <- matrix(0, nrow = length(h), ncol = length(s))
psi_mat_avg <- matrix(0, nrow = length(h), ncol = length(s))

rownames(pc_mat) <- h
colnames(pc_mat) <- s
rownames(pc_mat)[1] <- paste("H =",h[1])
colnames(pc_mat)[1] <- paste("Sig =",s[1])

rownames(mag_mat) <- rownames(pc_mat)
colnames(mag_mat) <- colnames(pc_mat)

rownames(know_mat) <- rownames(pc_mat)
colnames(know_mat) <- colnames(pc_mat)

rownames(llr_mat) <- rownames(pc_mat)
colnames(llr_mat) <- colnames(pc_mat)

rownames(psi_mat) <- rownames(pc_mat)
colnames(psi_mat) <- colnames(pc_mat)

ll <- function(this_x, this_mu, this_sigma) {
  0.5*(this_x - this_mu)^2/this_sigma^2 + 0.5*log(this_sigma^2)
}
# 
# xseq <- seq(-300,300)
# -2*(ll(x, 75, 140) - ll(x, -75, 140))

for(dummy in 1:20){
  print(dummy)
  for(k in 1:length(s)){
    for(j in 1:length(h)){
      mu <- 75;
      sigma <- s[k];
      hTrue <- h[j];
      hSubj <- hTrue;
      m <- 200;
      
      # Then, generate some random data to simulate the behavior of the eggs in each trial.
      r <- rnorm(m, mean = 0, sd = 1)
      switch <- ifelse(r > hTrue, 0,1)
      dist <- numeric()
      val <- 1;
      for(i in 1:m){
        if(switch[i] == 1){
          val <- ifelse(val == 1, 2, 1)
        }
        dist[i] <- val
      }
      
      x <- numeric()
      y <- numeric()
      for(i in 1:m){
        if(dist[i] == 1){
          x[i] <- rnorm(1, mu, sigma)
          y[i] <- rnorm(1, 0, sigma)
        }
        else{
          x[i] <- rnorm(1, -mu, sigma)
          y[i] <- rnorm(1, 0, sigma)
        }
      }
      
      # Now, based on the data above, the model will make a prediction of which chicken the egg came from.
      l <- numeric()
      psi <- numeric()
      llr <- numeric()
    
      l[1] <- 0
      for(n in 2:m){
        
       #bel <- sign(l[n - 1])*1
        bel <- ifelse(dist[i-1] == 1, 1, -1)
        
        psi[n] = bel + 
          log((1 - hSubj)/hSubj + exp(-bel)) - 
          log((1 - hSubj)/hSubj + exp( bel))
        
        llr[n] = -2*
          (
            ll(x[n], mu, sigma) -
              ll(x[n], -mu, sigma)
          )
        
        l[n] = psi[n] + llr[n];
      }
      
      # Percent correct
      pred <- ifelse(l > 0, 1, 2)
      pc <- sum(dist == pred)
      pc_mat[j,k] <- pc/m
      
      # Below is shown the number of times adding psi to the estimation did not affect the prediction.
      addbool <- sign(psi+llr) == sign(llr)
      #add_mat[j,k] <- sum(addbool[2:m])/(m-1)
      
      # Below is the number of times psi is greater in magnitude than llr
      greatbool <- abs(psi) > abs(llr)
      mag_mat[j,k] <- sum(greatbool[2:m])/(m-1)
      
      # Below is strength of knowledge from psi vs llr
      know <- abs(psi)
      know_mat[j,k] <- mean(know[2:m])
      
      llr_mat[j,k] <- mean(llr[2:m])
      psi_mat[j,k] <- mean(psi[2:m])
    }
  }
  pc_mat_avg <- pc_mat_avg + pc_mat
  mag_mat_avg <- mag_mat_avg + mag_mat
  know_mat_avg <- know_mat_avg + know_mat
  llr_mat_avg <- llr_mat_avg + llr_mat
  psi_mat_avg <- psi_mat_avg + psi_mat
}
pc_mat_avg <- pc_mat_avg/20
mag_mat_avg <- mag_mat_avg/20
know_mat_avg <- know_mat_avg/20

llr_mat_avg <- llr_mat_avg/20
psi_mat_avg <- psi_mat_avg/20

library(corrplot)
corrplot(mag_mat_avg, is.corr = FALSE, tl.srt = 90,title = "% psi magnitude > llr magnitude", mar=c(0,0,1,0))
corrplot(pc_mat_avg, is.corr = FALSE, cl.lim = c(0.5, 1), tl.srt = 90,title = "percent correct", mar=c(0,0,1,0))
corrplot(know_mat_avg, is.corr = FALSE, tl.srt = 90,title = "Avg psi magnitude", mar=c(0,0,1,0))
corrplot(llr_mat_avg, is.corr = FALSE,title = "Avg llr")
corrplot(psi_mat_avg, is.corr = FALSE,title = "Avg psi")
