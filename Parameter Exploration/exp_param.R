# First, set the parameters.
# h = .05, .2, .95
h <- seq(.05,.95,0.05)
s <- seq(30,200,10)
mu <- 75;
library(corrplot)


genSimData <- function(batch_num, m, feedback = F, pred = F){
  dummy <- array(data = NA, dim=c( m, 8,length(h),length(s), batch_num), 
                 dimnames = list(1:m, c("dist","pred", "crct", "x", "y", "l", "psi", "llr"), h, s, 1:batch_num))
  tick <- 1
  for(k in 1:length(s)){
    for(j in 1:length(h)){
      for(batch in 1:batch_num){
        print(paste0(tick,"/",length(s)*length(h)*batch_num))
        tick <- tick+1
        
        sigma <- s[k];
        hTrue <- h[j];
        hSubj <- hTrue;
        
        
        # Then, generate some random data to simulate the behavior of the eggs in each trial.
        r <- runif(m, min = 0, max = 1) # create a string of random values between 0 and 1
        switch <- ifelse(r > hTrue, 0,1) # select those values greater than hTrue to get switch points.
        dist <- numeric() # initialize dist, the vector that indicates the correct chicken
        val <- 1;
        for(i in 1:m){
          if(switch[i] == 1){ # If switch
            val <- ifelse(val == 1, 2, 1) # value changes
          }
          dist[i] <- val # dist stores the correct answer
        }
        
        # Generate the x,y position based on the correct chicken
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
        ll <- function(this_x, this_mu, this_sigma) {
          0.5*(this_x - this_mu)^2/this_sigma^2 + 0.5*log(this_sigma^2)
        }
        l[1] <- 0
        psi[1] <- 0
        llr[1] <- 0
        for(n in 2:m){
          
          if(feedback){
            l[n-1] <- sign(l[n-1])
          }
          
          psi[n] = l[n - 1] + 
            log((1 - hSubj)/hSubj + exp(-l[n - 1])) - 
            log((1 - hSubj)/hSubj + exp( l[n - 1]))
          
          llr[n] = -2*
            (
              ll(x[n], mu, sigma) -
                ll(x[n], -mu, sigma)
            )
          
          l[n] = psi[n] + llr[n];
        }
        pred <- ifelse(pred, ifelse(psi > 0, 1, 2), ifelse(l > 0, 1, 2))
        crct <- ifelse(dist == pred, 1, 0)
        
        dummy[ , , j, k, batch] <- cbind(dist, pred, crct, x, y, l, psi, llr)
        
      }
    }
  }
  return(dummy)
}

compareParameters <- function(dummy, type, n, pred = F){
  
  # Percent correct
  if(type == "pc"){
    my_mat <- apply(apply(dummy[ ,"crct", , , ],c(2,3,4), mean),c(1,2),mean)
  }
  
  # Below is shown the number of times adding \[Psi] to the estimation did not affect the prediction.
  if(type == "add"){
    mybool <- sign(dummy[ ,"psi", , , ]+dummy[ ,"llr", , , ]) != sign(dummy[ ,"llr", , , ])
    my_mat <- apply(apply(mybool,c(2,3,4), sum),c(1,2),mean)
  }
  
  # Below is the number of times psi is greater in magnitude than llr
  if(type == "mag"){
    mybool <- abs(dummy[ ,"psi", , , ]) > abs(dummy[ ,"llr", , , ])
    my_mat <- apply(apply(mybool,c(2,3,4), mean),c(1,2),mean)
  }
  
  # Below is strength of knowledge from psi vs llr
  if(type == "know"){
    mybool <- abs(dummy[ ,"psi", , , ])
    my_mat <- apply(apply(mybool,c(2,3,4), mean),c(1,2),mean)
  }
  
  corrplot(my_mat, is.corr = FALSE, tl.srt = 90,title = paste0(type, ", ntrials = ", n, ", pred = ", pred), mar=c(0,0,1,0), method = "number")
}
#corrplot(pc_mat_avg, is.corr = FALSE, cl.lim = c(0.5, 1), tl.srt = 90,title = "percent correct", mar=c(0,0,1,0))
#corrplot(know_mat_avg, is.corr = FALSE, tl.srt = 90,title = "Avg psi magnitude", mar=c(0,0,1,0))

dummy_data_50_pred <- genSimData(30, 50, pred = T)
dummy_data_100_pred <- genSimData(30, 100, pred = T)
dummy_data_150_pred <- genSimData(30, 200, pred = T)
dummy_data_200_pred <- genSimData(30, 200, pred = T)

dummy_data_50_est <- genSimData(30, 50, pred = F)
dummy_data_100_est <- genSimData(30, 100, pred = F)
dummy_data_150_est <- genSimData(30, 200, pred = F)
dummy_data_200_est <- genSimData(30, 200, pred = F)

compareParameters(dummy_data_150, "pc")
compareParameters(dummy_data_150, "add")

compareParameters(dummy_data_50_pred, "mag", n=50, pred = T)
compareParameters(dummy_data_100_pred, "mag", n=100, pred = T)
compareParameters(dummy_data_150_pred, "mag", n=150, pred = T)
compareParameters(dummy_data_200_pred, "mag", n=200, pred = T)
compareParameters(dummy_data_50_est, "mag", n=50)
compareParameters(dummy_data_100_est, "mag", n=100)
compareParameters(dummy_data_150_est, "mag", n=150)
compareParameters(dummy_data_200_est, "mag", n=200)


genSimData2 <- function(m, sigma, hTrue, hSubj, feedback = F){
  dummy <- array(data = NA, dim=c(m, 8), 
                 dimnames = list(1:m, c("dist","pred", "crct", "x", "y", "l", "psi", "llr"))
  )
  
  # Then, generate some random data to simulate the behavior of the eggs in each trial.
  r <- runif(m, min = 0, max = 1) # create a string of random values between 0 and 1
  switch <- ifelse(r > hTrue, 0,1) # select those values greater than hTrue to get switch points.
  dist <- numeric() # initialize dist, the vector that indicates the correct chicken
  val <- 1;
  for(i in 1:m){
    if(switch[i] == 1){ # If switch
      val <- ifelse(val == 1, 2, 1) # value changes
    }
    dist[i] <- val # dist stores the correct answer
  }
  
  # Generate the x,y position based on the correct chicken
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
  ll <- function(this_x, this_mu, this_sigma) {
    0.5*(this_x - this_mu)^2/this_sigma^2 + 0.5*log(this_sigma^2)
  }
  l[1] <- 0
  psi[1] <- 0
  llr[1] <- 0
  for(n in 2:m){
    
    if(feedback){
      l[n-1] <- sign(l[n-1])
    }
    
    psi[n] = l[n - 1] + 
      log((1 - hSubj)/hSubj + exp(-l[n - 1])) - 
      log((1 - hSubj)/hSubj + exp( l[n - 1]))
    
    llr[n] = -2*
      (
        ll(x[n], mu, sigma) -
          ll(x[n], -mu, sigma)
      )
    
    l[n] = psi[n] + llr[n];
  }
  pred <- ifelse(l > 0, 1, 2)
  crct <- ifelse(dist == pred, 1, 0)
  
  dummy <- cbind(dist, pred, crct, x, y, l, psi, llr)
  
  return(dummy)
}

m <- c(50,100,150)
sigma <- c(110,160)
hTrue <- c(0.05,0.95)
hSubj <- c(0.05,0.30,0.50,0.70,0.95)
tick <- 1
for(batch in 1:10){
  for(h in 1:3){
    for(i in 1:2){
      for(j in 1:2){
        for(k in 1:5){
          print(tick)
          tick <- tick+1
          dummy <- genSimData2(m[h], sigma[i], hTrue[j], hSubj[k])
        }
      }
    }
  }
}
