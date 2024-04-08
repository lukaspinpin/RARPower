library(foreach)
library(doParallel)
library(doRNG)
library(Rfast)

simulation <- function(N=200, dist="bern", K=1, par1 = c(0,0),  par2 = NULL, burnin=3, ar=NULL, nsim= 10^4, method="ER", measure="sd"){
  
  #Generate Result Data Frame
  colnames <-  c()
  for(i in 1:ncol(par1)){
    colnames <- c(colnames, paste("p1", i, sep = "", collapse = NULL))
  }
  if(!is.null(par2)){
    for(i in 1:ncol(par2)){
      colnames <- c(colnames, paste("p2", i, sep = "", collapse = NULL))
    }
  }
  for(j in 1:(K+1)){
    colnames <- c(colnames, paste( "n", j-1, sep = "", collapse = NULL))
    colnames <- c(colnames, paste("Var(n", j-1, ")",sep = "", collapse = NULL))
  }
  colnames <- c(colnames, "EMR", "Var_EMR", "Z_Type1", "BM_Type1", "Z_Power", "BM_Power", "Per_Superior", "Var_Superior")
  if(!is.null(ar)){
    colnames <- c(colnames, ar)
    colnames <- c(colnames, paste("BIAS(", ar, ")",sep = "", collapse = NULL))
  }
  result <- data.frame(matrix(ncol = length(colnames), nrow = max(nrow(par1),nrow(par2))))
  colnames(result) <- colnames
  
  cl <- parallel::makeCluster(1)  #increase according to the limitations of your cluster if you run it on the cluster
  doParallel::registerDoParallel(cl)
  
  #2-Armed Or Multi-armed
  ######  2  Armed  ##########################################################
  if(K==1){
    ########## Theoretical Variance Calculation ############################
    if(dist=="bern"){
      mean <- par1
      fun <- function(x){x*(1-x)}
      variance <-  data.frame(lapply(par1,fun))
      theta <- par1[,2]*(1-par1[,1]) + 0.5*(par1[,1]*par1[,2]+(1-par1[,1])*(1-par1[,2]))
      psi <- 0.5*theta +0.25
    }
    if(dist=="norm"){
      mean <- par1
      fun <- function(x){x^2}
      variance <-  data.frame(lapply(par2,fun))
      xn <- (mean[,1] - mean[,2])/sqrt(variance[,1]+variance[,2])
      theta <-  1 - dnorm(xn)
      psi <- 0.5*theta +0.25
    }
    if(dist=="expon"){
      fun <- function(x){1/x}
      mean <- data.frame(lapply(par1,fun))
      fun <- function(x){1/(x^2)}
      variance <-  data.frame(lapply(par1,fun))
      theta <-  par1[,1]/(par1[,1]+par1[,2])
      psi <- 0.5*theta +0.25
    }
    if(dist=="beta"){
      fun <- function(x){5/(x+5)}
      mean <- data.frame(lapply(par1,fun))
      fun <- function(x){5*x/((5+x)^2+(5+x+1))}
      variance <-  data.frame(lapply(par1,fun))
      theta <- c(0.5,0.5631895, 0.6317465, 0.704089, 0.7777655, 0.813794,0.8486965,
                 0.881401, 0.91132, 0.937633, 0.959517, 0.976518, 0.988456, 0.991853,
                  0.994524, 0.995621, 0.997357, 0.998557, 0.998991, 0.999326, 0.999584)
      psi <- 0.5*theta +0.25
      psi <- c(rep(0.5, length(theta)),theta[-1])
    }
    if(dist=="Lickert"){
      ## Variance needed to be estimated through simulations in extra Code
      mean <-  read.csv("Lickert_Mean.csv", header = TRUE)
      mean <- mean[,2:3]
      variance <-  read.csv("Lickert_Variance.csv", header = TRUE)
      variance <- variance[,2:3]
      theta <- c(0.5, 0.544, 0.593, 0.6446757, 0.696642, 0.7469708, 0.796875, 0.818026,
                  0.840608, 0.852698, 0.8653635, 0.878488, 0.892124, 0.906136, 0.92029,
                  0.93435, 0.9479362, 0.9606685, 0.9720625, 0.9817345, 0.9893505,
                  0.9947455, 0.998026, 0.999626)
      psi <- 0.5*theta +0.25
      psi <- c(rep(0.5, length(psi)),psi[-1])
    }
    ############## ER ######################################################
    if(method == "ER" ){
      result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
      result.ER <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD")) %dorng% { #%dorng% 
        
        result_part <-result_part_emp
        
        source('ER.R')
        p <- as.numeric(par1[i,])
        start <- length(p)
        result_part[1,1:start] <- p
        if(is.null(par2)){
          sim.ER <- sim_ER(N=N, dist=dist, K=1, par1 = p,  burnin= burnin, nsim=nsim, measure=measure)
        } else{
          sig <- as.numeric(par2[i,])
          start <- start+ length(sig)
          result_part[1,1:start] <- c(p, sig)
          sim.ER <- sim_ER(N=N, dist=dist, K=1, par1 = p, burnin= burnin, par2=sig, nsim=nsim)
        }
        
        
        result_part[1,(start+1)] <- mean(sim.ER[,5]); result_part[1,(start+2)] <- var(sim.ER[,5]) 
        result_part[1,(start+3)] <- mean(sim.ER[,6]); result_part[1,(start+4)] <- var(sim.ER[,6]) 
        
        result_part[1,(start+5)] <- mean(sim.ER[,1])
        result_part[1,(start+6)] <- var(sim.ER[,1])
        
        reject_Z <- sum(sim.ER[,2] < 0.05)
        reject_BM <- sum(sim.ER[,3] < 0.05)
        
        if(p[1]==p[2]){
          result_part[1,(start+7)] <- reject_Z/nsim
          result_part[1,(start+8)] <- reject_BM/nsim
        } else {
          result_part[1,(start+7)] <- NA
          result_part[1,(start+8)] <- NA
        }
        result_part[1,(start+9)] <- reject_Z/nsim
        result_part[1,(start+10)] <- reject_BM/nsim
        
        result_part[1,(start+11)] <- mean(sim.ER[,4])
        result_part[1,(start+12)] <- var(sim.ER[,4])
        
        result_part
        
      }
      result[,1:(length(colnames))] <- result.ER
    }
    ############## RPW ######################################################
    if(dist=="bern" && method == "RPW"){
      result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
      result.RPW <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD")) %dorng% { #%dorng% 
        source('RPW.R')
        result_part <-result_part_emp
        p <- as.numeric(par1[i,])
        start <- length(p)
        result_part[1,1:start] <- p
        sim.RPW <- sim_RPW(N=N, dist=dist, K=1, par1 = p,  nsim=nsim, burnin=burnin)
        result_part[1,(start+1)] <- mean(sim.RPW[,5]); result_part[1,(start+2)] <- var(sim.RPW[,5]) 
        result_part[1,(start+3)] <- mean(sim.RPW[,6]); result_part[1,(start+4)] <- var(sim.RPW[,6]) 
        
        result_part[1,(start+5)] <- mean(sim.RPW[,1])
        result_part[1,(start+6)] <- var(sim.RPW[,1])
        
        reject_Z <-  sum(sim.RPW[,2] < 0.05)
        reject_BM <- sum(sim.RPW[,3] < 0.05)
        
        if(p[1]==p[2]){
          result_part[1,(start+7)] <- reject_Z/nsim
          result_part[1,(start+8)] <- reject_BM/nsim
        } else {
          result_part[1,(start+7)] <- NA
          result_part[1,(start+8)] <- NA
        }
        result_part[1,(start+9)] <- reject_Z/nsim
        result_part[1,(start+10)] <- reject_BM/nsim
        
        result_part[1,(start+11)] <- mean(sim.RPW[,4])
        result_part[1,(start+12)] <- var(sim.RPW[,4])
        
        result_part

      }
      result[,1:(length(colnames))] <- result.RPW 
    }
    ############## Sequential ######################################################
    if(method == "SEQ_R_minF"){  ### change
      result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
      result.SEQ_R <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD")) %dorng% { #%dorng% 
        source('SEQ.R')
        result_part <-result_part_emp
        p <- as.numeric(par1[i,])
        start <- length(p)
        result_part[1,1:start] <- p
        sim.SEQ_R <- sim_SEQ(N=N, dist=dist, K=1, par1 = p,  nsim=nsim, measure=measure, ar=ar, burnin=burnin)
        result_part[1,(start+1)] <- mean(sim.SEQ_R[,5]); result_part[1,(start+2)] <- var(sim.SEQ_R[,5]) 
        result_part[1,(start+3)] <- mean(sim.SEQ_R[,6]); result_part[1,(start+4)] <- var(sim.SEQ_R[,6]) 
        
        result_part[1,(start+5)] <- mean(sim.SEQ_R[,1])
        result_part[1,(start+6)] <- var(sim.SEQ_R[,1])
        
        reject_Z <-  sum(sim.SEQ_R[,2] < 0.05)
        reject_BM <- sum(sim.SEQ_R[,3] < 0.05)
        
        if(p[1]==p[2]){
          result_part[1,(start+7)] <- reject_Z/nsim
          result_part[1,(start+8)] <- reject_BM/nsim
        } else {
          result_part[1,(start+7)] <- NA
          result_part[1,(start+8)] <- NA
        }
        result_part[1,(start+9)] <- reject_Z/nsim
        result_part[1,(start+10)] <- reject_BM/nsim
        
        result_part[1,(start+11)] <- mean(sim.SEQ_R[,4])
        result_part[1,(start+12)] <- var(sim.SEQ_R[,4])
        
        result_part[1,(start+13)] <- 1- sqrt(p[1])/(sqrt(p[1]) + sqrt(p[2]))

        result_part[1,(start+14)] <- mean(sim.SEQ_R[,7]) - result_part[1,(start+13)]
        result_part
      }
      result[,1:(length(colnames))] <- result.SEQ_R 
    }
    if(method == "SEQ_R_maxMR"){  ### change
      result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
      result.SEQ_R <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD")) %dorng% { #%dorng% 
        source('SEQ.R')
        result_part <-result_part_emp
        p <- as.numeric(par1[i,])
        mu <- as.numeric(mean[i,])
        var <- as.numeric(variance[i,])
        if(dist=="norm") {
          std <- as.numeric(par2[i,])
          start <- length(c(p,std))
          sim.SEQ_R <- sim_SEQ(N=N, dist=dist, K=1, par1 = p, par2 = std,  nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- c(p,std)
        } else {
          start <- length(p)
          sim.SEQ_R <- sim_SEQ(N=N, dist=dist, K=1, par1 = p,   nsim=nsim, ar=ar)
          result_part[1,1:start] <- p
        }

        result_part[1,(start+1)] <- mean(sim.SEQ_R[,5]); result_part[1,(start+2)] <- var(sim.SEQ_R[,5]) 
        result_part[1,(start+3)] <- mean(sim.SEQ_R[,6]); result_part[1,(start+4)] <- var(sim.SEQ_R[,6]) 
        
        result_part[1,(start+5)] <- mean(sim.SEQ_R[,1])
        result_part[1,(start+6)] <- var(sim.SEQ_R[,1])
        
        reject_Z <-  sum(sim.SEQ_R[,2] < 0.05)
        reject_BM <- sum(sim.SEQ_R[,3] < 0.05)
        
        
        if(p[1]==p[2]){
          result_part[1,(start+7)] <- reject_Z/nsim
          result_part[1,(start+8)] <- reject_BM/nsim
        } else {
          result_part[1,(start+7)] <- NA
          result_part[1,(start+8)] <- NA
        }
        result_part[1,(start+9)] <- reject_Z/nsim
        result_part[1,(start+10)] <- reject_BM/nsim
        
        result_part[1,(start+11)] <- mean(sim.SEQ_R[,4])
        result_part[1,(start+12)] <- var(sim.SEQ_R[,4])
        
        if(dist=="bern"){
          result_part[1,(start+13)] <- sqrt(1-p[1])/(sqrt(1-p[1]) + sqrt(1-p[2]))
        } else {
          if(mu[1]<=0 || mu[2]<=0){ 
            result_part[1,(start+13)] <- 0.5
          } else {
            result_part[1,(start+13)] <- sqrt(mu[2])*var[1]/(sqrt(mu[2])*var[1]+sqrt(mu[1])*var[2])
          }
        }
        result_part[1,(start+14)] <- mean(sim.SEQ_R[,7]) - result_part[1,(start+13)]
        result_part
        
      }
      result[,1:(length(colnames))] <- result.SEQ_R 
    }
    if(method == "SEQ_Neyman"){
      result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
      result.SEQ_Neyman <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD")) %dorng% { #%dorng% 
        source('SEQ.R')
        result_part <-result_part_emp
        p <- as.numeric(par1[i,])
        var <- as.numeric(variance[i,])
        if(dist=="norm") {
          std <- as.numeric(par2[i,])
          start <- length(c(p,std))
          sim.SEQ_Neyman <- sim_SEQ(N=N, dist=dist, K=1, par1 = p, par2 = std,  nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- c(p,std)
        } else {
          start <- length(p)
          sim.SEQ_Neyman <- sim_SEQ(N=N, dist=dist, K=1, par1 = p,  measure=measure, nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- p
        }
        result_part[1,(start+1)] <- mean(sim.SEQ_Neyman[,5]); result_part[1,(start+2)] <- var(sim.SEQ_Neyman[,5]) 
        result_part[1,(start+3)] <- mean(sim.SEQ_Neyman[,6]); result_part[1,(start+4)] <- var(sim.SEQ_Neyman[,6]) 
        
        result_part[1,(start+5)] <- mean(sim.SEQ_Neyman[,1])
        result_part[1,(start+6)] <- var(sim.SEQ_Neyman[,1])
        
        reject_Z <-  sum(sim.SEQ_Neyman[,2] < 0.05)
        reject_BM <- sum(sim.SEQ_Neyman[,3] < 0.05)
        
        if(p[1]==p[2]){
          result_part[1,(start+7)] <- reject_Z/nsim
          result_part[1,(start+8)] <- reject_BM/nsim
        } else {
          result_part[1,(start+7)] <- NA
          result_part[1,(start+8)] <- NA
        }
        result_part[1,(start+9)] <- reject_Z/nsim
        result_part[1,(start+10)] <- reject_BM/nsim
        
        result_part[1,(start+11)] <- mean(sim.SEQ_Neyman[,4])
        result_part[1,(start+12)] <- var(sim.SEQ_Neyman[,4])
        
        result_part[1,(start+13)] <- 1 - var[1]/(var[1]+var[2])
        result_part[1,(start+14)] <- mean(sim.SEQ_Neyman[,7]) - result_part[1,(start+13)]
        
        result_part
        
      }
      result[,1:(length(colnames))] <- result.SEQ_Neyman 
    }
    if(method == "SEQ_NPPower"){
      result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
      result.SEQ_NPPower <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD")) %dorng% { #%dorng% 
        source('SEQ.R')
        result_part <-result_part_emp
        p <- as.numeric(par1[i,])
        var <- as.numeric(variance[i,])
        if(dist=="norm") {
          std <- as.numeric(par2[i,])
          start <- length(c(p,std))
          sim.SEQ_NPPower <- sim_SEQ(N=N, dist=dist, K=1, par1 = p, par2 = std,  nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- c(p,std)
        } else {
          start <- length(p)
          sim.SEQ_NPPower <- sim_SEQ(N=N, dist=dist, K=1, par1 = p,   nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- p
        }
        result_part[1,(start+1)] <- mean(sim.SEQ_NPPower[,5]); result_part[1,(start+2)] <- var(sim.SEQ_NPPower[,5]) 
        result_part[1,(start+3)] <- mean(sim.SEQ_NPPower[,6]); result_part[1,(start+4)] <- var(sim.SEQ_NPPower[,6]) 
        
        result_part[1,(start+5)] <- mean(sim.SEQ_NPPower[,1])
        result_part[1,(start+6)] <- var(sim.SEQ_NPPower[,1])
        
        reject_Z <-  sum(sim.SEQ_NPPower[,2] < 0.05)
        reject_BM <- sum(sim.SEQ_NPPower[,3] < 0.05)
        
        if(p[1]==p[2]){
          result_part[1,(start+7)] <- reject_Z/nsim
          result_part[1,(start+8)] <- reject_BM/nsim
        } else {
          result_part[1,(start+7)] <- NA
          result_part[1,(start+8)] <- NA
        }
        result_part[1,(start+9)] <- reject_Z/nsim
        result_part[1,(start+10)] <- reject_BM/nsim
        
        result_part[1,(start+11)] <- mean(sim.SEQ_NPPower[,4])
        result_part[1,(start+12)] <- var(sim.SEQ_NPPower[,4])
        
        result_part[1,(start+13)] <- 1 - var[1]/(var[1]+var[2])
        result_part[1,(start+14)] <- mean(sim.SEQ_NPPower[,7]) - result_part[1,(start+13)]
        
        result_part
        
      }
      result[,1:(length(colnames))] <- result.SEQ_NPPower 
    }
    if(method == "SEQ_WMW"){
      result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
      result.SEQ_WMW <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD", "pseudorank")) %dorng% { #%dorng% 
        source('SEQ.R')
        result_part <-result_part_emp
        p <- as.numeric(par1[i,])
        var <- as.numeric(variance[i,])
        if(dist=="norm") {
          std <- as.numeric(par2[i,])
          start <- length(c(p,std))
          sim.SEQ_WMW <- sim_SEQ(N=N, dist=dist, K=1, par1 = p, par2 = std,  nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- c(p,std)
        } else {
          start <- length(p)
          sim.SEQ_WMW <- sim_SEQ(N=N, dist=dist, K=1, par1 = p,   nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- p
        }
        result_part[1,(start+1)] <- mean(sim.SEQ_WMW[,5]); result_part[1,(start+2)] <- var(sim.SEQ_WMW[,5]) 
        result_part[1,(start+3)] <- mean(sim.SEQ_WMW[,6]); result_part[1,(start+4)] <- var(sim.SEQ_WMW[,6]) 
        
        result_part[1,(start+5)] <- mean(sim.SEQ_WMW[,1])
        result_part[1,(start+6)] <- var(sim.SEQ_WMW[,1])
        
        reject_Z <-  sum(sim.SEQ_WMW[,2] < 0.05)
        reject_BM <- sum(sim.SEQ_WMW[,3] < 0.05)
        
        if(p[1]==p[2]){
          result_part[1,(start+7)] <- reject_Z/nsim
          result_part[1,(start+8)] <- reject_BM/nsim
        } else {
          result_part[1,(start+7)] <- NA
          result_part[1,(start+8)] <- NA
        }
        result_part[1,(start+9)] <- reject_Z/nsim
        result_part[1,(start+10)] <- reject_BM/nsim
        
        result_part[1,(start+11)] <- mean(sim.SEQ_WMW[,4])
        result_part[1,(start+12)] <- var(sim.SEQ_WMW[,4])
        
        result_part[1,(start+13)] <- as.numeric(psi[i])
        result_part[1,(start+14)] <- mean(sim.SEQ_WMW[,7]) - result_part[1,(start+13)]
        
        result_part
        
      }
      result[,1:(length(colnames))] <- result.SEQ_WMW 
    }
    if( method == "SEQ_WMW_log"){
      result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
      result.SEQ_WMW_log <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD", "pseudorank")) %dorng% { #%dorng% 
        source('SEQ.R')
        result_part <-result_part_emp
        p <- as.numeric(par1[i,])
        var <- as.numeric(variance[i,])
        if(dist=="norm") {
          std <- as.numeric(par2[i,])
          start <- length(c(p,std))
          sim.SEQ_WMW_log <- sim_SEQ(N=N, dist=dist, K=1, par1 = p, par2 = std,  nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- c(p,std)
        } else {
          start <- length(p)
          sim.SEQ_WMW_log <- sim_SEQ(N=N, dist=dist, K=1, par1 = p,   nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- p
        }
        result_part[1,(start+1)] <- mean(sim.SEQ_WMW_log[,5]); result_part[1,(start+2)] <- var(sim.SEQ_WMW_log[,5]) 
        result_part[1,(start+3)] <- mean(sim.SEQ_WMW_log[,6]); result_part[1,(start+4)] <- var(sim.SEQ_WMW_log[,6]) 
        
        result_part[1,(start+5)] <- mean(sim.SEQ_WMW_log[,1])
        result_part[1,(start+6)] <- var(sim.SEQ_WMW_log[,1])
        
        reject_Z <-  sum(sim.SEQ_WMW_log[,2] < 0.05)
        reject_BM <- sum(sim.SEQ_WMW_log[,3] < 0.05)
        
        if(p[1]==p[2]){
          result_part[1,(start+7)] <- reject_Z/nsim
          result_part[1,(start+8)] <- reject_BM/nsim
        } else {
          result_part[1,(start+7)] <- NA
          result_part[1,(start+8)] <- NA
        }
        result_part[1,(start+9)] <- reject_Z/nsim
        result_part[1,(start+10)] <- reject_BM/nsim
        
        result_part[1,(start+11)] <- mean(sim.SEQ_WMW_log[,4])
        result_part[1,(start+12)] <- var(sim.SEQ_WMW_log[,4])
        
        result_part[1,(start+13)] <- as.numeric(psi[i])^log(N-1)/(as.numeric(psi[i])^log(N-1)+(1-as.numeric(psi[i]))^log(N-1))
        result_part[1,(start+14)] <- mean(sim.SEQ_WMW_log[,7]) - result_part[1,(start+13)]
        
        result_part
        
      }
      result[,1:(length(colnames))] <- result.SEQ_WMW_log 
    }
    if( method == "SEQ_WMW_lin"){
      result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
      result.SEQ_WMW_lin <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD", "pseudorank")) %dorng% { #%dorng% 
        source('SEQ.R')
        result_part <-result_part_emp
        p <- as.numeric(par1[i,])
        var <- as.numeric(variance[i,])
        if(dist=="norm") {
          std <- as.numeric(par2[i,])
          start <- length(c(p,std))
          sim.SEQ_WMW_lin <- sim_SEQ(N=N, dist=dist, K=1, par1 = p, par2 = std,  nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- c(p,std)
        } else {
          start <- length(p)
          sim.SEQ_WMW_lin <- sim_SEQ(N=N, dist=dist, K=1, par1 = p,   nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- p
        }
        result_part[1,(start+1)] <- mean(sim.SEQ_WMW_lin[,5]); result_part[1,(start+2)] <- var(sim.SEQ_WMW_lin[,5]) 
        result_part[1,(start+3)] <- mean(sim.SEQ_WMW_lin[,6]); result_part[1,(start+4)] <- var(sim.SEQ_WMW_lin[,6]) 
        
        result_part[1,(start+5)] <- mean(sim.SEQ_WMW_lin[,1])
        result_part[1,(start+6)] <- var(sim.SEQ_WMW_lin[,1])
        
        reject_Z <-  sum(sim.SEQ_WMW_lin[,2] < 0.05)
        reject_BM <- sum(sim.SEQ_WMW_lin[,3] < 0.05)
        
        if(p[1]==p[2]){
          result_part[1,(start+7)] <- reject_Z/nsim
          result_part[1,(start+8)] <- reject_BM/nsim
        } else {
          result_part[1,(start+7)] <- NA
          result_part[1,(start+8)] <- NA
        }
        result_part[1,(start+9)] <- reject_Z/nsim
        result_part[1,(start+10)] <- reject_BM/nsim
        
        result_part[1,(start+11)] <- mean(sim.SEQ_WMW_lin[,4])
        result_part[1,(start+12)] <- var(sim.SEQ_WMW_lin[,4])
        
        Nmin <- min(100,N) #for computational reasons
        result_part[1,(start+13)] <- as.numeric(psi[i])^(Nmin-1)/(as.numeric(psi[i])^(Nmin-1)+(1-as.numeric(psi[i]))^(Nmin-1))
        result_part[1,(start+14)] <- mean(sim.SEQ_WMW_lin[,7]) - result_part[1,(start+13)]
        
        result_part
        
      }
      result[,1:(length(colnames))] <- result.SEQ_WMW_lin 
    }
    if( method == "ERADE"){
      result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
      result.ERADE <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD", "pseudorank")) %dorng% { #%dorng% 
        source('ERADE.R')
        result_part <-result_part_emp
        p <- as.numeric(par1[i,])
        mu <- as.numeric(mean[i,])
        var <- as.numeric(variance[i,])
        if(dist=="norm") {
          std <- as.numeric(par2[i,])
          start <- length(c(p,std))
          sim.ERADE <- sim_ERADE(N=N, dist=dist, K=1, par1 = p, par2 = std,  nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- c(p,std)
        } else {
          start <- length(p)
          sim.ERADE <- sim_ERADE(N=N, dist=dist, K=1, par1 = p,   measure= measure, nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- p
        }
        result_part[1,(start+1)] <- mean(sim.ERADE[,5]); result_part[1,(start+2)] <- var(sim.ERADE[,5]) 
        result_part[1,(start+3)] <- mean(sim.ERADE[,6]); result_part[1,(start+4)] <- var(sim.ERADE[,6]) 
        
        result_part[1,(start+5)] <- mean(sim.ERADE[,1])
        result_part[1,(start+6)] <- var(sim.ERADE[,1])
        
        reject_Z <-  sum(sim.ERADE[,2] < 0.05)
        reject_BM <- sum(sim.ERADE[,3] < 0.05)
        
        if(p[1]==p[2]){
          result_part[1,(start+7)] <- reject_Z/nsim
          result_part[1,(start+8)] <- reject_BM/nsim
        } else {
          result_part[1,(start+7)] <- NA
          result_part[1,(start+8)] <- NA
        }
        result_part[1,(start+9)] <- reject_Z/nsim
        result_part[1,(start+10)] <- reject_BM/nsim
        
        result_part[1,(start+11)] <- mean(sim.ERADE[,4])
        result_part[1,(start+12)] <- var(sim.ERADE[,4])
        
        if(ar=="R_minF"){
          result_part[1,(start+13)] <- 1- sqrt(p[1])/(sqrt(p[1]) + sqrt(p[2]))
        } 
        if(ar=="R_maxMR"){
          if(dist=="bern"){
            result_part[1,(start+13)] <- sqrt(1-p[1])/(sqrt(1-p[1]) + sqrt(1-p[2]))
          } else {
            if(mu[1]<=0 || mu[2]<=0){ 
              result_part[1,(start+13)] <- 0.5
            } else {
              result_part[1,(start+13)] <- sqrt(mu[2])*var[1]/(sqrt(mu[2])*var[1]+sqrt(mu[1])*var[2])
            }
          }
        }
        if(ar=="Neyman" || ar=="NPPower" ){
          result_part[1,(start+13)] <- 1 - var[1]/(var[1]+var[2])
        }
        if(ar=="WMW"){
          result_part[1,(start+13)] <- as.numeric(psi[i])
        }
        if(ar=="WMW_log"){
          result_part[1,(start+13)] <- as.numeric(psi[i])^log(N-1)/(as.numeric(psi[i])^log(N-1)+(1-as.numeric(psi[i]))^log(N-1))
        } 
        if(ar=="WMW_lin"){
          result_part[1,(start+13)] <- as.numeric(psi[i])^(N-1)/(as.numeric(psi[i])^(N-1)+(1-as.numeric(psi[i]))^(N-1))
        }
        result_part[1,(start+14)] <- mean(sim.ERADE[,7]) - result_part[1,(start+13)]
        
        result_part
        
      }
      result[,1:(length(colnames))] <- result.ERADE 
    }
    if( method == "DBCD"){
      result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
      result.DBCD <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD", "pseudorank")) %dorng% { #%dorng% 
        source('DBCD.R')
        result_part <-result_part_emp
        p <- as.numeric(par1[i,])
        mu <- as.numeric(mean[i,])
        var <- as.numeric(variance[i,])
        if(dist=="norm") {
          std <- as.numeric(par2[i,])
          start <- length(c(p,std))
          sim.DBCD <- sim_DBCD(N=N, dist=dist, K=1, par1 = p, par2 = std,  nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- c(p,std)
        } else {
          start <- length(p)
          sim.DBCD <- sim_DBCD(N=N, dist=dist, K=1, par1 = p,  measure= measure, nsim=nsim, ar=ar, burnin=burnin)
          result_part[1,1:start] <- p
        }
        result_part[1,(start+1)] <- mean(sim.DBCD[,5]); result_part[1,(start+2)] <- var(sim.DBCD[,5]) 
        result_part[1,(start+3)] <- mean(sim.DBCD[,6]); result_part[1,(start+4)] <- var(sim.DBCD[,6]) 
        
        result_part[1,(start+5)] <- mean(sim.DBCD[,1])
        result_part[1,(start+6)] <- var(sim.DBCD[,1])
        
        reject_Z <-  sum(sim.DBCD[,2] < 0.05)
        reject_BM <- sum(sim.DBCD[,3] < 0.05)
        
        if(p[1]==p[2]){
          result_part[1,(start+7)] <- reject_Z/nsim
          result_part[1,(start+8)] <- reject_BM/nsim
        } else {
          result_part[1,(start+7)] <- NA
          result_part[1,(start+8)] <- NA
        }
        result_part[1,(start+9)] <- reject_Z/nsim
        result_part[1,(start+10)] <- reject_BM/nsim
        
        result_part[1,(start+11)] <- mean(sim.DBCD[,4])
        result_part[1,(start+12)] <- var(sim.DBCD[,4])
        
        if(ar=="R_minF"){
          result_part[1,(start+13)] <- 1- sqrt(p[1])/(sqrt(p[1]) + sqrt(p[2]))
        } 
        if(ar=="R_maxMR"){
          if(dist=="bern"){
            result_part[1,(start+13)] <- sqrt(1-p[1])/(sqrt(1-p[1]) + sqrt(1-p[2]))
          } else {
            if(mu[1]<=0 || mu[2]<=0){ 
              result_part[1,(start+13)] <- 0.5
            } else {
              result_part[1,(start+13)] <- sqrt(mu[2])*var[1]/(sqrt(mu[2])*var[1]+sqrt(mu[1])*var[2])
            }
          }
        }
        if(ar=="Neyman" || ar=="NPPower" ){
          result_part[1,(start+13)] <- 1 - var[1]/(var[1]+var[2])
        }
        if(ar=="WMW"){
          result_part[1,(start+13)] <- as.numeric(psi[i])
        }
        if(ar=="WMW_log"){
          result_part[1,(start+13)] <- as.numeric(psi[i])^log(N-1)/(as.numeric(psi[i])^log(N-1)+(1-as.numeric(psi[i]))^log(N-1))
        } 
        if(ar=="WMW_lin"){
          result_part[1,(start+13)] <- as.numeric(psi[i])^(N-1)/(as.numeric(psi[i])^(N-1)+(1-as.numeric(psi[i]))^(N-1))
        }
        result_part[1,(start+14)] <- mean(sim.DBCD[,7]) - result_part[1,(start+13)]
        
        result_part
        
      }
      result[,1:(length(colnames))] <- result.DBCD 
    }
    if( method == "TS"){
      result_part_emp <- data.frame(matrix(ncol = length(colnames), nrow = 0))
      result.TS <- foreach(i = 1:length(par1[,1]), .combine='rbind', .packages=c("BSDA", "rankFD", "pseudorank", "rBeta2009")) %dorng% { #%dorng% 
        source('TSsimp.R')
        result_part <-result_part_emp
        p <- as.numeric(par1[i,])
        mu <- as.numeric(mean[i,])
        var <- as.numeric(variance[i,])
        if(dist=="norm") {
          std <- as.numeric(par2[i,])
          start <- length(c(p,std))
          sim.TS <- sim_TS(N=N, dist=dist, K=1, par1 = p, par2 = std,  nsim=nsim)
          result_part[1,1:start] <- c(p,std)
        } else {
          start <- length(p)
          sim.TS <- sim_TS(N=N, dist=dist, K=1, par1 = p,   nsim=nsim)
          result_part[1,1:start] <- p
        }
        result_part[1,(start+1)] <- mean(sim.TS[,5]); result_part[1,(start+2)] <- var(sim.TS[,5]) 
        result_part[1,(start+3)] <- mean(sim.TS[,6]); result_part[1,(start+4)] <- var(sim.TS[,6]) 
        
        result_part[1,(start+5)] <- mean(sim.TS[,1])
        result_part[1,(start+6)] <- var(sim.TS[,1])
        
        reject_Z <-  sum(sim.TS[,2] < 0.05)
        reject_BM <- sum(sim.TS[,3] < 0.05)
        
        if(p[1]==p[2]){
          result_part[1,(start+7)] <- reject_Z/nsim
          result_part[1,(start+8)] <- reject_BM/nsim
        } else {
          result_part[1,(start+7)] <- NA
          result_part[1,(start+8)] <- NA
        }
        result_part[1,(start+9)] <- reject_Z/nsim
        result_part[1,(start+10)] <- reject_BM/nsim
        
        result_part[1,(start+11)] <- mean(sim.TS[,4])
        result_part[1,(start+12)] <- var(sim.TS[,4])
        
        result_part[1,(start+13)] <- NaN
        result_part[1,(start+14)] <- NaN
        
        result_part
        
      }
      result[,1:(length(colnames))] <- result.TS 
    }
    
  } else {
    
  }
  
  
  #Print and save Results
  if(method=="ERADE" || method=="DBCD"){
    filename <-  paste(method, ar, N, dist,"All.csv",sep="_")
  } else {
    filename <-  paste(method, N, dist,"All.csv",sep="_")
  }
  
  write.csv(round(result,4), filename, row.names=TRUE)
  print(round(result,4))
  print(round(N*(1-result[7])))
  stopCluster(cl)
}

# Generate Parameter Data Frame
# Bernoulli
p1S <- seq(0,1,0.01) #0.05
p2S <- seq(0,1,0.01) #0.05
par1_Ber <- data.frame(sort(rep(p1S,max(length(p1S),length(p2S)))),rep(p2S,max(length(p1S),length(p2S))))
#par1_Ber <- data.frame(0.05,0.30)

#Normal
mu1 <- seq(0,4,0.1) #0.1
mu2 <- seq(0,4,0.1) #0.1
sig1 <- 1
sig2 <- 3
par1_Nor <- data.frame(sort(rep(mu1,max(length(mu1),length(mu2)))),rep(mu2,max(length(mu1),length(mu2))))
par2_Nor <- data.frame(rep(sig1,length(par1_Nor[,1])),rep(sig2,length(par1_Nor[,1])))

#beta
beta1 <- 100
beta2 <- c(100,90,80,70,60,55,50,45,40,35,30,25,20,18,16,15,13,11,10,9,8)
par1_beta <- data.frame(c(beta2,rep(beta1,length(beta2[-1]))), c(beta2,beta2[-1]))

#expon
exp1 <- 1
exp2 <- c(c(1,1.25,1.5,1.75,2,2.5),c(3:19), seq(20,100,4), 
            c(150,200,300,500,750,1000))
exp2 <-  1/exp2
par1_expon <- data.frame(c(exp2,rep(exp1,length(exp2[-1]))), c(exp2,exp2[-1]))

# Simulations
#Test
#simulation(N=100, dist="bern", K=1, par1 = par1_Ber , nsim=100,  burnin=3, method = "TS") #par1_Ber

#Scenarios
#N <- c(20, 50,60,80, 100, 150,200)
N <- c(50)
Distributions <- c("bern", "norm", "expon", "Lickert")
param_list <- expand.grid(N = N, Distributions = Distributions)
#N <- c(60)

#Distributions <- c("bern")
#param_list <- expand.grid(N = N, Distributions = Distributions)


for(i in 1:dim(param_list)[1]){
  nsim =10^4 #later 10^4!
  burn <- round(param_list[i,]$N/12)
  if(param_list[i,]$Distributions=="bern") {
   simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Ber, burnin=param_list[i,]$N/2, nsim=nsim, method = "ER", measure = "sd") #10^4
   # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Ber, nsim=nsim, method = "RPW")
 #   simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, burnin=burn, par1 = par1_Ber, nsim=nsim, method = "SEQ_R_minF", ar="R_minF", measure = "lrr")
   #simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, burnin=burn, par1 = par1_Ber, nsim=nsim, method = "SEQ_R_maxMR", ar="R_maxMR")
  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, burnin=burn, par1 = par1_Ber, nsim=nsim, method = "SEQ_Neyman", ar="Neyman", measure = "sd") #10^4
     simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Ber, nsim=nsim, method = "SEQ_NPPower", ar="NPPower") 
    # #simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Ber, nsim=nsim, method = "SEQ_WMW", ar="WMW") #10^4
    # #simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Ber, nsim=nsim, method = "SEQ_WMW_log", ar="WMW_log") #10^4
    # #simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Ber, nsim=nsim, method = "SEQ_WMW_lin", ar="WMW_lin") #10^4
 #   simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, burnin=burn, par1 = par1_Ber, nsim=nsim, method = "DBCD", ar="R_minF", measure = "lrr") #10^4
 #    simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, burnin=burn, par1 = par1_Ber, nsim=nsim, method = "DBCD", ar="R_maxMR") #10^4
   simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, burnin=burn, par1 = par1_Ber, nsim=nsim, method = "DBCD", ar="Neyman", measure = "sd") #10^4
     simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Ber, nsim=nsim, method = "DBCD", ar="NPPower")
    # # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Ber, nsim=nsim, method = "DBCD", ar="WMW") #10^4
    # # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Ber, nsim=nsim, method = "DBCD", ar="WMW_log") #10^4
    # # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Ber, nsim=nsim, method = "DBCD", ar="WMW_lin") #10^4
    # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, burnin=burn, par1 = par1_Ber, nsim=nsim, method = "ERADE", ar="R_minF", measure = "lrr") #10^4
  #   simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, burnin=burn, par1 = par1_Ber, nsim=nsim, method = "ERADE", ar="R_maxMR") #10^4
  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, burnin=burn, par1 = par1_Ber, nsim=nsim, method = "ERADE", ar="Neyman", measure = "sd") #10^4
     simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Ber, nsim=nsim, method = "ERADE", ar="NPPower")
    # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Ber, nsim=nsim, method = "ERADE", ar="WMW") #10^4
    # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Ber, nsim=nsim, method = "ERADE", ar="WMW_log") #10^4
    # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Ber, nsim=nsim, method = "ERADE", ar="WMW_lin") #10^4
   # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Ber, nsim=nsim, method = "TS") #10^4
  }
  if(param_list[i,]$Distributions=="norm") {
    simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "ER") #10^4
  #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "SEQ_R_maxMR", ar="R_maxMR")
    simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "SEQ_Neyman", ar="Neyman") #10^4
     simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "SEQ_NPPower", ar="NPPower")
 #   simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "SEQ_WMW", ar="WMW") #10^4
 #  #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "SEQ_WMW_log", ar="WMW_log") #10^4
 #  #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "SEQ_WMW_lin", ar="WMW_lin") #10^4
 #    # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "DBCD", ar="R_maxMR")
  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "DBCD", ar="Neyman")
     simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "DBCD", ar="NPPower")
 #    # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "DBCD", ar="WMW")
 #    # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "DBCD", ar="WMW_log")
 #    # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "DBCD", ar="WMW_lin")
 #    # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "ERADE", ar="R_maxMR")
 simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "ERADE", ar="Neyman")
     simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "ERADE", ar="NPPower")
 #    # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "ERADE", ar="WMW")
 #    # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "ERADE", ar="WMW_log")
 #    # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "ERADE", ar="WMW_lin")
  #   simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_Nor, par2 = par2_Nor, nsim=nsim, method = "TS")
  }
  if(param_list[i,]$Distributions=="expon") {
    simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "ER") #10^4
  # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "SEQ_R_maxMR", ar="R_maxMR")
   simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "SEQ_Neyman", ar="Neyman") #10^4
     simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "SEQ_NPPower", ar="NPPower") 
  # # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "SEQ_WMW", ar="WMW") #10^4
  # #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "SEQ_WMW_log", ar="WMW_log") #10^4
  # #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "SEQ_WMW_lin", ar="WMW_lin") #10^4
  #   # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "DBCD", ar="R_maxMR")
   simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "DBCD", ar="Neyman") #10^4
     simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "DBCD", ar="NPPower")
  #   # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "DBCD", ar="WMW") #10^4
  #   #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "DBCD", ar="WMW_log") #10^4
  #   #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "DBCD", ar="WMW_lin") #10^4
  #   # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "ERADE", ar="R_maxMR")
   simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "ERADE", ar="Neyman") #10^4
     simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "ERADE", ar="NPPower") 
  #   # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "ERADE", ar="WMW") #10^4
    #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "ERADE", ar="WMW_log") #10^4
    #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "ERADE", ar="WMW_lin") #10^4
    #simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_expon, nsim=nsim, method = "TS") #10^4
  }
  if(param_list[i,]$Distributions=="beta") {
  #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "ER") #10^4
   # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "SEQ_R_maxMR", ar="R_maxMR")
  #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "SEQ_Neyman", ar="Neyman") #10^4
  #   simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "SEQ_NPPower", ar="NPPower")
  # #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "SEQ_WMW", ar="WMW") #10^4
  # #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "SEQ_WMW_log", ar="WMW_log") #10^4
  # #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "SEQ_WMW_lin", ar="WMW_lin") #10^4
  #   # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "DBCD", ar="R_maxMR")
  #   #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "DBCD", ar="Neyman") #10^4
  #   simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "DBCD", ar="NPPower")
  #   #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "DBCD", ar="WMW") #10^4
  #   #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "DBCD", ar="WMW_log") #10^4
  #   #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "DBCD", ar="WMW_lin") #10^4
  #   # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "ERADE", ar="R_maxMR")
  #   #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "ERADE", ar="Neyman") #10^4
  #   simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "ERADE", ar="NPPower") 
    #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "ERADE", ar="WMW") #10^4
    #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "ERADE", ar="WMW_log") #10^4
    #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "ERADE", ar="WMW_lin") #10^4
   # simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1, par1 = par1_beta, nsim=nsim, method = "TS") #10^4
  }
  if(param_list[i,]$Distributions=="Lickert") {
     simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "ER") #10^4
  #   simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "SEQ_R_maxMR", ar="R_maxMR")
     simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "SEQ_Neyman", ar="Neyman") #10^4
      simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "SEQ_NPPower", ar="NPPower")
  #   simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "SEQ_WMW", ar="WMW") #10^4
  #   simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "SEQ_WMW_log", ar="WMW_log") #10^4
  #   simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "SEQ_WMW_lin", ar="WMW_lin") #10^4#
  #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "DBCD", ar="R_maxMR")
    simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "DBCD", ar="Neyman") #10^4
     simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "DBCD", ar="NPPower")
  #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "DBCD", ar="WMW") #10^4
  #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "DBCD", ar="WMW_log") #10^4
  #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "DBCD", ar="WMW_lin") #10^4
  #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "ERADE", ar="R_maxMR")
    simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "ERADE", ar="Neyman") #10^4
    simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "ERADE", ar="NPPower")
  #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "ERADE", ar="WMW") #10^4
  #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "ERADE", ar="WMW_log") #10^4
  #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "ERADE", ar="WMW_lin") #10^4
  #  simulation(N=param_list[i,]$N, dist=param_list[i,]$Distributions, K=1,par1 = par1_beta, nsim=nsim, method = "TS") #10^4
  }
}