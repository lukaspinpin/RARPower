###############################################################################
library(BSDA)
library(rankFD)
library(pseudorank)

source('Data_Generator.R')

sim_DBCD <-  function(N=200, dist="bern", K=1, par1 = c(0,0), par2 = NULL, burnin=3, measure="sd", ar="WMW",nsim=10^4, alpha=0.05, one.sided = FALSE){
  
  if(K==1){
    out = matrix(nrow = nsim, ncol = 7)
    for(i in 1:nsim){
      if(i==1){
        out[i,] = two_arm_DBCD(N=N, dist=dist,  par1 =par1, par2 = par2, measure=measure, burnin = burnin,ar=ar, first=FALSE)
      } else {
        out[i,] = two_arm_DBCD(N=N, dist=dist,  par1 =par1, par2 = par2, measure=measure, burnin = burnin,ar=ar)
      }
      
    }
  } else {
    out = matrix(nrow = nsim, ncol = 6+K)
  }
  return(out)
}



###############################################################################

two_arm_DBCD = function(N=200, dist="bern",  par1 = c(0,0), par2 = NULL, measure="sd", burnin=3, ar="WMW", gamma=2, Z_corrected=FALSE, first=FALSE){
  
  A = rep(NA, N)   # Vector of allocations
  X = rep(NA, N)   # Vector of responses
  est <- rep(NA, N) 
  
  #Burnin
  if(burnin>0){
    A[1:burnin] = 1
    X[1:burnin] = generate_data(N=burnin, dist=dist, parV=c(par1[1],par2[1]))
    A[(burnin+1):(2*burnin)] = 2
    X[(burnin+1):(2*burnin)] = generate_data(N=burnin, dist=dist, parV=c(par1[2],par2[2]))
    est[1:(2*burnin)] <- 0.5
  }
  
  for (i in (2*burnin+1):N){
    n1 = sum(A==1, na.rm = TRUE)
    n2 = sum(A==2, na.rm = TRUE)
    
    if(ar=="R_minF"){
      if(dist=="bern") {
        s1 = sum(X[A==1], na.rm = TRUE)
        p1.hat = (s1)/(n1)
        s2 = sum(X[A==2], na.rm = TRUE)
        p2.hat = (s2)/(n2)
        if(measure=="sd"){
          if((sqrt(p1.hat) + sqrt(p2.hat)) == 0){
            rho1.hat = 0.5
          } else {
            rho1.hat = sqrt(p1.hat)/(sqrt(p1.hat) + sqrt(p2.hat))
          }
        } else if(measure=="lrr"){
          if((sqrt(p1.hat)*(1-p2.hat) + sqrt(p2.hat)*(1-p1.hat)) == 0){
            rho1.hat = 0.5
          } else {
            rho1.hat = 1- sqrt(p2.hat)*(1-p1.hat) / (sqrt(p1.hat)*(1-p2.hat) + sqrt(p2.hat)*(1-p1.hat))
          }
        }
        p = 1-rho1.hat
      }
    }
    if(ar=="R_maxMR"){
      if(dist=="bern") {
        s1 = sum(X[A==1], na.rm = TRUE)
        p1.hat = (s1)/(n1)
        s2 = sum(X[A==2], na.rm = TRUE)
        p2.hat = (s2)/(n2)
        if((sqrt(1-p1.hat) + sqrt(1-p2.hat)) == 0){
          p = 0.5
        } else {
          p = sqrt(1-p1.hat)/(sqrt(1-p1.hat) + sqrt(1-p2.hat))
        }
      } else {
        mu1 <- mean(X[A==1], na.rm = TRUE)
        mu2 <- mean(X[A==2], na.rm = TRUE)
        sig1_N <- sd(X[A==1], na.rm = TRUE)
        sig2_N <- sd(X[A==2], na.rm = TRUE)
        if(mu1<=0 || mu2<=0 || sig1_N+sig2_N==0){
          p <- 0.5
        } else {
          p <- sqrt(mu1)*sig2_N/(sqrt(mu2)*sig1_N+sqrt(mu1)*sig2_N)
        }
      }
    }
    if(ar=="Neyman"){
      if(measure=="sd"){
        sig1_N <- sd(X[A==1], na.rm = TRUE)
        sig2_N <- sd(X[A==2], na.rm = TRUE)
        if((sig1_N+sig2_N)==0){
          p <- 0.5
        } else {
          p <- 1 - sig1_N/(sig1_N+sig2_N)
        }
      } else if(measure=="lrr"){
        s1 = sum(X[A==1], na.rm = TRUE)
        p1.hat = (s1)/(n1)
        s2 = sum(X[A==2], na.rm = TRUE)
        p2.hat = (s2)/(n2)
        if((sqrt((1-p1.hat)*p2.hat) + sqrt((1-p2.hat)*p1.hat))==0){
          p <- 0.5
        } else {
          p <- sqrt((1-p1.hat)*p2.hat) / (sqrt((1-p1.hat)*p2.hat) + sqrt((1-p2.hat)*p1.hat))
        }
      }
      est[i] <- p
    }
    if(ar=="NPPower"){
      x1 <- X[A==1]
      x2 <- X[A==2]
      x1 <- x1[!is.na(x1)]
      x2 <- x2[!is.na(x2)]
      midrank <-rank(c(x1,x2), ties.method = "average")
      midrankI1 <-rank(x1, ties.method = "average")
      midrankI2 <-rank(x2, ties.method = "average")
      r1 <- mean(midrank[1:n1])
      r2 <- mean(midrank[-(1:n1)])
      S1 <- sum((midrank[1:n1]-midrankI1-r1+(n1+1)/2)^2)/(n1-1)
      S2 <- sum((midrank[-(1:n1)]-midrankI2-r2+(n2+1)/2)^2)/(n2-1)
      
      tau0 <- sqrt(S1)/(N-n1)
      tau1 <- sqrt(S2)/(N-n2)
      if((tau0+tau1)==0){
        p <- 0.5
      } else {
        p <- tau1/(tau0+tau1)
      }
      est[i] <- p
    }
    if(ar=="WMW"){ 
      N_curr <- n1 +n2
      pseudorank <- pseudorank(X,as.factor(A), ties.method = "average")
      levels <- sort(unique(A), na.last=NA)
      d <- length(levels)
      psi <- rep(0,d)
      for(j in 1:d){
        rpsi <- mean(pseudorank[A[1:N_curr]==j])
        psi[j] <- (rpsi-0.5)/N_curr
      }
      prop <- psi/(d/2)
      p <- prop[2]
    }
    if(ar=="WMW_log"){
      N_curr <- n1 +n2
      pseudorank <- pseudorank(X,as.factor(A), ties.method = "average")
      levels <- sort(unique(A), na.last=NA)
      d <- length(levels)
      psi <- rep(0,d)
      for(j in 1:d){
        rpsi <- mean(pseudorank[A[1:N_curr]==j])
        psi[j] <- (rpsi-0.5)/N_curr
      }
      prop <- psi/(d/2)
      prop <- prop^log(N_curr)/sum(prop^log(N_curr))
      p <- prop[2]
    }
    if(ar=="WMW_lin"){
      N_curr <- n1 +n2
      pseudorank <- pseudorank(X,as.factor(A), ties.method = "average")
      levels <- sort(unique(A), na.last=NA)
      d <- length(levels)
      psi <- rep(0,d)
      for(j in 1:d){
        rpsi <- mean(pseudorank[A[1:N_curr]==j])
        psi[j] <- (rpsi-0.5)/N_curr
      }
      prop <- psi/(d/2)
      prop <- prop^N_curr/sum(prop^N_curr)
      p <- prop[2]
    }
    if(p==0){p <- 1/N}
    if(p==1){p <- 1-1/N}
    est[i] <- p
    rho.hat <- c(1-p,p)
    
    alloc.prop = c(n1/i, n2/i)
    
    phi1 = DBCD(rho.hat, alloc.prop, gamma)
    
    
    A[i] = rbinom(1, 1, phi1[2]) + 1
    X[i]  <- generate_data(N=1, dist=dist, parV=c(par1[A[i]],par2[A[i]])) 
  }
  
  
  if(first){
    filename <-  paste(ar, N, dist, par1[1], par1[2], "Trial1.pdf",sep="_")
    pdf(filename, height = 4, width = 8)
    # Set plot parameters
    plot_colors <- "#0072B2"  # blue color for line
    title <- "Temporal Behavior of Estimator"  # plot title
    x_label <- "Patient"  # x-axis label
    y_label <- "Estimator"  # y-axis label
    
    # Create plot with custom formatting
    plot(1:N, est, type="l", col=plot_colors, lwd=2, ylim=c(0,1), xlab=x_label, ylab=y_label, main=title)
    grid()  # add gridlines to plot
    box()  # add border to plot
    dev.off() 
  }
  
  # Response
  x0 <- X[A==1]
  x1 <- X[A==2]
  n0 <- length(x0)
  n1 <- length(x1)
  s0 <- sum(x0)
  s1 <- sum(x1)
  response <- (s0+s1)/N
  
  #Imbalance Measure
  sup <- superior(dist=dist,par1 = par1, par2 = par2)
  if(sup==2){ #2 means that no arm is theoretically supirior 
    #   imb <- NA
    per.sup <- NA
  } else {
    if(sup==0){ 
      #    low <- n1 ; lar <- n0
      per.sup <- n0/N
    }
    if(sup==1){ 
      #   low <- n0 ; lar <- n1
      per.sup <- n1/N
    }
  }
  
  ##### Z - Test & BM-Test
  # One-Point-Sample-Correction
  sigx <- sd(x0)
  sigy <- sd(x1)
  ##### Z - Test & BM-Test
  source('WaldTest.R')
  ##### Z - Test 
  if(dist=="bern"){
    Z_P <-wald.test.binary(x0,x1, measure = measure)
  } else {
    Z_P <- wald.test(x0,x1)
  }
  if(sigx==0 && sigy==0){
    BM_P <- Z_P
  } else {
    BM <- lawstat::brunner.munzel.test(x0,x1)
    if(BM$statistic == Inf){
      BM_P <- 0
    } else if(BM$statistic == -Inf) {
      BM_P <-  1
    } else {
      BM_P <- BM$p.value
    }
  }
  #estimator <- p
  n <- c(n0,n1)/N
  return(c(response, Z_P, BM_P, per.sup, n, p))
}

###############################################################################

DBCD = function(rho.hat, alloc.prop, gamma=2){
  
  num = rho.hat*(rho.hat/alloc.prop)^gamma
  den = sum(num)
  
  return(num/den)
  
}

###############################################################################

sim_DBCD(N=100, dist="bern", K=1, par1 = c(0.1,0.8), burnin=3, ar="WMW",nsim=10^2)
