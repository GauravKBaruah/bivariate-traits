
source("bivariate_plant_pollinator.R")

require(ggplot2)
require(tidyverse)

  params <- list()
  params$matrix <- matrix(c(1,1,1,1,1,1),nrow=3,ncol=2)
  A <- ncol(params$matrix)
  P <- nrow(params$matrix)
  params$pref_matrix <- matrix(1,nrow=P,ncol=A)
  
  #species-specific genetic variances.
  var1_A <- runif(A, 0.005, 0.009)
  var2_A <- runif(A, 0.02, 0.06)
  
  var1_P <- runif(P, 0.004,0.0092)
  var2_P <- runif(P, 0.01, 0.05)
  
  rho<- 0
  
  #genetic covariance: diagonal = variance; off-diagonal = covariance.
  params$G_A <- array(0,dim=c(A,2,2))
  params$G_P <- array(0,dim=c(P,2,2))
  params$E_A <- array(0,dim=c(A,2,2))
  params$E_P <- array(0,dim=c(P,2,2))
  for (i in seq_len(A)) {
    cov12 <- rho * sqrt(var1_A[i] * var2_A[i])
    #variance-covariance matrix
    params$G_A[i,,] <-matrix(c(var1_A[i], cov12,cov12,var2_A[i]),nrow = 2,byrow = TRUE)
    params$E_A[i,,] <- matrix(c(0.05,0,0,0.02),2,2)
  }
  for (i in seq_len(P)) {
    cov12p <- rho * sqrt(var1_P[i] * var2_P[i])
    params$G_P[i,,] <-  matrix(c(var1_P[i], cov12p,cov12p,var2_P[i]), nrow = 2,byrow = TRUE)
    params$E_P[i,,] <- matrix(c(0.05,0,0,0.02),2,2)
  }
  #G and E must be symmetric positive; G+E must
  
  params$w <- c(1.5,1.5)          # mutualistic widths for traits 1 and 2
  params$w_c <- c(0.1,0.15)      # competitive widths for traits 1 and 2
  params$alpha0 <- 1        # interspecific competition amplitude
  params$alpha_ii <- 1      # constant trait-independent self-regulation
  params$gi <- 1.5
  params$ki <- 0.25
  params$bw <- 5
  params$a <- 0.11
  params$Temp <- 24           
  params$H <- 0.2
  params$mut.strength <- 1
  params$trait_evolve <- 1    # 0 switches off evolution of both means
  params$points <- 6          # 6 per axis = 36 nodes per species gaussian quadrature
  params$degree_scale <- TRUE
  params$dganimals <- colSums(params$matrix>0)
  params$dgplants <- rowSums(params$matrix>0)
  
  params$SA <- A
  params$SP <- P
  params$S <- A+P
  params$evolve_G <- FALSE
  
  Na <- rep(0.5,A)
  Np <- rep(0.5,P)
  muA1 <- c(18,22) #temperature optimum  polls
  muP1 <- c(18,20,22) #temperature optimum of plants
  muA2 <- c(1,2) # secondary trait--- need to think about trait limits
  muP2 <- c(1,1.5,2) 
  ic <- c(Na,Np,muA1,muP1,muA2,muP2)
  
  tmax <- 1000
  raw <- deSolve::ode(func=eqs_type2_bivariate,y=ic,parms=params,
                      times=seq(0,tmax,length.out=101),
                      rtol=1e-9,atol=1e-12)
  OUT <- organize_results_bivariate(raw,pars=params) 
  
  sum((OUT %>% filter(time == tmax, type == "mu"))$v)

  plot_all(OUT)
  plot_snapshot(dat = OUT,pars = params,moment = 1000)
  
  
  



#treatment grid- many runs
fact<- expand.grid(temperature = c(0,5,10,15,25,30,34), #temperatuer env
                   rho = seq(-0.9,0,0.9),  #correlation between two traits
                   replicate=1:3)

params <- list()

for(r in 1:nrow(fact)){
  
  
  params$matrix <- matrix(c(1,1,0,1,1,1,0,1,0,1,0,0),nrow=4,ncol=3,byrow = T)
  A <- ncol(params$matrix)
  P <- nrow(params$matrix)
  params$pref_matrix <- matrix(1,nrow=P,ncol=A)
  
  #species-specific genetic variances.
  var1_A <- runif(A, 0.05, 0.15)
  var2_A <- runif(A, 0.002, 0.006)
  
  var1_P <- runif(P, 0.04,0.12)
  var2_P <- runif(P, 0.001, 0.005)
  
  rho <- fact$rho[r]
  # genetic covariance: diagonal = variance; off-diagonal = covariance.
  params$G_A <- array(0,dim=c(A,2,2))
  params$G_P <- array(0,dim=c(P,2,2))
  params$E_A <- array(0,dim=c(A,2,2))
  params$E_P <- array(0,dim=c(P,2,2))
  for (i in seq_len(A)) {
    cov12 <- rho * sqrt(var1_A[i] * var2_A[i])
    #variance-covariance matrix
    params$G_A[i,,] <-matrix(c(var1_A[i], cov12,cov12,var2_A[i]),nrow = 2,byrow = TRUE)
    params$E_A[i,,] <- matrix(c(0.05,0,0,0.02),2,2)
  }
  for (i in seq_len(P)) {
    cov12p <- rho * sqrt(var1_P[i] * var2_P[i])
    params$G_P[i,,] <-  matrix(c(var1_P[i], cov12p,cov12p,var2_P[i]), nrow = 2,byrow = TRUE)
    params$E_P[i,,] <- matrix(c(0.05,0,0,0.02),2,2)
  }
  #G and E must be symmetric positive; G+E must
  
  params$w <- c(1,1.5)          # mutualistic widths for traits 1 and 2
  params$w_c <- c(0.1,0.15)      # competitive widths for traits 1 and 2
  params$alpha0 <- 1        # interspecific competition amplitude
  params$alpha_ii <- 1      # constant trait-independent self-regulation
  params$gi <- 1.5
  params$ki <- 0.25
  params$bw <- 5
  params$a <- 0.1
  params$Temp <- fact$temperature[r]          
  params$H <- 0.2
  params$mut.strength <- 1
  params$trait_evolve <- 1    # 0 switches off evolution of both means
  params$points <- 6          # 6 per axis = 36 nodes per species gaussian quadrature
  params$degree_scale <- TRUE
  params$dganimals <- colSums(params$matrix>0)
  params$dgplants <- rowSums(params$matrix>0)
  
  params$SA <- A
  params$SP <- P
  params$S <- A+P
  params$evolve_G <- FALSE
  
  Na <- rep(0.5,A)
  Np <- rep(0.5,P)
  muA1 <- runif(A, min= fact$temperature[r] -5 , fact$temperature[r]+5)
  muP1 <- runif(P, min= fact$temperature[r] -5 , fact$temperature[r]+5) #temperature optimum of plants
  muA2 <- runif(A, min= 0.1, 1) # secondary trait--- need to think about trait limits
  muP2 <- runif(P, 0.1,1)
  ic <- c(Na,Np,muA1,muP1,muA2,muP2)
  
  tmax <- 5000
  raw <- deSolve::ode(func=eqs_type2_bivariate,y=ic,parms=params,
                      times=seq(0,tmax,length.out=101),
                      rtol=1e-9,atol=1e-12)
  OUT <- organize_results_bivariate(raw,pars=params) 
  plot_all(OUT)
  plot_snapshot(dat = OUT,pars = params,moment = 5000)
  
  
  survivingspecies <-  which((OUT %>% filter(time == tmax, type == "N"))$v > 1e-5)
  speciestraitlag<- (OUT %>% filter(time == tmax, type == "mu1"))$v - fact$temperature[r]
  
  fact$biomass[r] <- sum((OUT %>% filter(time == tmax, type == "N"))$v)
  fact$richness[r] <- survivingspecies
  fact$traitlag[r] <- mean(speciestraitlag[survivingspecies])
  #trait-lag
  print(r)

  
   
}