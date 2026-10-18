
rm(list=ls())
library(deSolve)
library(statmod)
#Network rows = plants; columns = animals.

#Smooth abundance cutoff
cutoff <- function(x) {
  ifelse(x < 1, (x > 0)*x^3*(10 + x*(-15 + 6*x)), 1)
}

#Gaussian quadrature for a correlated bivariate normal distribution.
#Start with independent standard normals u and v. Construct:
#z1 = mu1 + sigma1*u
#z2 = mu2 + sigma2*(rho*u + sqrt(1-rho^2)*v).
#This gives variances P11/P22 and covariance P12.
bivariate_quadrature <- function(mu1, mu2, P, points) {
  sigma1 <- sqrt(P[1,1])
  sigma2 <- sqrt(P[2,2])
  rho <- P[1,2]/(sigma1*sigma2)
  normal_rule <- statmod::gauss.quad.prob(points, dist="normal")
  nodes <- normal_rule$nodes
  weights <- normal_rule$weights
  z1 <- numeric(points^2)
  z2 <- numeric(points^2)
  weight <- numeric(points^2)
  index <- 0
  for (a in seq_len(points)) {
    for (b in seq_len(points)) {
      index <- index + 1
      z1[index] <- mu1 + sigma1*nodes[a]
      z2[index] <- mu2 + sigma2*(rho*nodes[a] + sqrt(1-rho^2)*nodes[b])
      weight[index] <- weights[a]*weights[b]
    }
  }
  data.frame(z1=z1, z2=z2, weight=weight)
}

# Type-II benefit for one animal species.
# Here G is a returned benefit name, not the genetic covariance matrix.
#Both focal and partner distributions are approximated by quadrature,
type_2_animals_bivariate <- function(i, Na, Np, muA, muP, Pa, Pp, pars) {
  animal_nodes <- bivariate_quadrature(muA[i,1],muA[i,2],Pa[i,,],pars$points)
  growth_benefit <- numeric(nrow(animal_nodes))
  trait1_benefit <- numeric(nrow(animal_nodes))
  trait2_benefit <- numeric(nrow(animal_nodes))
  degree <- 1
  if (pars$degree_scale) degree <- max(1,pars$dganimals[i])
  for (a in seq_len(nrow(animal_nodes))) {
    z1 <- animal_nodes$z1[a]
    z2 <- animal_nodes$z2[a]
    mutualistic_input <- 0
    for (k in seq_along(Np)) {
      plant_nodes <- bivariate_quadrature(muP[k,1],muP[k,2],Pp[k,,],pars$points)
      matching <- exp(-(z1-plant_nodes$z1)^2/pars$w[1]^2 -
                       (z2-plant_nodes$z2)^2/pars$w[2]^2)
      # Integral over both traits of partner species k.
      average_matching <- sum(matching*plant_nodes$weight)
      mutualistic_input <- mutualistic_input +
        Np[k]*pars$matrix[k,i]*pars$pref_matrix[k,i]*
        pars$mut.strength*average_matching
    }
    # degree scaling affects only the numerator
    benefit <- (mutualistic_input/degree)/(1+pars$H*mutualistic_input)
    growth_benefit[a] <- benefit*animal_nodes$weight[a]
    trait1_benefit[a] <- (z1-muA[i,1])*benefit*animal_nodes$weight[a]
    trait2_benefit[a] <- (z2-muA[i,2])*benefit*animal_nodes$weight[a]
  }
  list(G=sum(growth_benefit), B1=sum(trait1_benefit), B2=sum(trait2_benefit))
}

#type-II benefit for one plant species.
type_2_plants_bivariate <- function(i, Na, Np, muA, muP, Pa, Pp, pars) {
  plant_nodes <- bivariate_quadrature(muP[i,1],muP[i,2],Pp[i,,],pars$points)
  growth_benefit <- numeric(nrow(plant_nodes))
  trait1_benefit <- numeric(nrow(plant_nodes))
  trait2_benefit <- numeric(nrow(plant_nodes))
  degree <- 1
  if (pars$degree_scale) degree <- max(1,pars$dgplants[i])

  for (p in seq_len(nrow(plant_nodes))) {
    z1 <- plant_nodes$z1[p]
    z2 <- plant_nodes$z2[p]
    mutualistic_input <- 0
    for (k in seq_along(Na)) {
      animal_nodes <- bivariate_quadrature(muA[k,1],muA[k,2],Pa[k,,],pars$points)
      matching <- exp(-(z1-animal_nodes$z1)^2/pars$w[1]^2 -
                       (z2-animal_nodes$z2)^2/pars$w[2]^2)
      average_matching <- sum(matching*animal_nodes$weight)
      mutualistic_input <- mutualistic_input +
        Na[k]*pars$matrix[i,k]*pars$pref_matrix[i,k]*
        pars$mut.strength*average_matching
    }
    benefit <- (mutualistic_input/degree)/(1+pars$H*mutualistic_input)
    growth_benefit[p] <- benefit*plant_nodes$weight[p]
    trait1_benefit[p] <- (z1-muP[i,1])*benefit*plant_nodes$weight[p]
    trait2_benefit[p] <- (z2-muP[i,2])*benefit*plant_nodes$weight[p]
  }
  list(G=sum(growth_benefit), B1=sum(trait1_benefit), B2=sum(trait2_benefit))
}

# ODE 
#State order 
#c(Na, Np, muA1, muP1, muA2, muP2).
# For each species:
# dN/dt = N * mean growth
#dmu/dt = G %*% solve(P) %*% E[(z-mu)*r(z)].
eqs_type2_bivariate <- function(time, state, pars) {
  A <- ncol(pars$matrix)
  P <- nrow(pars$matrix)
  S <- A+P

  #Extract abundances and the two mean traits.
  Na <- pmax(0,state[seq_len(A)])
  Np <- pmax(0,state[A+seq_len(P)])
  muA1 <- state[S+seq_len(A)]
  muP1 <- state[S+A+seq_len(P)]
  muA2 <- state[2*S+seq_len(A)]
  muP2 <- state[2*S+A+seq_len(P)]
  muA <- cbind(muA1,muA2)
  muP <- cbind(muP1,muP2)

  #cenetic and environmental covariance matrices do not evolve.
  Ga <- pars$G_A
  Gp <- pars$G_P
  Pa <- Ga+pars$E_A
  Pp <- Gp+pars$E_P
   #intrinsic growth and its centered trait moments 
  # b(z1,T) = gi/s * exp(-(T-z1)^2/(2*s^2)) - ki.
  #s=bw-a*mu1 is held fixed when comparing individual phenotypes.
  Temp <- pars$Temp
  if (is.function(pars$Temp)) Temp <- pars$Temp(time)
  ba <- numeric(A)
  bp <- numeric(P)
  bar_ba <- matrix(0,A,2)
  bar_bp <- matrix(0,P,2)
  for (i in seq_len(A)) {
    thermal_width <- pars$bw-pars$a*muA1[i]
    total_width <- thermal_width^2+Pa[i,1,1]
    birth <- pars$gi/sqrt(total_width)*exp(-(Temp-muA1[i])^2/(2*total_width))
    ba[i] <- birth-pars$ki
    bar_ba[i,] <- Pa[i,,][,1]*birth*(Temp-muA1[i])/total_width
  }
  for (i in seq_len(P)) {
    thermal_width <- pars$bw-pars$a*muP1[i]
    total_width <- thermal_width^2+Pp[i,1,1]
    birth <- pars$gi/sqrt(total_width)*exp(-(Temp-muP1[i])^2/(2*total_width))
    bp[i] <- birth-pars$ki
    bar_bp[i,] <- Pp[i,,][,1]*birth*(Temp-muP1[i])/total_width
  }

  #Competition exact gaussian average -- gaussian*gaussian = gaussian
  #alpha_ij = sqrt(det(Wc)/det(Wc+Pi+Pj))*exp(-dm' inv(D) dm/2).
  #beta_ij = E[(z_i-mu_i)*alpha(z_i,z_j)] = -Pi inv(D) dm alpha_ij.
  #The competition term is subtracted from growth AND centered moments.
  Wc <- diag(pars$w_c^2/2) # original exp(-d1^2/wc1^2-d2^2/wc2^2)
  alphaA <- matrix(0,A,A)
  alphaP <- matrix(0,P,P)
  betaA1 <- matrix(0,A,A)
  betaA2 <- matrix(0,A,A)
  betaP1 <- matrix(0,P,P)
  betaP2 <- matrix(0,P,P)
  for (i in seq_len(A)) {
    for (j in seq_len(A)) {
      if (i==j) {
        alphaA[i,j] <- pars$alpha_ii
      } else {
        D <- Wc+Pa[i,,]+Pa[j,,]
        difference <- muA[i,]-muA[j,]
        alphaA[i,j] <- pars$alpha0*sqrt(det(Wc)/det(D))*
          exp(-as.numeric(t(difference)%*%solve(D)%*%difference)/2)
        beta <- -Pa[i,,]%*%solve(D)%*%difference*alphaA[i,j]
        betaA1[i,j] <- beta[1]
        betaA2[i,j] <- beta[2]
      }
    }
  }
  for (i in seq_len(P)) {
    for (j in seq_len(P)) {
      if (i==j) {
        alphaP[i,j] <- pars$alpha_ii
      } else {
        D <- Wc+Pp[i,,]+Pp[j,,]
        difference <- muP[i,]-muP[j,]
        alphaP[i,j] <- pars$alpha0*sqrt(det(Wc)/det(D))*
          exp(-as.numeric(t(difference)%*%solve(D)%*%difference)/2)
        beta <- -Pp[i,,]%*%solve(D)%*%difference*alphaP[i,j]
        betaP1[i,j] <- beta[1]
        betaP2[i,j] <- beta[2]
      }
    }
  }

  #Mutualism: numerical gaussian quadrature
  ai <- numeric(A)
  aj <- numeric(P)
  bi <- matrix(0,A,2)
  bj <- matrix(0,P,2)
  for (i in seq_len(A)) {
    benefit <- type_2_animals_bivariate(i,Na,Np,muA,muP,Pa,Pp,pars)
    ai[i] <- benefit$G
    bi[i,1] <- benefit$B1
    bi[i,2] <- benefit$B2
  }
  for (i in seq_len(P)) {
    benefit <- type_2_plants_bivariate(i,Na,Np,muA,muP,Pa,Pp,pars)
    aj[i] <- benefit$G
    bj[i,1] <- benefit$B1
    bj[i,2] <- benefit$B2
  }

  #abundance derivatives 
  dndt_a <- Na*(ba-as.numeric(alphaA%*%Na)+ai)*cutoff(Na/1e-6)
  dndt_p <- Np*(bp-as.numeric(alphaP%*%Np)+aj)*cutoff(Np/1e-6)

  #Mean-trait derivatives
  # S1/S2 below are centered fitness moments, NOT selection gradients.
  # Multiplying by G %*% solve(P) gives the inherited response.
  dudt_A <- matrix(0,A,2)
  dudt_P <- matrix(0,P,2)
  for (i in seq_len(A)) {
    S1 <- bar_ba[i,1]-sum(betaA1[i,]*Na)+bi[i,1]
    S2 <- bar_ba[i,2]-sum(betaA2[i,]*Na)+bi[i,2]
    selection_moments <- c(S1,S2)
    dudt_A[i,] <- pars$trait_evolve*Ga[i,,]%*%solve(Pa[i,,])%*%selection_moments
    if (Na[i]==0) dudt_A[i,] <- c(0,0)
  }
  for (i in seq_len(P)) {
    S1 <- bar_bp[i,1]-sum(betaP1[i,]*Np)+bj[i,1]
    S2 <- bar_bp[i,2]-sum(betaP2[i,]*Np)+bj[i,2]
    selection_moments <- c(S1,S2)
    dudt_P[i,] <- pars$trait_evolve*Gp[i,,]%*%solve(Pp[i,,])%*%selection_moments
    if (Np[i]==0) dudt_P[i,] <- c(0,0)
  }
  list(c(dndt_a,dndt_p,dudt_A[,1],dudt_P[,1],dudt_A[,2],dudt_P[,2]))
}

#results organizer
organize_results_bivariate <- function(sol,pars) {
  sol <- as.matrix(sol); S <- pars$S
  types <- c("N","mu1","mu2",if(pars$evolve_G)c("G11","G12","G22"))
  if(ncol(sol)!=1+S*length(types)) stop("Unexpected solution columns")
  identity <- data.frame(guild=c(rep("A",pars$SA),rep("P",pars$SP)),
                         species=c(seq_len(pars$SA),seq_len(pars$SP)))
  out <- lapply(seq_along(types),function(k) {
    data.frame(time=rep(sol[,1],each=S),
               guild=rep(identity$guild,nrow(sol)),
               species=rep(identity$species,nrow(sol)),type=types[k],
               v=as.vector(t(sol[,1+(k-1)*S+seq_len(S),drop=FALSE])))
  })
  ans <- do.call(rbind,out); rownames(ans) <- NULL
  for(nm in c("gamma","nestedness","C","theta","web.name"))
    if(length(pars[[nm]])==1) ans[[nm]] <- pars[[nm]]
  ans
}





# Population abundances
plot_density <- function(dat) {
  abundance <- dat[dat$type == "N", ]
  ggplot(abundance,aes(x = time, y = v, colour = factor(species))) +
    geom_line(linewidth = 0.8) +
    facet_wrap( ~ guild,scales = "free_y",labeller = as_labeller(c(A = "Pollinators", P = "Plants")) ) +
    labs(x = "Time",y = "Population density",colour = "Species") +
    theme_classic() +
    theme(legend.position = "none")
}


# Mean of trait 1
plot_trait1 <- function(dat) {
  trait1 <- dat[dat$type == "mu1", ]
  ggplot(trait1,aes(x = time, y = v, colour = factor(species))) +
    geom_line(linewidth = 0.8) +
    facet_wrap( ~ guild,scales = "free_y",labeller = as_labeller(c(A = "Pollinators", P = "Plants"))) +
    labs(x = "Time",y = "Mean trait 1",colour = "Species") +
    theme_classic() +
    theme(legend.position = "none")
}


# Mean of trait 2
plot_trait2 <- function(dat) {
  trait2 <- dat[dat$type == "mu2", ]
  ggplot(trait2,aes(x = time, y = v, colour = factor(species))) +
    geom_line(linewidth = 0.8) +facet_wrap( ~ guild,scales = "free_y",labeller = as_labeller(c(A = "Pollinators", P = "Plants"))) +
    labs(x = "Time",y = "Mean trait 2",colour = "Species") +
    theme_classic() +
    theme(legend.position = "none")
}


# Combine abundances and both trait means
plot_all <- function(dat) {
  cowplot::plot_grid(plot_density(dat), plot_trait1(dat),plot_trait2(dat),ncol = 1,
    align = "v",
    labels = c("Abundance", "Trait 1", "Trait 2"),
    label_size = 12)
}



plot_snapshot <- function(dat, pars, moment = 0,limits1 = NULL,limits2 = NULL,res = 100,levels = c(0.1, 0.5, 0.9)) {
  
  saved_times <- sort(unique(dat$time))
  snapshot_time <- saved_times[which.min(abs(saved_times - moment))]
  snap <- dat[dat$time == snapshot_time, ]
  species_table <- unique(snap[, c("guild", "species")])
  
  #
  species_table$N <- NA
  species_table$mu1 <- NA
  species_table$mu2 <- NA
  species_table$P11 <- NA
  species_table$P12 <- NA
  species_table$P22 <- NA
  
  for (i in seq_len(nrow(species_table))) {
    guild <- species_table$guild[i]
    species <- species_table$species[i]
    one_species <- snap[snap$guild == guild & snap$species == species,]
    species_table$N[i] <- one_species$v[one_species$type == "N"]
    species_table$mu1[i] <- one_species$v[ one_species$type == "mu1"]
    species_table$mu2[i] <- one_species$v[one_species$type == "mu2"]
    
    # Environmental and fixed genetic covariance matrices
    if (guild == "A") {
      G <- pars$G_A[species, , ]
      E <- pars$E_A[species, , ]
    } else {
      G <- pars$G_P[species, , ]
      E <- pars$E_P[species, , ]
    }
  
    P <- G + E
    
    species_table$P11[i] <- P[1, 1]
    species_table$P12[i] <- P[1, 2]
    species_table$P22[i] <- P[2, 2]
  }
  
  
  species_table <- species_table[species_table$N > 0, ]
  determinant <- with(species_table,P11 * P22 - P12^2)

  if (is.null(limits1)) {
    limits1 <- c(
      min(species_table$mu1 - 4 * sqrt(species_table$P11)),max(species_table$mu1 + 4 * sqrt(species_table$P11)))
  }
  
  if (is.null(limits2)) {
    limits2 <- c(min(species_table$mu2 - 4 * sqrt(species_table$P22)),
      max(species_table$mu2 + 4 * sqrt(species_table$P22)))
  }
  

  trait_grid <- expand.grid(
    z1 = seq(limits1[1], limits1[2], length.out = res),
    z2 = seq(limits2[1], limits2[2], length.out = res)
  )
  
  distributions <- vector("list", nrow(species_table))
  
  for (i in seq_len(nrow(species_table))) {
    
    mu1<- species_table$mu1[i]
    mu2 <- species_table$mu2[i]
     P11 <- species_table$P11[i]
    P12<- species_table$P12[i]
    P22<- species_table$P22[i]
    Delta <-P11 * P22 - P12^2
    
    z1 <-trait_grid$z1
    z2<- trait_grid$z2
    
    exponent <- -( P22 * (z1 - mu1)^2 - 2 *P12 * (z1 - mu1) *(z2 - mu2) +P11 *(z2 - mu2)^2) / (2 * Delta)
    
    # probability density p(z1, z2)
    probability_density <- exp(exponent)/(2 * pi * sqrt(Delta))
    
    one_distribution <- trait_grid
   one_distribution$guild <- species_table$guild[i]
    one_distribution$species <- species_table$species[i]
    
    #  density in trait space: N * p(z1, z2)
    one_distribution$density <-
      species_table$N[i] * probability_density
    one_distribution$relative_density <- exp(exponent)
    
    distributions[[i]] <- one_distribution
  }
  
  distributions <- do.call(rbind, distributions)
  
  ggplot2::ggplot(distributions,ggplot2::aes(x = z1,y = z2,group = factor(species),colour = factor(species))
  ) +ggplot2::geom_contour(ggplot2::aes(z = relative_density),
      breaks = levels,linewidth = 0.6) +
    ggplot2::geom_point(data = species_table,ggplot2::aes(x = mu1,
        y = mu2,colour = factor(species)),  inherit.aes = FALSE, size = 2) +
    ggplot2::facet_wrap( ~ guild,labeller = ggplot2::as_labeller(
        c(A = "Pollinators", P = "Plants"))) +
    ggplot2::labs(x = "Trait 1", y = "Trait 2",colour = "Species",title = paste("Trait distributions at time", snapshot_time) ) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none")
}

