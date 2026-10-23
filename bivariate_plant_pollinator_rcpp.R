# R interface to the C++ implementation of your attached fixed-covariance model.
# Keep this file and bivariate_plant_pollinator.cpp in the same directory.
# R must have Rcpp, statmod and deSolve installed, plus a working C++ compiler.
# Source this file ONCE before your simulation loop.

.bivariate_source_file <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
.bivariate_source_dir <- if (is.null(.bivariate_source_file)) getwd() else
  dirname(normalizePath(.bivariate_source_file))
Rcpp::sourceCpp(file.path(.bivariate_source_dir,"bivariate_plant_pollinator.cpp"))
rm(.bivariate_source_file,.bivariate_source_dir)

# Build one C++ model per simulation, after constructing all your parameters.
# Its fixed covariance transforms, quadrature offsets and competition inverses
# are calculated once, rather than at every derivative evaluation.
make_bivariate_rhs_cpp <- function(params) {

  params$matrix <- as.matrix(params$matrix)
  A <- ncol(params$matrix)
  P <- nrow(params$matrix)
  S <- A+P
  if (A<1 || P<1 || any(!is.finite(params$matrix)) || any(params$matrix<0))
    stop("matrix must be a finite, nonnegative plants-by-animals matrix")

  if (isTRUE(params$evolve_G))
    stop("This translates the attached model with FIXED G. Set evolve_G=FALSE.")

  # This model, like the supplied R script, uses scalar ecological parameters.
  for (name in c("gi","ki","bw","a","H","mut.strength",
                 "alpha0","alpha_ii","trait_evolve")) {
    value <- params[[name]]
    if (!is.numeric(value) || length(value)!=1 || !is.finite(value))
      stop(name," must be one finite number")
    params[[name]] <- as.double(value)
  }
  for (name in c("gi","ki","H","mut.strength","alpha0","alpha_ii","trait_evolve"))
    if (params[[name]]<0) stop(name," must be nonnegative")

  for (name in c("w","w_c")) {
    value <- params[[name]]
    if (!is.numeric(value) || length(value)!=2 || any(!is.finite(value)) || any(value<=0))
      stop(name," must contain two positive widths")
    params[[name]] <- as.double(value)
  }

  if (is.null(params$pref_matrix)) params$pref_matrix <- matrix(1,P,A)
  params$pref_matrix <- as.matrix(params$pref_matrix)
  if (!identical(dim(params$pref_matrix),c(P,A)) ||
      any(!is.finite(params$pref_matrix)) || any(params$pref_matrix<0))
    stop("pref_matrix must be finite, nonnegative, and match matrix dimensions")

  if (is.null(params$degree_scale)) params$degree_scale <- FALSE
  if (!is.logical(params$degree_scale) || length(params$degree_scale)!=1 ||
      is.na(params$degree_scale)) stop("degree_scale must be TRUE or FALSE")
  if (is.null(params$dganimals)) params$dganimals <- colSums(params$matrix>0)
  if (is.null(params$dgplants)) params$dgplants <- rowSums(params$matrix>0)
  if (params$degree_scale) {
    if (length(params$dganimals)!=A || length(params$dgplants)!=P ||
        any(!is.finite(c(params$dganimals,params$dgplants))) ||
        any(c(params$dganimals,params$dgplants)<0)) stop("Invalid degree vectors")
  }

  # Covariance arrays are [species, trait, trait]. G and E may be singular,
  # but their sum must be positive definite (including rho=+1/-1 treatments).
  for (guild in c("A","P")) {
    number <- if (guild=="A") A else P
    for (kind in c("G","E")) {
      name <- paste0(kind,"_",guild)
      covariance <- params[[name]]
      if (!is.numeric(covariance) ||
          !identical(dim(covariance),c(number,2L,2L)) ||
          any(!is.finite(covariance))) stop("Invalid covariance array: ",name)
      for (i in seq_len(number)) {
        one <- covariance[i,,]
        scale <- max(1,max(abs(one)))
        if (max(abs(one-t(one)))>1e-12*scale ||
            min(eigen(one,symmetric=TRUE,only.values=TRUE)$values)< -1e-12*scale)
          stop(name," must contain symmetric positive semidefinite matrices")
      }
    }
    total <- params[[paste0("G_",guild)]]+params[[paste0("E_",guild)]]
    for (i in seq_len(number))
      if (min(eigen(total[i,,],symmetric=TRUE,only.values=TRUE)$values)<=0)
        stop("G+E must be positive definite for every species")
  }

  if (is.null(params$points)) params$points <- 6L
  if (length(params$points)!=1 || !is.finite(params$points) ||
      params$points<1 || params$points!=as.integer(params$points))
    stop("points must be a positive integer")
  normal_rule <- statmod::gauss.quad.prob(params$points,dist="normal")

  # Copy all fixed values to C++. Each closure owns a separate model.
  pointer <- prepare_bivariate_cpp(params,normal_rule$nodes,normal_rule$weights)
  temperature <- params$Temp
  if (!is.function(temperature) &&
      (!is.numeric(temperature) || length(temperature)!=1 || !is.finite(temperature)))
    stop("Temp must be one finite number or a function(time)")

  function(time,state,parms) {
    Temp <- if (is.function(temperature)) temperature(time) else temperature
    list(bivariate_rhs_cpp(pointer,state,Temp))
  }
}

# Convenience function: automatically rebuild fixed calculations for each run.
# deSolve still controls integration and error tolerances. The numerical RHS
# calculations execute in C++; this is not a replacement ODE solver.
ode_bivariate_cpp <- function(y,times,parms,...) {
  rhs <- make_bivariate_rhs_cpp(parms)
  deSolve::ode(y=y,times=times,func=rhs,parms=parms,...)
}

# The results organizer and all plotting functions below are copied unchanged
# from your uploaded script. No rm(list=ls()) is executed when sourcing.
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

