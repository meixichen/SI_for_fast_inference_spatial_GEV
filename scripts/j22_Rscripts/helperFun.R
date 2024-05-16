multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Convert coordinates to long lat
convertCoords <- function(dt){
  # Variables for holding the coordinate system types 
  # Both ukgrid and irelandgrid
  ukgrid = "+init=epsg:27700"
  latlong = "+init=epsg:4326"
  irelandgrid = "+init=epsg:29903"
  # Set up for drawing on the map
  dt$ID <- c(1:dim(dt)[1])
  coordsNE <- cbind(Easting = as.numeric(as.character(dt$Easting)),
                    Northing = as.numeric(as.character(dt$Northing)))
  # Tells me which observations are from Ireland
  IRELANDBOOL <- nchar(dt$Station) > 5 & startsWith(dt$Station,prefix = "2")
  convertNEtoLL <- function(loc, Ireland, grid){
    GP_SP <- SpatialPointsDataFrame(loc[Ireland,], 
                                    data = data.frame(dt[Ireland,]$Station, 
                                                      dt[Ireland,]$ID), 
                                    proj4string = CRS(grid))
    # Convert from Eastings and Northings to Latitude and Longitude
    GP_SP_LL <- spTransform(GP_SP, CRS(latlong))
    # we also need to rename the columns
    colnames(GP_SP_LL@coords)[colnames(GP_SP_LL@coords) == "Easting"] <- "Longitude"
    colnames(GP_SP_LL@coords)[colnames(GP_SP_LL@coords) == "Northing"] <- "Latitude"
    coords <- data.table(long = GP_SP_LL@coords[,1], 
                         lat = GP_SP_LL@coords[,2], 
                         tempId = dt[Ireland,]$ID)
    return(coords)
  }
  
  if(sum(IRELANDBOOL) == 0){
    coords <- convertNEtoLL(coordsNE, !IRELANDBOOL, ukgrid)
  }else{
    coords <- rbind(convertNEtoLL(coordsNE, IRELANDBOOL, irelandgrid), 
                    convertNEtoLL(coordsNE, !IRELANDBOOL, ukgrid))
  }
  
  coords <- coords[order(coords$tempId),]
  coords <- coords[,1:2, with = F]
  coords <- as.matrix(coords)
  coords <- coords/10
  map <- maps::map('worldHires',
                   c('UK', 'Ireland', 'Isle of Man','Isle of Wight'),
                   xlim=c(-11,3), ylim=c(49,60.9), plot = FALSE)
  map_dt <- data.frame(long = map$x, lat = map$y)
  map_scaled <- data.table(long = map_dt$long/10 - min(coords[,1]),
                           lat = map_dt$lat/10 - min(coords[,2]))
  coords[,1] <- coords[,1]-min(coords[,1])
  coords[,2] <- coords[,2]-min(coords[,2])
  res <- list(coords = coords, map_scaled = map_scaled)
  return(res)
}


fit_model_psi <- function(data, desc){
  N <- dim(data)[1]
  # CJS submission: The following code is commented out because it is specific to the 
  # use case in Johannesson et al. (2022), e.g., specific to locations such as
  # Ireland
  #coords_obj <- convertCoords(data)
  #coords <- coords_obj$coords
  #map_scaled <- coords_obj$map_scaled
  coords <- cbind(data$long, data$lat)
  # CJS submission: modified mesh creation to match the mesh used in our paper
  mesh <- inla.mesh.2d(
    loc=coords,
    #offset = 0.08,
    #max.edge=0.07,
    ## discretization accuracy
    #cutoff=0.005,
    max.edge=2)
  A_mat <- inla.spde.make.A(mesh, loc=coords)
  N_s <- dim(A_mat)[2]
  prior.quantile.sd = 1; prior.quantile.range = 0.5
  spde <- 
    inla.spde2.pcmatern(mesh=mesh, prior.range=c(prior.quantile.range, 0.95), 
                        prior.sigma=c(prior.quantile.sd, 0.01))
  # Make the C matrix for generic0 random effect
  scale_psi <- 1/data$v_p
  stack_psi <- inla.stack(
    data=list(Y=data$psi), 
    A=list(A_mat,1),tag='dat',
    effect = list(list(s=1:N_s),
                  list(idx = 1:N, X_psi = desc)))
  formula_psi <- Y ~ -1 + X_psi +
    f(s, model=spde) + 
    f(idx,model="iid",hyper=list(theta=list(prior="pc.prec", param=c(.6,.01))))
  
  control.family = list(hyper = list(prec = list(
    initial = 0, fixed = TRUE)))
  
  model <- inla(formula_psi, family='gaussian', 
                control.compute=list(dic=TRUE),
                data=inla.stack.data(stack_psi), 
                control.predictor=list( A=inla.stack.A(stack_psi), compute = TRUE),
                control.family = control.family,
                scale = scale_psi)
  return(list(mdl = model, spde = spde, mesh = mesh))
} 
fit_model_tau <- function(data, desc){
  N <- dim(data)[1]
  #coords_obj <- convertCoords(data)
  #coords <- coords_obj$coords
  #map_scaled <- coords_obj$map_scaled
  coords <- cbind(data$long, data$lat)
  # CJS submission: modified mesh creation to match the mesh used in our paper
  mesh <- inla.mesh.2d(
    loc=coords,
    #offset = 0.08,
    #max.edge=0.07,
    ## discretization accuracy
    #cutoff=0.005,
    max.edge=2)
  A_mat <- inla.spde.make.A(mesh, loc=coords)
  N_s <- dim(A_mat)[2]
  prior.quantile.sd = 1; prior.quantile.range = 0.5
  spde <- 
    inla.spde2.pcmatern(mesh=mesh, prior.range=c(prior.quantile.range, 0.95), 
                        prior.sigma=c(prior.quantile.sd, 0.01))
  # Make the C matrix for generic0 random effect
  #sigma <- makeSigma_etay(data)
  #scale_psi <- 1/diag(sigma)[1:N]
  scale_tau <- 1/data$v_t
  stack_tau <- inla.stack(
    data=list(Y=data$tau), 
    A=list(A_mat,1),tag='dat',
    effect = list(list(s=1:N_s),
                  list(idx = 1:N, X_tau = desc)))
  formula_tau <- Y ~ -1 + X_tau +
    f(s, model=spde) + 
    f(idx,model="iid",hyper=list(theta=list(prior="pc.prec", param=c(.6,.01))))
  control.family = list(hyper = list(prec = list(
    initial = 0, fixed = TRUE)))
  
  model <- inla(formula_tau, family='gaussian', 
                control.compute=list(dic=TRUE),
                data=inla.stack.data(stack_tau), 
                control.predictor=list( A=inla.stack.A(stack_tau), compute = TRUE),
                control.family = control.family,
                scale = scale_tau)
  return(list(mdl = model, spde = spde, mesh = mesh))
}
fit_model_kappa <- function(data, desc){
  N <- dim(data)[1]
  # Make the C matrix for generic0 random effect
  #sigma <- makeSigma_etay(data)
  #scale_psi <- 1/diag(sigma)[1:N]
  scale_xi <- 1/data$v_k
  stack_xi <- inla.stack(
    data=list(Y = data$kappa), 
    A=list(1),tag='dat',
    effect = list(list(idx = 1:N, X_xi = desc)))
  formula_xi <- Y ~ -1 + X_xi + 
    f(idx,model="iid",hyper=list(theta=list(prior="pc.prec", param=c(.6,.01))))
  
  control.family = list(hyper = list(prec = list(
    initial = 0, fixed = TRUE)))
  
  model <- inla(formula_xi, family='gaussian', 
                data=inla.stack.data(stack_xi), 
                control.predictor=list( A=inla.stack.A(stack_xi), compute = TRUE),
                control.family = control.family,
                scale = scale_xi)
  return(list(mdl = model))
}
fit_model_gamma <- function(data, desc){
  N <- dim(data)[1]
  
  scale_tau <- 1/data$v_g
  stack_tau <- inla.stack(
    data=list(Y=data$gamma), 
    A=list(1),tag='dat',
    effect = list(list(idx = 1:N, X_tau = desc)))
  formula_tau <- Y ~ -1 + X_tau +
    f(idx,model="iid",hyper=list(theta=list(prior="pc.prec", param=c(.1,.05))))
  control.family = list(hyper = list(prec = list(
    initial = 0, fixed = TRUE)))
  
  model <- inla(formula_tau, family='gaussian', 
                data=inla.stack.data(stack_tau), 
                control.predictor=list( A=inla.stack.A(stack_tau), compute = TRUE),
                control.family = control.family,
                scale = scale_tau)
  return(list(mdl = model))
}

transform <- function(x, descriptor){
  if(descriptor == "BFIHOST"){return(x^2)}
  if(descriptor == "URBEXT"){return(log(x+1))}
  return(log(x))
}

makeQ_u <- function(s, rho, spde){
  Q <- inla.spde.precision(spde, theta=c(log(rho), log(s)))
  return(Q)
}

log_prior_logSp <- function(u,v){
  # Prob(sigma < sigma0) = alpha_1
  # Prob(range > range0) = alpha_2
  sigma0 = 1; range0 = 0.5
  alpha_1 = .95; alpha_2 = .01
  l1 <- -log(alpha_1)*range0
  l2 <- -log(alpha_2)/sigma0
  res <- log(l1) + log(l2) - v + u -l1*exp(-v)-l2*exp(u)
  return(res)
}

log_prior_logPrecision <- function(kappa){
  # Prob(sigma > u) = alpha
  u <- 0.6; alpha <- .05
  l <- -log(alpha)/u
  res <- log(l/2)-l*exp(-kappa/2)-kappa/2
  return(res)
}

log_prior_logPrecision2 <- function(kappa){
  # Prob(sigma > u) = alpha
  u <- 0.1; alpha <- .01
  l <- -log(alpha)/u
  res <- log(l/2)-l*exp(-kappa/2)-kappa/2
  return(res)
}

log_prior_logTheta <- function(kappa){
  sp_psi <- log_prior_logSp(kappa$u_psi, kappa$v_psi)
  sp_tau <- log_prior_logSp(kappa$u_tau, kappa$v_tau)
  pre_psi <- log_prior_logPrecision(kappa$kappa_psi)
  pre_tau <- log_prior_logPrecision(kappa$kappa_tau)
  pre_xi <- log_prior_logPrecision(kappa$kappa_xi)
  pre_gamma <- log_prior_logPrecision(kappa$kappa_gamma)
  res <- sp_psi + sp_tau + pre_psi + pre_tau + pre_xi + pre_gamma
  return(res)
}
log_prior_logTheta_noGamma <- function(kappa){
  sp_psi <- log_prior_logSp(kappa$u_psi, kappa$v_psi)
  sp_tau <- log_prior_logSp(kappa$u_tau, kappa$v_tau)
  pre_psi <- log_prior_logPrecision(kappa$kappa_psi)
  pre_tau <- log_prior_logPrecision(kappa$kappa_tau)
  pre_xi <- log_prior_logPrecision(kappa$kappa_xi)
  res <- sp_psi + sp_tau + pre_psi + pre_tau + pre_xi
  return(res)
}

