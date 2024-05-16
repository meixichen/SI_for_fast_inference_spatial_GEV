# Constants
xi00 <- 0
alp3 <- 0.8 # c_phi
alpha <- 4
beta <- 4
sig3 = (log(1 - (xi00+0.5)^alp3))*(1 - (xi00+0.5)^alp3)*(-(alp3)^(-1)*(xi00+0.5)^(-alp3+1)) # b_phi
b3 <- -sig3*(log(-log(1-0.5^alp3))) # a_phi
# Likelihood function to minimize
fn <- function(theta,y) {
  sig3 <- (log(1 - (xi00+0.5)^alp3))*(1 - (xi00+0.5)^alp3)*(-(alp3)^(-1)*(xi00+0.5)^(-alp3+1))
  b3 <- -sig3*(log(-log(1-0.5^alp3)))
  xitheta <- (1 - exp(-exp((theta[3]-b3)/sig3)))^(1/alp3) - 0.5
  -sum(dgev(y, loc = exp(theta[1]), scale = exp(theta[2] + theta[1]),
            shape = ((1-(exp(-exp((theta[3] - b3)/sig3))))^(1/alp3) - 0.5), log = TRUE)) -
    ((alpha - alp3)*log(xitheta + 0.5) + (beta-1)    *log(0.5 - xitheta) + (theta[3]-b3)/sig3 - exp((theta[3]-b3)/sig3)  )
}


####################### ML step ##############################################
dt_mles <- data.frame()
ms_start_t <- Sys.time() # START time
suppressWarnings(
  for(i in 1:n_st){
    station <- stations[i]
    currentS <- data_max %>% filter(Station == station)
    fitgev <- fgev(currentS$data)
    mu0 <- fitgev$estimate[1]
    sigma0 <- fitgev$estimate[2]
    xi0 <- fitgev$estimate[3]
    if (xi0 > 0) {
      xi0 <- min(xi0,0.45)
    } else {
      xi0 <- max(xi0,-0.45)
    }
    theta0 <- c(log(mu0),log(sigma0)-log(mu0),b3 + sig3*log(-log(1 - (xi0+0.5)^alp3)))
    GEV_fit <- nlm(fn, theta <- theta0, y <- currentS$data, hessian=TRUE)
    
    S_d <- try(solve(GEV_fit$hessian),silent=T)
    
    res <- data.frame(
      Station = station,
      psi = GEV_fit$estimate[1],
      tau = GEV_fit$estimate[2],
      kappa = GEV_fit$estimate[3],
      v_p = S_d[1,1],
      v_p_t = S_d[1,2],
      v_p_k = S_d[1,3],
      v_t = S_d[2,2],
      v_t_k = S_d[2,3],
      v_k = S_d[3,3]
    )
    dt_mles <- rbind(res, dt_mles)
  }
)

################################## Make covariates ready #####################################################################
# To make sure that the covariates are in the same order as ml estimates
# CJS submission: we modified the covariates here we don't have any
names <- c(
  "Int."
)
ncols <- length(names)
covariates <- matrix(c(
  rep(1,n_st) # Intercept
), ncol = 1, nrow = n_st
)
colnames(covariates) <- names

################################ Make X matrices for psi, tau and xi #########################################################
cov_names_psi <- c("Int.")
cov_names_tau <- c("Int.")
cov_names_kappa <- c("Int.")
X_psi <- as.matrix(covariates[, c("Int.")])
X_tau <- as.matrix(covariates[, c("Int.")])
X_kappa <- as.matrix(covariates[, c("Int.")])

############################### Make mesh and A matrices for spatial component ###############################################
coords <- cbind(locs_ms$long, locs_ms$lat)
# CJS submission: we modified the mesh to match the mesh we used in our main paper
mesh <- inla.mesh.2d(
  loc=coords,
  #offset = 0.08,
  #	cutoff = 0.005,
  max.edge = 2)
# Make A matrix - The same for psi and tau
A_mat <- inla.spde.make.A(mesh, loc=coords)
A_psi <- A_mat
A_tau <- A_mat

################# Fit INLA models, to have good starting proposal distribution ###############
# Add location to MLE data
dt_mles <- dt_mles %>%
  left_join(locs_ms %>% dplyr::select(Station,long,lat))
mdl_sep_psi <- fit_model_psi(data = dt_mles, desc = X_psi)
mdl_sep_tau <- fit_model_tau(data = dt_mles, desc = X_tau)
mdl_sep_xi <- fit_model_kappa(data = dt_mles, desc = X_kappa)

# The spde object is only used for getting the precision matrix for a given range and sd
# so priors do not matter here
d_spde <- inla.spde2.pcmatern(mesh, prior.range = c(.5, .5), prior.sigma = c(.5, .5))

# Make the Z matrix
Z <- bdiag(cbind(X_psi,A_psi), cbind(X_tau, A_tau), X_kappa)

N_psi <- dim(X_psi)[2]
N_tau <- dim(X_tau)[2]
N_xi <- dim(X_kappa)[2]
N_colsA <- dim(A_psi)[2]

################################################ Helper functions ###################################################
# Make the covariance matrix for the observations
makeSigma_etay <- function(dt){
  N_rows <- dim(dt)[1]
  Sigma_etay <- matrix(0, nrow = 3*N_rows, ncol = 3*N_rows)
  for(i in 1:N_rows){
    Sigma_etay[i,i] <- dt$v_p[i]
    Sigma_etay[i,i+N_rows] <- dt$v_p_t[i]
    Sigma_etay[i,i+2*N_rows] <- dt$v_p_k[i]
  }
  for(i in 1:N_rows){
    Sigma_etay[i+N_rows,i] <- dt$v_p_t[i]
    Sigma_etay[i+N_rows,i+N_rows] <- dt$v_t[i]
    Sigma_etay[i+N_rows,i+2*N_rows] <- dt$v_t_k[i]
  }
  for(i in 1:N_rows){
    Sigma_etay[i+2*N_rows,i] <- dt$v_p_k[i]
    Sigma_etay[i+2*N_rows,i+N_rows] <- dt$v_t_k[i]
    Sigma_etay[i+2*N_rows,i+2*N_rows] <- dt$v_k[i]
  }
  return(Matrix(Sigma_etay))
}
get_mode_mean_sd_for_hyperparameter <- function(mdl, fun, hyperparameter){
  x0 <- mdl$marginals.hyperpar[[hyperparameter]]
  E <- inla.emarginal(function(x) c(fun(x), fun(x)^2), x0)
  sd <- sqrt(E[2] - E[1]^2)
  mean <- E[1]
  
  mt_tmp <- inla.tmarginal(fun, x0)
  dat_tmp <- inla.smarginal(mt_tmp)
  mode <- dat_tmp$x[which(dat_tmp$y == max(dat_tmp$y))]
  return(list(sd = sd, mean = mean, mode = mode))
}
get_mode_mean_sd_for_hyperparameter <- function(mdl, fun, hyperparameter, calc_mode=T){
  x0 <- mdl$marginals.hyperpar[[hyperparameter]]
  E <- inla.emarginal(function(x) c(fun(x), fun(x)^2), x0)
  sd <- sqrt(E[2] - E[1]^2)
  mean <- E[1]
  
  if (calc_mode){
    mt_tmp <- inla.tmarginal(fun, x0)
    dat_tmp <- inla.smarginal(mt_tmp)
    mode <- dat_tmp$x[which(dat_tmp$y == max(dat_tmp$y))]
  } else{
    mode <- mean
  }
  return(list(sd = sd, mean = mean, mode = mode))
}

log_f_x_given_theta <- function(kappa, Q_x){
  res <- determinant(Q_x, logarithm = T)
  res <- .5*res$modulus[1]*res$sign
  return(res)
}
log_f_x_given_eta_hat_and_theta <- function(kappa, Q_x_given_etaHat, N, Z, eta_hat, Q_etay, K){
  mu_x_given_etaHat <- solve(a = Cholesky(Q_x_given_etaHat),b = b)
  determ <- determinant(Q_x_given_etaHat, logarithm = T)
  res <- -.5*t(mu_x_given_etaHat)%*%(Q_x_given_etaHat%*%mu_x_given_etaHat)
  res <- res + .5*determ$modulus[1]*determ$sign
  return(res)
}
calc_Q_x_given_etaHat <- function(Q_x, N, Q_etay){
  Q_x_given_etaHat <- Q_x
  Q_x_given_etaHat[1:(3*N), 1:(3*N)] <- Q_x[1:(3*N), 1:(3*N)] + Q_etay
  return(Q_x_given_etaHat)
}
calc_Q_x <- function(kappa = k_st, Q_beta_psi, Q_beta_tau, Q_beta_xi, N, Z, d_spde){
  
  Q_u_psi <- makeQ_u(s = exp(kappa$u_psi), rho = exp(kappa$v_psi), d_spde)
  Q_u_tau <- makeQ_u(s = exp(kappa$u_tau), rho = exp(kappa$v_tau), d_spde)
  
  Q_nu <- bdiag(Q_beta_psi,Q_u_psi,Q_beta_tau,Q_u_tau,Q_beta_xi)
  
  Q_epsilon <- bdiag(Diagonal(N,exp(kappa$kappa_psi)), Diagonal(N,exp(kappa$kappa_tau)),
                     Diagonal(N,exp(kappa$kappa_xi)))
  K <- dim(Q_nu)[1]
  Q_x[1:(3*N),  1:(3*N)] <- Q_epsilon
  Q_x[1:K + 3*N,1:(3*N)] <- -t(Z)%*%Q_epsilon
  Q_x[1:(3*N),  1:K + 3*N] <- - (Q_epsilon %*% Z)
  Q_x[1:K + 3*N, 1:K + 3*N] <- Q_nu + t(Z)%*%Q_epsilon%*%Z
  return(Q_x)
}

################################## Prior matrices for covariate coefficients ###################################

betaVarPrior <- 10000
K <- 2*N_colsA + N_psi +N_tau +N_xi
Q_beta_psi <- Diagonal(N_psi, 1/betaVarPrior)
Q_beta_tau <- Diagonal(N_tau, 1/betaVarPrior)
Q_beta_xi <- Diagonal(N_xi, 1/betaVarPrior)

################# get variance and mode for hyper parameters for proposal distribution ###############################

tau_psi_0 <- get_mode_mean_sd_for_hyperparameter(mdl_sep_psi$mdl, function(x) log(x), "Precision for idx")
tau_tau_0 <- get_mode_mean_sd_for_hyperparameter(mdl_sep_tau$mdl, function(x) log(x), "Precision for idx", F)
tau_xi_0 <- get_mode_mean_sd_for_hyperparameter(mdl_sep_xi$mdl, function(x) log(x), "Precision for idx")
rho_psi_0 <- get_mode_mean_sd_for_hyperparameter(mdl_sep_psi$mdl, function(x) log(x), "Range for s")
rho_tau_0 <- get_mode_mean_sd_for_hyperparameter(mdl_sep_tau$mdl, function(x) log(x), "Range for s")
sigma_psi_0 <- get_mode_mean_sd_for_hyperparameter(mdl_sep_psi$mdl, function(x) log(x), "Stdev for s")
sigma_tau_0 <- get_mode_mean_sd_for_hyperparameter(mdl_sep_tau$mdl, function(x) log(x), "Stdev for s")

# Number of stations
N <- nrow(dt_mles)
eta_hat <- matrix(c(dt_mles$psi, dt_mles$tau, dt_mles$kappa), nrow = 3 * N)

###### Make sigma_eta_y
sigma_eta_y <- makeSigma_etay(dt = dt_mles)
Q_etay <- solve(sigma_eta_y)

# Mean and covariance matrix of the proposal distribution in the first loop
kappa_0 <- c(tau_psi_0$mode, rho_psi_0$mode, sigma_psi_0$mode,
             tau_tau_0$mode, rho_tau_0$mode, sigma_tau_0$mode,
             tau_xi_0$mode) %>% as.matrix()

# This is another option for a covariance matrix for the proposal distribution
Sigma_kappa_0 <- diag(c(tau_psi_0$sd^2, rho_psi_0$sd^2, sigma_psi_0$sd^2,
                        tau_tau_0$sd^2, rho_tau_0$sd^2, sigma_tau_0$sd^2,
                        tau_xi_0$sd^2)*0.5
)
#################################### Model fitting - Hyperparameters (One chain) ##################
# Could also run more chains here parallel for example using mclapply from the parallel package
N_samples <- 1000 # Change this number to take more samples
Q_x <- Matrix(0,nrow = 3*N+K, ncol = 3*N+K)
B <- Matrix(0,nrow = 3*N, ncol = 3*N+K)
B[1:(3*N),1:(3*N)] <- diag(1,3*N)
b <- t(B)%*%(Q_etay%*%eta_hat)
kappa_k <- mvrnorm(mu = kappa_0, Sigma = Sigma_kappa_0)
kappa_mat <- matrix(NA, nrow = length(kappa_k), ncol = N_samples)
kappa_k_is_kappa_star <- F
accept <- 0
for(i in 1:N_samples){
  kappa_star <- mvrnorm(mu = kappa_k, Sigma = Sigma_kappa_0)
  k_st <- list(
    kappa_psi = kappa_star[1],
    v_psi = kappa_star[2],
    u_psi = kappa_star[3],
    kappa_tau = kappa_star[4],
    v_tau = kappa_star[5],
    u_tau = kappa_star[6],
    kappa_xi = kappa_star[7]
  )
  k_k <- list(
    kappa_psi = kappa_k[1],
    v_psi = kappa_k[2],
    u_psi = kappa_k[3],
    kappa_tau = kappa_k[4],
    v_tau = kappa_k[5],
    u_tau = kappa_k[6],
    kappa_xi = kappa_k[7]
  )
  if(i == 1){
    Q_x <- calc_Q_x(kappa = k_st, Q_beta_psi, Q_beta_tau, Q_beta_xi, N, Z, d_spde)
    Q_x_given_etaHat <- calc_Q_x_given_etaHat(Q_x = Q_x, N = N, Q_etay = Q_etay)
    l_x_etaHat_k_st <- log_f_x_given_eta_hat_and_theta(kappa = k_st, Q_x_given_etaHat = Q_x_given_etaHat,
                                                       N = N, Z = Z, eta_hat = eta_hat, Q_etay = Q_etay,K = K)
    l_x_k_st <- log_f_x_given_theta(kappa = k_st, Q_x)
    l_prior_k_st <- log_prior_logTheta_noGamma(k_st)
    
    Q_x <- calc_Q_x(kappa = k_k, Q_beta_psi, Q_beta_tau, Q_beta_xi, N, Z, d_spde)
    Q_x_given_etaHat <- calc_Q_x_given_etaHat(Q_x = Q_x, N = N, Q_etay = Q_etay)
    
    l_x_etaHat_k_k <- log_f_x_given_eta_hat_and_theta(kappa = k_k, Q_x_given_etaHat = Q_x_given_etaHat,
                                                      N = N, Z = Z, eta_hat = eta_hat, Q_etay = Q_etay,K = K)
    l_x_k_k <- log_f_x_given_theta(kappa = k_k, Q_x)
    l_prior_k_k <- log_prior_logTheta_noGamma(kappa = k_k)
  }else if(kappa_k_is_kappa_star){
    Q_x <- calc_Q_x(kappa = k_st, Q_beta_psi, Q_beta_tau, Q_beta_xi, N, Z, d_spde)
    Q_x_given_etaHat <- calc_Q_x_given_etaHat(Q_x = Q_x, N = N, Q_etay = Q_etay)
    
    l_x_etaHat_k_k <- l_x_etaHat_k_st
    l_x_k_k <- l_x_k_st
    l_prior_k_k <- l_prior_k_st
    
    l_x_etaHat_k_st <- log_f_x_given_eta_hat_and_theta(kappa = k_st, Q_x_given_etaHat = Q_x_given_etaHat,
                                                       N = N, Z = Z, eta_hat = eta_hat, Q_etay = Q_etay,K = K)
    l_x_k_st <- log_f_x_given_theta(kappa = k_st, Q_x)
    l_prior_k_st <- log_prior_logTheta_noGamma(k_st)
  }else{
    Q_x <- calc_Q_x(kappa = k_st, Q_beta_psi, Q_beta_tau, Q_beta_xi, N, Z, d_spde)
    Q_x_given_etaHat <- calc_Q_x_given_etaHat(Q_x = Q_x, N = N, Q_etay = Q_etay)
    
    l_x_etaHat_k_st <- log_f_x_given_eta_hat_and_theta(kappa = k_st, Q_x_given_etaHat = Q_x_given_etaHat,
                                                       N = N, Z = Z, eta_hat = eta_hat, Q_etay = Q_etay,K = K)
    l_x_k_st <- log_f_x_given_theta(kappa = k_st, Q_x)
    l_prior_k_st <- log_prior_logTheta_noGamma(k_st)
  }
  
  r <- l_prior_k_st-l_prior_k_k+l_x_k_st-l_x_k_k+l_x_etaHat_k_k-l_x_etaHat_k_st
  
  if(as.numeric(r) > log(runif(1,0,1))){
    kappa_k <- kappa_star
    kappa_k_is_kappa_star = T
    accept <- accept + 1
  }else{
    kappa_k <- kappa_k
    kappa_k_is_kappa_star = F
  }
  kappa_mat[,i] <- kappa_k
}
accept/N_samples

# Take first 20% as burn-in.
burnIn <- 0.2
chain <- kappa_mat[,seq(burnIn*N_samples,N_samples,1)]
sigma_psi <- sqrt(1/exp(chain[1,]))
sigma_tau <- sqrt(1/exp(chain[4,]))
sigma_xi <- sqrt(1/exp(chain[7,]))
r_psi <- exp(chain[2,])
s_psi <- exp(chain[3,])
r_tau <- exp(chain[5,])
s_tau <- exp(chain[6,])

#################################### Latent parameter posterior samples ####################################
N_x_loops <- dim(chain)[2]
x_mat <- matrix(NA, nrow = K+3*N, ncol = N_x_loops)
kappa_old <- rep(0,7)
Q_x <- Matrix(0,nrow = 3*N+K, ncol = 3*N+K)
B <- Matrix(0,nrow = 3*N, ncol = 3*N+K)
B[1:(3*N),1:(3*N)] <- Diagonal(3*N, 1)
b <- t(B)%*%(Q_etay%*%eta_hat)
for(i in 1:N_x_loops){
  k <- list(
    kappa_psi = chain[1,i],
    v_psi = chain[2,i],
    u_psi = chain[3,i],
    kappa_tau = chain[4,i],
    v_tau = chain[5,i],
    u_tau = chain[6,i],
    kappa_xi = chain[7,i]
  )
  Q_x <- calc_Q_x(kappa = k, Q_beta_psi = Q_beta_psi, Q_beta_tau = Q_beta_tau,
                  Q_beta_xi = Q_beta_xi, N = N, Z = Z, d_spde = d_spde)
  Q_x_given_etaHat <- calc_Q_x_given_etaHat(Q_x = Q_x,N = N, Q_etay = Q_etay)
  mu_x_given_etaHat <- solve(a = Cholesky(Q_x_given_etaHat),b = b)
  x <- inla.qsample(n = 1, mu = mu_x_given_etaHat, Q = Q_x_given_etaHat)
  x_mat[,i] <- x
}
ms_total_time <- difftime(Sys.time(), ms_start_t) # END time

############################## Posterior samples for psi, tau, kappa and gamma ########
eta <- x_mat[1:(3*N),]
psi_sam <- eta[1:N,]
tau_sam <- eta[1:N + N,]
phi_sam <- eta[1:N + 2*N,]
