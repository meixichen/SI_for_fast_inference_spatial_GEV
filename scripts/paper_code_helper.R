#' Simulate the small-scale data in the paper.
#'
#' @return Output to the global environment the following quantities:
#' simulated quantities including coordinates, parameter values,
#' observed extreme values, and 10% return levels.
simulate_data_small <- function() {
  # 1. Simulate location coordinates
  lon <- seq(0, 10, length.out = 20)
  lat <- seq(0, 10, length.out = 20)
  locs <- expand.grid(x = lon, y = lat)
  n_loc <- nrow(locs)

  # 2. Simulate log(b) from the log density surface of bivariate normal mixture
  G <- 2 # number of components
  pro <- c(0.4,0.6) # mixing proportion
  mu1 <- c(1, 0) # mean of the first component
  mu2 <- c(8, 7) # mean of the second component
  mean <- cbind(mu1, mu2)
  diag_entry <- c(0.5, 1)
  Sig1 <- diag(diag_entry[1], 2) # covariance matrix of the first component
  Sig2 <- diag(diag_entry[2], 2) # covariance matrix of the second component
  variance <- list(d=2, G=G, sigma=array(Sig1, Sig2), sigmasq=diag_entry, scale=diag_entry)
  parameters <- list(pro=pro, mean=mean, variance=variance)
  logb <- mclust::dens("VII", data = locs, logarithm = TRUE, parameters = parameters)
  logb <- (logb +40)/12# rescale log(b)
  b <- exp(logb)
  logb_mat <- matrix(logb, ncol=sqrt(n_loc))

  # 3. Simulate a from the log density surface of bivariate normal
  mu_a <- c(3,3)
  Sig_a <- diag(3,2)
  a <- mvtnorm::dmvnorm(x=locs, mean=mu_a, sigma=Sig_a, log=TRUE)
  a <- (a + 80)
  a_mat <- matrix(a, ncol=sqrt(n_loc))

  # 4. Simulate s
  mu_s <- c(4,4)
  Sig_s <- diag(2,2)
  logs <- mvtnorm::dmvnorm(x=locs, mean=mu_s, sigma=Sig_s, log=TRUE)
  logs <- (logs)/5
  s <- exp(logs)
  logs_mat <- matrix(logs, ncol=sqrt(n_loc))

  # 5. Simulate observed data
  data <- Map(rgev, n=sample(20:50, n_loc, replace=T),
              loc=a, scale=exp(logb), shape=exp(logs))

  # 6. 90% quantile of data at each location
  z_true <- unlist(Map(qgev, p=0.1, loc=a, scale=b,
                       shape=s, lower.tail=F))

  list(data = data, lon = lon, lat = lat, locs = locs,
       a = a, a_mat = a_mat,
       logb = logb, logb_mat = logb_mat,
       logs = logs, logs_mat = logs_mat,
       z_true = z_true)
}


#' Simulate the large-scale data in the paper.
#'
#' @return Output to the global environment the following quantities:
#' simulated quantities including coordinates, parameter values,
#' observed extreme values, and 10% return levels.
simulate_data_big <- function() {
  x <- y <- seq(0, 60, length=80)
  x <- x + rnorm(length(x), 0, 0.001)
  y <- y + rnorm(length(y), 0, 0.001)
  lon <- x
  lat <- y
  coor <- cbind(x, y)
  locs <- expand.grid(x,y)
  n_loc <- length(x)*length(y)
  # Simulate a
  ## nugget is a constant offset, sill is sigma, range is 2*ell
  a <- rgp(1, coor, cov.mod="powexp",
           mean=70, nugget=0, sill=8, range=20, smooth=1,
           grid=TRUE)
  a_mat <- matrix(a, ncol=sqrt(n_loc))

  # Simulate b
  logb <- a/10-4 +
    rgp(1, coor, cov.mod="powexp", mean=0, sill=0.01, range=20, smooth=1, grid=T)
  b <- exp(logb)
  logb_mat <- matrix(logb, ncol=sqrt(n_loc))

  # Simulate s
  logs <- -a/10+5  +
    rgp(1, coor, cov.mod="powexp",
        mean=0, nugget=0, sill=0.01, range=30, smooth=1,
        grid=TRUE)
  s <- exp(logs)
  logs_mat <- matrix(logs, ncol=sqrt(n_loc))

  # Simulate data
  a <- as.vector(a)
  logb <- as.vector(logb)
  b <- as.vector(b)
  logs <- as.vector(logs)
  s <- as.vector(s)
  data <- Map(rgev, n=sample(50, n_loc, replace=TRUE),
              loc=a, scale=b, shape=s)

  # Quantile
  z_true <- unlist(Map(qgev, p=0.1, loc=a, scale=b,
                       shape=s, lower.tail=F))

  list(data = data, lon = lon, lat = lat, locs = locs,
       a = a, a_mat = a_mat,
       logb = logb, logb_mat = logb_mat,
       logs = logs, logs_mat = logs_mat,
       z_true = z_true)
}


#' Get indices of quantities of interest from a matrix of samples.
#'
#' @param chain A matrix of samples with each column containing the samples
#' of a parameter.
#' @param names A vector of character strings, e.g, `c('a', 'b', 's')`, whose
#' indices in the MCMC chain we would like to obtain.
#' @param meshidxloc Optional vector of indices of the location of interest on the SPDE
#' mesh.
#' @return A list of vectors containing the indices of names in `names`.
get_ind_from_samples <- function(chain, names, meshidxloc=NULL) {
  chain_colnames <- colnames(chain)
  n_name <- length(names)
  out <- list()
  for (name in names) {
    name_ind <- which(chain_colnames == name)
    if (!is.null(meshidxloc)) {
      name_ind <- name_ind[meshidxloc]
    }
    out[[name]] <- name_ind
  }
  out
}


#' Get posterior mean of 10% return level.
#'
#' @param p_draws A matrix of GEV parameter samples.
#' @param u1_ind Vector of indices of a in the sample matrix.
#' @param u2_ind Same as above.
#' @param u3_ind Same as above.
#' @param maxsmooth If `TRUE`, use the parameterization in max-and-smooth.
#'
#' @return A vector of posterior draws of the 10% return levels at different
#' locations.
get_posterior_z10 <- function(p_draws, u1_ind, u2_ind, u3_ind,
                              maxsmooth=FALSE) {
  if(maxsmooth) {
    apply(p_draws, 1,
          function(all_draw) {
            psi <- all_draw[u1_ind]
            tau <- all_draw[u2_ind]
            phi <- all_draw[u3_ind]
            mapply(qgev, p=0.1, loc=exp(psi), scale=exp(psi+tau),
                   shape=phi2s(phi), lower.tail=F)
          })
  } else {
    apply(p_draws, 1,
          function(all_draw) {
            mapply(qgev, p=0.1, loc=all_draw[u1_ind],
                   scale=exp(all_draw[u2_ind]),
                   shape=exp(all_draw[u3_ind]), lower.tail=FALSE)
          })
  }
}

#' Get posterior mean and sd of the quantities of interest.
#'
#' @param fit A `spatialGEV_fit` object, as computed by [SpatialGEV::spatialGEV_fit()].
#' @param sam A `spatialGEV_sam` object, as computed by [SpatialGEV::spatialGEV_sample()].
#' @return A list of posterior means for a, b, s, and z10, as well as the
#' posterior samples of them.
spatialGEV_get_posteriors <- function(fit, sam) {
  meshidxloc <- fit$meshidxloc # extract location indices
  n_loc <- length(meshidxloc)

  # Mode estimates of the random effects
  # These coincide with the posterior means due to Normal approximation
  # a_mean <- fit$rep$par.random[meshidxloc]
  # logb_mean <- fit$rep$par.random[length(fit$rep$par.random)/3+meshidxloc]
  # logs_mean <- fit$rep$par.random[length(fit$rep$par.random)/3*2+meshidxloc]


  # Posterior draws of the random effects
  p_draws <- sam$parameter_draws

  if ("a"%in%fit$random) {
    a_labels <- paste0("a", 1:n_loc)
    a_ind <- which(colnames(p_draws)%in% a_labels) # a index
    a_mean <- apply(p_draws[,a_ind], 2, mean)
    a_sd <- apply(p_draws[,a_ind], 2, sd)
  } else {
    a_labels <- "a"
    a_ind <- which(colnames(p_draws)%in% a_labels)
    a_mean <- mean(p_draws[,a_ind])
    a_sd <- sd(p_draws[,a_ind])
  }
  if ("log_b"%in%fit$random) {
    logb_labels <- paste0("log_b", 1:n_loc)
    logb_ind <- which(colnames(p_draws)%in% logb_labels) # b index
    logb_mean <- apply(p_draws[,logb_ind], 2, mean)
    logb_sd <- apply(p_draws[,logb_ind], 2, sd)
  } else {
    logb_labels <- "log_b"
    logb_ind <- which(colnames(p_draws)%in% logb_labels)
    logb_mean <- mean(p_draws[,logb_ind])
    logb_sd <- sd(p_draws[,logb_ind])
  }
  if ("s"%in%fit$random) {
    s_labels <- paste0("s", 1:n_loc)
    s_ind <- which(colnames(p_draws)%in% s_labels) # s index
    s_mean <- apply(p_draws[,s_ind], 2, mean)
    s_sd <- apply(p_draws[,s_ind], 2, sd)
  } else {
    s_labels <- "s"
    s_ind <- which(colnames(p_draws)%in% s_labels)
    s_mean <- mean(p_draws[,s_ind])
    s_sd <- sd(p_draws[,s_ind])
  }

  # Get posterior draws of z_10 by transforming draws of a, b and s
  z_draws <- get_posterior_z10(p_draws, a_ind, logb_ind, s_ind)

  # Posterior mean an SD of z_10
  z_mean <- apply(z_draws, 1, mean)
  z_sd <- apply(z_draws, 1, sd)

  list(a_mean=a_mean, logb_mean=logb_mean, s_mean=s_mean, z_mean=z_mean,
       a_sd=a_sd, logb_sd=logb_sd, s_sd=s_sd, z_sd=z_sd,
       p_draws=p_draws, z_draws=z_draws)
}

#' Get posterior samples from the TMB report.
#'
#' @param report TMB report of the model given by [TMB::sdreport()].
#' @param meshidxloc Indices of the locations of interest on the mesh.
#' @param n_draw Number of posterior samples.
#'
#' @return Matrix of posterior draws of all quantities of interest.
get_posterior_from_tmb <- function(
  report,
  meshidxloc,
  n_draw=10000,
  u1_label="a",
  u2_label="logb",
  u3_label="logs",
  hyperparam_labels= c("beta_a","beta_b","beta_s",
                       "log_sigma_a","log_kappa_a",
                       "log_sigma_b","log_kappa_b",
                       "log_sigma_s","log_kappa_s")
) {
  mean_random <- report$par.random
  mean_fixed <- report$par.fixed
  jointPrec_mat <- report$jointPrecision
  C <- chol(jointPrec_mat)
  joint_cov <- backsolve(r = C,
                         x = backsolve(r = C, x = diag(nrow(jointPrec_mat)),
                                       transpose = TRUE)) # Cholesky decomp
  random_ind <- report$env$random #indices of random effects
  ind2rm <- setdiff(1:length(random_ind),
                    c(meshidxloc, # indices of a in the random effects vec
                      length(random_ind)/3 + meshidxloc, # indices of log_b in the random vec
                      length(random_ind)/3*2 + meshidxloc)) # indices of s in  the random vec
  mean_random <- mean_random[-ind2rm]
  joint_cov <- joint_cov[-ind2rm, -ind2rm]
  n_loc <- length(meshidxloc)
  par_names_random <- paste0(names(mean_random), 1:n_loc) # add location index
  par_names_fixed <- names(mean_fixed) # extract parameter names for the fixed effects
  joint_mean <- c(mean_random, mean_fixed)
  names(joint_mean) <- c(par_names_random, par_names_fixed) # modify parameter names
  fixed_ind <- (length(mean_random)+1):length(joint_mean) # indices of the fixed effects
  joint_post_draw <- mvtnorm::rmvnorm(n_draw, joint_mean, joint_cov)

  u1_ind <- paste0(u1_label, meshidxloc)
  u2_ind <- paste0(u2_label, meshidxloc)
  u3_ind <- paste0(u3_label, meshidxloc)
  colnames(joint_post_draw) <- c(u1_ind, u2_ind, u3_ind, hyperparam_labels)
  # Get return level posterior samples
  z_draws <- get_posterior_z10(joint_post_draw, u1_ind, u2_ind, u3_ind,
                               maxsmooth=(u1_label=="psi"))
  # Get posterior means
  u1 <- apply(joint_post_draw[,u1_ind], 2, mean)
  u2 <- apply(joint_post_draw[,u2_ind], 2, mean)
  u3 <- apply(joint_post_draw[,u3_ind], 2, mean)
  z_s <- apply(z_draws, 1, mean)

  list(joint_post_draw=joint_post_draw, z_draws=z_draws,
       u1=u1, u2=u2, u3=u3, z_s=z_s)
}

#' Process the MCMC samples from different chains.
#'
#' @param chains A list of lists. Each list contains `random`, `fixed`, `accept`, and `time`.
#' @param loc_ind Optional vector of indices of the location of interest on the SPDE mesh.
#' @return A list of z10 draws `z_draws` and posterior means of the quantities of interest `pos_means` as a `n_loc x 4` matrix.
get_random_posteriors_mcmc <- function(chains, loc_ind) {
  all_ind <- get_ind_from_samples(chains[[1]]$random,
                                  c("a", "log_b", "s"), loc_ind)
  lapply(chains,
         function(chain) {
           random_sam <- chain$random
           z_draws <- get_posterior_z10(random_sam,
                                        all_ind$a,
                                        all_ind$log_b,
                                        all_ind$s)
           z_mean <- apply(z_draws, 1, mean)
           a_mean <- apply(random_sam[, all_ind$a], 2, mean)
           logb_mean <- apply(random_sam[, all_ind$log_b], 2, mean)
           s_mean <- apply(random_sam[, all_ind$s], 2, mean)
           list(z_draws=z_draws,
                pos_means=cbind(a_mean, logb_mean,
                                s_mean, z_mean))
         })
}

#' Function to transform `theta = (psi, tau, phi)` to `(a, logb, s)`.
#'
#' @param theta A length 3 vector.
#' @return A length 3 vector on (a, logb, logs) scale.
j222paper_transform <- function(theta) {
  list2env(j22_constants(), envir = environment())
  a <- exp(theta[1])
  logb <- theta[1]+theta[2]
  s <- (1-(exp(-exp((theta[3] - b3)/sig3))))^(1/alp3) - 0.5
  logs <- log(s)
  c(a, logb, logs)
}

#' Function to transform `(a, logb, s)` to `(psi, tau, phi)`.
#' @param a
#' @param logb
#' @param s
#' @return A list of j22 parameters (psi, tau, phi)
paper2j22_transform <- function(a, logb, s) {
  list2env(j22_constants(), envir = environment())
  n_loc <- length(a)
  psi <- log(a)
  tau <- logb - psi
  phi <- b3 + sig3 * log(-log(1-(s+0.5)^alp3))
  list(psi=psi, tau=tau, phi=phi)
}

#' Transform `phi` in J22 paper to `s` in our paper.
#'
#' @param phi Numeric
#'
#' @return Shape parameter in the original scale.
phi2s <- function(phi) { # function to transform phi to shape
  list2env(j22_constants(), envir = environment())
  (1. - exp(-exp( (phi - b3) / sig3)))^(1./alp3) - .5
}

#' Calculate the constants for J22 GEV parameter basis.
#'
#' @return A list with elements `xi00`, `alpha`, `beta`, `alph3`, `b3`, `sig3`.
j22_constants <- function() {
  xi00 <- 0
  alp3 <- 0.8 # c_phi
  alpha <- 4
  beta <- 4
  sig3 <- (log(1 - (xi00+0.5)^alp3))*(1 - (xi00+0.5)^alp3)*
    (-(alp3)^(-1)*(xi00+0.5)^(-alp3+1)) # b_phi
  b3 <- -sig3*(log(-log(1-0.5^alp3))) # a_phi
  list(xi00 = xi00, alpha = alpha, beta = beta,
       alp3 = alp3, b3 = b3, sig3 = sig3)
}

#' Transform parameters `(a, b, s)` to `(a, log_b, log_s)`
#'
#' @param params Parameters `(a, b, s)`.
#' @return Parameters `(a, log_b, log_s)`.
evd2paper_transform <- function(params) {
  theta <- c(params[1], log(params[2]), log(params[3]))
  names(theta) <- c("a", "log_b", "log_s")
  theta
}

#' Jacobian of [evd2paper_transform()].
#'
#' @param params Parameters `(a, b, s)`.
#' @return Jacobian of `(a, log_b, log_s)` wrt `params`.
evd2paper_jacobian <- function(params) {
  diag(c(params[1], 1/params[2], 1/params[3]))
}


#' Max step using j22 method.
#'
#' @param y Vector of GEV observations.
#' @return List with elements `est` and `var`.
#'
#' @details `est` and `var` are in the j22 basis `(psi, tau, phi)`.
j22_maxstep <- function(y) {
  list2env(j22_constants(), envir = environment())
  # Likelihood function to minimize
  fn <- function(theta,y) {
    list2env(j22_constants(), envir = environment())
    xitheta <- (1 - exp(-exp((theta[3]-b3)/sig3)))^(1/alp3) - 0.5
    ans <- sum(evd::dgev(y, loc = exp(theta[1]), scale = exp(theta[2] + theta[1]),
                         shape = ((1-(exp(-exp((theta[3]-b3)/sig3))))^(1/alp3) - 0.5),
                         log = TRUE))
    ans <- ans + (alpha - alp3)*log(xitheta + 0.5) +
      (beta-1)*log(0.5 - xitheta) +
      (theta[3]-b3)/sig3 - exp((theta[3]-b3)/sig3)
    -ans
  }
  fitgev <- evd::fgev(y)
  # parameter transformations
  mu0 <- fitgev$estimate[1]
  sigma0 <- fitgev$estimate[2]
  xi0 <- fitgev$estimate[3]
  if (xi0 > 0) {
    xi0 <- min(xi0,0.45)
  } else {
    xi0 <- max(xi0,-0.45)
  }
  theta0 <- c(log(mu0),log(sigma0)-log(mu0),b3 + sig3*log(-log(1 - (xi0+0.5)^alp3)))
  GEV_fit <- nlm(fn, theta0, y = y, hessian = TRUE)
  theta <- GEV_fit$estimate
  S_d <- try(solve(GEV_fit$hessian),silent = TRUE)
  list(est = theta, var = S_d)
}

#' Max step using tmb method.
#'
#' @param y Vector of GEV observations.
#' @param s_prior The mean and standard deviation of the normal prior on `log_s`.
#' @return List with elements `est` and `var`.
#'
#' @details `est` and `var` are in the paper basis `(a, log_b, log_s)`.
tmb_maxstep <- function(y, s_prior) {
  # initialize optimization
  fitgev <- evd::fgev(yi)
  theta_init <- fitgev$estimate
  if(theta_init[3] <= 0) {
    # negative shape parameter:
    # pick smallest value such that support includes minimum value
    shape <- theta_init[2] / (theta_init[1] - .99 * min(yi))
    theta_init[3] <- shape
  }
  theta_init <- evd2paper_transform(theta_init)
  gev_adf <- TMB::MakeADFun(
    data = list(model = "model_gev",
                y = yi, reparam_s = 1,
                s_prior = s_prior),
    parameters = list(a = 0, log_b = 0, s = 0),
    DLL = "SpatialGEV_TMBExports",
    silent = TRUE
  )
  gev_fit <- tryCatch(
    # use evd initial value
    nlminb(
      ## start = ifelse(is.na(theta_init), 1, theta_init),
      start = theta_init,
      objective = gev_adf$fn,
      gradient = gev_adf$g,
      trace = 1
    ), error = function(e) {
      # start from (0, 0, 0)
      nlminb(
        start = c(0, 0, 0),
        objective = gev_adf$fn,
        gradient = gev_adf$g
      )
    })
  ## ev <- eigen(gev_adf$he(gev_fit$par))
  list(est = gev_fit$par,
       var = try(solve(gev_adf$he(gev_fit$par)), silent = TRUE))
}

#' Max step using evd method.
#'
#' @param y Vector of GEV observations.
#' @return List with elements `est` and `var`.
#'
#' @details `est` and `var` are in the paper basis `(a, log_b, log_s)`.
evd_maxstep <- function(y) {
  fitgev <- evd::fgev(yi)
  theta_est <- fitgev$estimate
  var_est <- fitgev$var.cov
  # convert to paper basis using delta method
  J <- evd2paper_jacobian(theta_est)
  var_est <- t(J) %*% var_est %*% J
  theta_est <- evd2paper_transform(theta_est)
  list(est = theta_est, var = var_est)
}

#' Max step for the GEV model.
#'
#' @param y Vector of GEV observations.
#' @param method Method used to perform optimization.  One of "tmb", "evd", or "j22".
#' @param s_prior For `method == "tmb"`, the mean and standard deviation of the normal prior on `log_s`.
#' @param basis Parameterization basis.  Either "paper" or "j22".
#'
#' @return List with elements `est` and `var` with dimension names `c("a, "log_b", "log_s")` if `basis == "paper"` and `c("psi", "tau", "phi")` if `basis == "j22"`.
#'
#' @details Minor modification to the code provided in Johanesson et al 2022.
maxsmooth_maxstep <- function(y,
                              method = c("tmb", "evd", "j22"),
                              s_prior =  c(0, 10),
                              basis = c("paper", "j22")) {
  basis <- match.arg(basis)
  method <- match.arg(method)
  out <- switch(method,
                tmb = tmb_maxstep(y, s_prior),
                evd = evd_maxstep(y),
                j22 = j22_maxstep(y))
  if(basis == "j22") {
    theta_names <- c("psi", "tau", "phi")
    if(method %in% c("tmb", "evd")) {
      # delta method
      J <- numDeriv::jacobian(paper2j22_transform,
                              out$est)
      out$var <- t(J) %*% out$var %*% J
      out$est <- paper2j22_transform(out$est)
    }
  } else if(basis == "paper") {
    theta_names <- c("a", "log_b", "log_s")
    if(method == "j22") {
      # delta method
      J <- numDeriv::jacobian(j222paper_transform,
                              out$est)
      out$var <- t(J) %*% out$var %*% J
      out$est <- j222paper_transform(out$est)
    }
  }
  names(out$est) <- theta_names
  colnames(out$var) <- theta_names
  rownames(out$var) <- theta_names
  out
}

#' Compare different methods for max step.
#'
#' @param x,y Names of method on axes.
#' @param type Which estimator to plot.
max_step_compare <- function(mle_set, var_set,
                             x, y, type = c("mle", "se")) {
  type <- match.arg(type)
  # create tibble for plotting
  plt <- lapply(1:length(mle_set), function(i) {
    method <- names(mle_set[i])
    mle <- mle_set[[i]]
    se <- t(apply(var_set[[i]], 3, diag))
    colnames(se) <- colnames(mle)
    bind_rows(list(
      as_tibble(mle) %>%
      mutate(method = !!method,
             loc_id = 1:n(),
             estimate = "mle",
             .before = a),
      as_tibble(se) %>%
      mutate(method = !!method,
             loc_id = 1:n(),
             estimate = "se",
             .before = a)
    ))
  }) %>% bind_rows()
  # plot itself
  plt %>%
    filter(estimate == type) %>%
    pivot_longer(cols = a:log_s,
                 names_to = "parameter",
                 values_to = "value") %>%
    pivot_wider(names_from = method, values_from = "value") %>%
    ggplot(aes(x = !!sym(x), y = !!sym(y))) +
    geom_point(aes(color = parameter)) +
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~ parameter, scales = "free")
}


#' Conditional posterior sampling from `p(u | theta, y)`.
#'
#' @param label "a" or "b" or "s".
#' @param dd Distance matrix.
#' @param obs_ind Indices of the label in the observation cov matrix.
#' @param theta_draws A vector of draws of theta in one iteration of sampling.
#' @return A vector of samples of `u`.
cond_post_sample <- function(label, dd, cov_mle, theta_draws) {
  cov_gp <- kernel_matern(dd,
                          exp(theta_draws[grep(paste0("log_sigma_", label),
                                               names(theta_draws))]),
                          exp(theta_draws[grep(paste0("log_kappa_", label),
                                               names(theta_draws))]))
  Sig_inv <- solve(cov_mle+cov_gp)
  sam_mu <- cov_gp%*%Sig_inv%*%dt_mles[-na_ind, label] +
    cov_mle%*%Sig_inv%*%rep(theta_draws[grep(paste0("beta_", label),
                                             names(theta_draws))], nrow(dd))
  sam_Sig <- cov_mle%*%Sig_inv %*% cov_gp
  tryCatch(mvtnorm::rmvnorm(n=1, mean=sam_mu, sigma=sam_Sig),
           error=function(e) NA)
}

#' Plot a quantity on a grid map.
#'
#' @param x Vector of longitude.
#' @param y Vector of latitude.
#' @param z Matrix of the quantity to plot.
#' @param title Character string.
#' @param x_lab Character string for x-axis label.
#' @param y_lab Character string for y-axis label.
#' @param cex Font size.
#'
#' @return A plot object.
grid_plot <- function(x, y, z,
                      title, x_lab="Longitude", y_lab="Latitude", cex=1.2) {
  image.plot(x=x, y=y, z=z,
             xlab=x_lab, ylab=y_lab, main=title,
             cex.lab=cex, cex.axis=cex, axis.args=list(cex.axis=cex))
}

#' Create customized ggplot theme.
#'
#' @param text_size Font size for text.
#'
#' @return ggplot object.
create_ggplot_theme <- function(text_size=13,
                                title_size=13) {
  theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size=text_size),
          plot.title = element_text(size = title_size))
}


#' Plot true vs estimated parameters as a scatterplot.
#'
#' @param true Vector of true values.
#' @param est Vector of estimated values.
#' @param axis_lim Vector of axis range.
#' @param title Character string.
#' @param x_lab X-axis label.
#' @param y_lab Y-axis label.
#' @param theme A ggplot theme to build on.
#'
#' @return A plot object.
true_vs_est_scatterplot <- function(true, est, axis_lim,
                                    title, x_lab, y_lab, theme) {
  ggplot() +
    geom_point(aes(x=true, y=est)) +
    xlim(axis_lim[1], axis_lim[2]) + ylim(axis_lim[1], axis_lim[2]) +
    geom_abline(linetype="dashed", linewidth=1, color="blue") +
    ylab(y_lab) + xlab(x_lab)+
    theme +
    ggtitle(title)
}

#' QQ-plot with diagonal line.
#'
#' @param zscore Vector of z-scores.
#' @param title Character string.
#' @param x_lab X-axis label.
#' @param y_lab Y-axis label.
#' @param theme A ggplot theme to build on.
#'
#' @return A plot object.
qq_plot <- function(zscore,
                    title, x_lab, y_lab, theme) {
  tibble(z = zscore) %>%
    ggplot(aes(sample = z)) +
    stat_qq_line(linewidth = 1, color = "blue") +
    stat_qq() +
    ylab(y_lab) + xlab(x_lab)+
    theme +
    ggtitle(title)
}

#' Plot values on a map.
#'
#' @param value Vector of values to plot.
#' @param zlim Range of the values to plot.
#' @param lon Vector of longitudes.
#' @param lat Vector of latitudes.
#' @param title Character string.
#'
#' @return A plot object.
map_plot <- function(value, zlim,
                     lon=locs$cell_lon, lat=locs$cell_lat, title="") {
  val <- fields::color.scale(value, col=viridisLite::viridis(10),
                             zlim=zlim)
  plot(lon, lat, col=val, axes=FALSE, pch=15,
       cex=1.2, xlab="", ylab="", main=title)
  maps::map("world", "Canada", add=TRUE)
  fields::image.plot(legend.only = TRUE,
                     zlim=zlim, col=viridisLite::viridis(10),
                     horizontal = FALSE, legend.shrink = 0.5,
                     cex=0.7, smallplot=c(0.85, 0.9, 0.5, 0.9))
  axis(1, at=seq(-140, -50, by=10))
  axis(2, at=seq(40, 80, by=10))
}


#' Plot the posterior predictive values against the observed values.
#'
#' @param obs Vector of observed quantiles at each location.
#' @param pred Vector of the posterior predictive quantiles at each location.
#'
#' @return A plot object.
pospred_vs_obs_plot <- function(obs, pred,
                                statistic_label="upper 10% quantile") {
  ggplot() + geom_point(aes(x=obs, y=pred)) +
    geom_abline(linetype="dashed", color="blue") +
    xlim(range(c(obs, pred))) + ylim(range(c(obs, pred))) +
    ylab(paste("Posterior", statistic_label, "of snowfall extremes")) +
    xlab(paste0(toupper(substr(statistic_label, 1, 1)),
                substr(statistic_label, 2, 1000000L),
                " of the observed snowfall extremes"))
}

#' Create a simple table with [kableExtra::kbl()].
#'
#' @param x Data frame for table.
#' @param caption Table caption.
#' @param digits Decimal places.
#'
#' @return The output of the call to [kableExtra::kbl()].
make_table <- function(x, caption, digits=3) {
  kableExtra::kbl(x, booktabs = TRUE, escape = FALSE,
                  caption = caption, digits=digits) %>%
    kable_styling(latex_options = c("hold_position"))
}
