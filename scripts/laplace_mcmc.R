require(TMB)
require(Matrix) # fast algorithms for (sparse) cholesky decompositions
require(rlang) # efficiently work inside environments

#' Get parameters of multivariate normal distribution in a TMB object.
#'
#' @param obj A TMB object created by [TMB:MakeADFun()].
#' @param select One of "all", "fixed", or "random".  See Details.
#' @return A list with two elements, `mean` and either `var` for `select == "fixed" or `prec` otherwise..
#'
#' @details The return value is computed as follows:
#'
#' For `select == "fixed"`:
#' ```
#' mean = obj$env$last.par[obj$env$lfixed()]
#' var = TMB::sdreport(obj)$cov.fixed
#' ```
#' For `select == "random"`:
#' ```
#' mean = obj$env$last.par[obj$env$lrandom()]
#' prec = obj$env$L.created.by.newton
#' ```
#' For `select == "all":
#' ```
#' mean = obj$env$last.par
#' prec = TMB::sdreport(obj, getJointPrecision=TRUE)$jointPrecision
#' ```
get_mvn <- function(obj, select = c("all", "fixed", "random")) {
  select <- match.arg(select)
  if(select == "all") {
    out <- list(
      mean = obj$env$last.par,
      prec = TMB::sdreport(obj, getJointPrecision=TRUE)$jointPr
    )
  } else if(select == "fixed") {
    out <- list(
      mean = rlang::with_env(obj$env, last.par[lfixed()]),
      var = TMB::sdreport(obj)$cov.fixed
    )
  } else if(select == "random") {
    out <- rlang::with_env(obj$env, {
      list(mean = last.par[lrandom()],
           prec = L.created.by.newton)
    })
  }
  out
}

#' Sample from a multivariate normal with sparse precision matrix.
#'
#' @param n Number of random draws.
#' @param mean Mean vector.
#' @param prec Sparse precision matrix, i.e., inheriting from [Matrix::sparseMatrix] or its Cholesky factor, i.e., inheriting from [Matrix::CHMfactor].
#'
#' @return A matrix with `n` rows, each of which is a draw from the corresponding normal distribution.
#'
#' @details If the matrix is provided in precision form, it is converted to Cholesky form using `Matrix::Cholesky(prec, super = TRUE)`.  Once it is of form [Matrix::CHMfactor], this function is essentially copied from local function `rmvnorm()` in function `MC()` defined in [TMB::MakeADFun()].
rmvn_prec <- function(n, mean, prec) {
  d <- ncol(prec) # number of mvn dimensions
  if(!is(prec, "CHMfactor")) {
    prec <- Matrix::Cholesky(prec, super = TRUE)
  }
  u <- matrix(rnorm(d*n),d,n)
  u <- Matrix::solve(prec,u,system="Lt")
  u <- Matrix::solve(prec,u,system="Pt")
  u <- t(as(u, "matrix") + mean)
}

#' Sample from the full Laplace approximation of a posterior distribution using MCMC
#'
#' @param obj Object returned by [TMB::MakeADFun()] which gives the negative Laplace log-posterior of the fixed-effect parameters via `obj$fn()`.
#' @param n_iter Number of MCMC iterations.
#' @param fixed_init Initial value of the fixed-effect parameters.
#' @param prop_sim Function taking arguments `prev` and `obj` which returns a draw `fixed_prop` from the proposal distribution.
#' @param prop_lpdf Either:
#' - A function taking arguments `new`, `prev`, and `obj` which returns the log-pdf of the proposal distribution; or
#' - `NULL`, which means that the proposal distribution is reversible and thus is not calculated in the Metropolis acceptance ratio.
#' @param random_ind Optional logical vector of random-effect parameters to return.  If `FALSE`, no random-effect sampling is done.  If missing, all random-effect parameters.
#' @param print_every Every how many iterations to print interactive output.  Default is `0` which means never.
#'
#' @return A list with elements:
#' \describe{
#'   \item{fixed}{Matrix of fixed-effect parameter draws.}
#'   \item{random}{Matrix of random-effect parameter draws.}
#'   \item{logpost}{Vector of fixed-effect log-posterior values.}
#'   \item{accept}{The acceptance rate of the MCMC.}
#' }
#'
#' @details Would be better to implement this as a class to store precomputed quantities.
#'
#' @note The external value of `obj` is modified as a result of this call.
laplace_mcmc <- function(obj,
                         n_iter,
                         fixed_init,
                         prop_sim,
                         prop_lpdf,
                         random_ind,
                         print_every = 0,
                         debug = FALSE) {
  # problem dimensions
  n_all <- length(obj$env$par)
  n_random <- sum(obj$env$lrandom())
  n_fixed <- n_all - n_random
  if(missing(random_ind)) {
    random_ind <- rep(TRUE, n_random)
  }
  is_reversible <- is.null(prop_lpdf)
  do_random <- any(random_ind)
  n_random <- n_random - sum(!random_ind)
  # one draw from p_Laplace(random | fixed, data).
  # intended to be called right after `obj$fn(fixed)`
  random_sim <- function() {
    mvn_par <- get_mvn(obj, select = "random")
    x <- rmvn_prec(1,
                   mean = mvn_par$mean, prec = mvn_par$prec)
    drop(x)[random_ind]
  }
  # calculate the log of p_Laplace(fixed | data)
  # returns -Inf if TMB returns any errors
  target_lp <- function(fixed) {
    nlp <- tryCatch(
      obj$fn(fixed)[1],
      error = function(e) Inf
    )
    if(is.na(nlp)) nlp <- Inf
    -nlp
  }
  # allocate memory
  fixed_out <- matrix(NA, n_iter, n_fixed)
  colnames(fixed_out) <- rlang::with_env(obj$env, {
    names(par)[lfixed()]
  })
  if(do_random) {
    random_out <- matrix(NA, n_iter, n_random)
    colnames(random_out) <- rlang::with_env(obj$env, {
      names(par)[lrandom()]
    })[random_ind]
  }
  logpost_out <- rep(NA, n_iter)
  n_accept <- 0 # number of accepted draws
  # initialize sampler
  fixed_curr <- fixed_init
  lp_curr <- target_lp(fixed_curr)
  if(!is.finite(lp_curr)) {
    stop("Error in call to `obj$fn(fixed_init)`.")
  }
  if(do_random) random_curr <- random_sim()
  if(is_reversible) lq_prop <- lq_curr <- 0
  if(debug) browser()
  # mcmc loop
  for(iter in 1:n_iter) {
    # generate proposal
    fixed_prop <- prop_sim(prev = fixed_curr, obj = obj)
    lp_prop <- target_lp(fixed_prop)
    if(!is_reversible) {
      lq_prop <- prop_lpdf(new = fixed_prop, prev = fixed_curr,
                           obj = obj)
      lq_curr <- prop_lpdf(new = fixed_curr, prev = fixed_prop,
                           obj = obj)
    }
    # accept/reject proposal
    log_acc <- (lp_prop - lq_prop) - (lp_curr - lq_curr)
    if(log(runif(1)) <= log_acc) {
      n_accept <- n_accept + 1
      fixed_curr <- fixed_prop
      if(do_random) random_curr <- random_sim()
      lp_curr <- lp_prop
      lq_curr <- lq_prop
    }
    # storage
    fixed_out[iter,] <- fixed_curr
    if(do_random) random_out[iter,] <- random_curr
    logpost_out[iter] <- lp_curr
    if((print_every > 0) && (iter %% print_every == 0)) {
      cat("Iteration", iter, "-- acc_rate:", signif(exp(log_acc), 2), "\n")
    }
  }
  list(fixed = fixed_out,
       random = random_out,
       accept = n_accept/n_iter)
}
