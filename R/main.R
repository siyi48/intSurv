#' Generate the matrices of basis functions by the method of sieves for
#' all nuisance functions
#'
#' This function generate polynomial splines that are needed to approximate
#' the nuisance functions \eqn{c_1(X)}, \eqn{c_0(X)}, and \eqn{\lambda(X)} by
#' the method of sieves.
#'
#' @param c1.mat the matrix that contains the continuous covariates which may
#' affect \eqn{c_1(X)}.
#' @param c1.discrete the matrix that contains the discrete covariates which may
#' affect \eqn{c_1(X)}. The default is `NULL`.
#' @param q.c1 the number of basis functions to approximate each continuous
#' covariate in \eqn{c_1(X)}. Here, it corresponds to the degree of freedom to
#' generate the polynomial basis.
#' @param c0.mat the matrix that contains the continuous covariates which may
#' affect \eqn{c_0(X)}.
#' @param c0.discrete the matrix that contains the discrete covariates which may
#' affect \eqn{c_0(X)}. The default is `NULL`.
#' @param q.c0 the number of basis functions to approximate each continuous
#' covariate in \eqn{c_0(X)}. Here, it corresponds to the degree of freedom to
#' generate the polynomial basis.
#' @param lambda.mat the matrix that contains the continuous covariates which
#' may affect \eqn{\lambda(X)}.
#' @param lambda.discrete the matrix that contains the discrete covariates
#' which may affect \eqn{\lambda(X)}.
#' @param q.lambda the number of basis functions to approximate each continuous
#' covariate in \eqn{\lambda(X)}. Here, it corresponds to the degree of freedom
#' to generate the polynomial basis.
#' @return a list of covariate matrices for the nuisance functions, including
#' \itemize{
#' \item `c1.sieve`: the covariate matrix for \eqn{c_1(X)} in the sieve
#' approximation
#' \item `c0.sieve`: the covariate matrix for \eqn{c_0(X)} in the sieve
#' approximation
#' \item `lambda.sieve`: the covariate matrix for \eqn{\lambda(X)} in the sieve
#' approximation
#' }
#' @export
sieve.mat <- function(c1.mat, c1.discrete = NULL, q.c1,
                      c0.mat, c0.discrete = NULL, q.c0,
                      lambda.mat, lambda.discrete = NULL, q.lambda){
  # create basis function for the nuisance function
  # (1) c1.sieve
  c1.ncol <- ncol(c1.mat)
  c1.sieve <- NULL
  for(k in 1:c1.ncol){
    basis.mat <- poly(c1.mat[,k], degree = q.c1)
    colnames(basis.mat) <- paste0("s1x", k, 1:q.c1)
    c1.sieve <- cbind(c1.sieve, basis.mat)
  }
  colnames(c1.discrete) <- paste0("s1x", (c1.ncol+1):(c1.ncol+ncol(c1.discrete)))
  c1.sieve <- cbind(c1.sieve, c1.discrete)

  # (2) c0.sieve
  c0.ncol <- ncol(c0.mat)
  c0.sieve <- NULL
  for(k in 1:c0.ncol){
    basis.mat <- poly(c0.mat[,k], degree = q.c0)
    colnames(basis.mat) <- paste0("s0x", k, 1:q.c0)
    c0.sieve <- cbind(c0.sieve, basis.mat)
  }
  colnames(c0.discrete) <- paste0("s0x", (c0.ncol+1):(c0.ncol+ncol(c0.discrete)))
  c0.sieve <- cbind(c0.sieve, c0.discrete)

  # (3) lambda.sieve
  lambda.ncol <- ncol(lambda.mat)
  lambda.sieve <- NULL
  for(k in 1:lambda.ncol){
    basis.mat <- poly(lambda.mat[,k], degree = q.lambda)
    colnames(basis.mat) <- paste0("s0ax", k, 1:q.lambda)
    lambda.sieve <- cbind(lambda.sieve, basis.mat)
  }
  colnames(lambda.discrete) <- paste0("s0ax",
                                      (lambda.ncol+1):(lambda.ncol+ncol(lambda.discrete)))
  lambda.sieve <- cbind(lambda.sieve, lambda.discrete)

  return(list(
    c1.sieve = c1.sieve,
    c0.sieve = c0.sieve,
    lambda.sieve = lambda.sieve
  ))
}

#' Integrative estimation of the heterogeneous treatment effect based on the
#' covariates from the basis approximation
#'
#' This function obtains the integrative estimator with its variance estimates
#' by minimizing the penalized log partial likelihood. The penalty function is
#' selected as the adaptive lasso.
#'
#' @param tau.mat the matrix that contains the covariate terms of \eqn{\tau(X)}.
#' @param c1.sieve the covariate matrix for \eqn{c_1(X)} in the sieve
#' approximation.
#' @param c0.sieve the covariate matrix for \eqn{c_0(X)} in the sieve
#' approximation.
#' @param lambda.sieve the covariate matrix for \eqn{\lambda(X)} in the sieve
#' approximation.
#' @param a the binary treatment assignment, where the active treatment group
#' should be encoded as `1`, and the control group should be encoded as `0`.
#' @param s the binary data source indicator, where the trial data should be
#' encoded as `1`, and the observational study should be encoded as `0`.
#' @param t the observed event time.
#' @param delta the event indicator.
#' @param nfolds the fold of cross-validation, the default value is `5`.
#' @return a list of estimators, including
#' \itemize{
#' \item `betahat`: the estimated treatment effect parameter
#' \item `ve.beta`: the variance estimate of the estimated treatment effect
#' \item `var.betahat`: the covariance matrix of the estimated treatment effect
#' parameter
#' \item `ate`: the estimated average treatment effect
#' \item `ve.ate`: the variance estimate of the average treatment effect
#' }
#' @import glmnet survival
#' @export
fitcox.int <- function(tau.mat, c1.sieve, c0.sieve, lambda.sieve,
                       a, s, t, delta, nfolds = 5){
  n <- nrow(tau.mat)
  c1.fitmat <- c1.sieve*s
  c0.fitmat <- c0.sieve*(1-s)
  lambda.fitmat <- lambda.sieve*a*(1-s)
  colnames(tau.mat) <- paste0("ax", 1:ncol(tau.mat))
  tau.fitmat <- tau.mat*a

  dat.fit <- cbind(c1.fitmat, c0.fitmat,
                   tau.fitmat, lambda.fitmat, t, delta)
  dat.fit <- data.frame(dat.fit)

  # initial value for adaptive lasso -- Cox model without penalization
  int.var <- paste(c(colnames(c1.fitmat), colnames(c0.fitmat),
                     colnames(tau.fitmat), colnames(lambda.fitmat),
                     "survival::strata(s)"),
                   collapse=" + ")
  formula.int <- reformulate(int.var, response = "survival::Surv(t, delta)")
  fit.int.ini <- survival::coxph(formula.int, data = dat.fit)
  par.int.ini <- fit.int.ini$coefficients
  theta.pen.ini <- par.int.ini[colnames(lambda.fitmat)]

  # fit adaptive lasso, with penalty on lambda(X)
  y.strata <- glmnet::stratifySurv(survival::Surv(t, delta), strata = s)
  xvar.mat <- dat.fit[,c(colnames(c1.fitmat), colnames(c0.fitmat),
                         colnames(tau.fitmat), colnames(lambda.fitmat))]
  xvar.mat <- data.matrix(xvar.mat)
  penalty.factor <- rep(0, length(par.int.ini)-1)
  penalty.factor[(ncol(c1.fitmat)+ncol(c0.fitmat)+ncol(tau.fitmat)+1):
                    (length(par.int.ini)-1)] <- 1/abs(theta.pen.ini)

  fit.cv <- suppressWarnings(glmnet::cv.glmnet(x = xvar.mat,
                              y = y.strata,
                              family = "cox",
                              penalty.factor = penalty.factor,
                              thresh = 1e-14,
                              nfolds = nfolds))
  par.int <- as.vector(coef(fit.cv, s = "lambda.min"))

  # After varaible selection, refit the Cox model
  int.allvar <- c(colnames(c1.fitmat), colnames(c0.fitmat),
                  colnames(tau.fitmat), colnames(lambda.fitmat))
  index0 <- which(par.int == 0)
  if(length(index0) == 0 || length(index0) == length(int.allvar)){
    int.var <- int.allvar
  }
  else{
    int.var <- int.allvar[-index0]
  }
  formula.int <- reformulate(c(int.var, "survival::strata(s)"),
                             response = "survival::Surv(t, delta)")
  fit.int <- survival::coxph(formula.int,
                             data = dat.fit)
  par.int.refit <- fit.int$coefficients
  coefse.int <- summary(fit.int)$coefficients[,3]
  betahat.int <- par.int.refit[colnames(tau.fitmat)]
  betase.int <- coefse.int[colnames(tau.fitmat)]
  var.betahat <- fit.int$var[(ncol(c1.fitmat) + ncol(c0.fitmat) + 1):
                               (ncol(c1.fitmat) + ncol(c0.fitmat) + ncol(tau.fitmat)),
                             (ncol(c1.fitmat) + ncol(c0.fitmat) + 1):
                               (ncol(c1.fitmat) + ncol(c0.fitmat) + ncol(tau.fitmat))]
  var.betahat.int <- var.betahat

  taux.int <- apply(tau.fitmat, 1, function(x) sum(x*betahat.int))
  ate.int <- sum((1 - s)*taux.int)/sum(1 - s)

  xbar.int <- matrix(colMeans(tau.fitmat[s==0,]), nrow = 1)
  xvar.int <- var(tau.fitmat[s==0,])/sum(1 - s)
  betamat <- matrix(betahat.int, nrow = 1)
  ateve.int <- as.numeric(xbar.int%*%var.betahat.int%*%t(xbar.int)) +
    as.numeric(betamat%*%xvar.int%*%t(betamat))

  return(list(
    # point estimate
    betahat = betahat.int,
    ve.beta = betase.int^2,
    var.betahat = var.betahat.int,
    ate = ate.int,
    ve.ate = ateve.int
  ))
}


#' Estimation of the heterogeneous treatment effect based on the covariates from
#' the basis approximation for the trial data
#'
#' This function obtains the trial estimator with its variance estimates by the
#' Cox model.
#'
#' @param tau.mat the matrix that contains the covariate terms of \eqn{\tau(X)}.
#' @param c1.sieve the covariate matrix for \eqn{c_1(X)} in the sieve
#' approximation.
#' @param a the binary treatment assignment, where the active treatment group
#' should be encoded as `1`, and the control group should be encoded as `0`.
#' @param s the binary data source indicator, where the trial data should be
#' encoded as `1`, and the observational study should be encoded as `0`.
#' @param t the observed event time.
#' @param delta the event indicator.
#' @return a list of estimators, including
#' \itemize{
#' \item `betahat`: the estimated treatment effect parameter
#' \item `ve.beta`: the variance estimate of the estimated treatment effect
#' \item `var.betahat`: the covariance matrix of the estimated treatment effect
#' parameter
#' \item `ate`: the estimated average treatment effect
#' \item `ve.ate`: the variance estimate of the average treatment effect
#' }
#' @import survival
#' @export
fitcox.rct <- function(tau.mat, c1.sieve, a, s = NULL, t, delta){
  if(is.null(s)){
    c1.fitmat <- c1.sieve
    colnames(tau.mat) <- paste0("ax", 1:ncol(tau.mat))
    tau.fitmat <- tau.mat*a
    t <- t
    a <- a
    delta <- delta
  }
  else{
    index.rct <- which(s == 1)
    c1.fitmat <- c1.sieve[index.rct,]
    colnames(tau.mat) <- paste0("ax", 1:ncol(tau.mat))
    tau.fitmat <- (tau.mat*a)[index.rct,]
    t <- t[index.rct]
    a <- a[index.rct]
    delta <- delta[index.rct]
  }
  n.rct <- nrow(c1.fitmat)
  dat.fit <- cbind(c1.fitmat, tau.fitmat, t, delta)
  dat.fit <- data.frame(dat.fit)

  # initial value for adaptive lasso -- Cox model without penalization
  rct.var <- paste(c(colnames(c1.fitmat), colnames(tau.fitmat)),
                   collapse=" + ")
  formula.rct <- reformulate(rct.var, response = "survival::Surv(t, delta)")
  fit.rct <- survival::coxph(formula.rct, data = dat.fit)
  par.rct <- fit.rct$coefficients
  betahat.rct <- par.rct[colnames(tau.fitmat)]
  coefse.rct <- summary(fit.rct)$coefficients[,3]
  betase.rct <- coefse.rct[colnames(tau.fitmat)]

  var.betahat <- fit.rct$var[(ncol(c1.fitmat) + 1):
                               (ncol(c1.fitmat) + ncol(tau.fitmat)),
                             (ncol(c1.fitmat) + 1):
                               (ncol(c1.fitmat) + ncol(tau.fitmat))]
  var.betahat.rct <- var.betahat

  taux.rct <- apply(tau.fitmat, 1, function(x) sum(x*betahat.rct))
  ate.rct <- mean(taux.rct)

  xbar.rct <- matrix(colMeans(tau.fitmat), nrow = 1)
  xvar.rct <- var(tau.fitmat)/n.rct
  betamat <- matrix(betahat.rct, nrow = 1)
  ateve.rct <- as.numeric(xbar.rct%*%var.betahat.rct%*%t(xbar.rct)) +
    as.numeric(betamat%*%xvar.rct%*%t(betamat))

  return(list(
    # point estimate
    betahat = betahat.rct,
    ve.beta = betase.rct^2,
    var.betahat = var.betahat.rct,
    ate = ate.rct,
    ve.ate = ateve.rct
  ))
}


#' Integrative estimation of the heterogeneous treatment effect with the use of
#' sieve approximation
#'
#' This function is the main function of obtaining the integrative estimator
#' with its variance estimates. The integrative estimator is obtained by
#' minimizing the penalized log partial likelihood, with the penalty function
#' selected as the adaptive lasso.
#'
#' @param tau.mat the matrix that contains the covariate terms of \eqn{\tau(X)}.
#' @param c1.mat the matrix that contains the continuous covariates which may
#' affect \eqn{c_1(X)}.
#' @param c1.discrete the matrix that contains the discrete covariates which may
#' affect \eqn{c_1(X)}. The default is `NULL`.
#' @param q.c1 the number of basis functions to approximate each continuous
#' covariate in \eqn{c_1(X)}. Here, it corresponds to the degree of freedom to
#' generate the polynomial basis.
#' @param c0.mat the matrix that contains the continuous covariates which may
#' affect \eqn{c_0(X)}.
#' @param c0.discrete the matrix that contains the discrete covariates which may
#' affect \eqn{c_0(X)}. The default is `NULL`.
#' @param q.c0 the number of basis functions to approximate each continuous
#' covariate in \eqn{c_0(X)}. Here, it corresponds to the degree of freedom to
#' generate the polynomial basis.
#' @param lambda.mat the matrix that contains the continuous covariates which
#' may affect \eqn{\lambda(X)}.
#' @param lambda.discrete the matrix that contains the discrete covariates
#' which may affect \eqn{\lambda(X)}.
#' @param q.lambda the number of basis functions to approximate each continuous
#' covariate in \eqn{\lambda(X)}. Here, it corresponds to the degree of freedom
#' to generate the polynomial basis.
#' @param a the binary treatment assignment, where the active treatment group
#' should be encoded as `1`, and the control group should be encoded as `0`.
#' @param s the binary data source indicator, where the trial data should be
#' encoded as `1`, and the observational study should be encoded as `0`.
#' @param t the observed event time.
#' @param delta the event indicator.
#' @param nfolds the fold of cross-validation, the default value is `5`.
#' @param type the type of estimator the user want to obtain. Available types
#' include
#' \itemize{
#' \item `int`: the integrative estimator that combines the trial data and the
#' observational study data;
#' \item `rct`: the trial estimator.
#' }
#' @return a list of estimators, including
#' \itemize{
#' \item `beta.est`: the estimated treatment effect parameter
#' \item `ve.beta`: the variance estimate of the estimated treatment effect
#' \item `cov.beta`: the covariance matrix of the estimated treatment effect
#' parameter
#' \item `ate.est`: the estimated average treatment effect
#' \item `ve.ate`: the variance estimate of the average treatment effect
#' }
#' @import glmnet survival
#' @export
#'
#' @examples
#' q.c1 <- 3
#' q.c0 <- 3
#' q.lambda <- 3
#' nfolds <- 5
#' c1.mat <- matrix(dat$x1, ncol = 1)
#' c1.discrete <- matrix(dat$x2, ncol = 1)
#' c0.mat <- matrix(dat$x1, ncol = 1)
#' c0.discrete <- matrix(dat$x2, ncol = 1)
#' lambda.mat <-  matrix(dat$x1, ncol = 1)
#' lambda.discrete <- cbind(1, dat$x2)
#' tau.mat <- cbind(1, dat$x1, dat$x1^2)
#' s <- dat$s
#' delta <- dat$delta
#' t <- dat$t
#' a <- dat$a
#' res.hte <- surv.hte(tau.mat = tau.mat, c1.mat = c1.mat,
#' c1.discrete = c1.discrete, q.c1 = q.c1,
#' c0.mat = c0.mat, c0.discrete = c0.discrete, q.c0 = q.c0,
#' lambda.mat = lambda.mat,
#' lambda.discrete = lambda.discrete, q.lambda = q.lambda,
#' a = a, s = s, t = t, delta = delta,
#' nfolds = nfolds, type = c("int", "rct"))
#' res.hte
surv.hte <- function(tau.mat, c1.mat, c1.discrete = NULL, q.c1,
                     c0.mat, c0.discrete = NULL, q.c0,
                     lambda.mat, lambda.discrete = NULL, q.lambda,
                     a, s, t, delta, nfolds = 5, type = "int"){
  xmat <- sieve.mat(c1.mat = c1.mat,
                    c1.discrete = c1.discrete, q.c1 = q.c1,
                    c0.mat = c0.mat,
                    c0.discrete = c0.discrete, q.c0 = q.c0,
                    lambda.mat = lambda.mat,
                    lambda.discrete = lambda.discrete, q.lambda = q.lambda)
  c1.sieve <- xmat$c1.sieve
  c0.sieve <- xmat$c0.sieve
  lambda.sieve <- xmat$lambda.sieve

  beta.est <- NULL
  ve.beta <- NULL
  cov.beta <- NULL
  name.list <- NULL
  ate.est <- NULL
  ve.ate <- NULL
  if("int" %in% type){
    res.int <- fitcox.int(tau.mat = tau.mat, c1.sieve = c1.sieve,
                          c0.sieve = c0.sieve, lambda.sieve = lambda.sieve,
                          a = a, s = s, t = t, delta = delta,
                          nfolds = nfolds)
    beta.est <- c(beta.est, list(res.int$betahat))
    ve.beta <- c(ve.beta, list(res.int$ve.beta))
    cov.beta <- c(cov.beta, list(res.int$var.betahat))
    name.list <- c(name.list, "int")
    ate.est <- c(ate.est, res.int$ate)
    names(ate.est)[length(ate.est)] <- "int"
    ve.ate <- c(ve.ate, res.int$ve.ate)
    names(ve.ate)[length(ve.ate)] <- "int"
  }

  if("rct" %in% type){
    res.rct <- fitcox.rct(tau.mat = tau.mat, c1.sieve = c1.sieve,
                          a = a, s = s, t = t, delta = delta)
    beta.est <- c(beta.est, list(res.rct$betahat))
    ve.beta <- c(ve.beta, list(res.rct$ve.beta))
    cov.beta <- c(cov.beta, list(res.rct$var.betahat))
    name.list <- c(name.list, "rct")
    ate.est <- c(ate.est, res.rct$ate)
    names(ate.est)[length(ate.est)] <- "rct"
    ve.ate <- c(ve.ate, res.rct$ve.ate)
    names(ve.ate)[length(ve.ate)] <- "rct"
  }
  names(beta.est) <- name.list
  names(ve.beta) <- name.list
  names(cov.beta) <- name.list

  return(list(beta.est = beta.est,
              ve.beta = ve.beta,
              cov.beta = cov.beta,
              ate.est = ate.est,
              ve.ate = ve.ate))
}

#' @importFrom stats coef poly reformulate var
NULL
