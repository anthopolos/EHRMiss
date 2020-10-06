#' @title Update latent class-specific regression coefficients in a probit model.
#'
#' @description Latent class-specific regression coefficients are updated based on a Bernoulli likelihood and probit link function. A normal prior distribution is assumed. The updating function can accomodate probit models that include a subject specific random effects.
#' @return A matrix of latent class-specific regression coefficients in which columns index the \code{K} latent classes.
update_betaLC_ProbitRE <- function(C, Z, tau, UObs, URe, prior.mu, prior.Sigma, subjectID){

  K <- length(table(C))
  n <- length(C)
  C_expand <- C[factor(subjectID)]

  values <- matrix(NA, nrow = ncol(UObs), ncol = K)

  taum <- matrix(tau, nrow = n, ncol = ncol(URe), byrow = TRUE)
  for (k in 1:K) {

    ind_sub <- which(C == k)
    ind_obs <- which(C_expand == k)
    subjectIDk <- subjectID[ind_obs]

    ### Posterior variance
    #Error of z is fixed to 1 so post_var does not get updated
    post_var <- solve(crossprod(as.matrix(UObs[ind_obs, ]), as.matrix(UObs[ind_obs, ])) + solve(prior.Sigma))

    ### Posterior mean calculation with and without RE
    if (!is.null(tau)) {

      ### Design matrix and subgroup of RE for class k
      URek <- URe[C_expand == k, ]
      tauk <- taum[ind_sub, ]

      sp_n_sub <- lapply(split(as.data.frame(URek), subjectIDk, drop = TRUE), as.matrix)
      convert_n_sub <- as.matrix(Matrix::bdiag(sp_n_sub))

      post_mean <- post_var %*% (t(as.matrix(UObs[ind_obs, ])) %*% (Z[ind_obs] - as.vector(convert_n_sub %*% c(t(tauk)))) + solve(prior.Sigma) %*% prior.mu)

    } else {

      post_mean <- post_var %*% (t(as.matrix(UObs[ind_obs, ])) %*% Z[ind_obs] + solve(prior.var) %*% prior.mean)

    }

    values[ , k] <- matrix(mnormt::rmnorm(1, mean = post_mean, varcov = post_var), nrow = ncol(UObs), ncol = 1)

  }

  return(values)

}
