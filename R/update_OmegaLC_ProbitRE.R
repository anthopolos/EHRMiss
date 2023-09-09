#' @title Update the latent class-specific variance-covariance matrices of the subject-level random effects in a probit model.
#'
#' @description The variance-covariance matrices of the random effects are updated conditional on class membership using an inverse-Wishart prior distribution. This updating function is based on a probit model (Bernoulli likelihood) and a multivariate normal distribution for the random effects.
#' @export
#' @return A \code{q} by \code{q} by \code{K} array of variance-covariance matrices for the random effects, where \code{q} is the number of random effects.
update_OmegaLC_ProbitRE <- function(C, tau, prior.scale, prior.df) {

  K <- length(table(C))
  n <- length(C)

  if (!is.null(dim(prior.scale))) {
    q <- dim(prior.scale)[1]
  } else {
    q <- 1
  }

  values <- array(NA, dim = c(q, q, K))
  taum <- matrix(tau, nrow = n, ncol = q, byrow = TRUE)

  for (k in 1:K) {

    ind_sub <- which(C == k)
    tauk <- taum[ind_sub, ]

    ### Posterior degrees of freedom
    df <- length(ind_sub)
    post_df <- df + prior.df

    ### Posterior scale matrix for class k
    # Prior mean is 0
    lik <- crossprod(tauk, tauk)
    post_scale <- lik + prior.scale # q x q

    values[ , , k] <- MCMCpack::riwish(v = post_df, S = post_scale)

  }

  return(values)

}
