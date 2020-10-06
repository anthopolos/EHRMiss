#' @title Update regression coefficients \code{delta} in multinomial probit model.
#' @description This uses a multivariate normal prior distribution to update the latent class-specific regression coefficients in the multinomial probit model.
#' @return An \code{m} by \code{K-1} matrix of latent class-specific regression coefficients in the multinomial probit model with \code{K=1} fixed to zero.
update_delta_MNP <- function(Z, delta, W, prior.mu, prior.Sigma){

  # Number of classes
  K <- ncol(delta) + 1

  # Update delta with proper conjugate, note the error distribution on latent variable is N(X\beta, 1) so posterior variance of delta does not get updated
  #See Chib 1993
  post_var <- solve(t(W) %*% W + solve(prior.Sigma))

  for (k in 1:(K - 1)) {

    post_mean <- post_var %*% (t(W) %*% Z[ , k] + solve(prior.Sigma) %*% prior.mu)
    delta[ , k] <- mnormt::rmnorm(1, mean = post_mean, varcov = post_var)

  }

  return(delta)

}
