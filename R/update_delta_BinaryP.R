#' @title Update regression coefficients \code{delta} in a binary probit model.
#'
#' @description This uses a normal prior distribution to update the regression coefficients in the probit model.
#' @return An \code{m} by 1 vector of regression coefficients in the binary probit model.
update_delta_BinaryP <- function(Z, delta, W, prior.mu, prior.Sigma){

  m <- ncol(as.matrix(W))
  K <- ncol(as.matrix(Z)) + 1

  # Update delta with proper conjugate, note the error distribution on latent variable is N(X\beta, 1) so posterior variance of delta does not get updated
  #See Chib 1993
  post_var <- solve(t(W) %*% W + solve(prior.Sigma))

  post_mean <- post_var %*% (t(W) %*% Z + solve(prior.Sigma) %*% prior.mu)

  delta <- matrix(mnormt::rmnorm(1, mean = post_mean, varcov = post_var), nrow = m, ncol = K - 1)

  return(delta)

}
