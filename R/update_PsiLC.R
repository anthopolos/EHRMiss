#' @title Update variance-covariance matrix \code{Psi} of the random effects for \code{Y1,\dots,YJ}.
#'
#' @description An inverse-Wishart prior distribution is used to update the variance-covariance matrix \code{Psi} of the random effects for each outcome variable. If only a random intercept is present, then the prior collapses to an inverse gamma distribution.
#'
#' @return A \code{J}-element list. Each element is latent class-specific variance-covariance matrix stored in a \code{q} by \code{q} by \code{K} array. For example, \code{Psi = list(Psi1 = array( , dim = c(q, q, K)), Psi2 = array( , dim = c(q, q, K)))}. If \code{q=1}, a scalar for the variance of the random intercept replaces the variance-covariance matrix.
update_PsiLC <- function(C, betaSub, bSub, Y, XSub, prior.scale, prior.df) {

  # Number of classes
  K <- length(table(C))

  #! Function restricts to the same XSub and q for each j
  #! Need to test to make sure it works with q = 1 versus q > 1
  n <- dim(XSub)[1]
  J <- ncol(Y)
  if (!is.null(dim(prior.scale))) {
    q <- dim(prior.scale)[1]
  } else {
    q <- 1
  }

  values <- as.list(1:J)

  ### Each Yj has random effects with its own latent class specific variance-covariance
  for (j in 1:J) {

    bSubjm <- matrix(bSub[[j]], nrow = n, ncol = q, byrow = TRUE)

    valuesTemp <- array(NA, dim = c(q, q, K))

    for (k in 1:K) {

      ind_sub <- which(C == k)
      bSubjmk <- bSubjm[ind_sub, ]
      XSubk <- as.matrix(XSub[ind_sub, ])

      # Posterior degrees of freedom
      df <- length(ind_sub)
      post_df <- df + prior.df

      # Posterior scale matrix for class k
      mu <- XSubk %*% betaSub[[j]][ , , k]
      lik <- crossprod(bSubjmk - mu, bSubjmk - mu)
      post_scale <- lik + prior.scale # q x q

      valuesTemp[ , , k] <- MCMCpack::riwish(v = post_df, S = post_scale)

    }

    values[[j]] <- valuesTemp
  }

  return(values)

}
