#' @title Update variance-covariance matrix \code{Sigma} of the \code{Y1,\dots,YJ}.
#'
#' @description An inverse-Wishart prior distribution is used to update the variance-covariance matrix \code{Sigma} of the longitudinal health outcomes \code{Y1,\dots,YJ}.
#'
#' @return A \code{J} (number of outcome variables) by \code{J} (number of outcome variables) by \code{K} (number of latent class) array of latent class-specific variance-covariance matrices.
update_SigmaLC <- function(C, betaObs, bSub, Y, XRe, XObs, prior.scale, prior.df, subjectID) {

  #! Restricts to the same XRe and XObs for each j

  # Background information
  n <- length(unique(subjectID))
  K <- length(table(C))
  J <- dim(Y)[2]
  q <- ncol(XRe)

  C_expand <- C[factor(subjectID)]

  values <- array(NA, dim = c(J, J, K))

  # Calculate latent class specific variance covariance matrices for the Y[ , j] (restricted to observations to subjects in class k)
  for (k in 1:K) {

    # Subject level information
    ind_sub <- which(C == k)

    # Observation level information
    ind_obs <- which(C_expand == k)
    Yk <- Y[ind_obs, ]
    XRek <- XRe[ind_obs, ]
    XObsk <- as.matrix(XObs[ind_obs, ])
    subjectIDk <- subjectID[ind_obs]

    Nk <- dim(Yk)[1]

    # Posterior degrees of freedom
    df <- Nk
    post_df <- df + prior.df

    ### Posterior scale matrix calculation
    mu <- matrix(NA, Nk, J)
    sp_XRe_sub <- lapply(split(as.data.frame(XRek), subjectIDk, drop = TRUE), as.matrix)
    convert_XRe_sub <- as.matrix(Matrix::bdiag(sp_XRe_sub))

    for (j in 1:J) {
      # Calculate latent class specific mu for each outcome j
      tempm <- matrix(bSub[[j]], nrow = n, ncol = q, byrow = TRUE)
      tempmk <- as.matrix(tempm[ind_sub, ])
      mu[ , j] <- as.vector(convert_XRe_sub %*% c(t(tempmk)) + XObsk %*% betaObs[[j]][ , k])
    }

    lik <- crossprod(Yk - mu, Yk - mu)
    post_scale <- lik + prior.scale

    values[ , , k] <- MCMCpack::riwish(v = post_df, S = post_scale)

  }

  return(values)

}
