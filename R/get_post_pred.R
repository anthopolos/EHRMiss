#' @title Get posterior predictive draws of the completed data at each iteration.
#'
#' @description The reference for the approach is from algorithm 3.7 on page 91 of Fruhwirth-Schnatter, S. (2006) Finite Mixture and Markov Switching Models. Springer Science & Business Media, New York. The latent classes are redrawn using the probabilities of latent class membership obtained from the latent class membership model.
#' @return Draws from the posterior predictive distribution at the current iteration.
get_post_pred <- function(J, priorPik, bSub, betaObs, Sigma, XRe, XObs, subjectIDY) {

  n <- dim(as.matrix(priorPik))[1]
  K <- ncol(as.matrix(priorPik))
  q <- length(bSub[[1]]) / n
  N <- length(subjectIDY)

  # Redraw C per algorithm 3.7
  Cdraw <- Hmisc::rMultinom(priorPik, 1)

  Cdraw_expand <- Cdraw[factor(subjectIDY)]

  store_Ydraw <- matrix(NA, nrow = N, ncol = J)

  for (k in 1:K) {

    # Observation level information
    ind_obs <- which(Cdraw_expand == k)
    subjectIDYk <- subjectIDY[ind_obs]
    XRek <- as.matrix(XRe[ind_obs, ])
    XObsk <- as.matrix(XObs[ind_obs, ])

    Nk <- length(subjectIDYk)

    # Subject level information
    ind_sub <- which(Cdraw == k)

    # N by n*q design matrix for random effects, to be multiplied by nq by 1 vector of stacked random effects
    sp_XRe_sub <- lapply(split(as.data.frame(XRek), subjectIDYk, drop = TRUE), as.matrix)
    convert_XRe_sub <- as.matrix(Matrix::bdiag(sp_XRe_sub))

    muObs <- matrix(NA, ncol = J, nrow = Nk)

    for (j in 1:J) {

      # Calculation for bj
      tempm <- matrix(bSub[[j]], nrow = n, ncol = q, byrow = TRUE)
      tempmk <- as.matrix(tempm[ind_sub, ])
      bSubj <- as.vector(convert_XRe_sub %*% c(t(tempmk)))

      muObs[ , j] <- XObsk %*% betaObs[[j]][ , k] + bSubj

    }

    # Generate Y using latent class specific variance-covariate of Y_js, making sure to place values based on ind_obs
    store_Ydraw[ind_obs, ] <- t(apply(muObs, 1, function(x) {mnormt::rmnorm(1, mean = x, varcov = Sigma[ , , k])}))

  }

  return(store_Ydraw)

}
