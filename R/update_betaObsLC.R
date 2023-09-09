#' @title Update regression coefficients \code{betaObs}.
#'
#' @description An multivariate normal prior distribution is used to update the regression coefficients \code{betaObs} corresponding to the covariates in \code{XObs}, for each outcome \code{Y1,\dots,YJ}.
#' @export
#' @return A \code{J} element list in which each element is an \code{s} by \code{K} matrix of latent class-specific regression coefficients.
update_betaObsLC <- function(C, betaObs, bSub, Sigma, Y, XRe, XObs, prior.mu, prior.Sigma, subjectID){

  ### Background information
  K <- length(table(C))
  J <- ncol(Y)
  s <- ncol(XObs)
  q <- ncol(XRe)
  n <- length(unique(subjectID))
  C_expand <- C[factor(subjectID)]

  bSubm <- lapply(bSub, function(x) matrix(x, nrow = n, ncol = q, byrow = TRUE))


  for (j in 1:J) {

    tempValues <- matrix(NA, nrow = s, ncol = K)

    for (k in 1:K) {

      # Subject level information
      ind_sub <- which(C == k)

      # Observation level information
      ind_obs <- which(C_expand == k)
      Yk <- Y[ind_obs, ]
      XRek <- XRe[ind_obs, ]
      XObsk <- XObs[ind_obs, ]
      subjectIDk <- subjectID[ind_obs]

      #! Function restricts to the same X and XRe for each j
      XX <- crossprod(XObsk, XObsk)


      # N by n*q design matrix for random effects, to be multiplied by nq by 1 vector of stacked random effects
      #sp_XRe_sub <- lapply(split(as.data.frame(XRek), subjectIDk, drop = TRUE), as.matrix)
      #convert_XRe_sub <- as.matrix(Matrix::bdiag(sp_XRe_sub))

      # Define the matrix of coefficients R per A.2 in Gelman BDA, page 580
      #R will be J x J
      R <- diag(1, nrow = J, ncol = J) - solve(diag(solve(Sigma[ , , k])) * diag(1, nrow = J, ncol = J)) %*% solve(Sigma[ , , k])

      ### Posterior variance for Y_j corresponding to beta_j regression coefficients to be updated
      # Conditional variance of Y_j
      condVarYj <- as.numeric(solve(solve(Sigma[ , , k])[j, j]))
      # Posterior variance of beta_j
      post_var <- solve(XX / condVarYj + solve(prior.Sigma))

      ### Posterior mean
      # Calculation for not j
      #Get vector of betaObs for class k for each Y_{-j} in a column
      betaObsNotj <- matrix(unlist(lapply(betaObs[-j], function(x){as.matrix(x)[ , k]})), ncol = J - 1, byrow = FALSE)


      #bSubNotj <- matrix(unlist(lapply(bSub[-j], function(x) {
        # Calculate contribution from bNotj subsetted to latent class k
       # tempm <- matrix(x, nrow = n, ncol = q, byrow = TRUE)
        #Retstrict to class k using ind_sub
        #tempmk <- as.matrix(tempm[ind_sub, ])
        #mu <- as.vector(convert_XRe_sub %*% c(t(tempmk)))
        #return(mu)
      #}
      #)), ncol = J - 1)
      #diffyNotj <- Yk[ , -j] - XObsk %*% betaObsNotj - bSubNotj  # N x J - 1

      # Calculation for bj
      #tempm <- matrix(bSub[[j]], nrow = n, ncol = q, byrow = TRUE)
      #tempmk <- as.matrix(tempm[ind_sub, ])
      #bSubj <- as.vector(convert_XRe_sub %*% c(t(tempmk)))

      #post_mean <- post_var %*% (t(XObsk) %*% (Yk[ , j] - bSubj - t(as.vector(R[j, -j]) %*% t(diffyNotj))) / condVarYj + t(prior.mu %*% solve(prior.Sigma)))


      ### Posterior mean as of April 26, 2022
      bSubNotj <- matrix(unlist(lapply(bSubm[-j], function(x) {

        # Calculate contribution from bNotj subsetted to latent class k
        reContj <- unlist(sapply(unique(subjectIDk), function(u) {
          selObs <- which(subjectIDk == u)
          XRekObs <- as.matrix(XRek)[selObs, ]

          xk <- as.matrix(x[u, ])
          mu <- as.matrix(XRekObs) %*% c(t(xk))
          return(mu)
        }))

      }

      )), ncol = J - 1)


      diffyNotj <- Yk[ , -j] - XObsk %*% betaObsNotj - bSubNotj  # N x J - 1


      # Calculation for bj
      #tempm <- matrix(bSub[[j]], nrow = n, ncol = q, byrow = TRUE)
      #tempmk <- as.matrix(tempm[ind_sub, ])
      #bSubj <- as.vector(convert_XRe_sub %*% c(t(tempmk)))

      bSubj <- unlist(sapply(unique(subjectIDk), function(u) {
        selObs <- which(subjectIDk == u)
        XRekObs <- as.matrix(XRek)[selObs, ]

        xk <- as.matrix(bSubm[[j]][u, ])
        mu <- as.matrix(XRekObs) %*% c(t(xk))
        return(mu)
      }))

      post_mean <- post_var %*% (t(XObsk) %*% (Yk[ , j] - bSubj - t(as.vector(R[j, -j]) %*% t(diffyNotj))) / condVarYj + t(prior.mu %*% solve(prior.Sigma)))


      tempValues[ , k] <- mnormt::rmnorm(1, mean = post_mean, varcov = post_var)

    }

    betaObs[[j]] <- tempValues

  }

  return(betaObs)

}
