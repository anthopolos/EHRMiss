#' @title Update subject level random effects \code{bSub}.
#'
#' @description A hierarchically-centered approach is used to update the random effects for each of \code{Y1,\dots,YJ}.
#' @export
#' @return A \code{J}-element list for each of \code{Y1,\dots,YJ}, with each element containing a vector of stacked subject level random effects, updated conditional on latent class membership.
update_bSubLC_fast <- function(C, bSub, betaObs, betaSub, Psi, Sigma, Y, XRe, XObs, XSub, subjectID) {

  ### Update random effects for each b[[j]], conditional on latent class
  #! Function assumes the same V, XObs, XSub for each j

  # Number of latent classes
  K <- length(table(C))

  # Number of outcomes
  J <- ncol(Y)

  # Unique number of subjects
  n <- length(unique(subjectID))

  # Number of random effects
  q <- ncol(XRe)

  C_expand <- C[factor(subjectID)]

  bSubm <- lapply(bSub, function(x) matrix(x, nrow = n, ncol = q, byrow = TRUE))

  values <- as.list(1:J)

  for (j in 1:J) {

    # For each Y_j, generate an n x q matrix of random effects conditional on latent class membership
    tempValues <- matrix(NA, nrow = n, ncol = q)

    for (k in 1:K) {

      # Subject level information
      ind_sub <- which(C == k)
      XSubk <- as.matrix(XSub[ind_sub, ])

      # Observation level information
      ind_obs <- which(C_expand == k)
      Yk <- Y[ind_obs, ]
      XRek <- XRe[ind_obs, ]
      XObsk <- as.matrix(XObs[ind_obs, ])
      subjectIDk <- subjectID[ind_obs]

      # Define the matrix of coefficients R per A.2 in Gelman BDA, page 580
      #R will be J x J
      R <- diag(1, nrow = J, ncol = J) - solve(diag(solve(Sigma[ , , k])) * diag(1, nrow = J, ncol = J)) %*% solve(Sigma[ , , k])

      # Gather parts outside of subject level "loop"
      # Posterior variance calculation pieces
      # Conditional variance of Y_j
      condVarYj <- as.numeric(solve(solve(Sigma[ , , k])[j, j]))
      Psijk <- Psi[[j]][ , , k]

      # Posterior mean calculation pieces
      # Contribution from not j
      betaObsNotj <- matrix(unlist(lapply(betaObs[-j], function(x){as.matrix(x)[ , k]})), ncol = J - 1, byrow = FALSE)

      #! Changed on April 26, 2022
      #bSubNotj <- matrix(unlist(lapply(bSub[-j], function(x) {
      #
        # Calculate contribution from bNotj subsetted to latent class k
      # tempm <- matrix(x, nrow = n, ncol = q, byrow = TRUE)
      # reContj <- unlist(sapply(unique(subjectIDk), function(u) {
      #   selObs <- which(subjectIDk == u)
      #   XRekObs <- as.matrix(XRek)[selObs, ]
      #   tempmk <- as.matrix(tempm[which(unique(subjectID) == u), ])
      #   mu <- as.matrix(XRekObs) %*% c(t(tempmk))
      #   return(mu)
      # }))
      #
      #}
      #
      #)), ncol = J - 1)

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

      #! Use of matrix(, ncol = J-1) has an implicit c(reContj) before placing in matrix

      diffyNotj <- Yk[ , -j] - XObsk %*% betaObsNotj - bSubNotj  # N x p - 1
      # condVarYj is a scalar
      likTempk <- (Yk[ , j] - XObsk %*% betaObs[[j]][ , k] - t(as.vector(R[j, -j]) %*% t(diffyNotj))) / condVarYj

      # Prior contribution
      ### Changed on April 26, 2022
      #priorTempk <- XSubk %*% betaSub[[j]][ , , k] %*% solve(Psijk)
      priorTempk <- XSub %*% betaSub[[j]][ , , k] %*% solve(Psijk)

      # Now calculate the updated random effects for subjects in class k using sapply
      tempValues[ind_sub, ] <- t(sapply(unique(subjectIDk), function(x) {

        selObs <- which(subjectIDk == x)

        # Posterior variance
        XRekObs <- as.matrix(XRek)[selObs, ]
        XReXRe <- crossprod(XRekObs, XRekObs)
        post_var <- solve(solve(Psijk) + XReXRe / condVarYj)

        # Posterior mean
        likTempkObs <- likTempk[selObs]

        ### Changed on April 26, 2022
        #priorkObs <- priorTempk[which(unique(subjectIDk) == x)]
        priorkObs <- priorTempk[x, ]
        likkObs <- t(XRekObs) %*% likTempkObs
        post_mean <- post_var %*% (likkObs + priorkObs)

        # Draw values for the ith subject
        values <- mnormt::rmnorm(1, mean = post_mean, varcov = post_var)
      }
      ))

    } # End of k loop

    ### Convert to a n*q stacked vector of random effects for each Y_j
    values[[j]] <- c(t(tempValues))

  }

  return(values) # nq x 1 stacked random effects for a J element list

}

