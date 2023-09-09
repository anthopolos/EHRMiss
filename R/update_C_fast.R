#' @title Update posterior latent class assignment \code{C} for each subject.
#'
#' @description This uses Bayes' theorem to update posterior latent class assignment. This is a ratio of the likelihood contribution of each subject to class \code{k} and the prior probability of belonging in class \code{k} to the likelihood contribution marginalized over latent classes. The function accomodates the case of \code{K} equal to 2 or greater than 2.
#' @export
#' @return A list of 2 elements. An \code{n} by 1 vector of latent class assignments called \code{C}. An \code{n} by \code{K} matrix of posterior probabilities of latent class assignment.
update_C_fast <- function(priorPik, betaObs, bSub, betaSub, Sigma, Psi, phi, tau, Omega, lambda, kappa, Theta,  Y, XRe, XObs, XSub, D, UObs, URe, M, Mvec, VObs, VRe, subjectIDY, subjectID, subjectIDM, modelVisit, modelResponse){

  ### Number of latent classes
  K <- ncol(as.matrix(priorPik))

  ### Number of subjects
  n <- length(unique(subjectIDY))

  ### Construct likelihood contribution for subject i to each class k
  J <- ncol(Y)
  N <- nrow(Y)
  q <- ncol(XRe)

  # For Y
  # Place holder for log likelihood contribution from Y_i, D_i, M_i
  # Factorization Pr(C_i = k \ | \ b_i, Y_i) \propto MVN_q(b_{1i}, \ | \ C, \Psi_1, betaSub_1) MVN_q(b_{2i}, \ | \ C, \Psi_2, betaSub_2) \prod_{t = 1}^{n_i} MVN_J(Y_{1i}, Y_{2i}, \ | \ b_{1i}, b_{2i}, \Sigma_2, C_i)
  # Log scale \log MVN_q(b_{1i}, \ | \ C, \Psi_1, betaSub_1) + \log MVN_q(b_{2i}, \ | \ C, \Psi_2, betaSub_2) + \sum_{t = 1}^{n_i} \log MVN_J(Y_{1i}, Y_{2i}, \ | \ b_{1i}, b_{2i}, \Sigma_2, C_i)

  likY <- matrix(NA, nrow = n, ncol = K)

  if (modelVisit == TRUE) {
    likD <- matrix(NA, nrow = n, ncol = K)
    # For D random effects piece
    #sp_URe_sub <- lapply(split(as.data.frame(URe), subjectID, drop = TRUE), as.matrix)
    #convert_URe_sub <- as.matrix(Matrix::bdiag(sp_URe_sub))

  }

  if (modelResponse == TRUE) {
    likM <- matrix(NA, nrow = n, ncol = K)
    # For M random effects piece
    #sp_VRe_sub <- lapply(split(as.data.frame(VRe), subjectIDM, drop = TRUE), as.matrix)
    #convert_VRe_sub <- as.matrix(Matrix::bdiag(sp_VRe_sub))

  }


  #! Changed on April 26, 2022
  bSubm <- lapply(bSub, function(x) matrix(x, nrow = n, ncol = q, byrow = TRUE))

  for (k in 1:K) {

    # Storage that changes by latent class
    muObs <- matrix(NA, nrow = N, ncol = J)
    muSub <- matrix(NA, nrow = n, ncol = J)
    likSub <- matrix(NA, nrow = n, ncol = J)

    # For f(Y | b, C) in llikObs, construct muObs = E[Y | b, C], which will be N x J
    #Y: N x J matrix of longitudinal outcomes for each outcome Y_j
    #muObs: N x J matrix of means for latent class k
    #Sigma: J x J variance-covariance specific to latent class k
    #llikObs: N x 1 vector

    for (j in 1:J) {

      # Random effects contribution to observation level
      #bSubjm <- matrix(bSub[[j]], nrow = n, ncol = q, byrow = TRUE)
      bSubjm <- bSubm[[j]]
      bSubjCont <- unlist(sapply(unique(subjectIDY), function(x) {
        selObs <- which(subjectIDY == x)
        XReObs <- as.matrix(XRe)[selObs, ]

        #! Changed on April 26, 2022
        #values <- as.matrix(XReObs) %*% c(t(as.matrix(bSubjm[which(unique(subjectIDY) == x), ])))

        values <- as.matrix(XReObs) %*% c(t(as.matrix(bSubjm[x, ])))
        return(values)
      }))

      # Observation level mean
      muObs[ , j] <- XObs %*% betaObs[[j]][ , k] + c(bSubjCont)

    }

    likObs <- sapply(1:nrow(Y), function(x) {
      mvtnorm::dmvnorm(Y[x, ], mean = muObs[x, ], sigma = Sigma[ , , k], log = FALSE)
    })
    #Conditional on random effects b, observations are independent, so can sum on log scale
    prodlikObs <- as.vector(tapply(likObs, subjectIDY, FUN = prod))

    # Random effects contribution
    #_{j} indexes Y_j
    # f(b_{1i}, b_{2i} | C_i, \Psi_1, \Psi_2) = f(b_{1i} | C_i, betaSub_1, \Psi_1) f(b_{2i}| C, betaSub_2, \Psi_2)
    for (j in 1:J) {
      # Subject level equation
      #bSubm <- matrix(bSub[[j]], nrow = n, ncol = q, byrow = TRUE)

      bSubjm <- bSubm[[j]]
      muSub <- XSub %*% betaSub[[j]][ , , k]
      #likSub[ , j] <- sapply(1:nrow(bSubm), function(x) {
      #  mvtnorm::dmvnorm(bSubm[x, ], mean = muSub[x, ], sigma = as.matrix(Psi[[j]][ , , k]), log = FALSE)
      likSub[ , j] <- sapply(1:nrow(bSubjm), function(x) {
          mvtnorm::dmvnorm(bSubjm[x, ], mean = muSub[x, ], sigma = as.matrix(Psi[[j]][ , , k]), log = FALSE)


      })
    }

    # Sum each random effects contribution for subject i since on log scale
    prodlikSub <- apply(likSub, 1, prod)
    # Log likelihood contribution of Y to class k
    likY[ , k] <- prodlikSub * prodlikObs

    # Log likelihood contribution of visit process D to class k
    if (modelVisit == TRUE) {

      # Observation level contribution
      taum <- matrix(tau, nrow = n, ncol = q, byrow = TRUE)

      # EDITS here
      # Yields T x n
      tauCont <- sapply(unique(subjectID), function(x) {
        selObs <- which(subjectID == x)
        UReObs <- as.matrix(URe)[selObs, ]
        #values <- as.matrix(UReObs) %*% c(t(as.matrix(taum[which(unique(subjectID) == x), ])))
        values <- as.matrix(UReObs) %*% c(t(as.matrix(taum[x, ])))
        return(values)
      })

      # Analogous to bSub area, c(tauCont) will yield n*T x 1
      muObs <- UObs %*% as.matrix(phi[ , k]) + c(tauCont)
      likObs <- dbinom(D, size = 1, prob = pnorm(muObs), log = FALSE)
      prodlikObs <- as.vector(tapply(likObs, subjectID, FUN = prod))

      # Subject level contribution
      muSub <- matrix(rep(0, n*q), nrow = n, ncol = q)
      likSub <- sapply(1:nrow(taum), function(x) {
        mvtnorm::dmvnorm(taum[x, ], mean = muSub[x, ], sigma = as.matrix(Omega[ , , k]), log = FALSE)
      })
      likD[ , k] <- prodlikObs * likSub

    }


    # Log likelihood contribution from subject i for the response model
    # Assume that f(M_{1i}, M_{1i},...M_{Jstar,i} | D = 1, C; ...) = f(M_{1i} | D = 1, C; ...)f(M_{2i} | D = 1, C; ...)...f(M_{Jstar,i} | D = 1, C; ...)
    if (modelResponse == TRUE) {

      #Jstar <- ncol(M)
      Jstar <- length(Mvec)

      ### Running product
      likMk <- 1

      for (j in 1:Jstar) {

        # For a given latent class, get the log likelihood from outcome j at the observation level. Given theta, we can sum on the log scale

        # Get the random effect contribution for j
        kappajm <- matrix(kappa[[j]], nrow = n, ncol = q, byrow = TRUE)
        kappajCont <- unlist(sapply(unique(subjectIDM), function(x) {
          selObs <- which(subjectIDM == x)
          VReObs <- as.matrix(VRe)[selObs, ]
          #values <- as.matrix(VReObs) %*% c(t(as.matrix(kappajm[which(unique(subjectIDM) == x), ])))
          values <- as.matrix(VReObs) %*% c(t(as.matrix(kappajm[x, ])))
          return(values)
        }))

        # Observation level likelihood
        muObsj <- VObs %*% lambda[[j]][ , k] + c(kappajCont)
        likObsj <- dbinom(M[ , Mvec[j]], size = 1, prob = pnorm(muObsj), log = FALSE)
        prodlikObsj <- as.vector(tapply(likObsj, subjectIDM, FUN = prod))

        # Subject level likelihood
        muSubj <- matrix(rep(0, n*q), nrow = n, ncol = q)
        likSubj <- sapply(1:nrow(kappajm), function(x) {
          mvtnorm::dmvnorm(kappajm[x, ], mean = muSubj[x, ], sigma = as.matrix(Theta[[j]][ , , k]), log = FALSE)
        })

        # Log likelihood contribution for response j
        likMj <- prodlikObsj * likSubj

        # Add above to running product over j responses in latent class k
        likMk <- likMk * likMj

      }

      likM[ , k] <- likMk

    }

  }


  # Numerator
  # priorPik has class K from undifferenced parameterization in position 1 so re-ordering here is unnecessary

  if (modelVisit == TRUE & modelResponse == TRUE) {
    pik_num <- priorPik * likY * likD * likM
  } else if (modelVisit == TRUE & modelResponse == FALSE) {
    pik_num <- priorPik * likY * likD
  } else if (modelVisit == FALSE & modelResponse == FALSE) {
    pik_num <- priorPik * likY
  }

  # Posterior probability of subject i belonging in each class
  pik <- t(apply(pik_num, 1, function(x) x / sum(x)))

  # Replace NA's with 1/K
  sel <- apply(pik, 1, function(x) sum(is.na(x)))
  pik[sel > 0, ] <- rep(rep(1 / K, times = K), times = sum(sel > 0))

  pik <- as.matrix(unname(pik))
  C <- Hmisc::rMultinom(pik, 1)

  list(C = C, pik = pik)

}
