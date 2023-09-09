#' @title Get posterior predictive checks including replicates from the joint model and a discrepancy measure computed using the completed data and the replicated completed data.
#'
#' @description This function generates posterior predictive draws of \code{Y1,\dots,YJ}, \code{D}, and \code{M1,\dots,MJ} based on redrawing all of the latent variables, including latent class memberships and random effects from each of the submodels for the longitudinal outcomes, visit process, and response process given a clinic visit. The posterior predictive draws are written to file. In addition, for posterior predictive checking based on the completed and replicated completed \code{Y1,\dots,YJ}, a discrepancy measure \code{T}, namely, the multivariate mean squared error, is calculated. References for using a discrepancy measure with completed and replicated completed data include: 1. Xu, D., Chatterjee, A. and Daniels, M. (2016) A Note on Posterior Predictive Checks to Assess Model Fit for Incomplete Data. Statistics in Medicine, 35, 5029-5039. 2. Gelman, A., Mechelen, I. V., Verbeke, G., Heitjan, D.F. and Meulders, M. (2005) Multiple Imputation for Model Checking: Completed-Data Plots with Missing and Latent Data. Biometrics, 61, 74-85.
#' @export
#' @return Draws of replicated longitudinal outcomes, visit process, and response process. Discrepancy measure computed using completed data and the replicated completed data.
get_model_checks <- function(K, C, priorPik, bSub, betaObs, betaSub, Psi, Sigma, Y, XRe, XObs, XSub, subjectIDY, Q, phi, Omega, D, URe, UObs, subjectID, L, lambda, Theta, M, VRe, VObs, subjectIDM, Mvec, modelVisit, modelResponse) {

  ### Posterior predictive draws based on redrawing all of the latent variables, including latent class and random effects
  #! x must be integer valued corresponding to rows of XSub

  # Longitudinal outcome model
  J <- ncol(as.matrix(Y))
  n <- dim(priorPik)[1]
  q <- ncol(as.matrix(XRe))
  N <- nrow(as.matrix(Y))

  # Complete data discrepancy measure TObs
  C_expandY <- C[factor(subjectIDY)]
  bSubm <- lapply(bSub, function(x) matrix(x, nrow = n, ncol = q, byrow = TRUE))

  # Replicated completed discrepany measure TRep based on redraw C per algorithm 3.7 (Fruwhirth)
  Cdraw <- Hmisc::rMultinom(priorPik, 1)
  Cdraw_expandY <- Cdraw[factor(subjectIDY)]

  # Redraw random effects for the longitudinal outcomes Y
  bSubdraw <- update_bSubLC_fast(C = Cdraw, bSub = bSub, betaObs = betaObs, betaSub = betaSub, Psi = Psi, Sigma = Sigma, Y = Y, XRe = XRe, XObs = XObs, XSub = XSub, subjectID = subjectIDY)
  bSubdrawm <- lapply(bSubdraw, function(x) matrix(x, nrow = n, ncol = q, byrow = TRUE))

  ### Store Ydraw for print to file
  store_Ydraw <- matrix(NA, nrow = N, ncol = J)


  ### Missing data processes
  if (modelVisit == TRUE) {
    # Redraw random effects for the visit process
    taudraw <- update_bLC_ProbitRE_fast(C = Cdraw, Omega = Omega, Z = Q, phi = phi, UObs = UObs, URe = URe, subjectID = subjectID)
    taudrawm <- matrix(taudraw, nrow = n, ncol = ncol(URe), byrow = TRUE)
    store_Ddraw <- matrix(NA, nrow = length(subjectID), ncol = 1)
  } else {
    store_Ddraw <- NULL
  }

  if (modelResponse == TRUE) {
    Jstar <- length(Mvec)

    kappadraw <- as.list(1:Jstar)
    for (j in 1:Jstar) {
      kappadraw[[j]] <- update_bLC_ProbitRE_fast(C = Cdraw, Omega = Theta[[j]], Z = L[ , j], phi = lambda[[j]], UObs = VObs, URe = VRe, subjectID = subjectIDM)
    }
    kappadrawm <- lapply(kappadraw, function(x) matrix(x, nrow = n, ncol = q, byrow = TRUE))

    store_Mdraw <- matrix(NA, nrow = length(subjectIDM), ncol = Jstar)

  } else {
    store_Mdraw <- NULL
  }

  ### Completed and replicated completed TObs, TRep to sum across latent classes
  TObs <- 0
  TRep <- 0


  for (k in 1:K) {

    ### Compute complete data discrepancy measure TObs
    ind_obs <- which(C_expandY == k)

    Yk <- Y[ind_obs, ]

    subjectIDYk <- subjectIDY[ind_obs]
    XRek <- as.matrix(XRe[ind_obs, ])
    XObsk <- as.matrix(XObs[ind_obs, ])

    muY <- matrix(NA, ncol = J, nrow = length(subjectIDYk))

    for (j in 1:J) {

      # Calculation for bj
      bSubj <- unlist(sapply(unique(subjectIDYk), function(u) {
        selObs <- which(subjectIDYk == u)
        XRekObs <- as.matrix(XRek)[selObs, ]

        xk <- as.matrix(bSubm[[j]][u, ])
        mu <- as.matrix(XRekObs) %*% c(t(xk))
        return(mu)
      }))

      #!tempm <- matrix(bSub[[j]], nrow = n, ncol = q, byrow = TRUE)
      #!tempmk <- as.matrix(tempm[ind_sub, ])
      #!bSubj <- as.vector(convert_XRe_sub %*% c(t(tempmk)))

      muY[ , j] <- XObsk %*% betaObs[[j]][ , k] + bSubj

    }

    # For latent class k, compute the discrepancy measure Mahalanobis distance under the observed and replicated data
    # 1 x J \times J x J \times J x 1 for distance for observation x
    diffObsk <- Yk - muY
    TObsk <- sapply(1:nrow(diffObsk), function(x) {
      t(as.matrix(diffObsk[x, ])) %*% solve(Sigma[ , , k]) %*% as.matrix(diffObsk[x, ])
    })

    # sum(TObsk) is the sum of the squared of the distances
    TObs <- TObs + sum(TObsk)




    ### Compute replicated completed data TRep
    ind_obs <- which(Cdraw_expandY == k)

    Yk <- Y[ind_obs, ]

    ### Redraw Y
    subjectIDYk <- subjectIDY[ind_obs]
    XRek <- as.matrix(XRe[ind_obs, ])
    XObsk <- as.matrix(XObs[ind_obs, ])

    # Subject level information
    #!ind_sub <- which(Cdraw == k)

    # N by n*q design matrix for random effects, to be multiplied by nq by 1 vector of stacked random effects
    #!sp_XRe_sub <- lapply(split(as.data.frame(XRek), subjectIDYk, drop = TRUE), as.matrix)
    #!convert_XRe_sub <- as.matrix(Matrix::bdiag(sp_XRe_sub))

    muY <- matrix(NA, ncol = J, nrow = length(subjectIDYk))

    for (j in 1:J) {

      # Calculation for bj
      #tempm <- matrix(bSub[[j]], nrow = n, ncol = q, byrow = TRUE)
      #tempmk <- as.matrix(tempm[ind_sub, ])
      #bSubj <- as.vector(convert_XRe_sub %*% c(t(tempmk)))
      # Calculation for bj
      bSubdrawj <- unlist(sapply(unique(subjectIDYk), function(u) {
        selObs <- which(subjectIDYk == u)
        XRekObs <- as.matrix(XRek)[selObs, ]

        xk <- as.matrix(bSubdrawm[[j]][u, ])
        mu <- as.matrix(XRekObs) %*% c(t(xk))
        return(mu)
      }))

      muY[ , j] <- XObsk %*% betaObs[[j]][ , k] + bSubdrawj

    }

    # Generate Y using latent class specific variance-covariate of Y_js, making sure to place values based on ind_obs
    Ydrawk <- t(apply(muY, 1, function(x) {mnormt::rmnorm(1, mean = x, varcov = Sigma[ , , k])}))
    diffRepk <- Ydrawk - muY
    TRepk <- sapply(1:nrow(diffRepk), function(x) {
      t(as.matrix(diffRepk[x, ])) %*% solve(Sigma[ , , k]) %*% as.matrix(diffRepk[x, ])
    })

    TRep <- TRep + sum(TRepk)

    store_Ydraw[ind_obs, ] <- Ydrawk


    ### Redraw D if modeling the visit process
    if (modelVisit == TRUE) {

      Cdraw_expand <- Cdraw[factor(subjectID)]

      ind_obs <- which(Cdraw_expand == k)
      URek <- URe[ind_obs, ]
      UObsk <- UObs[ind_obs, ]
      subjectIDk <- subjectID[ind_obs]

      Rek <- unlist(sapply(unique(subjectIDk), function(u) {
        selObs <- which(subjectIDk == u)
        URekObs <- as.matrix(URek)[selObs, ]

        xk <- as.matrix(taudrawm[u, ])
        mu <- as.matrix(URekObs) %*% c(t(xk))
        return(mu)
      }))

      lpk <- UObsk %*% phi[ , k] + as.matrix(c(Rek))

      store_Ddraw[ind_obs, 1] <- rbinom(length(lpk), size = 1, prob = pnorm(lpk))

    }

    ### Redraw M if modeling the visit process
    if (modelResponse == TRUE) {

      Cdraw_expand <- Cdraw[factor(subjectIDM)]

      ind_obs <- which(Cdraw_expand == k)
      VRek <- VRe[ind_obs, ]
      VObsk <- VObs[ind_obs, ]
      subjectIDMk <- subjectIDM[ind_obs]

      for (j in 1:Jstar) {
        Rek <- unlist(sapply(unique(subjectIDMk), function(u) {
          selObs <- which(subjectIDMk == u)
          VRekObs <- as.matrix(VRek)[selObs, ]

          xk <- as.matrix(kappadrawm[[j]][u, ])
          mu <- as.matrix(VRekObs) %*% c(t(xk))
          return(mu)
        }))

        lpk <- VObsk %*% lambda[[j]][ , k] + as.matrix(c(Rek))

        store_Mdraw[ind_obs, j] <- rbinom(length(lpk), size = 1, prob = pnorm(lpk))
      }


    }


  }

  store_T_completed <- c(TObs, TRep)

  list(store_T_completed = store_T_completed, store_Ydraw = store_Ydraw, store_Ddraw = store_Ddraw, store_Mdraw = store_Mdraw)

}
