#' @title Get discrepancy measure computed using the completed data and the replicated completed data.
#'
#' @description For a discrepancy measure for the model of \code{Y1,\dots,YJ}, the multivariate mean squared error is used. References for using a discrepancy measure with completed and replicated completed data include: 1. Xu, D., Chatterjee, A. and Daniels, M. (2016) A Note on Posterior Predictive Checks to Assess Model Fit for Incomplete Data. Statistics in Medicine, 35, 5029-5039. 2. Gelman, A., Mechelen, I. V., Verbeke, G., Heitjan, D.F. and Meulders, M. (2005) Multiple Imputation for Model Checking: Completed-Data Plots with Missing and Latent Data. Biometrics, 61, 74-85.
#' @return Discrepancy measure computed using completed data and the replicated completed data.
get_discrepancy_completed <- function(K, C, priorPik, bSub, betaObs, Sigma, Y, XRe, XObs, subjectIDY) {

  J <- ncol(as.matrix(Y))
  n <- dim(priorPik)[1]
  q <- ncol(as.matrix(XRe))
  N <- nrow(as.matrix(Y))

  # "observed" C expand
  C_expand <- C[factor(subjectIDY)]

  # Redraw C per algorithm 3.7
  Cdraw <- Hmisc::rMultinom(priorPik, 1)
  Cdraw_expand <- Cdraw[factor(subjectIDY)]

  ### TObs, Rep to sum across latent classes
  TObs <- 0
  TRep <- 0

  for (k in 1:K) {

    ### Calculation for "observed C"
    ind_obs <- which(C_expand == k)

    Yk <- Y[ind_obs, ]
    subjectIDYk <- subjectIDY[ind_obs]
    XRek <- as.matrix(XRe[ind_obs, ])
    XObsk <- as.matrix(XObs[ind_obs, ])

    Nk <- length(subjectIDYk)

    # Subject level information
    ind_sub <- which(C == k)

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

    # For latent class k, compute the discrepancy measure Mahalanobis distance under the observed and replicated data
    # 1 x J \times J x J \times J x 1 for distance for observation x
    diffObsk <- Yk - muObs
    TObsk <- sapply(1:nrow(diffObsk), function(x) {
      t(as.matrix(diffObsk[x, ])) %*% solve(Sigma[ , , k]) %*% as.matrix(diffObsk[x, ])
    })
    # sum(TObsk) is the sum of the squared of the distances
    TObs <- TObs + sum(TObsk)

    ### Calculation for replicated C and Y
    ind_obs <- which(Cdraw_expand == k)
    Yk <- Y[ind_obs, ]
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
    Ydrawk <- t(apply(muObs, 1, function(x) {mnormt::rmnorm(1, mean = x, varcov = Sigma[ , , k])}))

    diffRepk <- Ydrawk - muObs
    TRepk <- sapply(1:nrow(diffRepk), function(x) {
      t(as.matrix(diffRepk[x, ])) %*% solve(Sigma[ , , k]) %*% as.matrix(diffRepk[x, ])
    })

    TRep <- TRep + sum(TRepk)


  }

  store_T_completed <- c(TObs, TRep)

  list(store_T_completed = store_T_completed)

}
