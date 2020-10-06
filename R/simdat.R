#' @title Simulate data for \code{MVNYBinaryMiss}.
#'
#' @description This makes data generated from the shared parameter model in \code{MVNYBinaryMiss}.
#'
#' @param K A scalar for the assumed number of latent classes in the GMMs.
#' @param J A scalar for the number of longitudinal outcomes. \code{J} must be greater than 1.
#' @param data A data.frame with all analysis variables, except for the outcomes that will be generated, including \code{Y1,\dots,YJ} for the longitudinal health outcomes; \code{D} for the visit process; and \code{M1,\dots,MJ} for the response process for each health outcome given a clinic visit. The dataset must contain an integer variable for the patient ID of each longitudinal measurement called \code{subjectID}.
#' @param regf A list of formulas for the design matrices in each sub-model. The structure is detailed below:
#'
#' To improve Gibbs sampling properties, the multivariate model for \code{Y1,\dots,YJ} is updated based on a hierarchically-centered parameterisation. As such, the named formulas for the longitudinal health outcomes model are "YRe", which is the random effects design matrix that must include at least a column of 1's for a random intercept; "YObs", which is the observation-level fixed effects design matrix containing covariates that will not be associated with random effects; and, "YSub" which contains the formula for the random effects equations using subject-level covariates. Note that the formulas "YRe", "YObs", and "YSub" must be disjoint in the hierarchically-centered parameterisation.
#'
#' For the multinomial probit model of latent class membership, the named formula is "LatentClass", containing patient-level covariates that may be associated with probabilty of latent class membership.
#'
#' If a model for \code{D} is included, then \code{regf} contains named formulas "DObs" and "DRe" for the fixed effects design matrix and the random effects design matrix, respectively. "DRe" contains a subset of the columns of "DObs".
#'
#' If a model for any of \code{M1,\dots,MJ} is included, then \code{regf} contains named formulas "MObs" and "MRe" for the fixed effects design matrix and the random effects design matrix, respectively. "MRe" contains a subset of the columns of "MObs".
#'
#'Note that \code{Y1,\dots,YJ} are assumed to have the same design matrices, and the response processes \code{M1,\dots,MJ} are assumed to have the same design matrices.
#'
#' @param delta \code{K-1}-column matrix of true regression coefficients for the latent class membership model.
#' @param betaObs True regression coefficients for the design matrix specified by \code{regf} named formula "YObs". The object is a \code{J} element list in which each element is an \code{s} by \code{K} matrix, where \code{s} is the number of regression coefficients and \code{K} is the number of latent classes, where \code{s} is the number of regression coefficients and \code{K} is the number of latent classes.
#' @param betaSub True regression coefficients in the design matrix specified by \code{regf} named formula "YSub". The object is a \code{J} element list in which each element is an \code{p} by \code{q} by \code{K} array, where \code{p} is the number of regression coefficients, \code{q} is the number of random effects, and \code{K} is the number of latent classes.
#' @param Psi True values for the variance-covariance matrix of the subject-level random effects in the longitudinal health outcomes model. The object is a \code{J} element list in which each element is an \code{q} by \code{q} by \code{K} array, where \code{q} is the number of random effects, and \code{K} is the number of latent classes.
#' @param Sigma True variance-covariance matrix of the longitudinal health outcomes \code{Y1,\dots,YJ}. The object is a \code{J} by \code{J} by \code{K} array.
#' @param phi A \code{K}-column matrix of true regression coefficients in the design matrix specified by \code{regf} named formula "DObs".
#' @param Omega True variance-covariance matrix of the subject-level random effects in the visit process model. The object is an \code{q} by \code{q} by \code{K} array, where \code{q} is the number of random effects, and \code{K} is the number of latent classes.
#' @param responseMiss Logical. If \code{TRUE}, then response process missingness will be generated for the columns specified in \code{Mvec}.
#' @param lambda True regression coefficients for the design matrix specified by \code{regf} named formula "MObs". The object is a list of the same number of elements as \code{Mvec}. The order of the initial values must correspond to \code{Mvec}. If \code{responseMiss = FALSE}, then use \code{NULL}.
#' @param Theta True variance-covariance matrix of the subject-level random effects in the response process model. The object is a list of the same number of elements as \code{Mvec}. The order of the initial values must correspond to \code{Mvec}. Each element is a \code{q} by \code{q} by \code{K} array. If \code{responseMiss = FALSE}, then use \code{NULL}.
#' @param Mvec Integer-valued vector indicating which response process given a clinic visit \code{M1,\dots,MJ} will be modeled under a missing not at random mechanism. For example, if the response processes of \code{Y1} and \code{Y2}, which are \code{M1} and \code{M2}, are desired, then \code{Mvec} is \code{c(1,2)}.
#' @param full Logical. If \code{TRUE}, then the full data for \code{Y1,\dots,YJ} are returned, before introducing any missed visits or responses. If \code{FALSE}, the data with missed visits and missed responses are returned. In this case, the data may have fewer unique subjects than in the original \code{data} because some subjects may not have any observed clinci visits.
#' @export
#' @return A list with 5 elements. \code{data} returns the original data with \code{Y1,\dots,YJ} after introducing missing values; \code{YC1,\dots,YCJ} which do not have missing values; \code{D} which is a binary indicator for the response process; \code{M1,\dots,MJ} which is the response process corresponding to each longitudinal outcome; \code{PrV}, which is the probability of an observed clinic visit; and \code{PrM1,\dots,PrMJ}, which is the probability of a response given a clinic visit. \code{D} and \code{M1,\dots,MJ} are coded 1 if a visit (response) is observed, and 0 otherwise. \code{M1,\dots,MJ} will be \code{NA} when \code{D=0}. \code{C} is the true latent class membership for each subject. \code{ame} are the true population-averaged regression coefficients in the longitudinal health outcomes model. \code{weight} are the true latent class weights. \code{store_attributes} is a matrix with the attributes of the data generation.
simdat <- function(K, J, data, regf, delta, betaObs, betaSub, Psi, Sigma, phi, Omega, responseMiss, lambda, Theta, Mvec, full){

  ### Design matrix for latent class model
  W <- model.matrix(regf[["LatentClass"]], data = aggregate(data, by = list(data$subjectID), FUN = tail, n = 1))

  ### Generate latent class memberships
  lc <- function(K, delta, W){

    # Number of subjects
    n <- nrow(W)

    # Generate latent variable Z where covariance matrix is fixed for identifiability
    if (K > 2) {

      Z <- mnormt::rmnorm(n, varcov = diag(K - 1)) + W %*% delta
      # Categorical variable based on max value of Z
      C <- rep(NA, n)
      for (i in 1:n) {
        MAX <- max(Z[i, ]) < 0  #All Z[i, ] < 0
        if (MAX) {
          #if the max is less than 0, then set C = 0
          C[i] <- 0
        } else {
          #otherwise, set C equal to the index of the maximum
          C[i] <- which.max(Z[i, ])
        }
      }

      C <- C + 1

      # Get prior probabilities of class assignment
      #Construct original matrix R from Z based parameterization
      #We have been working with the differenced parameterization such that the coefficients represent a difference from baseline class K, where after a simple reparameterization, class K is switched to class 0 (the reference group moves from the Kth column to the first column)
      #If we assume the original parameterization is represented by R, then
      #Z_{i1} = R_{i1} - R_{iK}; Z_{i2} = R_{i2} - R_{iK}, and so forth
      #R_{iK} \sim N(0, 1) since \delta_K and \theta_K are fixed to 0
      #R is an n x K matrix, where the K^th column is that with \delta_k fixed to 0 (the K^th column is the reference column)

      priorPikReorder <- matrix(NA, nrow = n, ncol = K)

      R <- matrix(NA, nrow = n, ncol = K)
      R[ , K] <- rnorm(n, mean = 0, sd = 1)

      for (k in 1:(K - 1)) {
        R[ , k] <- Z[, k] + R[ , K]
      }

      # Covariance matrix
      VarCovP <- matrix(1, nrow = (K - 1), ncol = (K - 1))
      diag(VarCovP) <- 2

      # Values where distribution is evaluated
      vals <- matrix(0, nrow = n, ncol = (K - 1))
      for (k in 1:K) {
        mur_diffk <- R[ , -k] - R[ , k]
        priorPikReorder[ , k] <- mnormt::pmnorm(vals, mean = mur_diffk, varcov = VarCovP)
      }

      priorPikReorder <- as.matrix(unname(priorPikReorder))
      priorPik <- cbind(priorPikReorder[ , K], priorPikReorder[ , 1:(K - 1)])

    } else if (K == 2) {

      # delta is an m x 1 vector of regression coefficients
      Z <- mnormt::rmnorm(n, varcov = 1) + W %*% delta
      C <- ifelse(Z > 0, 1, 0)
      C <- C + 1

      ### Get prior class membership probabilities
      priorPik <- matrix(NA, nrow = n, ncol = K)
      priorPik[ , K] <- pnorm(Z)
      priorPik[ , K - 1] <- 1 - pnorm(Z)

    }

    list(C = C, priorPik = priorPik)

  }

  lc_res <- lc(K = K, delta = delta, W = W)
  C <- lc_res$C

  ### Background information for longitudinal data
  subjectID <- data[["subjectID"]]
  C_expand <- C[factor(subjectID)]
  N <- length(subjectID)
  n <- length(unique(subjectID))


  ### Generate longitudinal outcomes Y
  # XSub will include an intercept term and V will include an intercept term for the random intercept
  XSub <- model.matrix(regf[["YSub"]], data = aggregate(data, by = list(subjectID), FUN = tail, n = 1))
  XRe <- model.matrix(regf[["YRe"]], data = data)
  XObs <- model.matrix(regf[["YObs"]], data = data)
  p <- ncol(XSub)
  q <- ncol(XRe)

  Y <- matrix(NA, nrow = N, ncol = J)

  for (k in 1:K) {

    # Observation level information
    ind_obs <- which(C_expand == k)
    subjectIDk <- subjectID[ind_obs]
    XRek <- as.matrix(XRe[ind_obs, ])
    XObsk <- as.matrix(XObs[ind_obs, ])
    Nk <- length(ind_obs)

    # Subject level information
    ind_sub <- which(C == k)
    XSubk <- as.matrix(XSub[ind_sub, ])

    # Generate Y_j
    #A J element list with the latent class specific means of the random effects distribution (nk x q)
    mubSubk <- lapply(betaSub, function(x) XSubk %*% x[ , , k])
    #Storage linear predictor for latent class k
    lp <- matrix(NA, Nk, J)
    for (j in 1:J) {

      # Generate stacked random effects for subjects in class k, nk*q vector
      stackbSub <- c(t(apply(mubSubk[[j]], 1, function(x) {mnormt::rmnorm(1, mean = x, varcov = Psi[[j]][ , , k])})))
      sp_n_sub <- lapply(split(as.data.frame(XRek), subjectIDk, drop = TRUE), as.matrix)
      convert_n_sub <- as.matrix(Matrix::bdiag(sp_n_sub))

      # Calculate linear predictor
      lp[ , j] <- as.vector(convert_n_sub %*% stackbSub + XObsk %*% betaObs[[j]][ , k])

    }

    # Generate Y using latent class specific variance-covariate of Y_js, making sure to place values based on ind_obs
    Y[ind_obs, ] <- t(apply(lp, 1, function(x) {mnormt::rmnorm(1, mean = x, varcov = Sigma[ , , k])}))

  }

  #------------------------ Generate visit process
  # Design matrices for D, includes both observation and subject level covariates
  UObs <- model.matrix(regf[["DObs"]], data = data)
  URe <- model.matrix(regf[["DRe"]], data = data)
  q <- ncol(URe)

  D <- rep(NA, N)
  prV <- rep(NA, N)

  for (k in 1:K) {

    UObsk <- UObs[C_expand == k, ]
    URek <- URe[C_expand == k, ]
    subjectIDk <- subjectID[C_expand == k]
    Nk <- length(subjectIDk)
    nk <- length(unique(subjectIDk))

    stacktau <- c(t(mnormt::rmnorm(nk, mean = rep(0, q), varcov = as.matrix(Omega[, , k]))))
    sp_n_sub <- lapply(split(as.data.frame(URek), subjectIDk, drop = TRUE), as.matrix)
    convert_n_sub <- as.matrix(Matrix::bdiag(sp_n_sub))

    # Calculate linear predictor
    lp <- as.vector(convert_n_sub %*% stacktau + UObsk %*% phi[ , k])
    D[C_expand == k] <- rbinom(Nk, 1, prob = pnorm(lp))
    prV[C_expand == k] <- as.vector(pnorm(lp))

  }


  #--------------------- Generate response processes given a visit
  # If responseMiss == TRUE, then generate missingness in response given clinic visit according to MNAR process with choice of MAR missingness in addition

  if (responseMiss == TRUE) {

    # Design matrices for D, includes both observation and subject level covariates
    VObs <- model.matrix(regf[["MObs"]], data = data[D == 1, ])
    VRe <- model.matrix(regf[["MRe"]], data = data[D == 1, ])
    q <- ncol(VRe)

    subjectIDVisit <- subjectID[D == 1]
    NVisit <- length(subjectIDVisit)
    C_expandVisit <- C_expand[D == 1]

    # Jstar is the number of Y_j with missing responses conditional on having a visit, for which we will have a model M_j
    # Keep M the same dimension of N
    # Jstar is the number of Y_j that will have model for M_j
    Jstar <- length(Mvec)
    # the M_j for which there will be *no* MNAR response process
    notMvec <- which(!c(1:J) %in% Mvec)

    M <- matrix(NA, nrow = N, ncol = J)
    prM <- matrix(NA, nrow = N, ncol = J)

    # First handle MNAR response process M_j where j in Mvec
    for (j in 1:Jstar) {

      Mtemp <- rep(NA, NVisit)
      prMtemp <- rep(NA, NVisit)

      for (k in 1:K) {

        VObsk <- VObs[C_expandVisit == k, ]
        VRek <- VRe[C_expandVisit == k, ]
        subjectIDVisitk <- subjectIDVisit[C_expandVisit == k]
        NVisitk <- length(subjectIDVisitk)
        nVisitk <- length(unique(subjectIDVisitk))

        stackkappa <- c(t(mnormt::rmnorm(nVisitk, mean = rep(0, q), varcov = Theta[[j]][ , , k])))
        sp_n_sub <- lapply(split(as.data.frame(VRek), subjectIDVisitk, drop = TRUE), as.matrix)
        convert_n_sub <- as.matrix(Matrix::bdiag(sp_n_sub))

        # Calculate linear predictor
        lp <- as.vector(convert_n_sub %*% stackkappa + VObsk %*% lambda[[j]][ , k])
        Mtemp[C_expandVisit == k] <- rbinom(NVisitk, 1, prob = pnorm(lp))
        prMtemp[C_expandVisit == k] <- as.vector(pnorm(lp))

      }

      M[D == 1, Mvec[j]] <- Mtemp
      prM[D == 1, Mvec[j]] <- prMtemp

    } # End of Jstar loop

    # For columns not in Mvec, then we assume M_j = 1 for j not in Mvec (i.e. these Y_j are fully observed)
    for (j in 1:length(notMvec)) {
      M[D == 1, notMvec[j]] <- 1
    }


  } # End of responseMiss == TRUE

  #------------------------ Introduce missingness in Y
  ### Save a full data copy
  YC <- Y

  ### Visit process
  Y[D == 0, ] <- NA

  if (responseMiss == TRUE) {

    ### Response process given a visit
    #MNAR response process in Mvec
    for (j in 1:Jstar) {
      Mj <- M[ , Mvec[j]]
      Y[D == 1 & Mj == 0, Mvec[j]] <- NA
    }

    #For response processes in notMvec
    for (j in 1:length(notMvec)) {
      Mj <- M[ , notMvec[j]]
      Y[D == 1 & Mj == 0, notMvec[j]] <- NA
    }

  } else if (responseMiss == FALSE) {

    M <- matrix(rep(1, N*J), nrow = N, ncol = J)
    M[D == 0, ] <- NA
    prM <- matrix(rep(1, N*J), nrow = N, ncol = J)
    prM[D == 0] <- NA

  }

  #------------------------- Prepare the simulated dataset for return

  ### If not the full data analysis, return the dataset after dropping out subjects who were not observed at all over follow-up, i.e., D_{ij}==0 for j = 1,\dots,J

  if (full == FALSE) {

    #------------------------ Place Y, D, M in a data.frame
    colnames(Y) <- paste("Y", 1:ncol(Y), sep = "")
    colnames(YC) <- paste("YC", 1:ncol(YC), sep = "")
    colnames(M) <- paste("M", 1:J, sep = "")
    colnames(prM) <- paste("prM", 1:J, sep = "")

    data <- cbind(data, Y, YC, D, M, prV, prM)

    #----------------------- Here, complete data (i.e., the benchmark analysis) does not include subjects were not observed at all over follow-up
    ### Identify subjects who were not observed at all over follow-up because of 0 visits
    #! Note subjectID will be ordered below, if an integer this is not problematic

    countVisit <- as.vector(tapply(D, subjectID, function(x){sum(x)}))
    cat("Number of subjects to be removed due to subjects with 0 visits:", sum(countVisit == 0), "\n")

    ### Adjust the true model parameters and complete data to account for the subjects who were not observed at all over follow-up
    # First, restrict priorPik and C to subjects with at least 1 visit
    priorPik <- lc_res$priorPik[countVisit > 0, ]
    C <- C[countVisit > 0]

    ### Second restrict data frame to observations from subjects with at least one visit
    # Remove observations from subjects with 0 visits from the data, including the complete data and the data with missingness
    countVisit_expand <- countVisit[factor(subjectID)]
    data <- subset(data, subset = countVisit_expand > 0)
    cat("Number of observations removed due to subjects with 0 visits:", sum(countVisit_expand == 0), "\n")

  } else if (full == TRUE) {

    ### If full == TRUE, then we return the dataset that includes records for subjects who did not have any follow-up visits
    colnames(Y) <- paste("Y", 1:ncol(Y), sep = "")
    colnames(YC) <- paste("YC", 1:ncol(YC), sep = "")
    colnames(M) <- paste("M", 1:J, sep = "")
    colnames(prM) <- paste("prM", 1:J, sep = "")

    data <- cbind(data, Y, YC, D, M, prV, prM)

    ### Latent class membership probabilities
    priorPik <- lc_res$priorPik

  }

  #----------------- True AME and class weights based on returned dataset

  ### Calculate average marginal parameters using true latent class membership weights
  ame <- as.list(1:J)
  for (j in 1:J) {
    betaObsj <- betaObs[[j]]
    # betaSub[[j]] organized with columns indexing latent class; first p rows are for first random effects equation, second p are for second and so forth
    betaSubj <- matrix(sapply(betaSub[[j]], abind::abind), ncol = K, nrow = p*q)
    ame[[j]] <- c(apply(priorPik %*% t(betaSubj), 2, mean), apply(priorPik %*% t(betaObsj), 2, mean))

  }

  ### Save true weights based on priorPik
  weight <- as.vector(apply(priorPik, 2, mean))

  #----------------------- Summarize attributes of simulated data
  # Number of subjects, Number of observations, proportion in each latent class, overall proportion missed visits, overall proportion missed response in M_1, etc. proportion missed visit by latent class, proportion missed response in M_1 by latent class, and so forth
  store_attributes <- rbind(c(sum(as.vector(table(C))), rep(NA, K - length(sum(as.vector(table(C)))))),
                            c(dim(data)[1], rep(NA, K - length(dim(data)[1]))),
                            as.vector(table(C)) / sum(as.vector(table(C))),
                            c(1 - mean(data$D), rep(NA, K - length(1 - mean(data$D)))))

  #Each M_j
  temp <- NULL
  for (j in 1:J) {
    Mj <- data[ , paste("M", j, sep = "")]
    temp <- rbind(temp, c(1 - mean(Mj[data$D == 1]), rep(NA, K - length(1 - mean(Mj[data$D == 1])))))
  }
  store_attributes <- rbind(store_attributes, temp)


  # Missingness by latent class
  #Visit
  store_attributes <- rbind(store_attributes, as.vector(table(data$D, C[factor(data$subjectID)])[1, ]) / as.vector(apply(as.matrix(table(data$D, C[factor(data$subjectID)])), 2, sum)))

  #Response
  temp <- NULL
  for (j in 1:J) {
    Mj <- paste("M", j, sep = "")
    dataTemp <- data[ , c("D", Mj)]
    dataTemp$C <- C[factor(data$subjectID)]
    dataTemp <- subset(dataTemp, subset = D == 1)
    temp <- rbind(temp, as.vector(table(dataTemp[[Mj]], dataTemp$C)[1, ]) / as.vector(apply(as.matrix(table(dataTemp[[Mj]], dataTemp$C)), 2, sum)))

  }

  store_attributes <- rbind(store_attributes, temp)

  #----------------- Return information
  list(data = data, C = C, ame = ame, weight = weight, store_attributes = store_attributes)


}
