#' @title Get the likelihood after integrating out latent class membership and the random effects.
#'
#' @description The likelihood after integrating out latent class membership and random effects is used in the calculation of DIC3. The DIC3 is the modified DIC proposed for the latent variable setting. The reference for the DIC3 is in Celeux, G., Forbes, F., Robert, C. P., & Titterington, D. M. (2006). Deviance Information Criteria for Missing Data Models. Bayesian Analysis, 1(4), 651-674. This likelihood is also used the calculation of BIC1 and BIC2. When \code{modelVisit = TRUE} or \code{modelResponse = TRUE}, numerical integration is used to estimate the marginal likelihoods. An option to estimate these densities using a Monte Carlo approach is available, but the option is turned off \code{monteCarlo = FALSE} because of long computation time.
#' @export
#' @return Observed data likelihood contribution of each of \code{n} subjects at the current iteration.
get_observed_lik_sub <- function(priorPik, betaObs, betaSub, Sigma, Psi, Y, XRe, XObs, XSub, subjectIDY, phi, Omega, D, URe, UObs, subjectID, lambda, Theta, M, VRe, VObs, subjectIDM, Mvec, modelVisit, modelResponse, monteCarlo = FALSE) {

  J <- ncol(Y)
  n <- length(unique(subjectIDY))
  K <- dim(priorPik)[2]
  q <- ncol(as.matrix(XRe))
  s <- ncol(as.matrix(XObs))
  p <- ncol(as.matrix(XSub))

  ### Longitudinal model for Y
  # Use a K loop for the likelihood contribution by latent class
  # Storage for likelihood contribution of each subject to each class
  llik_y_pik <- matrix(NA, nrow = n, ncol = K)

  for (k in 1:K) {

    Sigmak <- Sigma[ , , k]

    ### Get each subject's observed data likelihood contribution, i.e., after integrating out the random effects
    # The variance-covariance matrix for subject i will be a J x ni \times J x ni matrix
    # Along the block diagonal are the variance covariance matrices as in univariate setting (ni \times ni)

    llik_y_pikTemp <- sapply(unique(subjectIDY), function(x) {

      selObs <- which(subjectIDY == x)

      ni <- length(selObs)
      # Y measurements
      Yi <- matrix(Y[selObs, ], ncol = J, nrow = ni)
      # Random effects design matrix
      XRei <- matrix(XRe[selObs, ], ncol = q, nrow = ni)
      # Observation level design matrix
      XObsi <- matrix(XObs[selObs, ], ncol = s, nrow = ni)
      # Expand subject level matrix
      XSubi <- matrix(rep(XSub[x, ], ni), ncol = p, nrow = ni, byrow = TRUE)

      # Storage for muki
      muki <- matrix(NA, ncol = J, nrow = ni)

      Rki <- as.list(1:J)
      for (j in 1:J) {

        # Transform this to save by class then outcome
        Psikj <- Psi[[j]][ , , k]

        Rki[[j]] <- (XRei %*% Psikj %*% t(XRei)) + diag(Sigmak[j, j], nrow = ni, ncol = ni)
        #! Will need to consider what to do if more than 1 subject level covariate wrt interaction with random effects
        muki[ , j] <- XObsi %*% betaObs[[j]][ , k] + XSubi %*% betaSub[[j]][ , , k]

      }

      ### Takes a list of matrices and makes a block diagonal matrix
      Rki <- as.matrix(Matrix::bdiag(Rki))

      for (j in 1:(J - 1)) {

        for (l in (j + 1):J) {

          sigmajl <- matrix(rep(Sigmak[j, l], ni*ni), nrow = ni, ncol = ni)
          selRow <- ((j - 1)*ni + 1):(j*ni)
          selCol <- ((l - 1)*ni + 1):(l*ni)
          Rki[selRow, selCol] <- Rki[selCol, selRow] <- sigmajl

        }

      }

      mvtnorm::dmvnorm(c(Yi), mean = c(muki), sigma = Rki, log = FALSE)

    })

    llik_y_pik[ , k] <- llik_y_pikTemp

  } # End of k loop

  ### Visit process model

  if (modelVisit == TRUE) {

    # Storage for the likelihood contribution of each subject to latent class k
    llik_d_pik <- matrix(NA, nrow = n, ncol = K)

    if (monteCarlo == TRUE) {

      # Number of iterations to integrate over REs in probit setting
      B <- 100

      # To integrate over the random effects, at MCMC iteration l, draw repeat this process 500 times
      #draw n random effects
      #compute density

      llik_d_pikTemp <- matrix(NA, nrow = B, ncol = length(subjectID))

      for (k in 1:K) {

        for (l in 1:B) {

          stacktauDraw <- c(t(mnormt::rmnorm(n, mean = rep(0, q), varcov = as.matrix(Omega[, , k]))))

          sp_n_sub <- lapply(split(as.data.frame(URe), subjectID, drop = TRUE), as.matrix)
          convert_n_sub <- as.matrix(Matrix::bdiag(sp_n_sub))

          # Calculate linear predictor
          lp <- as.vector(convert_n_sub %*% stacktauDraw + UObs %*% phi[ , k])

          # Store for latent class K, B \times N
          llik_d_pikTemp[l , ] <- dbinom(D, size = 1, prob = pnorm(lp))

        }

        # Each subjects contribution to each latent class
        #Evaluate \int_\tau \prod_{j=1}^J f(d_{ij} \, | \, \tau_i) f(\tau_i) \partial \tau_i
        #For each subject, we need to compute the approximation \frac{1}{B} \sum_{b=1}^B \prod_{j=1}^J f(d_{ij} \, | \, \phi_k^l, \omega_k^l,\tau_i^{l,b})
        #Thus, llik_d_pikTemp is B by N, then becomes n by B using the first apply, then the row wise average should be for each subject
        llik_d_pik[ , k] <- apply(apply(llik_d_pikTemp, 1, function(x) {tapply(x, subjectID, prod)}), 1, mean)

      } # End of K loop

    } # End of monteCarlo = TRUE

    if (monteCarlo == FALSE) {

      # For each latent class, we use numerical integration to integrate over the random effects
      #! currently only supports q = 1

      for (k in 1:K) {

        store_num_int <- sapply(1:length(subjectID), function(r) {

          # Function for numberical integration over random effects in probit regression
          get_num_int <- function(x) {

            pnorm(c(UObs[r, ] %*% phi[ , k]) + x)^D[r] * (1 - pnorm(c(UObs[r, ] %*% phi[ , k]) + x))^(1 - D[r]) * (1 / sqrt(2 * pi * Omega[ , , k])) * exp(-x^2 / (2 * Omega[ , , k]))

          }

          vals <- seq(-5, 5, length = 1000)
          dist <- vals[-1] - vals[-length(vals)]

          sum(get_num_int(x = vals[-1]) * dist)

        })

        llik_d_pik[ , k] <- tapply(store_num_int , subjectID, prod)

      } # End of the K loop

    } # end of MonteCarlo = FALSE

  } # Close of model visit = TRUE



  ### Response process model
  #! Only supports q = 1 for the time being
  if (modelResponse == TRUE) {

    Jstar <- length(Mvec)
    # Storage for the latent class contributions of each subject
    llik_m_pik <- matrix(NA, nrow = n, ncol = K)

    if (monteCarlo == TRUE) {

      for (k in 1:K) {

        ### Running product over j in Jstar for the contribution to class k
        likMk <- 1

        for (j in 1:Jstar) {

          llik_m_pikTemp <- matrix(NA, nrow = B, ncol = length(subjectIDM))

          for (l in 1:B) {

            stackkappaDraw <- c(t(mnormt::rmnorm(n, mean = rep(0, q), varcov = as.matrix(Theta[[j]][, , k]))))

            sp_n_sub <- lapply(split(as.data.frame(VRe), subjectIDM, drop = TRUE), as.matrix)
            convert_n_sub <- as.matrix(Matrix::bdiag(sp_n_sub))

            # Calculate linear predictor
            lp <- as.vector(convert_n_sub %*% stackkappaDraw + VObs %*% lambda[[j]][ , k])

            # Store for latent class K, B \times N
            llik_m_pikTemp[l , ] <- dbinom(M[ , Mvec[j]], size = 1, prob = pnorm(lp))
          } # End fo B loop

          # M_j's contribution to latent class k
          likMk <- likMk * apply(apply(llik_m_pikTemp, 1, function(x) {tapply(x, subjectIDM, prod)}), 1, mean)

        } # End of the Jstar loop

        llik_m_pik[ , k] <- likMk

      } # End of K loop

    } # End of Monte Carlo == TRUE

    if (monteCarlo == FALSE) {

      for (k in 1:K) {

        ### Running product over j in Jstar for the contribution to class k
        likMk <- 1

        for (j in 1:Jstar) {

          store_num_int <- sapply(1:length(subjectIDM), function(r) {

            # Function for numberical integration over random effects in probit regression
            get_num_int <- function(x) {

              pnorm(c(VObs[r, ] %*% lambda[[j]][ , k]) + x)^M[ , Mvec[j]][r] * (1 - pnorm(c(VObs[r, ] %*% lambda[[j]][ , k]) + x))^(1 - M[ , Mvec[j]][r]) * (1 / sqrt(2 * pi * Theta[[j]][ , , k])) * exp(-x^2 / (2 * Theta[[j]][ , , k]))

            }

            vals <- seq(-5, 5, length = 1000)
            dist <- vals[-1] - vals[-length(vals)]

            sum(get_num_int(x = vals[-1]) * dist)

          })

          likMk <- likMk * tapply(store_num_int , subjectIDM, prod)

        }  # End of Jstar loop

        llik_m_pik[ , k] <- likMk

      } # End of K loop

    } # End of Monte Carlo == FALSE

  } # End of modelResponse = TRUE


  # Compute \sum_{i = 1}^n log \sum_{k = 1}^K \pi_k f(Y | \beta_k) f(D | \phi_k) f(M | \lambda_k)
  if (modelVisit == TRUE & modelResponse == TRUE) {

    # n-length vector of likelihood contributions of each subject
    store_llikTemp <- apply(priorPik * llik_y_pik * llik_d_pik * llik_m_pik, 1, sum)

  } else if (modelVisit == TRUE & modelResponse == FALSE) {

    # n-length vector of likelihood contributions of each subject
    store_llikTemp <- apply(priorPik * llik_y_pik * llik_d_pik, 1, sum)

  } else if (modelVisit == FALSE & modelResponse == FALSE) {

    # n-length vector of likelihood contribution of each subject
    store_llikTemp <- apply(priorPik * llik_y_pik, 1, sum)

  }

  ### Log likelihood estimate
  return(store_llikTemp)

}
