#' @title Update prior probabilities for the latent class membership model.
#'
#' @description This function calculates the multinomial or binary probabilities based on the assumption of independence among latent classes.
#' @export
#' @return An \code{K} matrix of prior probabilities.
update_priorPik <- function(Z) {

  K <- ncol(as.matrix(Z)) + 1
  n <- nrow(as.matrix(Z))

  if (K > 2) {

    ### To obtain the prior probabilities, use differenced Z.
    # K = 3 as an example
    # z^*_{i1},...,z^*_{i3}
    # Difference with respect to K=3: z_{i1} = z^*_{i1} - z^*_{i3}, z_{i2} = z^*_{i2} - z^*_{i3},  z_{i3} = z^*_{i3} - z^*_{i3}
    # If we take Pr(z_{i1} > z_{i3}, z_{i1} > z_{i2}), then we obtain:
    # Pr(z_{i1} - z_{i3} > 0, z_{i1} - z_{i2} > 0) = Pr([z_{i1}^* - z_{i3}^* - (z_{i3}^* - z_{i3}^*) > 0], [z_{i1}^* - z_{i3}^* - (z_{i2}^* - z_{i3}^*) > 0])
    # = Pr([z_{i1}^* - z_{i3}^* > 0], [z_{i1}^* - z_{i2}^* > 0]) #This is what we want for the multivariate normal probability
    # = Pr([z_{i3}^* - z_{i1}^* < 0], [z_{i2}^* - z_{i1}^* < 0])

    # Note: In the parametrization, below makes the reference class in column 1
    # C_i = 1 if max(z_{i1},z_{i2}) \leq 0; C_i = k if max(z_{i1},z_{i2}) = z_{ik} \geq 0

    # Variance covariance for differenced variable
    # Var([z^*_{i1} - z^*_{i3}]) = 2
    # Cov(z_{i1}^* - z_{i3}^*, z_{i2}^* - z_{i3}^*) = 1

    # Base covariance matrix
    VarCovP <- matrix(1, nrow = (K - 1), ncol = (K - 1))
    diag(VarCovP) <- 2

    # priorPik
    priorPik <- matrix(NA, nrow = n, ncol = K)

    # Values where distribution is evaluated
    mu_err <- matrix(0, nrow = n, ncol = (K - 1))
    val_diff <- matrix(cbind(Z, 0), nrow = n, ncol = K)

    # mu_diffK will be n by K-1
    for (k in 1:K) {

      val_diffk <- val_diff[ , -k] - val_diff[ , k]

      # k always takes position 1
      #if (k == 1) {
      #  if (!is.null(tau2)) {VarCovP <- VarCovP + diag(tau2, nrow = (K - 1), ncol = (K - 1))}
      #  if (!is.null(xi2)) {VarCovP <- VarCovP + diag(xi2, nrow = (K - 1), ncol = (K - 1))}
      #  if (!is.null(gamma2)) {VarCovP <- VarCovP + diag(gamma2, nrow = (K - 1), ncol = (K - 1))}
      #}
      #if (k > 1) {
      #  if (!is.null(tau2)) {VarCovP <- VarCovP + diag(c(tau2[(k - 1)], tau2[-(k - 1)]), nrow = (K - 1), ncol = (K - 1))}
      #  if (!is.null(xi2)) {VarCovP <- VarCovP + diag(c(xi2[(k - 1)], xi2[-(k - 1)]), nrow = (K - 1), ncol = (K - 1))}
      # if (!is.null(gamma2)) {VarCovP <- VarCovP + diag(c(gamma2[(k - 1)], gamma2[-(k - 1)]), nrow = (K - 1), ncol = (K - 1))}
      #}

      priorPik[ , k] <- mnormt::pmnorm(((-1) * val_diffk), mean = mu_err, varcov = VarCovP)

    }

    priorPik <- cbind(priorPik[ , K], priorPik[ , 1:(K - 1)])

  } else if (K ==2) {

    priorPik <- matrix(NA, nrow = n, ncol = K)
    priorPik[ , K] <- pnorm(Z, mean = 0, sd = 1)
    priorPik[ , 1] <- 1 - priorPik[ , K]

  }

  return(priorPik)

}
