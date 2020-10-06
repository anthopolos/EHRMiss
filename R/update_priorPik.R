#' @title Update prior probabilities for the latent class membership model.
#'
#' @description This function calculates the multinomial or binary probabilities based on the assumption of independence among latent classes.
#' @return An \code{K} matrix of prior probabilities.
update_priorPik <- function(Z) {

  K <- ncol(as.matrix(Z)) + 1
  n <- nrow(as.matrix(Z))

  if (K > 2) {
    # Prior probability of class assignment
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
    # Prior probability of class assignment
    priorPik <- matrix(NA, nrow = n, ncol = K)

    priorPik[ , K] <- pnorm(Z)
    priorPik[ , 1] <- 1 - priorPik[ , K]

  }

  return(priorPik)

}
