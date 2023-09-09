#' @title Get the effective sample size for the BIC calculation.
#'
#' @description This calculation of the effective sample size is used in BIC2. The reference for the approach is in Jones, R. H. (2011). Bayesian information criterion for longitudinal and clustered data. Statistics in Medicine, 30(25), 3050-3056. I extended the approach to the setting of multiple longitudinal health outcomes.
#' @export
#' @return Effective sample size at the current iteration.
get_ESS <- function(priorPik, Psi, Sigma, Y, XRe, subjectIDY) {

  q <- ncol(as.matrix(XRe))
  K <- ncol(priorPik)
  J <- ncol(Y)

  subjectIDYSub <- unique(subjectIDY)

  ### What is the sample size for each subject based on marginalizing over latent class and random effects
  ESS_sub <- sapply(unique(subjectIDY), function(x) {

    ni <- length(which(subjectIDY == x))
    # Y measurements
    Yi <- matrix(Y[which(subjectIDY == x), ], ncol = J, nrow = ni)
    # Random effects design matrix
    XRei <- matrix(XRe[which(subjectIDY == x), ], ncol = q, nrow = ni)

    # Prior probability vector for x^th subject
    priorPiki <- matrix(priorPik[which(subjectIDYSub == x), ], nrow = K, ncol = 1)

    margRi <- as.list(1:J)

    for (j in 1:J) {

      # Marginalize between subject variance-covariance computation
      Psiji <- Reduce("+", sapply(1:K, function(u) {
        priorPiki[u] * Psi[[j]][ , , u]
      }, simplify = FALSE))

      Sigmaji <- sum(Sigma[j, j, ] * priorPiki)

      margRi[[j]] <- (XRei %*% Psiji %*% t(XRei)) + diag(Sigmaji, nrow = ni, ncol = ni)

    }

    # Takes a list of matrices and makes a block diagonal matrix
    margRi <- as.matrix(Matrix::bdiag(margRi))

    # Marginalized covariance between Y_j and Y_j'
    for (j in 1:(J - 1)) {

      for (l in (j + 1):J) {

        sigmajli <- sum(Sigma[j, l, ] * priorPiki)
        sigmajli_expand <- matrix(rep(sigmajli, ni*ni), nrow = ni, ncol = ni)

        selRow <- ((j - 1)*ni + 1):(j*ni)
        selCol <- ((l - 1)*ni + 1):(l*ni)

        margRi[selRow, selCol] <- margRi[selCol, selRow] <- sigmajli_expand

      }

    }

    sum(solve(cov2cor(margRi)))

  })

  sum(ESS_sub)

}
