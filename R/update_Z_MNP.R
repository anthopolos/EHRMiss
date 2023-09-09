#' @title Update latent variable \code{Z} in multinomial probit model.
#'
#' @description This uses a truncated normal distribution to update latent variable Z. This updating function must be used with \code{K} greater than 2.
#' @export
#' @return An \code{n} by \code{K - 1} matrix of latent variable \code{Z}.
update_Z_MNP <- function(Z, C, delta, W){

  # Number of subjects
  n <- length(C)

  # Number of classes
  K <- ncol(delta) + 1

  # Transform C to multinomial binary matrix
  CMN <- matrix(0, nrow = n, ncol = K)
  for (i in 1:n) {
    for (k in 1:K) {
      CMN[i , C[i]] <- 1
    }
  }

  ### Lower, upper matrices for K - 1 columns
  lower <- upper <- matrix(NA, nrow = n, ncol = (K - 1))

  for (i in 1:n) {
    for (k in 2:K) {
      if (CMN[i, k] == 1) {
        lower[i, (k - 1)] <- max(c(0, Z[i, -(k - 1)]))
        upper[i, (k - 1)] <- Inf
      } else if (CMN[i, k] == 0) {
        lower[i, (k - 1)] <- -Inf
        upper[i, (k - 1)] <- max(c(0, Z[i, -(k - 1)]))
      }
    }
  }

  for (k in 1:(K - 1)) {

    muk <- W %*% delta[ , k]

    Z[ , k] <- truncnorm::rtruncnorm(n = n, mean = muk, sd = 1, a = lower[ , k], b = upper[ , k])

  }

  return(Z)

}
