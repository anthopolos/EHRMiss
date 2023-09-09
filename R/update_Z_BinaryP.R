#' @title Update latent variable \code{Z} in a binary probit model.
#'
#' @description This uses a truncated normal distribution to update latent variable Z. This updating function must be used with \code{K} equals 2.
#' @export
#' @return An \code{n} by 1 vector latent variable \code{Z}.
update_Z_BinaryP <- function(Z, C, delta, W){

  n <- length(C)

  values <- rep(NA, n)

  mu <- W %*% delta

  values[C == 1] <- truncnorm::rtruncnorm(n = n, mean = mu, sd = 1, a = -Inf, b = 0)[C == 1]
  values[C == 2] <- truncnorm::rtruncnorm(n = n, mean = mu, sd = 1, a = 0, b = Inf)[C == 2]

  return(values)

}
