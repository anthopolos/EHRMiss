#' @title Update latent normal draws in a probit model.
#'
#' @description The latent normal draws for a probit model with or without a subject specific random effects are updated.
#' @export
#' @return An observation-level vector of latent normal draws.
update_ZLC_ProbitRE <- function(C, phi, tau, D, UObs, URe, subjectID){

  K <- length(table(C))
  n <- length(C)
  C_expand <- C[factor(subjectID)]

  taum <- matrix(tau, nrow = n, ncol = ncol(URe), byrow = TRUE)
  values <- rep(NA, length(D))

  for (k in 1:K) {

    ### Subject level
    ind_sub <- which(C == k)

    ### Observation level
    ind_obs <- which(C_expand == k)

    Dk <- D[ind_obs]
    UObsk <- UObs[C_expand == k, ]
    subjectIDk <- subjectID[C_expand == k]

    if (!is.null(tau)) {

      ### Design matrix and subgroup of RE for class k
      URek <- URe[C_expand == k, ]
      tauk <- taum[ind_sub, ]

      sp_n_sub <- lapply(split(as.data.frame(URek), subjectIDk, drop = TRUE), as.matrix)
      convert_n_sub <- as.matrix(Matrix::bdiag(sp_n_sub))

      mu <- UObsk %*% phi[ , k] + as.vector(convert_n_sub %*% c(t(tauk)))

    } else {

      mu <- UObsk %*% phi[ , k]

    }


    values[C_expand == k & D == 0] <- truncnorm::rtruncnorm(n = length(ind_obs), a = -Inf, b = 0, mean = mu, sd = 1)[Dk == 0]
    values[C_expand == k & D == 1] <- truncnorm::rtruncnorm(n = length(ind_obs), a = 0, b = Inf, mean = mu, sd = 1)[Dk == 1]

  }

  return(values)

}
