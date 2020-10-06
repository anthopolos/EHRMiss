#' @title Update random effects conditional on latent class membership in a probit model.
#'
#' @description The random effects in a probit model are updated conditional on class membership. A normal prior distribution is assumed.
#' @return A stacked \code{n}\code{q} by 1 vector of random effects for each subject.
update_bLC_ProbitRE_fast <- function(C, Omega, Z, phi, UObs, URe, subjectID) {

  K <- length(table(C))
  n <- length(C)
  C_expand <- C[factor(subjectID)]

  tempValues <- matrix(NA, nrow = n, ncol = ncol(URe))

  for (k in 1:K) {

    ### Subject level information
    ind_sub <- which(C == k)

    ### Observation level information
    ind_obs <- which(C_expand == k)
    URek <- URe[ind_obs, ]
    UObsk <- UObs[ind_obs, ]
    subjectIDk <- subjectID[ind_obs]

    ### Parts for Posterior variance of tau
    Omegak <- Omega[ , , k]

    ### Parts for posterior mean of tau, prior mean is 0 so cancels out
    likTempk <- Z[ind_obs] - UObsk %*% phi[ , k]

    ### Update tau for each subject
    tempValues[ind_sub, ] <- t(sapply(unique(subjectIDk), function(x) {

      selObs <- which(subjectIDk == x)

      # Posterior variance
      URekObs <- as.matrix(URek)[selObs, ]
      UU <- crossprod(URekObs, URekObs)
      post_var <- solve(solve(Omegak) + UU)

      # Posterior mean
      likTempkObs <- likTempk[selObs]
      lik <- t(URekObs) %*% (likTempkObs)
      post_mean <- post_var %*% lik

      # Draw values for the ith subject
      values <- mnormt::rmnorm(1, mean = post_mean, varcov = post_var)

    }
    ))

  } # End of k loop

  ### Convert to a n*q stacked vector of random effects for each Y_j
  values <- c(t(tempValues))

  return(values) # nq x 1 stacked random effects for a J element list

}
