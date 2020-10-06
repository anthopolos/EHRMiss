#' @title Update latent class specific regression coefficients \code{betaSub}.
#'
#' @description An multivariate normal prior distribution is used to update the regression coefficients \code{betaSub} corresponding to the covariates in \code{XSub},  for each outcome \code{Y1,\dots,YJ}.
#'
#' @return A \code{J} element list in which each element is a a \code{p} (number of covariates) by \code{q} (number of random effects) by \code{K} (number of latent classes) array of latent class-specific regression coefficients associated with \code{XSub}.
update_betaSubLC <- function(C, bSub, Psi, Y, XSub, prior.mu, prior.Sigma) {

  K <- length(table(C))

  ### Assumes the same XSub and q for each j
  J <- ncol(Y)
  n <- nrow(XSub)
  p <- ncol(XSub)
  q <- length(bSub[[1]]) / n

  values <- as.list(1:J)

  for (j in 1:J) {

    bSubjm <- matrix(bSub[[j]], nrow = n, ncol = q, byrow = TRUE)

    ### Place holder for latent class specific regression coefficients for outcome j
    valuesTemp <- array(NA, dim = c(p, q, K))

    for (k in 1:K) {

      ind_sub <- which(C == k)
      bSubjmk <- as.matrix(bSubjm[ind_sub, ])
      XSubk <- as.matrix(XSub[ind_sub, ])
      XX <- crossprod(XSubk, XSubk)

      Psik <- as.vector(diag(as.matrix(Psi[[j]][ , , k])))

      marginal_posterior <- sapply(1:length(Psik), FUN = function(s){
        # Posterior variance for beta in the q^{th} random effect equation
        post_var <- solve(solve(prior.Sigma) + XX / Psik[s])

        # Posterior mean
        #Likelihood contribution
        lik <- (t(XSubk) %*% (bSubjmk[ , s])) / Psik[s]   # p x n \times n x 1
        #Prior contribution
        prior <- solve(prior.Sigma) %*% prior.mu # p x p \times p x 1, contribution of prior

        post_mean <- post_var %*% (as.matrix(lik + prior))

        betaSubj <- mnormt::rmnorm(1, mean = post_mean, varcov = post_var)

        return(betaSubj)

      }, simplify = FALSE)

      valuesTemp[ , , k] <- matrix(unlist(marginal_posterior), nrow = p, ncol = q, byrow = FALSE)

    }

    values[[j]] <- valuesTemp

  }

  return(values)

}
