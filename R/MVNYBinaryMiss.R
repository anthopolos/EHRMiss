#' @title Bayesian shared parameter model for missing data in EHRs.
#'
#' @description The model is a Bayesian shared parameter model to handle missing data in \code{J} longitudinal health outcomes collected in electronic health records (EHRs). The sub-models for the longitudinal health outcomes, visit process, and response process given a clinic visit are growth mixture models (GMMs). Longitudinal health outcomes \code{Y1,\dots,YJ} are assumed to follow a multivariate normal distribution. The visit process \code{D} is defined as a binary indicator equal to 1 if patient has a clinic visit and 0 otherwise. To correspond to \code{Y1,\dots,YJ}, \code{M1,\dots,MJ} are binary indicators equal to 1 if a response is observed given a clinic visit and 0 otherwise. The visit process and response process given a visit are modeled using a probit link function. A multinomial probit model is used to model the probability of a subject belonging to a latent class. Conditional on a patient's latent class membership \code{C}, the longitudinal health outcomes, visit process, and response process given a clinic visit are assumed to be independent.
#'
#' @param K A scalar for the assumed number of latent classes in the GMMs.
#' @param J A scalar for the number of longitudinal outcomes. \code{J} must be greater than 1.
#' @param data A data.frame with all analysis variables. The data.frame must include variables named \code{Y1,\dots,YJ} for the longitudinal health outcomes; \code{D} for the visit process; \code{M1,\dots,MJ} for the response process for each health outcome given a clinic visit; and a consecutive, integer-valued variable for the patient ID of each longitudinal measurement called \code{subjectID}. Note that the data.frame must contain only the \code{Y1,\dots,YJ} and \code{M1,\dots,MJ} to be used in model fitting.
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
#' @param imputeResponse Logical. If \code{TRUE}, then missing values indicated by the response process given a clinic visit \code{M1,\dots,MJ} are imputed. If \code{FALSE}, then no imputation is conducted. In this case, the analysis dataset will be restricted to records with all \code{Y1,\dots,YJ} observed, and models for the visit process and response process given a clinic visit are not allowed.
#' @param Mvec Integer-valued vector indicating which response process given a clinic visit \code{M1,\dots,MJ} will be modeled under a missing not at random mechanism. For example, if the response processes of \code{Y1} and \code{Y2}, which are \code{M1} and \code{M2}, are desired, then \code{Mvec} is \code{c(1,2)}.
#' @param modelVisit Logical. If \code{TRUE}, then a model for the visit process \code{D} is implemented.
#' @param modelResponse Logical. If \code{TRUE}, then a model for at least one of the response processes given a clinic visit \code{M1,\dots,MJ}, as indicated by \code{Mvec}, is implemented.
#' @param priors A nested list of prior parameters. All regression coefficients are assigned a normal prior distribution. All variance-covariances are assigned an inverse Wishart prior distribution. If any parameters are not needed in the model, the element in the list can be set to NULL. The prior distributions are not allowed to vary by latent class. Details are provided in the \code{Details} section.
#' @param inits A list of initial values. If any parameters are not needed in the model, the element in the list can be set to NULL. The initial values are allowed to vary by latent class. Details are provided in the \code{Details} section.
#' @param n.samples Number of MCMC iterations.
#' @param burn Number of MCMC samples to discard as a burn-in.
#' @param monitor Logical. If \code{TRUE}, then after the burn-in period finishes, trace plots of the first four regression coefficients associated with \code{regf} formula named "YSub" are displayed and are periodically updated according to the \code{update} parameter.
#' @param update A scalar for the interval at which to print the iteration number and running size of the latent classes. If monitor is TRUE, then trace plots will also be updated at this interval.
#' @param modelComparison Logical. If \code{TRUE}, then model information criteria are returned, including the BIC and DIC3. In addition the log pseudo-marginal likelihood (LPML) statistic is computed. Model comparison is allowed only under \code{imputeResponse} is \code{TRUE}. Note that under an MNAR visit process or response process process given a clinic visit, these model comparison statistics may not have the desired interpretation. See Daniels and Hogan's Missing Data in Longitudinal Studies (2008).
#' @param lastIterations A scalar for the number of last iterations for which to carry out model comparison (if modelComparison equals \code{TRUE}), model checking, and write out draws of imputed longitudinal outcomes \code{Y1,\dots,YJ}, latent class membership \code{C}, and random effects \code{b} from the longitudinal outcomes submodel. If \code{NULL}, then this will be equal to \code{n.samples - burn}. This option was included because these additional computations can markedly increase model fitting time.
#' @param sims Logical. If \code{TRUE}, then the average class membership probabilities are returned with MCMC samples, and components of an applied analysis are not executed including model comparison, posterior predictive draws of the longitudinal health outcomes, posterior predictive checking, and written text files of imputations of the longitudinal health outcomes.
#'
#' @details Details are given on the specification of \code{priors} and \code{inits} and on the model comparison analysis.
#' \describe{
#' \item{A.	The nested list \code{priors} must contain the information on the prior distributions in the order as given below.}{
#' \enumerate{
#' \item Multivariate normal prior distribution on the regression coefficients in the latent class membership model. The list entry is of the form, for example, \code{list(rep(0, m), diag(10, m))}, where \code{m} is the number of regression coefficients. \code{rep(0, m)} assigns a prior mean of 0 to all regression coefficients. A diagonal variance-covariance matrix with variance 10 is assigned.
#'
#' \item Multivariate normal prior distribution on the regression coefficients in the design matrix specified by \code{regf} named formula "YObs". The list entry is specified in the same way as in the latent class membership model.
#'
#' \item Multivariate normal prior distribution on the regression coefficients in the design matrix specified by \code{regf} named formula "YSub". The list entry is specified in the same way as in the latent class membership model.
#'
#' \item Inverse-Wishart prior distribution on the variance-covariance matrix of the subject-level random effects in the longitudinal health outcomes model. The form is, for example, \code{list(diag(0.25, q), q + 2)}, where \code{q} is the number of random effects. The first element is the prior scale matrix of dimension \code{q} by \code{q}, and the second element is the prior degrees of freedom.
#'
#' \item Inverse-Wishart prior distribution on the variance-covariance matrix of the longitudinal health outcomes \code{Y1,\dots,YJ}. The form is, for example, \code{list(diag(0.25, R), R + 2)}, where \code{R} is the number of longitudinal health outcomes. The first element is the prior scale matrix of dimension \code{R} by \code{R}, and the second element is the prior degrees of freedom.
#'
#' \item If \code{modelVisit = TRUE}, a multivariate normal prior distribution on the regression coefficients in the design matrix specified by \code{regf} named formula "DObs". The list entry follows the same format as in the latent class membership model. If \code{modelVisit = FALSE}, then use \code{list(NULL)}.
#'
#' \item If \code{modelVisit = TRUE}, an inverse-Wishart prior distribution on the variance-covariance matrix of the subject-level random effects in the visit process model. The form follows the same format as in the longitudinal health outcomes model. The number of random effects must be the same as in the longitudinal health outcomes model. If \code{modelVisit = FALSE}, then use \code{list(NULL)}.
#'
#' \item If \code{modelResponse = TRUE}, a multivariate normal prior distribution on the regression coefficients in the design matrix specified by \code{regf} named formula "MObs". The list entry follows the same format as in the latent class membership model. If \code{modelResponse = FALSE}, then use \code{list(NULL)}.
#'
#' \item If \code{modelResponse = TRUE}, an inverse-Wishart prior distribution on the variance-covariance matrix of the subject-level random effects in the response process model. The form follows the same format as in the longitudinal health outcomes model. The number of random effects must be the same as in the longitudinal health outcomes model. If \code{modelResponse = FALSE}, then use \code{list(NULL)}.
#' }
#'
#' An example of the completed list object for specification of the prior distributions is: \code{priors <- list(list(rep(0, m), diag(1, m)), list(rep(0, s), diag(100, s)), list(rep(0, p), Sigma = diag(10000, p)), list(1, 1), list(diag(c(0.5, 0.6), J), (J + 2)), list(rep(0, f), diag(10, f)), list(1, 0.1), list(rep(0, e), diag(10, e)), list(1, 1))}. \code{m} is the number of covariates in the latent class membbership model; \code{s} is the number of covariates in the "YObs" design matrix; \code{p} is the number of covariates in the "YSub" design matrix; \code{f} is the number of covariates in the "DObs" design matrix; and, \code{e} is the number of covariates in the "MObs" design matrix.}
#'
#' \item{B.	The list \code{inits} must contain the information on initial values for each of the parameters, in the same order as in \code{priors}.}{
#' \enumerate{
#' \item Initial values for the regression coefficients in the latent class membership model. The list entry is of the form, for example, \code{matrix(rep(0, m * (K - 1)), nrow = m, ncol = (K - 1))}, where \code{m} is the number of regression coefficients, \code{K} is the assumed number of latent classes.
#'
#' \item Initial values for the regression coefficients in the design matrix specified by \code{regf} named formula "YObs". The object is a \code{J} element list in which each element is an \code{s} by \code{K} matrix, where \code{s} is the number of regression coefficients and \code{K} is the number of latent classes. An example for \code{J=2} is \code{list(matrix(rnorm(s * K), ncol = K, nrow = s), matrix(rnorm(s * K), ncol = K, nrow = s))}.
#'
#' \item Initial values for the regression coefficients in the design matrix specified by \code{regf} named formula "YSub". The object is a \code{J} element list in which each element is an \code{p} by \code{q} by \code{K} array, where \code{p} is the number of regression coefficients, \code{q} is the number of random effects, and \code{K} is the number of latent classes. An example for \code{J=2} is \code{list(array(rnorm(p*q*K), dim = c(p, q, K)), array(rnorm(p*q*K), dim = c(p, q, K)))}.
#'
#' \item Initial values for the variance-covariance matrix of the subject-level random effects in the longitudinal health outcomes model. The object is a \code{J} element list in which each element is an \code{q} by \code{q} by \code{K} array, where \code{q} is the number of random effects, and \code{K} is the number of latent classes. For an example, see the initial values for the regression coefficients in the design matrix specified by \code{regf} named formula "YSub".
#'
#' \item Initial values on the variance-covariance matrix of the longitudinal health outcomes \code{Y1,\dots,YJ}. The object is a \code{J} by \code{J} by \code{K} array.
#'
#' \item If \code{modelVisit = TRUE}, initial values on the regression coefficients in the design matrix specified by \code{regf} named formula "DObs". An example is \code{list(matrix(rnorm(f*K), ncol=K, nrow=f))}. If \code{modelVisit = FALSE}, then use \code{list(NULL)}.
#'
#' \item If \code{modelVisit = TRUE}, initial values on the variance-covariance matrix of the subject-level random effects in the visit process model. The object is an \code{q} by \code{q} by \code{K} array, where \code{q} is the number of random effects, and \code{K} is the number of latent classes. If \code{modelVisit = FALSE}, then use \code{list(NULL)}.
#'
#' \item If \code{modelResponse = TRUE}, initial values on the regression coefficients in the design matrix specified by \code{regf} named formula "MObs". The object is a list of the same number of elements as \code{Mvec}. The order of the initial values must correspond to \code{Mvec}. For example, if \code{Mvec = c(2)}, then \code{list(matrix(rnorm(e*K), ncol=K, nrow=e))}. If \code{modelResponse = FALSE}, then use \code{list(NULL)}.
#'
#' \item If \code{modelResponse = TRUE}, initial values on the variance-covariance matrix of the subject-level random effects in the response process model. The object is a list of the same number of elements as \code{Mvec}. The order of the initial values must correspond to \code{Mvec}. An example for \code{Mvec = c(2)} is \code{list(array(diag(0.5, q), dim = c(q, q, K))}. If \code{modelResponse = FALSE}, then use \code{list(NULL)}.}}}
#'
#' @export
#' @return Model summaries are printed to the console, including posterior means and 95\% credible intervals, posterior latent class assignment, and model comparison statistics. A label switching diagnostic using Stephen's method from the \code{label.switching} package in \code{R} is printed.
#' A list of matrices of stored MCMC samples (\code{n.samples-burn}) for unknown parameters, including:
#' \enumerate{
#' \item \code{store_delta}: Regression coefficients from the latent class membership model.
#' \item \code{store_pi}: Posterior probabilities of latent class membership for each subject and latent class.
#' \item \code{store_priorPi}: Probabilities of latent class membership for each subject and latent class from the multinomial probit model of latent class membership.
#' \item \code{store_betaObs}: Regression coefficients from the design matrix from formula "YObs" in the longitudinal health outcomes model.
#' \item \code{store_betaSub}: Regression coefficients from the design matrix from formula "YSub" in the longitudinal health outcomes model.
#' \item \code{store_Psi}: Variance-covariance of the subject-level random effects in the longitudinal health outcomes model.
#' \item \code{store_Sigma}: Variance-covariance of the  longitudinal health outcomes model.
#' \item \code{store_phi}: Regression coefficients from the design matrix from formula "DObs" in the visit process model. If \code{modelVisit = FALSE}, this will be \code{NULL}.
#' \item \code{store_Omega}: Variance-covariance of the subject-level random effects in the visit process model. If \code{modelVisit = FALSE}, this will be \code{NULL}.
#' \item \code{store_lambda}: Regression coefficients from the design matrix from formula "MObs" in the response process model. If \code{modelResponse = FALSE}, this will be \code{NULL}.
#' \item \code{store_Theta}: Variance-covariance of the subject-level random effects in the response process model. If \code{modelResponse = FALSE}, this will be \code{NULL}.
#' \item \code{store_ame}: Population-averaged regression coefficients in the longitudinal health outcomes model obtained by averaging over the latent class membership probabilities.
#' \item \code{store_class_weight}: Latent class membership probabilities averaged from the subject-level to the latent class-level. If \code{sims = FALSE}, these will not be returned.
#' }
#'
#' If \code{sims = FALSE}, then after burn-in, samples for draws from the posterior predictive distribution, the discrepancy measure for posterior predictive checking, and imputations for the missing longitudinal health outcomes are written to separate comma-separated text files in the working directory. In addition, text files are printed for saved draws of the random effects from the longitudinal outcomes submodel \code{store_bSub.txt} and the discrete latent class membership variable \code{store_C.txt}. Saved random effects can be used for re-constructing patient-specific curves. Saved latent class membership draws may be useful for assessing the key conditional independence assumption.
#'
#'Draws from the posterior predictive distribution are stored in \code{store_Ydraw.txt}. These are samples of the replicated completed longitudinal health outcomes. See Gelman, A., Mechelen, I. V., Verbeke, G., Heitjan, D.F. and Meulders, M. (2005) Multiple Imputation for Model Checking: Completed-Data Plots with Missing and Latent Data. Biometrics, 61, 74-85. This file is of dimension \code{n.samples - burn} by \code{N}, where \code{N} is the number of observations given a clinic visit (i.e., \code{length(Y)}). In addition posterior predictive draws of the visit process \code{D}, and the response process given a clinic visit \code{M1,\dots,MJ} are stored in \code{store_Ddraw.txt} and \code{store_Mdraw.txt}. Under MNAR, these posterior predictive draws can be used for model checking based on the observed and replicated observed data. See the following references: 1. Xu, Dandan, Arkendu Chatterjee, and Michael Daniels. "A note on posterior predictive checks to assess model fit for incomplete data." Statistics in medicine 35.27 (2016): 5029-5039.; and 2. Daniels, Michael J., Arkendu S. Chatterjee, and Chenguang Wang. "Bayesian model selection for incomplete data using the posterior predictive distribution." Biometrics 68.4 (2012): 1055-1063.
#'
#'Samples of a measure of discrepancy are stored in \code{store_T_completed.txt}. This file is of dimension \code{n.samples - burn} by \code{2}, where the first column is the discrepancy measure computed using the completed data, and the second column is the discrepancy measure using the replicated completed data. The discrepancy measure is the multivariate mean square error. See \code{?get_discrepancy_plot}.
#'
#' Samples of the imputations for \code{Y1,\dots,YJ} are written to \code{store_miss_Y1.txt,\dots,\code{store_miss_YJ.txt}}. If a longitudinal health outcome has no missing values, no file is generated. Note that imputations will not be generated for patient-time windows during which none of the longitudinal outcomes \code{Y1,\dots,YJ} are observed.
#'
#' @examples
#'
#'
#' #------------------------------- Load data from the EHRMiss package
#' data(growth)
#' names(growth)
#' dim(growth)
#'
#' ### Construct consecutive integer valued subjectID
#' growth$subjectID <- seq(1, length(unique(growth$subjectID)), by = 1)[factor(growth$subjectID)]
#'
#' #------------------------------ Model details
#' ### Number of outcomes
#' J <- 2
#' ### Number of latent classes
#' K <- 2

#' #------------------------------- Specify the formulas for the design matrices and put the
#' #formulas in regf
#' regf <- list(LatentClass = ~ 1 + birthweight,
#'             YRe = ~ 1,
#'             YObs = ~ -1 + time,
#'             YSub = ~ 1,
#'             DObs = ~ 1 + time,
#'             DRe = ~ 1,
#'             MObs = ~ 1 + time,
#'             MRe = ~ 1)
#'
#' #----------------------- MCMC preparation
#' m <- length(all.vars(regf[["LatentClass"]])) + 1
#' p <- length(all.vars(regf[["YSub"]])) + 1
#' s <- length(all.vars(regf[["YObs"]]))
#' f <- length(all.vars(regf[["DObs"]])) + 1
#' e <- length(all.vars(regf[["MObs"]])) + 1
#'
#' ### Number of random effects, assumed the same for #'all models
#' q <- length(all.vars(regf[["YRe"]])) + 1
#'
#' ### Prior distributions
#' priors <- list(list(rep(0, m), diag(1, m)), list(rep(0, s), diag(100, s)),
#'               list(rep(0, p), diag(10000, p)),
#'               list(1, 1),
#'               list(diag(c(0.5, 0.6), J), (J + 2)),
#'               list(rep(0, f), diag(100, f)),
#'               list(1, 1),
#'               list(rep(0, e), diag(100, e)),
#'               list(scale = 1, df = 1))
#'
#' ### Initial values
#' inits <- list(matrix(rep(0, m * (K - 1)), nrow = m, ncol = (K - 1)),
#'              list(matrix(rnorm(s * K), nrow = s, ncol = K),
#'              matrix(rnorm(s * K), nrow = s, ncol = K)),
#'              list(array(rnorm(p * q * K), dim = c(p, q, K)),
#'              array(rnorm(p * q * K), dim = c(p, q, K))),
#'              list(array(rep(0.4, K), dim = c(q, q, K)),
#'              array(rep(0.4, K), dim = c(q, q, K))),
#'              array(c(1, 0, 0, 1, 0.5, 0, 0, 0.5), dim = c(J, J, K)),
#'              matrix(rnorm(f * K), nrow = f, ncol = K),
#'              array(rep(0.5, K), dim = c(q, q, K)),
#'              list(matrix(rnorm(e * K), ncol = K)),
#'              list(array(rep(0.5, K), dim = c(q, q, K))))

#' #-------------------------- Fit model with MNAR visit process, MNAR response process for Y2
#'
#' n.samples <- 100
#' burn <- 50
#' update <- 10
#' monitor <- TRUE
#'
#'
#' ### MNAR Visit process
#' ### MNAR response process for Y2
#' res <- MVNYBinaryMiss(K = K, J = J, data = growth,
#'  regf = regf, imputeResponse = TRUE,
#' Mvec = 2, modelVisit = TRUE, modelResponse = TRUE, priors =  priors, inits = inits,
#' n.samples = n.samples, burn = burn, monitor = monitor, update = update,
#' modelComparison = TRUE, lastIterations = 10, sims = FALSE)
#' ### Get Bayesian posterior predictive p-value, see ?get_discrepancy_plot
#' # Set working directory to where discrepancy samples are written
#' store_T_completed <- read.table("store_T_completed.txt", header = FALSE, sep = ",")
#' get_discrepancy_plot(store_T_completed)
MVNYBinaryMiss <- function(K, J, data, regf, imputeResponse, Mvec, modelVisit, modelResponse, priors, inits, n.samples, burn, monitor, update, modelComparison, lastIterations, sims = FALSE) {


  ### Imputations are conducted only for Y | D = 1 & M = 0
  if (imputeResponse == TRUE) {

    #------------------------ For modeling choices based on imputing Y for missing response
    # Model visit process or not with total observations N
    # Restrict to observed visits for Y and M so modeling Y | D = 1 and M | D = 1, total observations NVisit
    # Model response or not, total observations NVisit

    ### Modeling components for model with visit process model
    if (modelVisit == TRUE) {

      # Subject ID for the visit process model
      subjectID <- data$subjectID
      # D for the visit process model
      D <- data$D
      # Observation level design matrix for D
      UObs <- model.matrix(regf[["DObs"]], data = data)
      # Random effects design matrix design matrix for D
      URe <- model.matrix(regf[["DRe"]], data = data)

    } else {

      subjectID <- NULL
      D <- NULL
      UObs <- NULL
      URe <- NULL

    }

    ### Restrict to restrict to D = 1 for Y
    data <- subset(data, subset = D == 1)
    cat("Number of obs. after restricting to observed visits:", dim(data)[1], "\n")

    # Identify columns with outcome Y
    y_col_index <- grepl("^Y[1-9]$", colnames(data), perl = TRUE)

    #! Change on 1/22
    ### Restrict to observed visits with at least 1 Y_{jit} observed
    #data$exc <- sapply(1:nrow(data), function(x) {
    #  exc <- 1 * (sum(is.na(data[x , paste("Y", 1:J, sep = "")])) == J)
    #  return(exc)
    #}
    #)

    data$exc <- sapply(1:nrow(data), function(x) {
      exc <- 1 * (sum(is.na(data[x , y_col_index])) == J)
      return(exc)
    }
    )
    data <- subset(data, exc == 0)

    # Subject ID for Y is NVisit by 1
    subjectIDY <- data$subjectID

    # Multivariate Y, NVisit by J
    #! Change on Jan 24, 2023
    #Y <- data.matrix(data[ , paste("Y", 1:J, sep = "")])
    Y <- data.matrix(data[ , y_col_index])

    # Random effects design matrix for Y, includes at least intercept
    XRe <- model.matrix(regf[["YRe"]], data = data)

    # Observation level design matrix for Y, no intercept. Includes time-invariant and fully observed subject level covariates (e.g., race) that will not enter the random effects equations
    XObs <- model.matrix(regf[["YObs"]], data = data)

    # Subject level design matrix for Y, including intercept and subject level covariates that enter into random effects equations
    XSub <- model.matrix(regf[["YSub"]], data = aggregate(data, by = list(subjectIDY), FUN = tail, n = 1))


    ### Components for imputation of response process process missingness given D = 1
    # Response indicator generated for each Y_j
    # M is NVisit x J with 0 as missing and 1 as observed
    #Response variables in the data must have names indexing Y (e.g., if Y2 is missing in response, then name of M is M2)
    #! Change on Jan 24, 2023
    #M <- data.matrix(data[ , paste("M", 1:J, sep = "")])

    # Identify columns with response measurement M
    m_col_index <- grepl("^M[1-9]$", colnames(data), perl = TRUE)
    M <- data.matrix(data[ , m_col_index])

    ### Modeling components if we are modeling the response processes given D = 1
    if (modelResponse == TRUE) {

      # Subject ID for M | D
      Mvec <- Mvec
      subjectIDM <- data$subjectID
      # Observation level design matrix for M | D = 1
      VObs <- model.matrix(regf[["MObs"]], data = data)
      # Random effects design matrix design matrix for M | D = 1
      VRe <- model.matrix(regf[["MRe"]], data = data)

    } else {
      # For update C defaults
      Mvec <- NULL
      subjectIDM <- NULL
      VObs <- NULL
      VRe <- NULL
    }

    ### Modeling components for latent class membership model
    W <- model.matrix(regf[["LatentClass"]], data = aggregate(data, by = list(subjectIDY), FUN = tail, n = 1))

  }


  #------------------------ For modeling choices based on no imputations
  # For Y, restrict to observed records for Y so modeling Y | D = 1, M = 1, analysis  based on NComplete (complete cases based on observed Y1 *and* Y2)
  if (imputeResponse == FALSE) {

    # M is N x J indicator matrix
    # M[ , j] = 1 if observed; 0 if missing given D = 1; and NA if D = 0
    #! Change on Jan 24, 2023
    #M <- data.matrix(data[ , paste("M", 1:J, sep = "")])
    # Identify columns with response measurement M
    m_col_index <- grepl("^M[1-9]$", colnames(data), perl = TRUE)
    M <- data.matrix(data[ , m_col_index])

    # Logical: TRUE if observed responses on all Y_j for subject i
    allObs <- apply(M, 1, function(x) sum(x, na.rm = TRUE) == J)

    ### For Y, restrict dataset to D = 1 & M_j = 1 all j
    # Restrict dataset to only observed response for all Y_j per subject i and observation window t
    data <- subset(data, allObs == TRUE)
    cat("Number of obs. after restricting to complete pairs:", dim(data)[1], "\n")


    ### Modeling components for model Y | M = 1 (implies D = 1)
    # Subject ID for Y is NCompletePairs by 1
    # Re-number so subscript is not out of bounds
    subjectIDY <- seq(1, length(unique(data$subjectID)), by = 1)[factor(data$subjectID)]

    # Multivariate Y, NCompletePairs by J
    #! Change on Jan 24, 2023
    #Y <- data.matrix(data[ , paste("Y", 1:J, sep = "")])
    # Identify columns with outcome Y
    y_col_index <- grepl("^Y[1-9]$", colnames(data), perl = TRUE)
    Y <- data.matrix(data[ , y_col_index])

    # Random effects design matrix for Y, includes at least intercept
    XRe <- model.matrix(regf[["YRe"]], data = data)

    # Observation level design matrix for Y, no intercept. Includes time-invariant and fully observed subject level covariates (e.g., race) that will not enter the random effects equations
    XObs <- model.matrix(regf[["YObs"]], data = data)

    # Subject level design matrix for Y, including intercept and subject level covariates that enter into random effects equations
    XSub <- model.matrix(regf[["YSub"]], data = aggregate(data, by = list(subjectIDY), FUN = tail, n = 1))

    ### Modeling components for latent class membership model
    W <- model.matrix(regf[["LatentClass"]], data = aggregate(data, by = list(subjectIDY), FUN = tail, n = 1))

    # Update C function defaults
    subjectID <- NULL
    D <- NULL
    UObs <- NULL
    URe <- NULL

    subjectIDM <- NULL
    VObs <- NULL
    VRe <- NULL
    Mvec <- NULL

  }


  #---------------------------- Data checks
  if (length(unique(subjectIDY)) != dim(W)[1]) {
    stop(paste("No. unique subjects in longitudinal outcomes model does not equal no. subjects in latent class membership model"))
  }

  if (modelVisit & length(unique(subjectID)) != dim(W)[1]) {
    stop(paste("No. unique subjects in visit process model does not equal no. subjects in latent class membership model"))
  }

  if (modelResponse & length(unique(subjectIDM)) != dim(W)[1]) {
    stop(paste("No. unique subjects in response process model does not equal no. subjects in latent class membership model"))
  }

  #-----------------------------------  Background information
  ### Number of unique subjects
  n <- nrow(W)
  ### Number of latent classes
  K <- K
  ### Number of outcomes
  J <- ncol(Y)

  ### Number of covariates in each design matrix
  m <- ncol(W)
  s <- ncol(XObs)
  p <- ncol(XSub)
  ### q will be the number of random effects in the longitudinal outcomes model
  q <- ncol(XRe)
  if (modelVisit == TRUE) f <- ncol(UObs)
  if (modelResponse == TRUE) e <- ncol(VObs)


  #-------------------------------------------- Prior distributions
  ### Multinomial probit model for latent class assignment
  prior_mu_delta <- priors[[1]][[1]]
  prior_Sigma_delta <- priors[[1]][[2]]

  ### Multivariate outcomes model
  prior_mu_betaObs <- priors[[2]][[1]]
  prior_Sigma_betaObs <- priors[[2]][[2]]

  prior_mu_betaSub <- priors[[3]][[1]]
  prior_Sigma_betaSub <- priors[[3]][[2]]

  prior_scale_Psi <- priors[[4]][[1]]
  prior_df_Psi <- priors[[4]][[2]]

  prior_scale_Sigma <- priors[[5]][[1]]
  prior_df_Sigma <- priors[[5]][[2]]

  ### D visit process model
  if (modelVisit == TRUE) {
    prior_mu_phi <- priors[[6]][[1]]
    prior_Sigma_phi <- priors[[6]][[2]]
    prior_scale_Omega <- priors[[7]][[1]]
    prior_df_Omega <- priors[[7]][[2]]
  }

  ### M response process given a visit model
  if (modelResponse == TRUE) {
    prior_mu_lambda <- priors[[8]][[1]]
    prior_Sigma_lambda <- priors[[8]][[2]]
    prior_scale_Iota <- priors[[9]][[1]]
    prior_df_Iota <- priors[[9]][[2]]
  }


  #---------------------------------------------- Initialization
  ### Initialize C based on number of subjects in X above
  delta <- inits[[1]]
  C <- sample(1:K, size = n, prob = rep(1/K, K), replace = TRUE)
  #! Change on 4/19/2022
  if (K > 2) {
    #Z <- mnormt::rmnorm(n, varcov = diag(K - 1)) + W %*% delta
    # K-1 by K-1 covariance matrix
    VarCovP <- matrix(1, nrow = (K - 1), ncol = (K - 1))
    diag(VarCovP) <- 2
    err <- do.call("rbind", sapply(1:n, function(x) {res <- mnormt::rmnorm(1, mean = rep(0, (K-1)), varcov = VarCovP)
    return(res) }, simplify = FALSE))

    Z <- err + W %*% delta

  } else if (K == 2) {
    #Z <- mnormt::rmnorm(n, varcov = 1) + W %*% delta
    err <- mnormt::rmnorm(n, varcov = diag(K - 1))
    Z <- err + W %*% delta

  }

  ### Longitudinal multivariate Y
  betaObs <- inits[[2]]
  betaSub <- inits[[3]]
  Psi <- inits[[4]]
  Sigma <- inits[[5]]

  ### Inits for bSub for each j
  bSub <- as.list(1:J)
  for (j in 1:J) {
    bSubm <- matrix(NA, nrow = n, ncol = q)
    for (k in 1:K) {
      ind_sub <- which(C == k)
      mubSub <- as.matrix(XSub[ind_sub, ]) %*% betaSub[[j]][ , , k]
      bSubm[ind_sub, ] <- t(apply(mubSub, 1, function(x) {
        mnormt::rmnorm(1, mean = x, varcov = as.matrix(Psi[[j]][ , , k]))
      }))
    }

    bSub[[j]] <- c(t(bSubm))
  }

  if (imputeResponse == TRUE) {

    C_expand <- C[factor(subjectIDY)]

    for (j in 1:J) {
      for (k in 1:K) {
        obs <- M[ , j] == 1
        miss <- M[ , j] == 0
        Y[C_expand == k & miss , j] <- mean(Y[C_expand == k & obs , j])
      }
    }

  }


  ### Visit process D
  if (modelVisit == TRUE) {

    phi <- inits[[6]]
    Q <- rnorm(dim(UObs)[1])
    Omega <- inits[[7]]
    tau <- matrix(NA, nrow = n, ncol = ncol(URe))

    for (k in 1:K) {
      ind_sub <- which(C == k)
      tau[ind_sub, ] <- mnormt::rmnorm(length(ind_sub), mean = rep(0, ncol(URe)), varcov = Omega[ , , k])
    }

    tau <- c(t(tau))

  } else {

    phi <-  NULL
    tau <- NULL

  }


  ### Response process M | D = 1
  if (modelResponse == TRUE) {

    # Mvec indicates which columns of M will be modeled
    Jstar <- length(Mvec)
    lambda <- inits[[8]]
    L <- matrix(rnorm(length(subjectIDM) * Jstar), nrow = dim(VObs)[1], ncol = Jstar)
    Theta <- inits[[9]]
    kappa <- as.list(1:Jstar)

    # theta is a J element list with nq times 1 stacked random effects for each Y_j with response process missingness
    for (j in 1:Jstar) {
      kappaTemp <- matrix(NA, nrow = n, ncol = ncol(VRe))
      for (k in 1:K) {
        ind_sub <- which(C == k)
        kappaTemp[ind_sub, ] <- mnormt::rmnorm(length(ind_sub), mean = rep(0, ncol(VRe)), varcov = Theta[[j]][ , , k])
      }
      kappa[[j]] <- c(t(kappaTemp))
    }
  } else {
    lambda <-  NULL
    kappa <- NULL
  }

  #---------------------------------------------------------- Storage of MCMC samples

  ### Relevant to model comparison
  #! Changed 3/17/2023 to new.burn to accomodate saving only a lastIterations of samples
  if (!is.null(lastIterations)) {
    new.burn <- n.samples - lastIterations
  } else if (is.null(lastIterations)) {
    new.burn <- burn
  }



  ### Multinomial class membership model
  store_delta <- matrix(NA, nrow = (n.samples - burn), ncol = ((K - 1)*m), dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(rep(paste("Class", rep(2:K), sep = ""), each = length(colnames(W))), colnames(W), sep = "_")))

  store_pi <- matrix(NA, nrow = (n.samples - burn), ncol = (K*n), dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(paste("Class", rep(1:K, each = n), sep = ""), "Subject", 1:n, sep = "_")))

  store_priorPi <- matrix(NA, nrow = (n.samples - burn), ncol = (K*n), dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(paste("Class", rep(1:K, each = n), sep = ""), "Subject", 1:n, sep = "_")))

  ### Longitudinal outcomes
  #Stored by first outcome, first latent class, each variable, first outcome, second latent class, each variable, and so forth
  store_betaObs <- matrix(NA, nrow = n.samples - burn, ncol = length(c(unlist(betaObs))), dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(rep(paste("Y", rep(1:J), sep = ""), each = s*K), rep(rep(paste("Class", rep(1:K), sep = ""), each = s), times = J), rep(colnames(XObs), times = K*J), sep = "_")))

  store_betaSub <- matrix(NA, nrow = n.samples - burn, ncol = p*q*K*J, dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(rep(paste("Y", rep(1:J), sep = ""), each = p*q*K), rep(rep(paste("Class", rep(1:K), sep = ""), each = p*q), times = J), rep(rep(paste("RE", rep(1:q), sep = ""), each = p), times = J*K), rep(colnames(XSub), times = q*K*J), sep = "_")))

  store_Sigma <- matrix(NA, nrow = (n.samples - burn), ncol = J*J*K, dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(rep(paste("Class", rep(1:K), sep = ""), each = J*J), rep(paste("Sigma", rep(1:J, times = J), rep(1:J, each = J), sep = ""), times = K), sep = "_")))

  store_Psi <- matrix(NA, nrow = (n.samples - burn), ncol = q*q*K*J, dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(rep(paste("Y", rep(1:J), sep = ""), each = q*q*K), rep(rep(paste("Class", rep(1:K), sep = ""), each = q*q), times = J), rep(rep(paste("Psi", rep(1:q, times = q), rep(1:q, each = q), sep = ""), times = K), times = J), sep = "_")))

  store_ame <- matrix(NA, nrow = n.samples - burn, ncol = (p*q + s)*J, dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste("AME", rep(paste("Y", rep(1:J), sep = ""), each = (p*q + s)), rep(c(paste(rep(paste("RE", rep(1:q), sep = ""), each = p), rep(colnames(XSub), times = q), sep = "_"), colnames(XObs)), times = J), sep = "_")))

  ### Visit process D model
  if (modelVisit == TRUE) {
    store_phi <- matrix(NA, nrow = (n.samples - burn), ncol = (K*f), dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"),            paste("D", paste("Class", rep(1:K, each = f), sep = ""), colnames(UObs), sep = "_")))
    store_Omega <- matrix(NA, nrow = (n.samples - burn), ncol = q*q*K, dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste("D", rep(paste("Class", rep(1:K), sep = ""), each = q*q), rep(paste("Omega", rep(1:q, times = q), rep(1:q, each = q), sep = ""), times = K), sep = "_")))
  }

  ### Response process M model
  if (modelResponse == TRUE) {
    store_lambda <- matrix(NA, nrow = n.samples - burn, ncol = length(c(unlist(lambda))), dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(rep(paste("M", rep(Mvec), sep = ""), each = e*K), rep(rep(paste("Class", rep(1:K), sep = ""), each = e), times = Jstar), rep(colnames(VObs), times = K*Jstar), sep = "_")))
    store_Theta <- matrix(NA, nrow = (n.samples - burn), ncol = q*q*K*Jstar, dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(rep(paste("M", Mvec, sep = ""), each = q*q*K), rep(rep(paste("Class", rep(1:K), sep = ""), each = q*q), times = Jstar), rep(rep(paste("Theta", rep(1:q, times = q), rep(1:q, each = q), sep = ""), times = K), times = Jstar), sep = "_")))

  }

  ### Model comparison based on information criteria storage
  #! Changed 3/17/2023 to new.burn to accomodate saving only a lastIterations of samples
  if (sims == FALSE & imputeResponse == TRUE & modelComparison == TRUE) {

    #Effective sample size
    store_ESS <- rep(NA, n.samples - new.burn - 1)
    store_observed_lik_sub <- matrix(NA, nrow = (n.samples - new.burn - 1), ncol = n)
    store_observed_llik <- rep(NA, n.samples - new.burn - 1)
  }

  ### For data simulation statistics
  if (sims == TRUE) {
    store_class_weight <- matrix(NA, nrow = n.samples - burn, ncol = K, dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste("Weight", paste("Class", rep(1:K), sep = ""), sep = "_")))
  }

  #----------------------------------- Start time for loop
  ptm <- proc.time()

  for (i in 1:n.samples) {

    if (K > 2) {
      ### Update parameters in multinomial probit model
      Z <- update_Z_MNP(Z = Z, C = C, delta = delta, W = W)
      delta <- update_delta_MNP(Z = Z, delta = delta, W = W, prior.mu = prior_mu_delta, prior.Sigma = prior_Sigma_delta)

    } else {
      ### Update parameters in binary probit model
      Z <- update_Z_BinaryP(Z = Z, C = C, delta = delta, W = W)
      delta <- update_delta_BinaryP(Z = Z, delta = delta, W = W, prior.mu = prior_mu_delta, prior.Sigma = prior_Sigma_delta)
    }

    # Update prior probability of class membership based on Z
    priorPik <- update_priorPik(Z = Z)

    ### Update parameters in Y model
    # Update ... fast for Sigma and betaObs is slower than the functions below
    Sigma <- update_SigmaLC(C = C, betaObs = betaObs, bSub = bSub, Y = Y, XRe = XRe,  XObs = XObs, prior.scale = prior_scale_Sigma, prior.df = prior_df_Sigma, subjectID = subjectIDY)

    betaObs <- update_betaObsLC(C = C, betaObs = betaObs, bSub = bSub, Sigma = Sigma, Y = Y, XRe = XRe, XObs = XObs, prior.mu = prior_mu_betaObs, prior.Sigma = prior_Sigma_betaObs, subjectID = subjectIDY)

    bSub <- update_bSubLC_fast(C = C, bSub = bSub, betaObs = betaObs, betaSub = betaSub, Psi = Psi, Sigma = Sigma, Y = Y, XRe = XRe, XObs = XObs, XSub = XSub, subjectID = subjectIDY)

    Psi <- update_PsiLC(C = C, betaSub = betaSub, bSub = bSub, Y = Y, XSub = XSub, prior.scale = prior_scale_Psi, prior.df = prior_df_Psi)

    betaSub <- update_betaSubLC(C = C, bSub = bSub, Psi = Psi, Y = Y, XSub = XSub, prior.mu = prior_mu_betaSub, prior.Sigma = prior_Sigma_betaSub)

    ### Imputation of Y for Y_j with missing response | D = 1
    if (imputeResponse == TRUE) {

      C_expand <- C[factor(subjectIDY)]
      muY <- matrix(NA, nrow = nrow(Y), ncol = J)

      # Calculate muY conditional on class membership for each Y_j
      for (k in 1:K) {
        ind_sub <- which(C == k)
        ind_obs <- which(C_expand == k)

        XRek <- as.matrix(XRe[ind_obs, ])
        subjectIDYk <- subjectIDY[ind_obs]

        for (j in 1:J) {

          # Calculation for random effects of each subject
          sp_XRe_sub <- lapply(split(as.data.frame(XRek), subjectIDYk, drop = TRUE), as.matrix)
          convert_XRe_sub <- as.matrix(Matrix::bdiag(sp_XRe_sub))
          tempm <- matrix(bSub[[j]], nrow = n, ncol = q, byrow = TRUE)
          tempmk <- as.matrix(tempm[ind_sub, ])
          bSubj <- as.vector(convert_XRe_sub %*% c(t(tempmk)))
          # For each subject and observation time, the mu_{it} will use the latent class specific regression coefficients and the random effects according to posterior latent class assignment
          muY[C_expand == k, j] <- as.matrix(XObs[C_expand == k, ]) %*% betaObs[[j]][ , k] + bSubj
        }
      }

      for (t in 1:nrow(Y)) {

        obs <- M[t, ] == 1
        miss <- M[t, ] == 0
        k <- C_expand[t]

        SigmaObs <- solve(Sigma[obs, obs, k])

        muMiss <- muY[t, miss] + Sigma[miss, obs, k] %*% SigmaObs %*% (Y[t, obs] - muY[t, obs])
        SigmaMiss <- Sigma[miss, miss, k] - Sigma[miss, obs, k] %*% SigmaObs %*% Sigma[obs, miss, k]

        if (any(miss)) {
          Y[t, miss] <- mnormt::rmnorm(1, mean = muMiss, varcov = SigmaMiss)
        } else {
          Y[t, ] <- Y[t, ]
        }

      }

    }


    ### Update parameters in D model
    if (modelVisit == TRUE) {

      Q <- update_ZLC_ProbitRE(C = C, phi = phi, tau = tau, D = D, UObs = UObs, URe = URe, subjectID = subjectID)
      phi <- update_betaLC_ProbitRE(C = C, Z = Q, tau = tau, UObs = UObs, URe = URe, prior.mu = prior_mu_phi, prior.Sigma = prior_Sigma_phi, subjectID = subjectID)
      tau <- update_bLC_ProbitRE_fast(C = C, Omega = Omega, Z = Q, phi = phi, UObs = UObs, URe = URe, subjectID = subjectID)
      Omega <- update_OmegaLC_ProbitRE(C = C, tau = tau, prior.scale = prior_scale_Omega, prior.df = prior_df_Omega)

    }

    ### Update parameters in M | D == 1 model
    #M[ , Mvec[j]] will give the response in M that needs to be modeled 1,\dots,Jstar
    if (modelResponse == TRUE) {
      # Assume that f(M_1, M_2,...M_Jstar | D = 1, C; ...) = f(M_1 | D = 1, C; ...)f(M_2 | D = 1, C; ...)...f(M_Jstar | D = 1, C; ...)
      for (j in 1:Jstar) {
        L[ , j] <- update_ZLC_ProbitRE(C = C, phi = lambda[[j]], tau = kappa[[j]], D = M[ , Mvec[j]], UObs = VObs, URe = VRe, subjectID = subjectIDM)
        lambda[[j]] <- update_betaLC_ProbitRE(C = C, Z = L[ , j], tau = kappa[[j]], UObs = VObs, URe = VRe, prior.mu = prior_mu_lambda, prior.Sigma = prior_Sigma_lambda, subjectID = subjectIDM)
        kappa[[j]] <- update_bLC_ProbitRE_fast(C = C, Omega = Theta[[j]], Z = L[ , j], phi = lambda[[j]], UObs = VObs, URe = VRe, subjectID = subjectIDM)
        Theta[[j]] <- update_OmegaLC_ProbitRE(C = C, tau = kappa[[j]], prior.scale = prior_scale_Iota, prior.df = prior_df_Iota)
      }
    }


    ### For real data analysis, save model comparison and model checking pieces

    if (sims == FALSE) {

      #! Changed 3/17/2023 to new.burn to accomodate saving only a lastIterations of samples
      if (i > (new.burn + 1)) {

        # Model comparison
        if (imputeResponse == TRUE & modelComparison == TRUE) {

          # Effective sample size
          store_ESS[i - new.burn - 1] <- get_ESS(priorPik = priorPik, Psi = Psi, Sigma = Sigma, Y = Y, XRe = XRe, subjectIDY = subjectIDY)

          ### Observed data log likelihood for BIC calculation
          #store_observed_lik_sub: n.samples-burn \times n storage
          store_observed_lik_sub[i - new.burn - 1, ] <- get_observed_lik_sub(priorPik = priorPik, betaObs = betaObs, betaSub = betaSub, Sigma = Sigma, Psi = Psi, Y = Y, XRe = XRe, XObs = XObs, XSub = XSub, subjectIDY = subjectIDY, phi = phi, Omega = Omega, D = D, URe = URe, UObs = UObs, subjectID = subjectID, lambda = lambda, Theta = Theta, M = M, VRe = VRe, VObs = VObs, subjectIDM = subjectIDM, Mvec = Mvec, modelVisit = modelVisit, modelResponse = modelResponse, monteCarlo = FALSE)
          #write.table(t(c(store_observed_lik_sub[i - burn - 1, ])), file = "store_observed_lik_sub.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)

          #store_observed_llik: n.samples-burn length storage
          store_observed_llik[i - new.burn - 1] <- sum(log(store_observed_lik_sub[i - new.burn - 1, ]))
          #write.table(store_observed_llik[i - burn - 1], file = "store_observed_llik.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)

        } # End of model comparison computations

        ### Discrepancy measure using completed datasets per Gelman et al. 2005
        ### Replicates of Y, D, M to allow posterior predictive checks based on replicated observed data per Xu et al. 2016
        # discrepancy completed is based on computing T^rep using replicated completed data (C, b, Y_obs, Y_mis) and comparing to T^obs calculated using the completed data Y_obs, Y_mis. This hold for the MNAR MNAR, MNAR MAR, and MAR MAR analyses. For MCAR analysis, this will simply be replicates of the observed data Y | D = 1, M = 1.
        #Calculation for T_{obs} is based on current iteration of C so this step must come before redrawing C
        #Calculation for T_{rep} is based on redrawing C using priorPik, redrawing b, and then redrawing Y
        model_checks <- get_model_checks(K = K, C = C, priorPik = priorPik, bSub = bSub, betaObs = betaObs, betaSub = betaSub, Psi = Psi, Sigma = Sigma, Y = Y, XRe = XRe, XObs = XObs, XSub = XSub, subjectIDY = subjectIDY, Q = Q, phi = phi, Omega = Omega, D = D, URe = URe, UObs = UObs, subjectID = subjectID, L = L, lambda = lambda, Theta = Theta, M = M, VRe = VRe, VObs = VObs, subjectIDM = subjectIDM, Mvec = Mvec, modelVisit = modelVisit, modelResponse = modelResponse)

        # Summary discrepancy measures based on replicated completed data
        write.table(t(c(model_checks[["store_T_completed"]])), file = "store_T_completed.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)

        write.table(t(c(model_checks[["store_Ydraw"]])), file = "store_Ydraw.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)

        if (modelVisit == TRUE) {
          write.table(t(c(model_checks[["store_Ddraw"]])), file = "store_Ddraw.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        }

        if (modelResponse == TRUE) {
          write.table(t(c(model_checks[["store_Mdraw"]])), file = "store_Mdraw.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        }


      } # End of i > new.burn

    } # End of sims


    ### Update latent class membership
    postClass <- update_C_fast(priorPik = priorPik, betaObs = betaObs, bSub = bSub, betaSub = betaSub, Sigma = Sigma, Psi = Psi, phi = phi, tau = tau, Omega = Omega, lambda = lambda, kappa = kappa, Theta = Theta,  Y = Y, XRe = XRe, XObs = XObs, XSub = XSub, D = D, UObs = UObs, URe = URe, M = M, Mvec = Mvec, VObs = VObs, VRe = VRe, subjectIDY = subjectIDY, subjectID = subjectID, subjectIDM = subjectIDM, modelVisit = modelVisit, modelResponse = modelResponse)

    C <- postClass[["C"]]
    pik <- postClass[["pik"]]

    ### Track iteration
    if (i %% update == 0) {
      cat("Iteration:", i, "\n")
      cat("Class size:", table(C), "\n")
    }

    ### Store samples
    if (i > burn) {

      ### Store samples after burn-in
      # Multinomial model
      store_delta[i - burn, ] <- c(delta)
      store_priorPi[i - burn, ] <- c(priorPik)
      store_pi[i - burn, ] <- c(pik)
      # Longitudinal model
      store_betaObs[i - burn, ] <- c(unlist(betaObs))
      store_betaSub[i - burn, ] <- c(unlist(betaSub))
      store_Psi[i - burn, ] <- c(unlist(Psi))
      store_Sigma[i - burn, ] <- c(Sigma)
      # AME calculation for each Y_j and storage
      ame <- as.list(1:J)
      for (j in 1:J) {
        betaObsj <- betaObs[[j]]
        # Organize betaSub[[j]] with columns indexing latent class; first p rows are for first random effects equation, second p are for second and so forth
        betaSubj <- matrix(sapply(betaSub[[j]], abind::abind), ncol = K, nrow = p*q)
        ame[[j]] <- c(apply(priorPik %*% t(betaSubj), 2, mean), apply(priorPik %*% t(betaObsj), 2, mean))

      }
      store_ame[i - burn, ] <- c(unlist(ame))

      ### For visit process

      if (modelVisit == TRUE) {
        store_phi[i - burn, ] <- c(phi)
        store_Omega[i - burn, ] <- c(Omega)
      } else {
        store_phi <- NULL
        store_Omega <- NULL
      }

      if (modelResponse == TRUE) {
        store_lambda[i - burn, ] <- c(unlist(lambda))
        store_Theta[i - burn, ] <- c(unlist(Theta))
      } else {
        store_lambda <- NULL
        store_Theta <- NULL
      }

      ### For data simulation statistics
      if (sims == TRUE) {
        store_class_weight[i - burn, ] <- as.vector(apply(priorPik, 2, mean))
      }


      if (monitor == TRUE & (i - burn) %% update == 0) {
        par(mfrow = c(2, 2))
        # Monitor selected parameters
        plot(store_betaSub[1:(i - burn), 1], type = "l")
        plot(store_betaSub[1:(i - burn), 2], type = "l")
        plot(store_betaSub[1:(i - burn), 3], type = "l")
        plot(store_betaSub[1:(i - burn), 4], type = "l")
      }

    }


    ### Changed as of 3/17/2023
    if (i > new.burn) {

      ### Data analysis post analysis
      if (sims == FALSE) {

        if (imputeResponse == TRUE) {

          ### Save imputations for inspection and to form complete data (not replicated though)
          for (j in 1:J) {

            if (any(M[ , j] == 0)) {
              fileName <- paste(paste("store", "miss", paste("Y", j, sep = ""), sep = "_"), ".txt", sep = "")
              write.table(t(c(Y[M[ , j] == 0, j])), file = fileName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
            }
          }
        }

        #### Save draws of C for posterior predictive checks post-estimation
        write.table(t(c(C)), file = "store_C.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)

        ### Save draws of random effects for posterior predictive checks post-estimation
        write.table(t(c(unlist(bSub))), file = "store_bSub.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)


      }



    } # End of i > new.burn


  } # End MCMC iteration loop


  ### Register total time
  tottime <- proc.time() - ptm
  cat("Total minutes elapsed:", tottime/60, "\n")

  ### Summary of model output

  # Samples to return in a named list
  if (sims == TRUE) {

    S <- list(store_delta = store_delta, store_pi = store_pi, store_priorPi = store_priorPi, store_betaObs = store_betaObs, store_betaSub = store_betaSub, store_Psi = store_Psi, store_Sigma = store_Sigma, store_phi = store_phi, store_Omega = store_Omega, store_lambda = store_lambda, store_Theta = store_Theta, store_ame = store_ame, store_class_weight = store_class_weight)

  } else if (sims == FALSE) {

    S <- list(store_delta = store_delta, store_pi = store_pi, store_priorPi = store_priorPi, store_betaObs = store_betaObs, store_betaSub = store_betaSub, store_Psi = store_Psi, store_Sigma = store_Sigma, store_phi = store_phi, store_Omega = store_Omega, store_lambda = store_lambda, store_Theta = store_Theta, store_ame = store_ame)

  }

  # If sims == FALSE, print to screen

  if (sims == FALSE) {

    # Model comparison under MAR or MNAR assumption, using information criteria
    if (imputeResponse == TRUE & modelComparison == TRUE) {

      ### DIC3 using predictive density estimate
      #Sum of the log of the predictive density estimate per subject obtained by averaging over MCMC iterations
      Dbar_3 <- -2 * mean(store_observed_llik)
      Dtilde_3 <-  -2 * sum(log(apply(store_observed_lik_sub, 2, mean)))
      pD_3 <- Dbar_3 - Dtilde_3
      DIC_3 <- Dbar_3 + pD_3

      ### BIC
      #Dimension of the Bayesian multivariate GMM:
      #J*K*(p + s) number of fixed effects in XObs and Sub for each latent class and outcome
      #(K-1)*m number of fixed effects in w
      #K*(J*(J + 1))/2 number of elements in the variance covairance matrix of Y
      #Psi J*K * (q * (q + 1) / 2) for each latent class and outcome
      dK <-  (J * K * (p + s)) + ((K - 1) * m) + (K * (J * (J + 1)) / 2) + (J * K * (q * (q + 1) / 2))

      if (modelVisit == TRUE) {
        dK <- dK + (K * f) + K * (q * (q + 1) / 2)
      }

      if (modelResponse == TRUE) {
        dK <- dK + (K * e) + K * (q * (q + 1) / 2)
      }

      #Average ESS
      ESS_bar <- mean(store_ESS)

      ### BIC 1 is calculated using the log likelihood after integrating out the random effects and length of subjectIDY
      # Approximate the maximum likelihood estimator by taking the parameter values that maximize the log of the observed data likelihood over MCMC samples. This is per Fruhwirth-Schnatter, S., & Pyne, S. (2010). Bayesian Inference for Finite Mixtures of Univariate and Multivariate Skew-Normal and Skew-t Distributions. Biostatistics, 11(2), 317-336. https://doi.org/10.1093/biostatistics/kxp062
      # Following SAS (Jones 2011) use the number of subjects as the sample size
      BIC_1 <- -2 * max(store_observed_llik) + dK * log(length(unique(subjectIDY)))

      ### BIC 2 is calculated using the log likelihood after integrating out the random effects and using the effective sample size
      BIC_2 <- -2 * max(store_observed_llik) + dK * log(ESS_bar)

      ### Add calculations for conditional predictive ordinate (see Gelfand and Dey 1994) and log pseudo marginal likleihood
      # n-length vector of CPOs
      CPO <- apply(store_observed_lik_sub, 2, function(x) {
        1 / mean(1 / x)
      })
      # Largest LPML indicates the better fitting model
      LPML <- sum(log(CPO))

      ### Table of model comparison
      table_model_comparison <- matrix(c(BIC_1, BIC_2, DIC_3, LPML), nrow = 1)
      colnames(table_model_comparison) <- c("BIC1", "BIC2", "DIC3", "LPML")
      rownames(table_model_comparison) <- c("value")

      ### Table of effective number of parameters
      table_pD <- matrix(c(Dbar_3, Dtilde_3, pD_3), nrow = 3, ncol = 1, byrow = FALSE)
      colnames(table_pD) <- c("DIC3")
      rownames(table_pD) <- c("Dbar", "Dtilde", "pD")

    } # End of model comparison


    # Posterior means and 95% Credible intervals after burn in
    objs <- c("store_delta", "store_betaObs", "store_Sigma", "store_betaSub", "store_Psi", "store_ame")
    if (modelVisit == TRUE) {
      objs <- c(objs, "store_phi", "store_Omega")
    }
    if (modelResponse == TRUE) {
      objs <- c(objs, "store_lambda", "store_Theta")
    }

    temp <- lapply(S[objs], function(X){
      apply(X, 2, function(x){c(round(mean(x), digits = 4), round(quantile(x, probs = c(0.025, 0.975)), digits = 4))})
    })
    table_posterior_summary <- t(do.call("cbind", temp))
    colnames(table_posterior_summary) <- c("Post. Mean", "2.5 %", "97.5%")


    # Posterior classification
    postPi <- matrix(c(apply(S[["store_pi"]], 2, mean)), nrow = n, ncol = K, byrow = FALSE)
    postPredC <- apply(postPi, 1, which.max)

    #No. subjects in class k with post prob > 0.95, .90, .8
    sub_gr95 <- sub_gr90 <- sub_gr80 <-  rep(NA, K)
    for (k in 1:K) {
      temp <- postPi[postPredC == k, k]
      sub_gr95[k] <- sum(temp >= .95)
      sub_gr90[k] <- sum(temp >= .90)
      sub_gr80[k] <- sum(temp >= .80)
    }

    #Mean and median posterior probability by latent class
    meanPostPi <- medianPostPi <-  rep(NA, K)
    for (k in 1:K) {
      temp <- postPi[postPredC == k, k]
      meanPostPi[k] <- round(mean(temp), digits = 2)
      medianPostPi[k] <- round(median(temp), digits = 2)
    }

    table_latent_class_unit_assignment <- rbind(table(postPredC), sub_gr95, sub_gr90, sub_gr80, meanPostPi, medianPostPi)
    rownames(table_latent_class_unit_assignment) <- c("Predicted class size", "No. subjects with probability at least 0.95",
                                                      "No. subjects with probability at least 0.90",
                                                      "No. subjects with probability at least 0.80",
                                                      "Mean probability", "Median probability")
    colnames(table_latent_class_unit_assignment) <- c(paste("Class", 1:K))


    # Evaluation of label switching using Stephen's method
    temp <- array(S[["store_pi"]], dim = c(dim(store_pi)[1], n, K))
    perm_res <- label.switching::stephens(temp)

    #Number of rows not equal to the identity
    identity <- matrix(rep(c(1:K)), nrow = dim(store_pi)[1], ncol = K, byrow = TRUE)
    if (isTRUE(all.equal(perm_res[[1]], identity))) {
     cat("No evidence of label switching problem using Stephen's method from the label.switching package", "\n", "\n")
    } else {
     cat("Warning: Stephen's method from the label.switching package has detected label switching", "\n", "\n")
    }

    # Print output to screen
    cat("Background information", "\n")
    cat("Number of subjects:", n, "\n")
    cat("Number of observations:", length(subjectIDY), "\n")
    cat("Number of latent classes:", K, "\n", "\n")

    cat("Posterior latent class assignment:", "\n")
    print(table_latent_class_unit_assignment)

    cat("\n", "Reference class in latent class membership model:", 1, "\n")
    cat("Posterior means and 95% credible intervals:", "\n")
    print(table_posterior_summary)
    cat("\n", "Footnotes for posterior means and 95% credible intervals:", "\n")
    cat("Elements of variance-covariances are indexed by row and then column.",   "\n")
    cat("RE indexes the random effects equations. If there is only a random intercept, then this will be RE1.", "\n")
    cat("AME indicates the population-averaged regression coefficients.",   "\n")

    cat("\n", "Model comparison statistics:", "\n")
    if (imputeResponse == TRUE & modelComparison == TRUE) {
      print(table_model_comparison)
      cat("\n", "Footnotes for model comparison:", "\n")
      cat("BIC1:", "Computed using number of observations equal to the number of observed clinic visits",   "\n")
      cat("BIC2:", "Computed using number of observations equal to the effective sample size from the longitudinal health outcomes model",   "\n")
      cat("DIC calculation details", "\n")
      print(table_pD)
    } else {
      cat("Not available", "\n")
    }

  } # End summary if sims = FALSE


  ### Return list of MCMC samples
  return(S)

}

