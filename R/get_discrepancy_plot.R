#' @title Get discrepancy plot using completed and replicated completed data.
#'
#' @description Overall adequacy of the fitted model can be examine using Bayesian posterior predictive p-values. At each MCMC iteration, a discrepancy measure is computed using the completed and replicated completed data. The Bayesian predictive p-value denotes the probability that the discrepancy measure under the replicated data is greater than that under the observed data. A p-value near 0.5 indicates adequate model fit, while a p-value outside the range of 0.05 and 0.95 is considered to suggest a lack of model fit.
#'
#' The discrepancy measure used in this function is a multivariate mean squared error, computed as in:
#'Daniels, M.J. and Hogan, J.W. (2008) Missing Data in Longitudinal Studies: Strategies for Bayesian Modeling and Sensitivity Analysis. Chapman and Hall/CRC, Boca Raton, p. 67.
#'
#' A reference for using the completed data versus replicated completed data is in:
#' Gelman, A., Mechelen, I. V., Verbeke, G., Heitjan, D.F. and Meulders, M. (2005) Multiple Imputation for Model Checking: Completed-Data Plots with Missing and Latent Data. Biometrics, 61, 74-85.
#' @param discrepancySamples Discrepancy samples written to \code{store_T_completed.txt} in \code{MVNYBinaryMiss}, when \code{sims = FALSE}.
#' @export
#' @return A scatter plot of the discrepancy measure calculated using the completed data versus the replicated completed data. The plot is annotated with the posterior predictive p-value.
get_discrepancy_plot <- function(discrepancySamples) {

  theme_plot <- ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold", size = 20),  axis.title.y = ggplot2::element_text(size = 10, face = "bold"), axis.text.y = ggplot2::element_text(size = 8, face = "bold"), axis.title.x = ggplot2::element_text(size = 10, face = "bold"), axis.text.x = ggplot2::element_text(size = 8, face = "bold"), legend.title = ggplot2::element_text(colour = "black", size = 10, face = "bold"), legend.text = ggplot2::element_text(colour = "black", size = 8, face = "bold"), legend.key.size = ggplot2::unit(1, "cm"), strip.text.x = ggplot2::element_text(size = 10, face = "bold"))

  # Column names
  colnames(discrepancySamples) <- c("TObs", "TRep")

  pvalue <- with(discrepancySamples, round(sum(TRep > TObs) / dim(discrepancySamples)[1], digits = 2))

  # Prepare text
  grob <- grid::grobTree(grid::textGrob(paste("Bayesian p-value", "=", pvalue, sep = " "), x = 0.1,  y = 0.95, hjust = 0, gp = grid::gpar(col = "black", fontsize = 10, fontface = "bold")))

  # Make plot of observed T versus replicated T
  p <- ggplot2::ggplot(discrepancySamples, ggplot2::aes(x = TObs, y = TRep)) + ggplot2::geom_point(shape = 1)
  p <- p + ggplot2::geom_abline(intercept = 0, slope = 1, color = "red", size = 0.5) + ggplot2::annotation_custom(grob) + ggplot2::scale_y_continuous(name = "Replicated T") + ggplot2::scale_x_continuous(name = "Completed T") + theme_plot

  print(p)

}
