#' Calculate Nagelkerke pseudo R-squared for a linear model
#'
#' This function computes the Nagelkerke pseudo R-squared for a linear model (`lm` object)
#' by comparing the log-likelihood of the fitted model to that of the null model (intercept-only).
#'
#' @param model An object of class \code{lm}.
#'
#' @return A numeric value representing the Nagelkerke pseudo R-squared.
#'
#' @details
#' Nagelkerke R² is a normalized version of the likelihood ratio R², scaled to have a maximum of 1.
#' It is commonly used for generalized linear models but can also be applied to linear models.
#' 
#' The formula used is: 
#' \deqn{R^2_{Nagelkerke} = \frac{1 - \exp\left(\frac{2}{n} (LL_{null} - LL_{model})\right)}{1 - \exp\left(\frac{2}{n} LL_{null}\right)}}
#' where \eqn{LL_{model}} is the log-likelihood of the fitted model and \eqn{LL_{null}} is the log-likelihood of the null model.
#'
#' @references Nagelkerke, N. J. D. (1991). A note on a general definition of the coefficient of determination. \emph{Biometrika}, 78(3), 691–692.
#'
#' @examples
#' data(mtcars)
#' fit <- lm(mpg ~ wt + hp, data = mtcars)
#' nagelkerke_r2_lm(fit)
#'
#' @export
nagelkerke_r2_lm <- function(model) {
  if (!inherits(model, "lm")) {
    stop("Input must be a linear model object of class 'lm'.")
  }
  
  # Extract response variable from model frame
  mf <- model.frame(model)
  y <- model.response(mf)
  
  # Number of observations
  n <- length(y)
  
  # Fit null model (intercept only) with same data
  null_model <- lm(y ~ 1)
  
  # Compute log-likelihoods
  ll_model <- logLik(model)
  ll_null <- logLik(null_model)
  
  # Nagelkerke pseudo-R² calculation
  r2_nagelkerke <- (1 - exp((2 / n) * (ll_null - ll_model))) / 
                   (1 - exp((2 / n) * as.numeric(ll_null)))
  
  return(as.numeric(r2_nagelkerke))
}
