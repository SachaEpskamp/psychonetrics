\name{find_penalized_lambda}
\alias{find_penalized_lambda}
\alias{compute_lambda_max}
\alias{compute_penalized_ebic}
\alias{apply_beta_min}
\title{
Automated Lambda Selection for Penalized Estimation
}
\description{
Selects the optimal penalty strength (lambda) for penalized maximum likelihood (PML) or penalized full-information maximum likelihood (PFIML) models via EBIC-based grid search. A log-spaced grid of lambda values is searched from \code{lambda_max} (the smallest lambda for which all penalized parameters are exactly zero, derived from KKT conditions) down to \code{lambda_max * lambda_min_ratio}. Each lambda is optimized using ISTA (Iterative Shrinkage-Thresholding Algorithm) with warm starts from the previous solution, and the model with the best EBIC is selected. After selection, a \code{beta_min} threshold is applied to set near-zero penalized parameters exactly to zero.

This function is called automatically by \code{\link{runmodel}} when the model contains auto-select penalty parameters (i.e., \code{penalty_lambda = NA} in the parameter table, which is the default when using \code{estimator = "PML"} or \code{estimator = "PFIML"}). It can also be called directly for more control over the search.
}
\usage{
find_penalized_lambda(model, criterion = "ebic.5", nLambda = 50,
                      lambda_min_ratio = 1e-4, beta_min,
                      verbose)
}
\arguments{
  \item{model}{
A \code{psychonetrics} model with \code{estimator = "PML"} or \code{estimator = "PFIML"}. The model should contain parameters with \code{penalty_lambda = NA} (auto-select), which is the default when specifying \code{estimator = "PML"} or \code{"PFIML"} in model constructors such as \code{\link{ggm}}, \code{\link{lvm}}, etc.
}
  \item{criterion}{
Character string specifying the information criterion used to select the best lambda. One of:
\describe{
  \item{\code{"bic"}}{BIC (equivalent to EBIC with gamma = 0)}
  \item{\code{"ebic.25"}}{EBIC with gamma = 0.25}
  \item{\code{"ebic.5"}}{EBIC with gamma = 0.5 (default; recommended for network models)}
  \item{\code{"ebic.75"}}{EBIC with gamma = 0.75}
  \item{\code{"ebic1"}}{EBIC with gamma = 1.0}
}
Higher gamma values penalize model complexity more strongly and tend to produce sparser solutions. Gamma = 0.5 is widely recommended for network models (Foygel & Drton, 2010).
}
  \item{nLambda}{
Integer, the number of lambda values in the search grid. Default is 50. The grid is log-spaced between \code{lambda_max} and \code{lambda_max * lambda_min_ratio}.
}
  \item{lambda_min_ratio}{
Numeric, the ratio of the smallest to the largest lambda in the grid. Default is \code{1e-4}, so that \code{lambda_min = lambda_max * 1e-4}.
}
  \item{beta_min}{
Numeric threshold for zeroing out small penalized parameter estimates. Parameters with absolute value below \code{beta_min} are set to zero, both during EBIC evaluation (for parameter counting) and in the final model. Default is \code{sqrt(log(p) / n)} where \code{p} is the number of penalized parameters and \code{n} is the total sample size.
}
  \item{verbose}{
Logical, should progress messages be printed? Defaults to the model's \code{verbose} setting.
}
}
\details{
The search proceeds as follows:

\enumerate{
  \item \strong{Compute lambda_max}: The analytical upper bound is derived from KKT (Karush-Kuhn-Tucker) conditions at the null-penalized model (all penalized parameters fixed at zero). Specifically, \code{lambda_max = max(|grad_j|) / alpha} where \code{grad_j} are the gradient elements for penalized parameters at the null solution and \code{alpha} is the elastic net mixing parameter.

  \item \strong{Build lambda grid}: A log-spaced sequence of \code{nLambda} values from \code{lambda_max} down to \code{lambda_max * lambda_min_ratio}.

  \item \strong{Grid search with warm starts}: For each lambda (from largest to smallest), the model is optimized using ISTA. The solution from the previous lambda serves as the starting point (warm start), which substantially reduces computation time. The EBIC (Extended Bayesian Information Criterion) is computed using the unpenalized log-likelihood (Convention A), counting only parameters with \code{|estimate| >= beta_min} as non-zero.

  \item \strong{Beta-min thresholding}: Penalized parameters with \code{|estimate| < beta_min} are set exactly to zero in the final model.
}

The EBIC is computed as:
\deqn{EBIC = -2 \cdot LL + k \cdot \log(n) + 4 \cdot k \cdot \gamma \cdot \log(p_v)}{EBIC = -2*LL + k*log(n) + 4*k*gamma*log(nVar)}
where \eqn{LL} is the unpenalized log-likelihood, \eqn{k} is the effective number of parameters (unpenalized free parameters plus non-zero penalized parameters), \eqn{n} is the sample size, \eqn{\gamma}{gamma} is the EBIC tuning parameter, and \eqn{p_v}{nVar} is the number of observed variables.

The returned model stores metadata about the search in \code{model@optim$lambda_search}, a list containing: \code{criterion}, \code{gamma}, \code{best_lambda}, \code{best_criterion}, \code{beta_min}, \code{lambda_max}, \code{nLambda}, and \code{n_thresholded}.
}
\value{
An object of class \code{psychonetrics} (\link{psychonetrics-class}) with:
\itemize{
  \item The selected lambda assigned to all auto-select penalty parameters
  \item Parameters optimized at the best lambda
  \item Near-zero penalized parameters set to exactly zero (via \code{beta_min})
  \item Search metadata stored in \code{model@optim$lambda_search}
}

After running \code{find_penalized_lambda}, use \code{\link{refit}} to obtain valid standard errors and fit indices via post-selection inference.
}
\references{
Foygel, R., & Drton, M. (2010). Extended Bayesian Information Criteria for Gaussian Graphical Models. In \emph{Advances in Neural Information Processing Systems} (Vol. 23).
}
\author{
Sacha Epskamp
}
\seealso{
\code{\link{penalize}} for manually setting penalty parameters,
\code{\link{runmodel}} for running models (calls \code{find_penalized_lambda} automatically when auto-select penalties are present),
\code{\link{refit}} for post-selection inference after penalized estimation.
}
\examples{
\donttest{
# Load bfi data from psych package:
library("psychTools")
library("dplyr")
data(bfi)

ConsData <- bfi \%>\%
  select(A1:A5) \%>\%
  na.omit

vars <- names(ConsData)

# Automatic lambda selection (called implicitly by runmodel):
mod <- ggm(ConsData, vars = vars, estimator = "PML")
mod <- mod \%>\% runmodel  # automatically searches for best lambda

# Check lambda search results:
mod@optim$lambda_search

# Post-selection refit for standard errors:
mod_refit <- mod \%>\% refit
mod_refit \%>\% parameters

# Direct call with custom criterion:
mod2 <- ggm(ConsData, vars = vars, estimator = "PML")
mod2 <- find_penalized_lambda(mod2, criterion = "ebic1",
                              nLambda = 100)
}
}
