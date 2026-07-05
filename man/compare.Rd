\name{compare}
\alias{compare}
\alias{print.psychonetrics_compare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model comparison
}
\description{
This function will print a table comparing multiple models on chi-square, AIC and BIC. When all compared models are fitted with a robust estimator, a Satorra-Bentler-family scaled chi-square difference test is added.
}
\usage{
compare(\dots, scaled.test.method = c("auto", "satorra.bentler.2001",
        "satorra.bentler.2010", "satorra.2000"))

\method{print}{psychonetrics_compare}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
Any number of \code{psychonetrics} models. Can be named to change the rownames of the output.
}
\item{scaled.test.method}{
The scaled chi-square difference test to use when all compared models carry a robust scaling factor (see Details). \code{"auto"} (the default) selects the appropriate test from the models' scaled-test type, matching \pkg{lavaan}'s \code{lavTestLRT}: the exact Satorra (2000) scaled-and-shifted test for the scaled-shifted (\code{MLMV}, \code{WLSMV}) and mean-and-variance-adjusted (\code{MLMVS}) estimators, and the Satorra-Bentler (2001) mean-scaled test for \code{MLM} / \code{MLR}. The remaining options force a specific test: \code{"satorra.bentler.2001"} (a mean-scaled test), \code{"satorra.bentler.2010"} (a strictly positive mean-scaled variant), or \code{"satorra.2000"} (an exact, scaled-and-shifted test, appropriate for the scaled-shifted estimators such as \code{MLMV} and \code{WLSMV}).
}
\item{x}{Output of the \code{compare} function.}
}
\details{
When every compared model is fitted with a robust estimator (\code{MLM}, \code{MLMV}, \code{MLMVS} or \code{MLR}, or one of the least-squares estimators \code{DWLS}, \code{WLS} and \code{ULS}), two additional columns are added: \code{Chisq_diff_scaled} and \code{p_value_scaled}. These contain the Satorra-Bentler-family scaled chi-square difference test for each adjacent (nested) pair of models, sorted by degrees of freedom. The ordinary (unscaled) \code{Chisq_diff}, \code{DF_diff} and \code{p_value} columns are always present and unchanged. When the models are not fitted with a robust estimator (e.g. ordinary \code{ML} or \code{FIML}), the scaled columns are omitted entirely.

For an adjacent nested pair with the more constrained model \eqn{M_0} (degrees of freedom \eqn{r_0}, unscaled chi-square \eqn{T_0}, scaling factor \eqn{c_0}) nested in the less constrained model \eqn{M_1} (\eqn{r_1}, \eqn{T_1}, \eqn{c_1}), and \eqn{m = r_0 - r_1}:
\itemize{
\item \code{"satorra.bentler.2001"} computes \eqn{c_d = (r_0 c_0 - r_1 c_1)/m} and \eqn{T_d = (T_0 - T_1)/c_d} (Satorra & Bentler, 2001).
\item \code{"satorra.bentler.2010"} replaces \eqn{c_1} by the scaling factor of \eqn{M_1} evaluated at \eqn{M_0}'s estimates, guaranteeing a positive statistic (Satorra & Bentler, 2010).
\item \code{"satorra.2000"} computes an exact, scaled-and-shifted statistic from the asymptotic covariance of the sample statistics (Satorra, 2000; Asparouhov & Muthen, 2010); this is the appropriate difference test for the scaled-shifted estimators \code{MLMV} and \code{WLSMV}.
}
The statistic is set to \code{NA} (with a warning) when the implied scaling factor is non-positive or a required building block cannot be computed.
}
\value{
A data frame with chi-square values, degrees of freedoms, RMSEAs, AICs, and BICs, and---when all models are fitted with a robust estimator---scaled chi-square difference columns.
}
\references{
Satorra, A. (2000). Scaled and adjusted restricted tests in multi-sample analysis of moment structures. In R. D. H. Heijmans, D. S. G. Pollock, & A. Satorra (Eds.), \emph{Innovations in multivariate statistical analysis} (pp. 233-247). Kluwer Academic Publishers.

Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square test statistic for moment structure analysis. \emph{Psychometrika}, 66(4), 507-514.

Satorra, A., & Bentler, P. M. (2010). Ensuring positiveness of the scaled difference chi-square test statistic. \emph{Psychometrika}, 75(2), 243-248.

Asparouhov, T., & Muthen, B. (2010). \emph{Simple second order chi-square correction}. Mplus technical appendix.
}
\author{
Sacha Epskamp
}
\examples{
\dontrun{
library("dplyr")

# Load data:
data(StarWars)

# Originally a partial measurement invariance model in
# https://psychonetrics.org/files/Epskampetal2017.pdf:
mod1 <- lvm(StarWars, lambda = matrix(1, 10), vars = paste0("Q", 1:10),
            identification = "variance", latents = "Force",
            estimator = "MLM") \%>\% runmodel

# A model fixing two prequel loadings to zero:
mod2 <- mod1 \%>\% fixpar("lambda", row = "Q4", col = "Force") \%>\%
  fixpar("lambda", row = "Q5", col = "Force") \%>\% runmodel

# Compare with the scaled (Satorra-Bentler 2001) difference test:
compare(full = mod1, restricted = mod2)
}
}
