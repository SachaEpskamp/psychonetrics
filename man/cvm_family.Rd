\name{cvm}
\alias{cvm}
\alias{cnm}

\title{
Composite variable family of psychonetrics models
}
\description{
This is the family of models that models the data as a structural equation model (SEM) with both latent variables (common factors) and composites, following the factor- and composite-based SEM framework of Schamberger, Schuberth, Henseler and Rosseel (2025). A latent variable is measured by its indicators with error (\code{lambda}), whereas a composite is an exact weighted sum of its indicators (\code{weights}). Latent variables and composites (jointly termed constructs) can be related to each other through regressions (\code{beta}) and (residual) covariances or a network. As in the \code{\link{lvm}} family, the variance--covariance structures in the model can be modeled in different ways, including as Gaussian graphical models (Epskamp, Rhemtulla and Borsboom, 2017). The \code{latent} argument sets the model for the constructs, the \code{residual} argument sets the model for the residuals of the indicators of latent variables, and the \code{composite} argument sets the model for the indicators of composites (the matrices with the suffix \code{_c}; see Details). The wrapper \code{cnm()} sets \code{composite = "ggm"} for the composite network model (CNM), which models a network of indicators within every composite. This family is new in version 0.18 and experimental: arguments and defaults may still change.
}
\usage{
cvm(data, lambda, weights, latent = c("cov", "ggm"),
                   residual = c("cov", "chol", "prec", "ggm", "cor"),
                   composite = c("cov", "chol", "prec", "ggm", "cor"),
                   sigma_zeta = "full", omega_zeta = "full",
                   delta_zeta = "full", sigma_epsilon = "diag",
                   kappa_epsilon = "diag", omega_epsilon = "zero",
                   lowertri_epsilon = "diag", delta_epsilon = "diag",
                   rho_epsilon = "zero", SD_epsilon = "diag",
                   sigma_c = "full", kappa_c = "full", omega_c = "full",
                   lowertri_c = "full", delta_c = "full", rho_c = "full",
                   SD_c = "full", beta = "zero", nu, nu_eta,
                   identify = TRUE, identification = c("loadings",
                   "variance"), bad_marker_crit = 0.1, vars, latents, groups,
                   groupvar, covs,
                   means, nobs, missing = "auto", equal = "none",
                   baseline_saturated = TRUE, estimator = "default",
                   optimizer, storedata = FALSE, WLS.W, covtype =
                   c("choose", "ML", "UB"), standardize = c("none", "z",
                   "quantile"), sampleStats, verbose = FALSE,
                   bootstrap = FALSE, boot_sub, boot_resample,
                   penalty_lambda = NA, penalty_alpha = 1,
                   penalize_matrices)

cnm(\dots)
}
\arguments{
  \item{data}{
A data frame encoding the data used in the analysis. Can be missing if \code{covs} and \code{nobs} are supplied.
}
  \item{lambda}{
A model matrix encoding the factor loading structure of the latent variables. Each row indicates an observed variable and each column a construct. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constraints. Columns of composites must be zero. Can be omitted if the model only contains composites. For multiple groups, this argument can be an array with each slice encoding such a matrix.
}
  \item{weights}{
A model matrix encoding which observed variables form which composite. Each row indicates an observed variable and each column a construct. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constraints. Columns of latent variables must be zero. Every construct must have elements in either \code{lambda} (a latent variable) or \code{weights} (a composite), and a variable can only be an indicator of a single composite and cannot also load on a latent variable. Can be omitted if the model only contains latent variables (in which case the model equals a model from the \code{\link{lvm}} family). A composite with a single indicator represents that observed variable itself, which allows observed variables to be used in the structural part of the model.
}
  \item{latent}{
The type of model used for the (residual) variances and covariances of the constructs: \code{"cov"} (default) to model a variance--covariance matrix (\code{sigma_zeta}), or \code{"ggm"} to model a Gaussian graphical model (\code{omega_zeta} and \code{delta_zeta}). See Details for the restrictions on composites.
}
  \item{residual}{
The type of model used for the residuals of the indicators of latent variables: \code{"cov"} (default), \code{"chol"}, \code{"prec"}, \code{"ggm"} or \code{"cor"}, as in the \code{\link{lvm}} family.
}
  \item{composite}{
The type of model used for the variances and covariances of the indicators of composites (the \code{_c} matrices): \code{"cov"} (default) to model a variance--covariance matrix (\code{sigma_c}), \code{"chol"} to model a Cholesky decomposition (\code{lowertri_c}), \code{"prec"} to model a precision matrix (\code{kappa_c}), \code{"ggm"} to model a Gaussian graphical model (\code{omega_c} and \code{delta_c}), and \code{"cor"} to model a correlation matrix with separate standard deviations (\code{rho_c} and \code{SD_c}).
}
  \item{sigma_zeta}{
Only used when \code{latent = "cov"}. Either \code{"full"} to estimate every element freely, \code{"diag"} to only include diagonal elements, or a matrix of the dimensions construct x construct with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. The diagonal elements of composites are not estimated but computed (see Details). For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega_zeta}{
Only used when \code{latent = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"zero"} to set all elements to zero, or a matrix of the dimensions construct x construct with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_zeta}{
Only used when \code{latent = "ggm"}. Either \code{"diag"} or \code{"zero"}, or a matrix of the dimensions construct x construct with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. The elements of composites are not estimated but computed (see Details). For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{sigma_epsilon}{
Only used when \code{residual = "cov"}. Either \code{"full"} to estimate every element freely, \code{"diag"} to only include diagonal elements, or a matrix with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. This matrix, as all other \code{_epsilon} matrices, has one row and column for every observed variable that is not an indicator of a composite, in the order in which these appear in the data. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_epsilon}{
Only used when \code{residual = "prec"}. Either \code{"full"}, \code{"diag"}, or a matrix as described for \code{sigma_epsilon}.
}
  \item{omega_epsilon}{
Only used when \code{residual = "ggm"}. Either \code{"full"}, \code{"zero"}, or a matrix as described for \code{sigma_epsilon}.
}
  \item{lowertri_epsilon}{
Only used when \code{residual = "chol"}. Either \code{"full"}, \code{"diag"}, or a matrix as described for \code{sigma_epsilon}.
}
  \item{delta_epsilon}{
Only used when \code{residual = "ggm"}. Either \code{"diag"}, \code{"zero"}, or a matrix as described for \code{sigma_epsilon}.
}
  \item{rho_epsilon}{
Only used when \code{residual = "cor"}. Either \code{"full"}, \code{"zero"}, or a matrix as described for \code{sigma_epsilon}.
}
  \item{SD_epsilon}{
Only used when \code{residual = "cor"}. Either \code{"diag"}, \code{"zero"}, or a matrix as described for \code{sigma_epsilon}.
}
  \item{sigma_c}{
Only used when \code{composite = "cov"}. The variances and covariances of the indicators of composites. Either \code{"full"} (default) to freely estimate all variances and all covariances between indicators of the same composite, \code{"diag"} to only include the variances, or a matrix with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. This matrix, as all other \code{_c} matrices, has one row and column for every indicator of a composite, in the order in which these appear in the data. Elements connecting indicators of different composites are always fixed to zero (see Details). For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_c}{
Only used when \code{composite = "prec"}. Either \code{"full"}, \code{"diag"}, or a matrix as described for \code{sigma_c}.
}
  \item{omega_c}{
Only used when \code{composite = "ggm"}. The partial correlation network of the indicators within every composite. Either \code{"full"} (default) to estimate all edges between indicators of the same composite, \code{"zero"} to set all elements to zero, or a matrix as described for \code{sigma_c}.
}
  \item{lowertri_c}{
Only used when \code{composite = "chol"}. Either \code{"full"}, \code{"diag"}, or a matrix as described for \code{sigma_c}.
}
  \item{delta_c}{
Only used when \code{composite = "ggm"}. Either \code{"full"} or \code{"diag"} (both estimate all scaling parameters), \code{"zero"}, or a matrix as described for \code{sigma_c}.
}
  \item{rho_c}{
Only used when \code{composite = "cor"}. Either \code{"full"}, \code{"zero"}, or a matrix as described for \code{sigma_c}.
}
  \item{SD_c}{
Only used when \code{composite = "cor"}. Either \code{"full"} or \code{"diag"} (both estimate all standard deviations), \code{"zero"}, or a matrix as described for \code{sigma_c}.
}
  \item{beta}{
A model matrix encoding the structural relations between the constructs (latent variables and composites), with element \code{[i, j]} encoding the effect of construct \code{j} on construct \code{i}. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{nu}{
Optional vector encoding the intercepts of the indicators of latent variables and the means of the indicators of composites. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free elements, and higher integers to indicate equality constraints. For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
}
  \item{nu_eta}{
Optional vector encoding the intercepts of the constructs. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free intercepts, and higher integers to indicate equality constraints. The intercepts of composites are not estimated but computed (see Details). For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
}
  \item{identify}{
Logical, should the model be automatically identified?
}
  \item{identification}{
Type of identification used for the latent variables. \code{"loadings"} to fix the first factor loadings to 1, and \code{"variance"} to fix the (residual) variance or scaling of the latent variables to 1. Composites are always scaled by fixing the weight of one indicator to 1 (see \code{bad_marker_crit} and Details).
}
  \item{bad_marker_crit}{
A composite is scaled by fixing the weight of its first indicator to 1, unless that indicator is a bad choice: if the starting value of the standardized weight of the first indicator is smaller than \code{bad_marker_crit} times the largest standardized starting weight of the composite, the indicator with the largest starting weight is used to scale the composite instead, and a message is shown. Set to 0 to always use the first indicator. See Details.
}
  \item{vars}{
An optional character vector encoding the variables used in the analysis. Must equal names of the dataset in \code{data}.
}
  \item{latents}{
An optional character vector with names of the constructs (latent variables and composites).
}
  \item{groups}{
Deprecated. Use \code{groupvar} instead. An optional string indicating the name of the group variable in \code{data}.
}
  \item{groupvar}{
An optional string indicating the name of the group variable in \code{data}.
}
  \item{covs}{
A sample variance--covariance matrix, or a list/array of such matrices for multiple groups. Make sure \code{covtype} argument is set correctly to the type of covariances used.
}
  \item{means}{
A vector of sample means, or a list/matrix containing such vectors for multiple groups.
}
  \item{nobs}{
The number of observations used in \code{covs} and \code{means}, or a vector of such numbers of observations for multiple groups.
}
  \item{missing}{
How should missingness be handled in computing the sample covariances and number of observations when \code{data} is used. Can be \code{"auto"} (default) for automatic detection, \code{"listwise"} for listwise deletion, or \code{"pairwise"} for pairwise deletion. When \code{"auto"}, the function checks for missing data and switches ML to FIML.
}
  \item{equal}{
A character vector indicating which matrices should be constrained equal across groups.
}
  \item{baseline_saturated}{
A logical indicating if the baseline and saturated model should be included. Mostly used internally and NOT Recommended to be used manually.
}
  \item{estimator}{
The estimator to be used. Currently implemented are \code{"ML"} for maximum likelihood estimation, \code{"FIML"} for full-information maximum likelihood estimation, \code{"PML"} for penalized maximum likelihood estimation, \code{"PFIML"} for penalized full-information maximum likelihood estimation, \code{"ULS"} for unweighted least squares estimation, \code{"WLS"} for weighted least squares estimation, and \code{"DWLS"} for diagonally weighted least squares estimation. Defaults to \code{"ML"}. If missing data is detected, \code{"ML"} is automatically switched to \code{"FIML"} and \code{"PML"} to \code{"PFIML"}.
}
  \item{optimizer}{
The optimizer to be used. Can be one of \code{"nlminb"} (the default R \code{nlminb} function), \code{"ucminf"} (from the \code{optimr} package), \code{"nloptr_TNEWTON"} (preconditioned truncated Newton via \code{nloptr}), and \code{"LBFGS++"} (pure C++ L-BFGS-B). Defaults to \code{"nlminb"}.
}
  \item{storedata}{
Logical, should the raw data be stored? Needed for bootstrapping (see \code{\link{bootstrap}}).
}
  \item{WLS.W}{
The weights matrix used in WLS estimation (experimental)
}
  \item{covtype}{
If 'covs' is used, this is the type of covariance (maximum likelihood or unbiased) the input covariance matrix represents. Set to \code{"ML"} for maximum likelihood estimates (denominator n) and \code{"UB"} to unbiased estimates (denominator n-1). The default will try to find the type used, by investigating which is most likely to result from integer valued datasets.
}
  \item{standardize}{
Which standardization method should be used? \code{"none"} (default) for no standardization, \code{"z"} for z-scores, and \code{"quantile"} for a non-parametric transformation to the quantiles of the marginal standard normal distribution.
}
  \item{sampleStats}{
An optional sample statistics object. Mostly used internally.
}
  \item{verbose}{
Logical, should progress be printed to the console?
}
  \item{bootstrap}{
Should the data be bootstrapped? If \code{TRUE} the data are resampled and a bootstrap sample is created. These must be aggregated using \code{\link{aggregate_bootstraps}}! Can be \code{TRUE} or \code{FALSE}. Can also be \code{"nonparametric"} (which sets \code{boot_sub = 1} and \code{boot_resample = TRUE}) or \code{"case"} (which sets \code{boot_sub = 0.75} and \code{boot_resample = FALSE}).
}
  \item{boot_sub}{
Proportion of cases to be subsampled (\code{round(boot_sub * N)}).
}
  \item{boot_resample}{
Logical, should the bootstrap be with replacement (\code{TRUE}) or without replacement (\code{FALSE})
}
  \item{penalty_lambda}{
Numeric penalty strength for penalized ML estimation (PML/PFIML). \code{NA} (default) triggers automatic selection via EBIC-based grid search when a penalized estimator is used; set to a specific numeric value to use a fixed penalty strength (0 = no penalty). See \code{\link{find_penalized_lambda}} and \code{\link{penalize}}.
}
  \item{penalty_alpha}{
Elastic net mixing parameter: 1 = LASSO (default), 0 = ridge.
}
  \item{penalize_matrices}{
Character vector of matrix names to penalize. If missing, the network structures in the model are penalized (\code{omega_c} or \code{kappa_c}, \code{omega_epsilon} or \code{kappa_epsilon}, and \code{omega_zeta}); the weights are never penalized by default. See Details.
}
  \item{\dots}{
Arguments sent to \code{cvm}
}
}
\details{
\strong{The model.} Every observed variable is either an indicator of latent variables or an indicator of one composite. For the indicators of latent variables \eqn{\boldsymbol{y}_f}{y_f}, the indicators of composites \eqn{\boldsymbol{y}_c}{y_c}, and the constructs \eqn{\boldsymbol{\eta}}{eta} (containing both the latent variables \eqn{\boldsymbol{\eta}_f}{eta_f} and the composites \eqn{\boldsymbol{\eta}_c}{eta_c}), the model is:
\deqn{\boldsymbol{y}_f = \boldsymbol{\nu}_f + \boldsymbol{\Lambda} \boldsymbol{\eta}_f + \boldsymbol{\varepsilon}}{y_f = nu_f + lambda * eta_f + epsilon}
\deqn{\boldsymbol{\eta}_c = \boldsymbol{W}^{\top} \boldsymbol{y}_c}{eta_c = t(weights) * y_c}
\deqn{\boldsymbol{\eta} = \boldsymbol{\nu}_\eta + \boldsymbol{B} \boldsymbol{\eta} + \boldsymbol{\zeta}}{eta = nu_eta + beta * eta + zeta}
in which \eqn{\boldsymbol{\Lambda}}{lambda} contains the factor loadings (\code{lambda}), \eqn{\boldsymbol{W}}{weights} the composite weights (\code{weights}) and \eqn{\boldsymbol{B}}{beta} the regressions between the constructs (\code{beta}). In contrast to a latent variable, a composite has no measurement error: it is fully determined by its indicators. The composite model does not explain the covariances between the indicators of a composite. Instead, the composite acts as the single weighted sum through which its indicators relate to all other variables in the model (Dijkstra, 2017; Schuberth, Henseler and Dijkstra, 2018; Henseler, 2021).

\strong{The \code{_c} matrices: covariances of the composite indicators.} The matrices with the suffix \code{_c} (\code{sigma_c}, \code{kappa_c}, \code{omega_c}, \code{delta_c}, \code{lowertri_c}, \code{rho_c} and \code{SD_c}) model the variance--covariance matrix of the \emph{indicators of composites} (termed composite indicators by Bollen and Bauldry, 2011), in the same way in which the suffix \code{_zeta} refers to the constructs and the suffix \code{_epsilon} refers to the residuals of the indicators of latent variables. Note that the suffix \code{_c} thus does not refer to the relations \emph{between} composites, which are modeled in \code{beta} and the \code{_zeta} matrices. The \code{_c} matrices have one row and column for every composite indicator and are block-diagonal, with one block per composite: variances and covariances of indicators of the same composite are parameters, which are freely estimated by default (\code{"full"}), whereas covariances between indicators of different composites are implied by the relations between the composites and are therefore always fixed to zero in these matrices. With \code{composite = "ggm"} (or the wrapper \code{cnm()}), \code{omega_c} encodes a partial correlation network of the indicators within every composite. Elements of this network can be constrained to zero, for example using \code{\link{prune}} and \code{\link{stepup}}; these functions never add edges between indicators of different composites. A model in which these matrices are saturated within every composite (or diagonal) can be transformed between the types with \code{\link{transmod}}: for example, \code{transmod(x, composite = "ggm")} turns a \code{cvm} model into the equivalent \code{cnm} model, and \code{transmod(x, composite = "cov")} does the reverse. The types \code{residual} and \code{latent} can be transformed in the same way (\code{latent} only between \code{"cov"} and \code{"ggm"}).

\strong{Computed quantities.} The variance of a composite follows from its weights and the covariances of its indicators. As a result, the (residual) variances of composites in \code{sigma_zeta} (or the scaling of composites in \code{delta_zeta} when \code{latent = "ggm"}) and the intercepts of composites in \code{nu_eta} are not free parameters: they are computed from the other parameters. These elements are shown in the parameter table and in \code{\link{getmatrix}} without standard errors. The elements of \code{nu} are intercepts for the indicators of latent variables and means for the composite indicators. In addition, \code{\link{getmatrix}} returns the following derived matrices: \code{lambda} also contains the composite loadings (the regression coefficients of the composite indicators on their composite) in the columns of composites, \code{sigma_eta} and \code{mu_eta} contain the model-implied variance--covariance matrix and means of the constructs, and \code{sigma_epsilon_c} contains the implied (singular) residual variance--covariance matrix of the composite indicators.

\strong{Identification and scaling.} A composite is scaled by fixing the weight of one of its indicators to 1, by default the first. The weights of the other indicators are then relative to this indicator, which makes the first indicator a poor choice when its own weight is (near) zero: the other weights then become extremely large and are not well defined. Therefore, as in \pkg{lavaan}, another indicator is used automatically when the starting value of the weight of the first indicator is near zero relative to the other weights (controlled with \code{bad_marker_crit}), in which case a message is shown. To scale a composite with a specific indicator, free the fixed weight with \code{\link{freepar}} and fix the weight of the other indicator to 1 with \code{\link{fixpar}}; the same can be used to fix all weights of a composite (for example to form a sum score). The scaling does not change the fit of the model, but it does change the weights, the effects on and of the composite, and their standard errors. A composite with freely estimated weights must be related to at least one other construct through \code{beta}, \code{sigma_zeta} or \code{omega_zeta}: the weights of an isolated composite are not identified (Dijkstra, 2017), and a warning is given when such a model is evaluated. Because the weights are estimated from the relations of a composite with the other constructs, they tend to have larger standard errors than factor loadings, in particular when these relations are weak or the sample is small. As for any structural equation model, models that only differ in the direction of effects between constructs can be equivalent.

\strong{Latent network models.} With \code{latent = "ggm"}, the (residual) variance--covariance matrix of the constructs is modeled as a Gaussian graphical model. Because the scaling of a composite is computed, a composite can then either be part of the latent network (\code{omega_zeta}) or have incoming effects in \code{beta}, but not both. A composite without incoming effects can thus be connected in \code{omega_zeta}, while a composite with incoming effects only has a (computed) residual variance. This is handled in the identification of the model: for a composite without incoming effects the effects on that composite, and for a composite with incoming effects its edges in the latent network, are marked such that modification indices and model search functions skip them. A model that violates the restriction cannot be evaluated and gives an error. If effects on composites are added or removed by hand (\code{\link{freepar}}, \code{\link{fixpar}}), use \code{\link{identify}} to update which parameters are marked.

\strong{Estimation.} Estimation is implemented in both R and C++ (the default), which can be compared using \code{\link{usecpp}}. Starting values for the weights are obtained with partial least squares (mode B), which is consistent for the composite model (Dijkstra, 2017). The intercepts of the indicators of latent variables are started such that the implied means equal the sample means (a latent variable that is predicted by composites has a nonzero mean). The optimizer and its settings are the same as for all other model families. Composite weights are only weakly identified, which makes the fit function flat in the direction of the weights: the optimizer can stop slightly before the optimum, such that the chi-square is typically accurate to about two decimals and estimates can differ in the second or third decimal from software with stricter convergence criteria. In 1500 simulated data sets with n = 250, the chi-square was more than 0.01 too high in 18 (at most 0.24; up to 0.9 in models with a saturated network, and at most 0.03 with n = 1000) and estimates were within 0.1 standard error of the optimum in 99\% of the data sets; in one weakly identified data set the optimizer stopped on a flat ridge, with an estimate 2.3 standard errors from the optimum and a very large reported standard error. Very large standard errors of weights or effects are therefore a reason to estimate the model again with stricter criteria. These can be set with \code{\link{setoptimizer}}, for example \code{setoptimizer(x, "nlminb", optim.args = list(control = list(rel.tol = 1e-10, eval.max = 20000, iter.max = 10000)))}. After estimation a warning is given when the optimizer does not report convergence, when the indicator that scales a composite turns out to have a negligible weight, when a composite with estimated weights is not related to any other construct (its weights are then not identified), or when the computed residual variance of a composite is not positive (its predictors then explain more than its variance). Standard errors are based on the expected Fisher information. Ordinal data and correlation input are currently not supported for this family.

\strong{Penalized estimation.} With \code{estimator = "PML"} (or \code{"PFIML"} with missing data) the network structures in the model are estimated with a LASSO or elastic net penalty, for example the networks of indicators within composites in \code{cnm()}, after which \code{\link{refit}} re-estimates the selected structure without the penalty to obtain standard errors and fit measures. The penalized estimates are obtained with a first-order (proximal gradient) optimizer, which does not use the expected Hessian: for models from this family it can need very many iterations, and it stops less close to the optimum in the weights than the unpenalized estimators do. Use \code{\link{refit}} for the final estimates.

\strong{Relation to lavaan.} Models with \code{latent = "cov"}, \code{residual = "cov"} and \code{composite = "cov"} can be estimated with \pkg{lavaan} (Rosseel, 2012; version 0.7 or later) using the \code{<~} operator, and can be converted to and from \pkg{lavaan} with \code{\link{tolavaan}} and \code{\link{fromlavaan}}. By default \pkg{lavaan} fixes the variances and covariances of the composite indicators to their sample values, whereas they are estimated here (the option \code{composites.cov = "free"} in \pkg{lavaan}): this gives the same estimates, standard errors of weights, loadings and regressions, chi-square and degrees of freedom with complete data, but a different number of parameters. Network structures (\code{"ggm"} and \code{"prec"}) have no \pkg{lavaan} counterpart.
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\references{
Bollen, K. A., & Bauldry, S. (2011). Three Cs in measurement models: Causal indicators, composite indicators, and covariates. \emph{Psychological Methods, 16}(3), 265--284. \doi{10.1037/a0024448}

Dijkstra, T. K. (2017). A perfect match between a model and a mode. In H. Latan & R. Noonan (Eds.), \emph{Partial least squares path modeling} (pp. 55--80). Springer.

Epskamp, S., Rhemtulla, M., & Borsboom, D. (2017). Generalized network psychometrics: Combining network and latent variable models. \emph{Psychometrika, 82}(4), 904--927.

Henseler, J. (2021). \emph{Composite-based structural equation modeling: Analyzing latent and emergent variables}. Guilford Press.

Rosseel, Y. (2012). lavaan: An R package for structural equation modeling. \emph{Journal of Statistical Software, 48}(2), 1--36. \doi{10.18637/jss.v048.i02}

Schamberger, T., Schuberth, F., Henseler, J., & Rosseel, Y. (2025). Factor- and composite-based structural equation modeling: A new approach to incorporate composites in the traditional SEM framework. \emph{arXiv preprint} arXiv:2508.06112. \doi{10.48550/arXiv.2508.06112}

Schuberth, F., Henseler, J., & Dijkstra, T. K. (2018). Confirmatory composite analysis. \emph{Frontiers in Psychology, 9}, 2541. \doi{10.3389/fpsyg.2018.02541}
}
\author{
Sacha Epskamp
}
\seealso{
\code{\link{lvm}}, \code{\link{varcov}}, \code{\link{fromlavaan}}, \code{\link{tolavaan}}, \code{\link{penalize}}
}
\examples{
library("dplyr")

# Load data:
data(StarWars)

# Variables to use: Q1 (being a Star Wars fan), three questions on the original
# trilogy (Q5 - Q7) and three questions on the sequel trilogy (Q8 - Q10):
vars <- c("Q1", paste0("Q", 5:10))

# Weights matrix: two composites, and Q1 as a composite with a single indicator
# (this represents the observed variable itself):
weights <- matrix(0, 7, 3)
weights[1, 1] <- 1
weights[2:4, 2] <- 1
weights[5:7, 3] <- 1

# Structural model: both composites predict being a fan (and are correlated):
beta <- matrix(0, 3, 3)
beta[1, 2:3] <- 1

# Form the model:
mod <- cvm(StarWars, weights = weights, beta = beta, vars = vars,
           latents = c("Fan", "Originals", "Sequels"))

# Run the model:
mod <- mod \%>\% runmodel

# Inspect fit and parameters:
mod \%>\% fit
mod \%>\% parameters

# Estimated weights and the composite loadings:
getmatrix(mod, "weights")
getmatrix(mod, "lambda")

\donttest{
### Composite network model ###
# Model the indicators within every composite as a network:
mod_cnm <- cnm(StarWars, weights = weights, beta = beta, vars = vars,
           latents = c("Fan", "Originals", "Sequels")) \%>\% runmodel

# Remove non-significant edges:
mod_cnm_pruned <- mod_cnm \%>\% prune(alpha = 0.05)

# Compare the models:
compare(full = mod_cnm, pruned = mod_cnm_pruned)

# Network of indicators per composite:
getmatrix(mod_cnm_pruned, "omega_c")

# The cvm model above and the (unpruned) cnm model are equivalent. A fitted
# model can be transformed from one to the other without estimating it again:
mod_cnm2 <- transmod(mod, composite = "ggm", keep_computed = TRUE)
mod_cnm2 \%>\% parameters

### A latent variable and a composite ###
# The sequel questions as a latent variable, predicted by a composite of the
# questions on the original trilogy:
vars2 <- paste0("Q", 5:10)
lambda2 <- matrix(0, 6, 2)
lambda2[4:6, 2] <- 1
weights2 <- matrix(0, 6, 2)
weights2[1:3, 1] <- 1
beta2 <- matrix(0, 2, 2)
beta2[2, 1] <- 1

mod2 <- cvm(StarWars, lambda = lambda2, weights = weights2, beta = beta2,
            vars = vars2, latents = c("Originals", "Sequels")) \%>\% runmodel
mod2 \%>\% parameters

### A network of constructs ###
# The three constructs of the first example in a latent network:
mod3 <- cvm(StarWars, weights = weights, vars = vars, latent = "ggm",
            latents = c("Fan", "Originals", "Sequels")) \%>\% runmodel
getmatrix(mod3, "omega_zeta")
}
}
