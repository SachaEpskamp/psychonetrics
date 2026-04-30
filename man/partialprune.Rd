\name{partialprune}
\alias{partialprune}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Partial pruning of multi-group models
}
\description{
This function will search for a multi-group model with equality constrains on some but not all parameters. This is called partial pruning (Epskamp, Isvoranu, & Cheung, 2020; Haslbeck, 2020).  The algorithm is as follows: 1. remove all parameters not significant at alpha in all groups (without equality constrains), 2. create a union model with all parameters included in any group included in all groups and constrained equal. 3. Stepwise free equality constrains of the parameter that features the largest sum of modification indices until BIC can no longer be improved. 4. Select and return (by default) the best model according to BIC (original model, pruned model, union model and partially pruned model). 
}
\usage{
partialprune(x, alpha = 0.01, matrices, verbose,
             combinefun = c("unionmodel", "intersectionmodel", "identity"),
             return = c("partialprune","best","union_equal","prune"),
             criterion = "bic",  best = c("lowest","highest"),
             final_prune = c("saturated","partialprune"),
             useMIs = c("joint","simple"),
             release_model = c("pruned","saturated"),
             identify = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{alpha}{
Significance level to use.
}
  \item{matrices}{
Vector of strings indicating which matrices should be pruned. Will default to network structures.
}
  \item{verbose}{
Logical, should messages be printed?
}
  \item{combinefun}{
Function used to combine models of different groups.
}
  \item{return}{
What model to return? \code{"best"} for best fitting model (according to BIC), \code{"partialprune"} for the partialpruned model, \code{"union_equal"} for the union model with equality constraints, and \code{"prune"} for the originally pruned model without equality constraints.
}
\item{best}{
Should the lowest or the highest index of \code{criterion} be used to select the final model?
}
  \item{criterion}{
What criterion to use for the model selection in the last step? Defaults to \code{"bic"} for BIC selection.
}
\item{final_prune}{
Should the last prune step be based on removing edges not significant in the last model in the partialprune algorithm or the first model (saturated model) in the algorithm? Defaults to \code{"saturated"}. Set to \code{"partialprune"} to mimic psychonetrics < 0.13.1 behavior.
}
\item{useMIs}{
How to rank candidate equality constraints in the stepwise release loop. The default \code{"joint"} (introduced in 0.15.4) uses the joint multivariate score / Lagrange-multiplier test (one statistic per equality-constrained parameter, df = G - 1; see \code{\link{equalityScoreTest}}). The legacy \code{"simple"} option sums the per-group univariate modification indices instead, which reproduces the behavior of psychonetrics <= 0.15.3. For two-group models the two options are equivalent; they differ for G >= 3.
}
\item{release_model}{
Specifies the alternative model used in the BIC comparison when deciding whether to release an equality constraint in the stepwise loop. With the default \code{"pruned"} (new in 0.15.5), the just-released parameter is immediately fixed to zero in any group where the initial per-group prune (Step 1) had already removed it, BEFORE refitting and comparing BIC. The equality-constrained model is therefore tested against the asymmetric structure suggested by the initial per-group prune, rather than against a fully saturated per-group alternative. This improves specificity at low sample sizes by penalizing the alternative less heavily for spurious cross-group differences. Set to \code{"saturated"} to reproduce the psychonetrics 0.15.4 behavior in which the released parameter is freed in all groups with no per-group zero pattern applied until the final prune step.
}
\item{identify}{
Logical, passed to the internal \code{\link{prune}}, \code{\link{unionmodel}}, \code{\link{intersectionmodel}}, \code{\link{groupequal}}, and \code{\link{groupfree}} calls. When \code{FALSE}, \code{\link{identify}} is not invoked on the intermediate models, so user-fixed parameters (e.g. \code{fixpar("beta", 1, 1, value = 1)} in a multi-group Ising model) are not re-freed when equality constraints are introduced on \code{matrices}. Defaults to \code{TRUE}, which preserves the historical behavior.
}
  \item{\dots}{
Arguments sent to \code{\link{prune}}.
}
}

\references{
Epskamp, S., Isvoranu, A. M., & Cheung, M. (2020). Meta-analytic gaussian network aggregation. PsyArXiv preprint. DOI:10.31234/osf.io/236w8.

Haslbeck, J. (2020). Estimating Group Differences in Network Models using Moderation Analysis. PsyArXiv preprint. DOI:10.31234/osf.io/926pv.
}

\examples{
\donttest{
# Load bfi data from psych package:
library("psychTools")
data(bfi)

# Also load dplyr for the pipe operator:
library("dplyr")

# Let's take the extroversion items, and gender:
ExData <- bfi \%>\%
  select(E1:E5, gender) \%>\% 
  na.omit # Let's remove missingness (otherwise use Estimator = "FIML)

# Define variables:
vars <- names(ExData)[1:5]

# Saturated estimation:
mod_saturated <- ggm(ExData, 
                     vars = vars,
                     groups = "gender")

# Partial prune model:
mod_partial <- mod_saturated %>% 
  partialprune

# Obtain the networks:
getmatrix(mod_partial, "omega")

# Differences:
getmatrix(mod_partial, "omega")[[1]] -
  getmatrix(mod_partial, "omega")[[2]]

# Difference detected in edge 4 - 5
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
