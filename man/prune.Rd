\name{prune}
\alias{prune}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stepdown model search by pruning non-significant parameters.
}
\description{
This function will (recursively) remove parameters that are not significant and refit the model. 
}
\usage{
prune(x, alpha = 0.01,  adjust = c("none", "holm", "hochberg", "hommel", "bonferroni", 
"BH", "BY", "fdr"), matrices, runmodel = TRUE, recursive = TRUE, verbose = TRUE, 
log = TRUE, identify = TRUE, startreduce = 1, limit = Inf, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{alpha}{
Significance level to use.
}
  \item{adjust}{
p-value adjustment method to use. See \code{p.adjust}.
}
  \item{matrices}{
Vector of strings indicating which matrices should be pruned. Will default to network structures.
}
  \item{runmodel}{
Logical, should the model be evaluated after pruning?
}
  \item{recursive}{
Logical, should the pruning process be repeated?
}
  \item{verbose}{
Logical, should messages be printed?
}
  \item{log}{
Logical, should the log be updated?
}
  \item{identify}{
Logical, should models be identified automatically?
}
  \item{startreduce}{
A numeric value indicating a factor with which the starting values should be reduced. Can be useful when encountering numeric problems.
}
  \item{limit}{
The maximum number of parameters to be pruned.
}
  \item{\dots}{
Arguments sent to \link{\code{runmodel}}
}
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\author{
Sacha Epskamp
}