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
partialprune(x, alpha = 0.01, matrices, verbose, combinefun = unionmodel, 
             return = c("best","partialprune","union_equal","prune"),
               criterion = "bic",  best = c("lowest","highest"), ...)
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
What model to retur? \code{"best"} for best fitting model (according to BIC, \code{"partialprune"} for the partialpruned model, \code{"union_equal"} for the union model with equality constraints, and \code{"prune"} for the originally pruned model without equality constraints.)
}
\item{best}{
Should the lowest or the highest index of \code{criterion} be used to select the final model?
}
  \item{criterion}{
What criterion to use for the model selection in the last step? Defaults to \code{"bic"} for BIC selection.
}
  \item{\dots}{
Arguments sent to \code{\link{prune}}.
}
}

\references{
Epskamp, S., Isvoranu, A. M., & Cheung, M. (2020). Meta-analytic gaussian network aggregation. PsyArXiv preprint. DOI:10.31234/osf.io/236w8.

Haslbeck, J. (2020). Estimating Group Differences in Network Models using Moderation Analysis. PsyArXiv preprint. DOI:10.31234/osf.io/926pv.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
