\name{getmatrix}
\alias{getmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract an estimated matrix
}
\description{
This function will extract an estimated matrix, and will either return a single matrix for single group models or a list of such matrices for multiple group models. 
}
\usage{
getmatrix(x, matrix, group, threshold = FALSE, alpha = 0.01,
           adjust = c("none", "holm", "hochberg", "hommel",
           "bonferroni", "BH", "BY", "fdr"), mode = c("tested",
           "all"), diag = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{matrix}{
String indicating the matrix to be extracted.
}
  \item{group}{
Integer indicating the group for the matrix to be extracted.
}
\item{threshold}{Logical. Should the matrix be thresholded (non-significant values set to zero? Can also be a value with an absolute threshold below wich parameters are set to zero.)}
  \item{alpha}{
Significance level to use.
}
  \item{adjust}{
p-value adjustment method to use. See \code{p.adjust}.
}
\item{mode}{
Mode for adjusting for multiple comparisons. Should all parameters be considered as the total number of tests or only the tested parameters (parameters of interest)?
}
\item{diag}{Set to FALSE to set diagonal elements to zero.}
}
\value{
A matrix of parameter estimates, of a list of such matrices for multiple group models.
}
\author{
Sacha Epskamp
}
\examples{# Load bfi data from psych package:
library("psychTools")
data(bfi)

# Also load dplyr for the pipe operator:
library("dplyr")

# Let's take the agreeableness items, and gender:
ConsData <- bfi \%>\% 
  select(A1:A5, gender) \%>\% 
  na.omit # Let's remove missingness (otherwise use Estimator = "FIML)

# Define variables:
vars <- names(ConsData)[1:5]

# Let's fit a full GGM:
mod <- ggm(ConsData, vars = vars, omega = "full")

# Run model:
mod <- mod \%>\% runmodel

# Obtain network:
mod \%>\% getmatrix("omega")
}