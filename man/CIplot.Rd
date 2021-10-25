\name{CIplot}
\alias{CIplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Analytic Confidence Intervals
}
\description{
Function to plot analytic confidence intervals (CI) of matrix elements estimated in psychonetrics.
}
\usage{
CIplot(x, matrices, alpha_ci = 0.05, 
      alpha_color = c(0.05, 0.01, 0.001, 1e-04), 
      labels, labels2, labelstart, print = TRUE, 
      major_break = 0.2, minor_break = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{matrices}{
Vector of strings indicating the matrices to plot CIs for
}
  \item{alpha_ci}{
The alpha level used for the CIs
}
  \item{alpha_color}{
A vector of alphas used for coloring the CIs
}
  \item{labels}{
The labels for the variables associated with the rows of a matrix.
}
  \item{labels2}{
The labels for the variables associated with the columns of a matrix. Defaults to the value of \code{labels} for square matrices.
}
  \item{labelstart}{
The value to determine if labels are printed to the right or to the left of the CI
}
  \item{print}{
Logical, should the plots also be printed? Only works when one matrix is used in 'matrices'
}
  \item{major_break}{
Numeric indicating the step size between major breaks
}
  \item{minor_break}{
Numeric indicating the step size between minor breaks
}
}
\value{
A single ggplot2 object, or a list of ggplot2 objects for each matrix requested.
}
\author{
Sacha Epskamp
}
\examples{
### Example from ?ggm ###
# Load bfi data from psych package:
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

# Let's fit an empty GGM:
mod0 <- ggm(ConsData, vars = vars)

# Run the model:
mod0 <- mod0 \%>\% runmodel

# Labels:
labels <- c(
  "indifferent to the feelings of others",
  "inquire about others' well-being",
  "comfort others",
  "love children",
  "make people feel at ease")

# Plot the CIs:
CIplot(mod0,  "omega", labels = labels, labelstart = 0.2)

### Example from ?gvar ###
library("dplyr")
library("graphicalVAR")

beta <- matrix(c(
  0,0.5,
  0.5,0
),2,2,byrow=TRUE)
kappa <- diag(2)
simData <- graphicalVARsim(50, beta, kappa)

# Form model:
model <- gvar(simData)

# Evaluate model:
model <- model \%>\% runmodel

# Plot the CIs:
CIplot(model,  "beta")

}
