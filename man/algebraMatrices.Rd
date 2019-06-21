\name{duplicationMatrix}
\alias{duplicationMatrix}
\alias{eliminationMatrix}
\alias{diagonalizationMatrix}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model matrices used in derivatives
}
\description{
These matrices are used in the analytic gradients
}
\usage{
duplicationMatrix(n, diag = TRUE)

eliminationMatrix(n, diag = TRUE)

diagonalizationMatrix(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of rows and columns in the original matrix
}
  \item{diag}{
Logical indicating if the diagonal should be included (set to FALSE for derivative of vech(x))
}
}
\author{
Sacha Epskamp
}
\value{
A sparse matrix
}
\examples{
# Duplication matrix for 10 variables:
duplicationMatrix(10)

# Elimination matrix for 10 variables:
eliminationMatrix(10)

# Diagonailzation matrix for 10 variables:
diagonalizationMatrix(10)
}