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
getmatrix(x, matrix, group)
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
}
\author{
Sacha Epskamp
}