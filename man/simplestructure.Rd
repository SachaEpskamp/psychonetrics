\name{simplestructure}
\alias{simplestructure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate factor loadings matrix with simple structure
}
\description{
This function generates the input for \code{lambda} arguments in latent variable models using a simple structure. The input is a vector with an element for each variable indicating the factor the variable loads on.
}
\usage{
simplestructure(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector indicating which factor each indicator loads on.
}
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
