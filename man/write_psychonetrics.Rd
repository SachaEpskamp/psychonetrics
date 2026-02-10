\name{write_psychonetrics}
\alias{write_psychonetrics}
\title{
Write comprehensive model output to a text file
}
\description{
Writes a comprehensive plain-text output file containing all model information, similar to Mplus or LISREL output. The file includes model specification, sample information, parameter estimates, fit measures, model matrices, modification indices, and the model logbook. This file can be shared as supplementary material in publications.
}
\usage{
write_psychonetrics(x, file = "psychonetrics_output.txt",
                    matrices = TRUE, MIs = TRUE, logbook = TRUE)
}
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{file}{
Character string specifying the path to the output file. Defaults to \code{"psychonetrics_output.txt"}.
}
  \item{matrices}{
Logical. Should the full estimated model matrices be included in the output? Defaults to \code{TRUE}. Set to \code{FALSE} for more compact output.
}
  \item{MIs}{
Logical. Should modification indices be included? Only included if modification indices have been computed via \code{\link{addMIs}}. Defaults to \code{TRUE}.
}
  \item{logbook}{
Logical. Should the model logbook be included? Defaults to \code{TRUE}.
}
}
\value{
Invisibly returns the file path of the written output.
}
\author{
Sacha Epskamp
}
\examples{
\donttest{
library("dplyr")

# Load data:
data(StarWars)

# Simple CFA:
Lambda <- matrix(1, 4)
mod <- lvm(StarWars, lambda = Lambda, vars = c("Q1","Q5","Q6","Q7"),
           identification = "variance", latents = "Originals")
mod <- mod \%>\% runmodel \%>\% addfit

# Write output:
write_psychonetrics(mod, file = tempfile(fileext = ".txt"))
}
}
