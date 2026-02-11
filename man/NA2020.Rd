\name{NA2020}
\alias{NA2020}
\docType{data}
\title{
Network Analysis 2020 Self-Report Data
}
\description{
A subset of self-report data collected in 2020 by Adela-Maria Isvoranu as part of a graduate-level course on network analysis at the University of Amsterdam. Students in the course collected the data from a convenience sample. The full dataset contains responses to 87 items covering topics such as sleep, well-being, and social functioning. This subset contains 8 items used in Chapter 6 of the textbook \emph{Network Psychometrics with R} (Epskamp, Isvoranu, & Haslbeck, 2022) to demonstrate Gaussian Graphical Model estimation with psychonetrics. All items are rated on a 1--7 Likert scale. The dataset contains some missing values, making it suitable for demonstrating FIML estimation.
}
\usage{data("NA2020")}
\format{
  A data frame with 501 observations on 8 variables.
  \describe{
    \item{\code{regular_sleep}}{I try to keep a regular sleep pattern (Q10)}
    \item{\code{worried_sleep}}{I am worried about my current sleeping behavior (Q13)}
    \item{\code{sleep_interfere}}{My sleep interferes with my daily functioning (Q14)}
    \item{\code{happy_health}}{I am happy with my physical health (Q68)}
    \item{\code{optimistic_future}}{I feel optimistic about the future (Q70)}
    \item{\code{very_happy}}{I am very happy (Q75)}
    \item{\code{feel_alone}}{I often feel alone (Q77)}
    \item{\code{happy_love_life}}{I am happy with my love life (Q80)}
  }
}
\source{
The full dataset is available at \url{https://osf.io/45n6d/}. Data collected by Adela-Maria Isvoranu as part of a graduate course on network analysis at the University of Amsterdam (2020).
}
\references{
Epskamp, S., Haslbeck, J. M. B., Isvoranu, A. M., & Van Borkulo, C. D. (2022). Pairwise Markov random fields. In A. M. Isvoranu, S. Epskamp, L. J. Waldorp, & D. Borsboom (Eds.), Network psychometrics with R: A guide for behavioral and social scientists (pp. 95–118). Routledge.
}
\examples{
data(NA2020)

\donttest{
# Estimate a GGM using FIML (as in Chapter 6):
library(dplyr)
mod <- ggm(NA2020, estimator = "FIML") \%>\% runmodel

# Prune non-significant edges:
mod_pruned <- mod \%>\% prune(alpha = 0.05)

# Inspect parameters:
mod_pruned \%>\% parameters
}
}
\keyword{datasets}
