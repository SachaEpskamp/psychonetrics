\name{update_baseline_saturated}
\alias{update_baseline}
\alias{update_saturated}
\title{
Set a custom baseline or saturated reference model
}
\description{
Convenience setters that overwrite the \code{$baseline} or \code{$saturated} slot of a \code{psychonetrics} model's \code{@baseline_saturated} list. Together with the \code{baseline_saturated = FALSE} argument exposed by every model constructor, they make it straightforward to attach hand-built reference models before calling \code{\link{runmodel}}.
}
\usage{
update_baseline(x, baseline)
update_saturated(x, saturated)
}
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{baseline}{
A \code{psychonetrics} model to use as the baseline reference.
}
  \item{saturated}{
A \code{psychonetrics} model to use as the saturated reference.
}
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class}) with the relevant slot overwritten.
}
\details{
Pass \code{baseline_saturated = FALSE} to the model constructor to skip auto-construction of the default reference models, attach your own with these setters, and run the target. \code{\link{runmodel}} will run any attached baseline / saturated that is not yet \code{@computed}.
}
\examples{
\dontrun{
mod <- modelfun(..., baseline_saturated = FALSE)
mod <- update_baseline(mod, mod_baseline)
mod <- update_saturated(mod, mod_saturated)
mod <- runmodel(mod)
}
}
\author{
Sacha Epskamp
}
