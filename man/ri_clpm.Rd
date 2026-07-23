\name{ri_clpm}
\alias{ri_clpm}
\alias{ri_clpm_stationary}
\alias{ri_clpm_search}
\alias{print.ri_clpm_search}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random intercept cross-lagged panel models
}
\description{
The \code{ri_clpm} function sets up a random intercept cross-lagged panel model (RI-CLPM) as a special case of the latent variable model (\code{\link{lvm}}) framework. The \code{ri_clpm_stationary} function adds stationarity (equality over time) constraints to such a model, and \code{ri_clpm_search} runs a sequential search that adds these constraints one at a time, retaining each only if it does not worsen fit.
}
\usage{
ri_clpm(data, vars,
        datatype = c("auto", "wide", "long"),
        idvar, beepvar,
        standardize = c("none", "z", "quantile", "z_per_wave"),
        lambda,
        type = c("cor", "cov", "chol", "prec", "ggm"),
        verbose = FALSE, \dots)

ri_clpm_stationary(x,
        stationary = c("intercepts", "contemporaneous",
                       "innovation", "temporal"))

ri_clpm_search(x,
        criterion = c("BIC", "AIC", "Chisq"),
        alpha = 0.05,
        include_panelvar = TRUE,
        verbose = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  A data frame encoding the data used in the analysis. Can be missing if \code{covs} and \code{nobs} are supplied (via \code{\dots}).
}

  \item{vars}{
  Required argument. For wide-format data (the default), a \emph{matrix} with each row indicating a variable and each column indicating a measurement (wave); each cell must contain the name of the variable in \code{data} corresponding to that variable at that wave. The row names of this matrix are used as variable names. For long-format data, a character vector of variable names (see \code{datatype}, \code{idvar} and \code{beepvar}).
  }

  \item{datatype}{
  The data format. \code{"auto"} (default) detects \code{"wide"} when \code{vars} is a matrix and \code{"long"} when \code{vars} is a character vector.
  }

  \item{idvar}{
  Only used for long-format data: the name of the variable encoding the subject/cluster ID.
  }

  \item{beepvar}{
  Optional, only used for long-format data: the name of the variable encoding the measurement occasion (wave). If omitted, observations are assumed to be ordered within each subject.
  }

  \item{standardize}{
  Standardization to apply to the data before estimation. \code{"none"} (default) applies none, \code{"z"} standardizes each variable to mean 0 and SD 1 pooling across all waves, \code{"z_per_wave"} standardizes each variable separately at each wave, and \code{"quantile"} applies a (pooled) nonparametric quantile transformation per variable.
  }

  \item{lambda}{
  A model matrix encoding the factor loading structure. Latent variables are not yet supported for the RI-CLPM; supplying \code{lambda} currently raises an error. When omitted, an observed-variable model is created.
  }

  \item{type}{
  The parameterization used for the latent (innovation and between-person random intercept) covariance structure: \code{"cor"} (correlations with separate standard deviations; the default), \code{"cov"} (covariances), \code{"chol"} (Cholesky), \code{"prec"} (precision), or \code{"ggm"} (Gaussian graphical model / partial correlations). Because the latent covariance matrix is block-diagonal over time, this parameterization is applied to both the within-person (innovation) blocks and the between-person (random intercept) block. The default \code{"cor"} parameterization separates the innovation correlations from the innovation standard deviations, so that \code{ri_clpm_search} can test equality of the contemporaneous relations and equality of the innovation variances as two distinct, nested steps.
  }

  \item{verbose}{
  Logical, should progress be printed to the console?
  }

  \item{x}{
  A \code{psychonetrics} RI-CLPM model (as returned by \code{ri_clpm}).
  }

  \item{stationary}{
  The part of the model to constrain to be stationary (equal over time): \code{"intercepts"} (also frees the random-intercept means), \code{"temporal"} (the cross-lagged/autoregressive effects), \code{"contemporaneous"} (the innovation covariances/network, waves 2 and later), or \code{"innovation"} (the innovation variances, waves 2 and later).
  }

  \item{criterion}{
  The criterion used by \code{ri_clpm_search} to decide whether to retain each constraint: \code{"BIC"} (default) or \code{"AIC"} retain the constraint when the criterion does not increase; \code{"Chisq"} retains it when the likelihood-ratio test is not significant at level \code{alpha}.
  }

  \item{alpha}{
  Significance level used when \code{criterion = "Chisq"}.
  }

  \item{include_panelvar}{
  Logical. If \code{TRUE} (default) and all stationarity constraints are retained, \code{ri_clpm_search} additionally tests the panelVAR model, which treats the first wave as the stationary distribution of the process (wave 1 endogenous).
  }

  \item{\dots}{
  Arguments sent to \code{\link{lvm}} (for \code{ri_clpm}) or to \code{\link{runmodel}} (for \code{ri_clpm_search}).
  }

}
\details{
The RI-CLPM decomposes each observed score into a stable between-person random intercept and a within-person deviation (innovation). The within-person deviations follow a lag-1 (cross-lagged) structure over time. In this implementation the innovations and the random intercepts are modeled as latent variables under the \code{\link{lvm}} framework, with the latent covariance matrix being block-diagonal over time, so that the within-person (contemporaneous) and between-person covariance structures can each be modeled as a covariance matrix, a Cholesky decomposition, a precision matrix, or a Gaussian graphical model (GGM).

\code{ri_clpm_search} formalizes the typical analysis workflow as a sequence akin to measurement-invariance testing. Starting from the unconstrained model it adds, in order, stationarity of the temporal effects, the contemporaneous relations (the innovation correlations / GGM), the intercepts and the innovation variances, and finally the panel(G)VAR model (which additionally treats the first wave as the stationary distribution of the process). Each more constrained model is compared to the currently selected one and the constraint is kept only if it does not worsen fit according to \code{criterion}; the search stops at the first rejection. With the default \code{type = "cor"} parameterization the contemporaneous (correlation) step and the innovation-variance (standard deviation) step are cleanly separated and nested; the innovation variances are never equated with the exogenous wave-1 variances nor with the between-person random-intercept variances.

Structural missing waves (\code{NA} entries in the \code{vars} design matrix, i.e. a variable not measured at a particular wave) are not yet supported. Incomplete cases (missing data for some subjects) are handled automatically through full-information maximum likelihood.
}
\value{
\code{ri_clpm} and \code{ri_clpm_stationary} return a \code{psychonetrics} model object. \code{ri_clpm_search} returns a list of class \code{"ri_clpm_search"} with elements \code{selected} (the chosen \code{psychonetrics} model), \code{selected_name}, \code{models} (all fitted models), \code{comparison} (a \code{compare} table), \code{path} (a data frame describing each step and decision), \code{criterion} and \code{alpha}.
}
\author{
Sacha Epskamp
}
\seealso{
\code{\link{lvm}}, \code{\link{dlvm1}}, \code{\link{panelvar}}, \code{\link{panelgvar}}
}
\examples{
# Simulate a small panel data set (3 variables, 4 waves):
set.seed(1)
nPerson <- 300; nVar <- 3; nWave <- 4

# Between-person random intercepts:
RI <- matrix(rnorm(nPerson * nVar), nPerson, nVar)

# Lag-1 temporal structure:
beta <- diag(0.3, nVar); beta[2, 1] <- 0.2; beta[3, 2] <- 0.2

# Within-person process:
within <- array(0, c(nPerson, nVar, nWave))
within[, , 1] <- matrix(rnorm(nPerson * nVar), nPerson, nVar)
for (w in 2:nWave){
  within[, , w] <- within[, , w - 1] \%*\% t(beta) +
                   matrix(rnorm(nPerson * nVar, sd = 0.8), nPerson, nVar)
}

# Observed scores = within-person deviation + random intercept:
data <- do.call(cbind, lapply(1:nWave, function(w) within[, , w] + RI))
colnames(data) <- as.vector(outer(paste0("V", 1:nVar), 1:nWave, paste, sep = "_"))
data <- as.data.frame(data)

# Design matrix (rows = variables, columns = waves):
design <- matrix(colnames(data), nVar, nWave)
rownames(design) <- paste0("V", 1:nVar)

# Form the unconstrained RI-CLPM with GGM innovation/between structure:
mod <- ri_clpm(data, design, type = "ggm")

\donttest{
# Run the unconstrained model:
mod <- runmodel(mod)

# Add stationarity constraints step by step:
mod_temp <- runmodel(ri_clpm_stationary(mod,      "temporal"))
mod_con  <- runmodel(ri_clpm_stationary(mod_temp, "contemporaneous"))
compare(base = mod, temporal = mod_temp, contemporaneous = mod_con)

# Or run the whole sequential search automatically:
res <- ri_clpm_search(mod, criterion = "BIC")
res
res$selected
}
}
