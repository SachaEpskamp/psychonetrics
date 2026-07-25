\name{ri_clpm}
\alias{ri_clpm}
\alias{ri_clpm_stationary}
\alias{ri_clpm_search}
\alias{ri_clpn}
\alias{ri_clpn_stationary}
\alias{ri_clpn_search}
\alias{print.ri_clpm_search}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random intercept cross-lagged panel models and networks
}
\description{
The \code{ri_clpm} function sets up a random intercept cross-lagged panel model (RI-CLPM) as a special case of the latent variable model (\code{\link{lvm}}) framework. The \code{ri_clpm_stationary} function adds stationarity (equality over time) constraints to such a model, and \code{ri_clpm_search} runs a guarded decision-tree search over these constraints. The functions \code{ri_clpn}, \code{ri_clpn_stationary} and \code{ri_clpn_search} are the network-parameterization counterparts (see Details).
}
\usage{
ri_clpm(data, vars,
        datatype = c("auto", "wide", "long"),
        idvar, beepvar,
        standardize = c("none", "z", "quantile", "z_per_wave"),
        lambda,
        type = c("cor", "cov", "chol", "prec", "ggm"),
        verbose = FALSE, \dots)

ri_clpn(data, vars,
        type = c("ggm", "cor", "cov", "chol", "prec"), \dots)

ri_clpm_stationary(x,
        stationary = c("intercepts", "contemporaneous",
                       "innovation", "temporal"))

ri_clpn_stationary(\dots)

ri_clpm_search(x,
        criterion = c("BIC", "AIC", "Chisq"),
        alpha = 0.05,
        include_panelvar = TRUE,
        verbose = TRUE,
        on_rejection = c("release", "stop"), \dots)

ri_clpn_search(\dots)
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
  The parameterization used for the latent (innovation and between-person random intercept) covariance structure: \code{"cor"} (correlations with separate standard deviations; the default of \code{ri_clpm}), \code{"cov"} (covariances), \code{"chol"} (Cholesky), \code{"prec"} (precision), or \code{"ggm"} (Gaussian graphical model / partial correlations; the default of \code{ri_clpn}). Because the latent covariance matrix is block-diagonal over time, this parameterization is applied to both the within-person (innovation) blocks and the between-person (random intercept) block. The default \code{"cor"} parameterization separates the innovation correlations from the innovation standard deviations, so that \code{ri_clpm_search} can test equality of the contemporaneous relations and equality of the innovation variances as two distinct, nested steps.
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
  The criterion used by \code{ri_clpm_search} to decide whether to retain each constraint: \code{"BIC"} (default) or \code{"AIC"} retain the constraint when the criterion does not increase; \code{"Chisq"} retains it when the likelihood-ratio test is not significant at level \code{alpha}. The default is \code{"BIC"} (unchanged since the introduction of \code{ri_clpm_search}).
  }

  \item{alpha}{
  Significance level used when \code{criterion = "Chisq"}.
  }

  \item{include_panelvar}{
  Logical. If \code{TRUE} (default) and the covariance-stationarity constraints are retained, \code{ri_clpm_search} additionally tests the panel(G)VAR model, which treats the first wave as the stationary distribution of the process (wave 1 endogenous). See the \code{"wave1_endogenous"} step in the Details section.
  }

  \item{on_rejection}{
  What \code{ri_clpm_search} does when a constraint is rejected. \code{"release"} (default) runs the guarded decision tree: only a rejected \emph{temporal} constraint stops the search, and a rejection of any other constraint releases that constraint (it is not imposed) but the search continues with the remaining steps. \code{"stop"} reproduces the pre-0.17.2 behavior of stopping at the first rejection, and is provided for reproducibility of analyses run with earlier versions.
  }

  \item{\dots}{
  Arguments sent to \code{\link{lvm}} (for \code{ri_clpm}), to \code{ri_clpm} (for \code{ri_clpn}), to \code{ri_clpm_stationary} (for \code{ri_clpn_stationary}), or to \code{\link{runmodel}} (for \code{ri_clpm_search}) / \code{ri_clpm_search} (for \code{ri_clpn_search}).
  }

}
\details{
The RI-CLPM decomposes each observed score into a stable between-person random intercept and a within-person deviation (innovation). The within-person deviations follow a lag-1 (cross-lagged) structure over time. In this implementation the innovations and the random intercepts are modeled as latent variables under the \code{\link{lvm}} framework, with the latent covariance matrix being block-diagonal over time, so that the within-person (contemporaneous) and between-person covariance structures can each be modeled as a covariance matrix, a Cholesky decomposition, a precision matrix, or a Gaussian graphical model (GGM).

\bold{The RI-CLPN.} When the innovation and between-person structures are parameterized as Gaussian graphical models (\code{type = "ggm"}), the RI-CLPM becomes the \emph{random-intercept cross-lagged panel network model} (RI-CLPN): the contemporaneous and the between-person relations are then partial-correlation networks rather than covariance matrices, while the temporal (cross-lagged) structure is the directed lag-1 network in either case. The function \code{ri_clpn} is a thin wrapper around \code{ri_clpm} that uses \code{type = "ggm"} as its default (the \code{type} argument can still be overridden), mirroring how \code{\link{panelgvar}} wraps \code{\link{panelvar}}. Because \code{ri_clpm_stationary} and \code{ri_clpm_search} operate on an already specified model object, no parameterization is involved there, and \code{ri_clpn_stationary} and \code{ri_clpn_search} are pure aliases of those functions, provided so that a network-parameterized analysis can be written entirely in RI-CLPN notation.

\bold{The stationarity search.} \code{ri_clpm_search} formalizes the typical analysis workflow as a \emph{guarded decision tree} over the four stationarity constraints, closed off by the panel(G)VAR model. The steps are visited in the order \code{"temporal"}, \code{"contemporaneous"}, \code{"intercepts"}, \code{"innovation"}, \code{"wave1_endogenous"}, and each step is tested against the current working model (the base model plus all constraints retained so far). A constraint is retained if it does not worsen fit according to \code{criterion}: with \code{criterion = "Chisq"} when the likelihood-ratio test has \eqn{p \ge} \code{alpha}, and with \code{"AIC"} or \code{"BIC"} when that criterion does not increase. A candidate model that did not converge, or that produced a negative or non-finite chi-square or a non-finite criterion, is never retained; such a fit counts as a rejection of its constraint and the reason is reported in the verdict table.

The tree treats the steps differently according to the targets of inference:

\itemize{
\item \code{"temporal"} is the only hard \bold{gate}. The temporal (cross-lagged) effects are the primary target of inference in a cross-lagged panel analysis. If they are not stationary there is no single lag-1 temporal network to report, and the target of inference changes to a \emph{set} of per-transition temporal networks (one for each wave-to-wave step). The search therefore stops entirely, the unconstrained base model is selected, and the remaining steps are recorded as \code{"not tested"}, since each of them would only be interpretable conditional on a temporal structure that does not hold. When this happens, check whether the waves are equally spaced in time: unequal measurement intervals produce wave-varying lag-1 effects even when the underlying process is perfectly stationary, and so mimic temporal non-stationarity.

\item \code{"contemporaneous"}, \code{"intercepts"} and \code{"innovation"} are nuisance- or secondary-structure questions, so a rejection is \bold{released} rather than fatal: the constraint is not imposed (that set of parameters stays wave-specific) but the search continues to the next step. Each of these constraints is therefore tested exactly once, conditional only on the constraints that survived their own test. This reaches endpoints that a strict chain cannot reach, for example wave-varying contemporaneous relations combined with stationary intercepts and stationary innovation variances.

\item \code{"wave1_endogenous"} tests the panel(G)VAR model, which additionally treats the first wave as a draw from the stationary distribution of the within-person process. It is gated on \emph{covariance} stationarity and therefore requires the temporal, contemporaneous \emph{and} innovation constraints to all have been retained: a within-person process whose innovation covariance structure varies over waves has no stationary distribution, so there is nothing for an endogenous first wave to be drawn from and the step is undefined. In that case the step is reported as \code{"not applicable"} and the search finishes. When the step does apply, the working model is compared to \code{panelvar(vars = ..., sampleStats = ..., within_latent = type, between_latent = type, baseline = "none")}; if the intercept constraint was released, that panel model is fit with \code{stationary_means = FALSE}, so that both models carry a saturated (wave-specific) mean structure and the comparison isolates the wave-1 block. The panel model is the \emph{constrained} special case here: it has \eqn{m(m+1)/2} fewer free parameters, and hence \eqn{m(m+1)/2} \emph{more} degrees of freedom, than the working RI-CLPM, namely the freed exogenous wave-1 within-person covariance block (\eqn{m} = number of variables). The likelihood-ratio test is thus the chi-square of the panel model minus the chi-square of the working model on \eqn{m(m+1)/2} degrees of freedom. This degrees-of-freedom gap is verified at run time and a warning is issued when it differs, which indicates that the constraint bookkeeping of the working model was disturbed upstream (for example by \code{\link{prune}} or \code{\link{fixpar}}).
}

Setting \code{on_rejection = "stop"} reproduces the pre-0.17.2 behavior exactly: the same tree is walked but the search breaks at the first rejection, so no constraint is ever released and the panel(G)VAR step is only reached when all four constraints are retained.

With the default \code{type = "cor"} parameterization the contemporaneous (correlation) step and the innovation-variance (standard deviation) step are cleanly separated and nested; the innovation variances are never equated with the exogenous wave-1 variances nor with the between-person random-intercept variances.

\bold{Extracting and plotting the networks.} Because the model is set up in the \code{\link{lvm}} framework, its matrices span all waves at once and the temporal, contemporaneous and between-person networks are blocks of those matrices. Those blocks can be addressed directly, using the same matrix names as \code{\link{panelvar}} / \code{\link{panelgvar}} together with a \code{time} argument:

\preformatted{getmatrix(x, "beta", time = 2)               # temporal network, wave 2 -> 3
getmatrix(x, "PDC",  time = 2)               # ... as partial directed correlations
getmatrix(x, "omega_zeta_within", time = 3)  # contemporaneous network at wave 3
getmatrix(x, "omega_zeta_between")           # between-person network
getmatrix(x, "nu_between")                   # random-intercept means
CIplot(x, "omega_zeta_within", time = 3)     # confidence intervals of that block}

Blocks carry the variable names of the design matrix, so they can be passed straight to \code{qgraph}. Omitting \code{time} returns one block per wave, with waves that were equated by \code{ri_clpm_stationary} collapsed into a single entry. See \code{\link{getmatrix}} and \code{\link{CIplot}} for the details.

Structural missing waves (\code{NA} entries in the \code{vars} design matrix, i.e. a variable not measured at a particular wave) are not yet supported. Incomplete cases (missing data for some subjects) are handled automatically through full-information maximum likelihood.
}
\value{
\code{ri_clpm}, \code{ri_clpn} and \code{ri_clpm_stationary} return a \code{psychonetrics} model object. \code{ri_clpm_search} returns a list of class \code{"ri_clpm_search"} with elements \code{selected} (the chosen \code{psychonetrics} model), \code{selected_name} (a compact description of the selected model, for example \code{"temporal+contemporaneous+innovation (means free)"} or \code{"panelVAR (means free)"}), \code{models} (all fitted models), \code{comparison} (a \code{\link{compare}} table; note that its rows are sorted by degrees of freedom and that consecutive rows need not be nested when constraints were released), \code{path} (a data frame describing each fitted model and decision), \code{verdicts}, \code{criterion}, \code{alpha} and \code{on_rejection}.

The \code{verdicts} element is a data frame with one row per step and the columns \code{step} (\code{"temporal"}, \code{"contemporaneous"}, \code{"intercepts"}, \code{"innovation"}, \code{"wave1_endogenous"}), \code{tested} (logical), \code{decision} (\code{"retained"}, \code{"released"}, \code{"rejected (gate)"}, \code{"not tested"} or \code{"not applicable"}), \code{Chisq_diff}, \code{DF_diff}, \code{p_value}, \code{dAIC}, \code{dBIC} (each relative to the working model that the step was tested against) and \code{note} (a short explanation, e.g. why a step was not applicable or why a fit was discarded).
}
\author{
Sacha Epskamp
}
\seealso{
\code{\link{lvm}}, \code{\link{dlvm1}}, \code{\link{panelvar}}, \code{\link{panelgvar}}, \code{\link{getmatrix}}, \code{\link{CIplot}}
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

# Form the unconstrained RI-CLPN (= RI-CLPM with GGM innovation/between
# structures). Equivalent to ri_clpm(data, design, type = "ggm"):
mod <- ri_clpn(data, design)

\donttest{
# Run the unconstrained model:
mod <- runmodel(mod)

# Add stationarity constraints step by step:
mod_temp <- runmodel(ri_clpm_stationary(mod,      "temporal"))
mod_con  <- runmodel(ri_clpm_stationary(mod_temp, "contemporaneous"))
compare(base = mod, temporal = mod_temp, contemporaneous = mod_con)

# Or run the whole guarded decision-tree search automatically:
res <- ri_clpn_search(mod, criterion = "BIC")
res
res$verdicts
res$selected

# Extract the networks as nVar x nVar blocks (rather than the full
# latent-sized matrices of the underlying lvm model):
getmatrix(mod, "omega_zeta_between")           # between-person network
getmatrix(mod, "beta", time = 1)               # temporal network, T1 -> T2
getmatrix(mod, "PDC", time = 1)                # ... as partial directed correlations
getmatrix(mod, "omega_zeta_within", time = 2)  # contemporaneous network at wave 2

# The pre-0.17.2 stop-at-first-rejection chain:
ri_clpn_search(mod, criterion = "BIC", on_rejection = "stop")
}
}
