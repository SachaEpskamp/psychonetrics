\name{lvm}
\alias{lvm}
\alias{lnm}
\alias{rnm}
\alias{lrnm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Continuous latent variable family of psychonetrics models
}
\description{
This is the family of models that models the data as a structural equation model (SEM), allowing the latent and residual variance-covariance matrices to be further modeled as networks. The \code{latent} and \code{residual} arguments can be used to define what latent and residual models are used respectively: \code{"cov"} (default) models a variance-covariance matrix directly, \code{"chol"} models a Cholesky decomposition, \code{"prec"} models a precision matrix, and \code{"ggm"} models a Gaussian graphical model (Epskamp, Rhemtulla and Borsboom, 2017). The wrapper \code{lnm()} sets \code{latent = "ggm"} for the latent network model (LNM), the wrapper  \code{rnm()} sets \code{residual = "ggm"} for the residual network model (RNM), and the wrapper  \code{lrnm()} combines the LNM and RNM.
}
\usage{
lvm(data, lambda, latent = c("cov", "chol", "prec",
                    "ggm"), residual = c("cov", "chol", "prec", "ggm"),
                    sigma_zeta = "full", kappa_zeta = "full", omega_zeta =
                    "full", lowertri_zeta = "full", delta_zeta = "full",
                    sigma_epsilon = "empty", kappa_epsilon = "empty",
                    omega_epsilon = "empty", lowertri_epsilon = "empty",
                    delta_epsilon = "empty", beta = "empty", nu, nu_eta,
                    identify = TRUE, identification = c("loadings",
                    "variance"), vars, latents, groups, covs, means, nobs,
                    missing = "listwise", equal = "none",
                    baseline_saturated = TRUE, estimator = "ML", optimizer
                    = "default", storedata = FALSE, WLS.V, covtype =
                    c("choose", "ML", "UB"), sampleStats)

lnm(\dots)
rnm(\dots)
lrnm(\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame encoding the data used in the analysis. Can be missing if \code{covs} and \code{nobs} are supplied.
}
  \item{lambda}{
A model matrix encoding the factor loading structure. Each row indicates an indicator and each column a latent. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{latent}{
The type of latent model used. See description.
}
  \item{residual}{
The type of residual model used. See description.
}
  \item{sigma_zeta}{
Only used when \code{latent = "cov"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_zeta}{
Only used when \code{latent = "prec"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega_zeta}{
Only used when \code{latent = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_zeta}{
Only used when \code{latent = "chol"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_zeta}{
Only used when \code{latent = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{sigma_epsilon}{
Only used when \code{residual = "cov"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_epsilon}{
Only used when \code{residual = "prec"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega_epsilon}{
Only used when \code{residual = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_epsilon}{
Only used when \code{residual = "chol"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_epsilon}{
Only used when \code{residual = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{beta}{
A model matrix encoding the structural relations between latent variables. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{nu}{
Optional vector encoding the intercepts of the observed variables. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free intercepts, and higher integers to indicate equality constrains. For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
}
  \item{nu_eta}{
Optional vector encoding the intercepts of the latent variables. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free intercepts, and higher integers to indicate equality constrains. For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
}
  \item{identify}{
Logical, should the model be automatically identified?
}
  \item{identification}{
Type of identification used. \code{"loadings"} to fix the first factor loadings to 1, and \code{"variance"} to fix the diagonal of the latent variable model matrix (sigma_zeta, lowertri_zeta, delta_zeta or kappa_zeta) to 1.
}
  \item{vars}{
An optional character vector encoding the variables used in the analysis. Must equal names of the dataset in \code{data}.
}
  \item{latents}{
An optional character vector with names of the latent variables.
}
  \item{groups}{
An optional string indicating the name of the group variable in \code{data}.
}
  \item{covs}{
A sample variance--covariance matrix, or a list/array of such matrices for multiple groups. Make sure \code{covtype} argument is set correctly to the type of covariances used. 
}
  \item{means}{
A vector of sample means, or a list/matrix containing such vectors for multiple groups.
}
  \item{nobs}{
The number of observations used in \code{covs} and \code{means}, or a vector of such numbers of observations for multiple groups.
}
  \item{missing}{
How should missingness be handled in computing the sample covariances and number of observations when \code{data} is used. Can be \code{"listwise"} for listwise deletion, or \code{"pairwise"} for pairwise deletion.
}
  \item{equal}{
A character vector indicating which matrices should be constrained equal across groups. 
}
  \item{baseline_saturated}{
A logical indicating if the baseline and saturated model should be included. Mostly used internally and NOT Recommended to be used manually.
}
  \item{estimator}{
The estimator to be used. Currently implemented are \code{"ML"} for maximum likelihood estimation, \code{"FIML"} for full-information maximum likelihood estimation, \code{"ULS"} for unweighted least squares estimation, \code{"WLS"} for weighted least squares estimation, and \code{"DWLS"} for diagonally weighted least squares estimation.
}
  \item{optimizer}{
The optimizer to be used. Usually either \code{"nlminb"} (with box constrains) or \code{"ucminf"} (ignoring box constrains), but any optimizer supported by \code{optimr} can be used.
}
  \item{storedata}{
Logical, should the raw data be stored? Needed for bootstrapping (see \code{\link{bootstrap}}).
}
  \item{WLS.V}{
The weights matrix used in WLS estimation (experimental)
}
  \item{sampleStats}{
An optional sample statistics object. Mostly used internally. 
}
\item{covtype}{
  If 'covs' is used, this is the type of covariance (maximum likelihood or unbiased) the input covariance matrix represents. Set to \code{"ML"} for maximum likelihood estimates (denominator n) and \code{"UB"} to unbiased estimates (denominator n-1). The default will try to find the type used, by investigating which is most likely to result from integer valued datasets.
}
\item{\dots}{
Arguments sent to \code{varcov}
}
}
\details{
The model used in this family is:

\eqn{\mathrm{var}( \boldsymbol{y} ) = \boldsymbol{\Lambda} (\boldsymbol{I} - \boldsymbol{B})^{-1} \boldsymbol{\Sigma}_{\zeta}  (\boldsymbol{I} - \boldsymbol{B})^{-1\top}  \boldsymbol{\Lambda}^{\top} +  \boldsymbol{\Sigma}_{\varepsilon} }{var(y) = lambda * (I - beta)^(-1) * sigma_zeta * (I - beta)^(-1') * lambda' + sigma_epsilon}

\eqn{\mathcal{E}( \boldsymbol{y} ) = \boldsymbol{\nu} +  \boldsymbol{\Lambda}  (\boldsymbol{I} - \boldsymbol{B})^{-1} \boldsymbol{\nu}_eta}{E(y) = nu + lambda * (I - beta)^(-1) * nu_eta}

in which the latent covariance matrix can further be modeled in three ways. With \code{latent = "chol"} as Cholesky decomposition:

\eqn{\boldsymbol{\Sigma}_{\zeta} = \boldsymbol{L}_{\zeta}\boldsymbol{L}_{\zeta}}{sigma_zeta = lowertri_zeta * lowertri_zeta},

with \code{latent = "prec"} as Precision matrix:

\eqn{\boldsymbol{\Sigma}_{\zeta} = \boldsymbol{K}_{\zeta}^{-1}}{sigma_zeta = kappa_zeta^(-1)},

and finally with \code{latent = "ggm"} as Gaussian graphical model:

\eqn{\boldsymbol{\Sigma}_{\zeta} = \boldsymbol{\Delta}_{\zeta}(\boldsymbol{I} - \boldsymbol{\Omega}_{\zeta})^(-1) \boldsymbol{\Delta}_{\zeta}}{sigma_zeta = delta_zeta * (I - omega_zeta)^(-1) * delta_zeta}.

Likewise, the residual covariance matrix can also further be modeled in three ways. With \code{residual = "chol"} as Cholesky decomposition:

\eqn{\boldsymbol{\Sigma}_{\varepsilon} = \boldsymbol{L}_{\varepsilon}\boldsymbol{L}_{\varepsilon}}{sigma_epsilon = lowertri_epsilon * lowertri_epsilon},

with \code{latent = "prec"} as Precision matrix:

\eqn{\boldsymbol{\Sigma}_{\varepsilon} = \boldsymbol{K}_{\varepsilon}^{-1}}{sigma_epsilon = kappa_epsilon^(-1)},

and finally with \code{latent = "ggm"} as Gaussian graphical model:

\eqn{\boldsymbol{\Sigma}_{\varepsilon} = \boldsymbol{\Delta}_{\varepsilon}(\boldsymbol{I} - \boldsymbol{\Omega}_{\varepsilon})^(-1) \boldsymbol{\Delta}_{\varepsilon}}{sigma_epsilon = delta_epsilon * (I - omega_epsilon)^(-1) * delta_epsilon}.
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\references{
Epskamp, S., Rhemtulla, M., & Borsboom, D. (2017). Generalized network psychometrics: Combining network and latent variable models. Psychometrika, 82(4), 904-927.
}
\author{
Sacha Epskamp
}

\examples{
library("dplyr")

### Confirmatory Factor Analysis ###

# Example also shown in https://youtu.be/Hdu5z-fwuk8

# Load data:
data(StarWars)

# Originals only:
Lambda <- matrix(1,4)

# Model:
mod0 <- lvm(StarWars, lambda = Lambda, vars = c("Q1","Q5","Q6","Q7"), 
            identification = "variance", latents = "Originals")
            
# Run model:
mod0 <- mod0 \%>\% runmodel

# Evaluate fit:
mod0 \%>\% fit

\donttest{
# Full analysis
# Factor loadings matrix:
Lambda <- matrix(0, 10, 3)
Lambda[1:4,1] <- 1
Lambda[c(1,5:7),2] <- 1
Lambda[c(1,8:10),3] <- 1

# Observed variables:
obsvars <- paste0("Q",1:10)

# Latents:
latents <- c("Prequels","Original","Sequels")

# Make model:
mod1 <- lvm(StarWars, lambda = Lambda, vars = obsvars, 
            identification = "variance", latents = latents)

# Run model:
mod1 <- mod1 \%>\% runmodel

# Look at fit:
mod1

# Look at parameter estimates:
mod1 \%>\% parameters

# Look at modification indices:
mod1 \%>\% MIs

# Add and refit:
mod2 <- mod1 \%>\% freepar("sigma_epsilon","Q10","Q4") \%>\% runmodel

# Compare:
compare(original = mod1, adjusted = mod2)

# Fit measures:
mod2 \%>\% fit

### Path diagrams ###
# semPlot is not (yet) supported by default, but can be used as follows:
# Load packages:
library("semPlot")

# Estimates:
lambdaEst <- getmatrix(mod2, "lambda")
psiEst <- getmatrix(mod2, "sigma_zeta")
thetaEst <- getmatrix(mod2, "sigma_epsilon")

# LISREL Model: LY = Lambda (lambda-y), TE = Theta (theta-epsilon), PS = Psi
mod <- lisrelModel(LY =  lambdaEst, PS = psiEst, TE = thetaEst)

# Plot with semPlot:
semPaths(mod, "std", "est", as.expression = "nodes")


# We can make this nicer (set whatLabels = "none" to hide labels):
semPaths(mod,

# this argument controls what the color of edges represent. In this case, 
# standardized parameters:
    what = "std", 
    
# This argument controls what the edge labels represent. In this case, parameter 
# estimates:
    whatLabels = "est", 
    
# This argument draws the node and edge labels as mathematical exprssions:    
    as.expression = "nodes", 
  
# This will plot residuals as arrows, closer to what we use in class:
    style = "lisrel",
    
# This makes the residuals larger:
    residScale = 10, 
    
# qgraph colorblind friendly theme:
    theme = "colorblind",
    
# tree layout options are "tree", "tree2", and "tree3":
    layout = "tree2", 

# This makes the latent covariances connect at a cardinal center point:
    cardinal = "lat cov",

# Changes curve into rounded straight lines:
    curvePivot = TRUE, 
    
# Size of manifest variables:
    sizeMan = 4, 
    
# Size of latent varibales:
    sizeLat = 10, 
    
# Size of edge labels:
    edge.label.cex = 1,
    
# Sets the margins:
    mar = c(9,1,8,1), 
    
# Prevents re-ordering of ovbserved variables:
    reorder = FALSE, 
    
# Width of the plot:
    width = 8, 
    
# Height of plot:
    height = 5, 

# Colors according to latents:
    groups = "latents",
    
# Pastel colors:
    pastel = TRUE, 
    
# Disable borders:
    borders = FALSE 
    )
    
# Use arguments filetype = "pdf" and filename = "semPlotExample1" to store PDF

### Latent Network Modeling ###

# Latent network model:
lnm <- lvm(StarWars, lambda = Lambda, vars = obsvars,
           latents = latents, identification = "variance",
           latent = "ggm")

# Run model:
lnm <- lnm \%>\% runmodel

# Look at parameters:
lnm \%>\% parameters

# Remove non-sig latent edge:
lnm <- lnm \%>\% prune(alpha = 0.05)

# Compare to the original CFA model:
compare(cfa = mod1, lnm = lnm)

# Plot network:
library("qgraph")
qgraph(lnm@modelmatrices[[1]]$omega_zeta, labels = latents,
       theme = "colorblind", vsize = 10)

# A wrapper for the latent network model is the lnm function:
lnm2 <- lnm(StarWars, lambda = Lambda, vars = obsvars,
            latents = latents, identification = "variance")
lnm2 <- lnm2 \%>\% runmodel \%>\% prune(alpha = 0.05)
compare(lnm, lnm2) # Is the same as the model before.

# I could also estimate a "residual network model", which adds partial correlations to 
# the residual level:
# This can be done using lvm(..., residal = "ggm") or with rnm(...)
rnm <- rnm(StarWars, lambda = Lambda, vars = obsvars,
           latents = latents, identification = "variance")
# Stepup search:
rnm <- rnm \%>\% stepup

# It will estimate the same model (with link Q10 - Q4) as above. In the case of only one 
# partial correlation, There is no difference between residual covariances (SEM) or 
# residual partial correlations (RNM).


# For more information on latent and residual network models, see:
# Epskamp, S., Rhemtulla, M.T., & Borsboom, D. Generalized Network Psychometrics: 
# Combining Network and Latent Variable Models 
# (2017). Psychometrika. doi:10.1007/s11336-017-9557-x

### Gaussian graphical models ###

# All psychonetrics functions (e.g., lvm, lnm, rnm...) allow input via a covariance 
# matrix, with the "covs" and "nobs" arguments.
# The following fits a baseline GGM network with no edges:
S <- (nrow(StarWars) - 1)/ (nrow(StarWars)) * cov(StarWars[,1:10])
ggmmod <- ggm(covs = S, nobs = nrow(StarWars), omega = "empty")

# Run model with stepup search and pruning:
ggmmod <- ggmmod \%>\% stepup(greedy=TRUE) \%>\% prune

# Fit measures:
ggmmod \%>\% fit

# Plot network:
nodeNames <- c(
"I am a huge Star Wars\\nfan! (star what?)",
"I would trust this person\\nwith my democracy.",
"I enjoyed the story of\\nAnakin's early life.",
"The special effects in\\nthis scene are awful (Battle of\\nGeonosis).",
"I would trust this person\\nwith my life.",
"I found Darth Vader's big\\nreveal in 'Empire' one of the greatest
moments in movie history.",
"The special effects in\\nthis scene are amazing (Death Star\\nExplosion).",
"If possible, I would\\ndefinitely buy this\\ndroid.",
"The story in the Star\\nWars sequels is an improvement to\\nthe previous movies.",
"The special effects in\\nthis scene are marvellous (Starkiller\\nBase Firing)."
)
library("qgraph")
qgraph(as.matrix(ggmmod@modelmatrices[[1]]$omega), nodeNames = nodeNames, 
    legend.cex = 0.25,  theme = "colorblind", layout = "spring")

# We can actually compare this model statistically (note they are not nested) to the 
# latent variable model:
compare(original_cfa = mod1, adjusted_cfa = mod2, exploratory_ggm = ggmmod)


### Meausrement invariance ###
# Let's say we are interested in seeing if people >= 30 like the original trilogy better 
# than people < 30.
# First we can make a grouping variable:
StarWars$agegroup <- ifelse(StarWars$Q12 < 30, "young", "less young")

# Let's look at the distribution:
table(StarWars$agegroup) # Pretty even...

# Observed variables:
obsvars <- paste0("Q",1:10)

# Let's look at the mean scores:
StarWars \%>\% group_by(agegroup) \%>\% summarize_each_(funs(mean),vars = obsvars)
# Less young people seem to score higher on prequel questions and lower on other 
# questions

# Factor loadings matrix:
Lambda <- matrix(0, 10, 3)
Lambda[1:4,1] <- 1
Lambda[c(1,5:7),2] <- 1
Lambda[c(1,8:10),3] <- 1

# Residual covariances:
Theta <- diag(1, 10)
Theta[4,10] <- Theta[10,4] <- 1

# Latents:
latents <- c("Prequels","Original","Sequels")

# Make model:
mod_configural <- lvm(StarWars, lambda = Lambda, vars = obsvars,
            latents = latents, sigma_epsilon = Theta,
            identification = "variance",
            groups =  "agegroup")

# Run model:
mod_configural <- mod_configural \%>\% runmodel

# Look at fit:
mod_configural
mod_configural \%>\% fit

# Looks good, let's try weak invariance:
mod_weak <- mod_configural \%>\% groupequal("lambda") \%>\% runmodel

# Compare models:
compare(configural = mod_configural, weak = mod_weak)

# weak invariance can be accepted, let's try strong:
mod_strong <- mod_weak \%>\% groupequal("nu") \%>\% runmodel
# Means are automatically identified

# Compare models:
compare(configural = mod_configural, weak = mod_weak, strong = mod_strong)

# Questionable p-value and AIC difference, but ok BIC difference. This is quite good, but 
# let's take a look. I have not yet implemented LM tests for equality constrains, but we 
# can look at something called "equality-free" MIs:
mod_strong \%>\% MIs(matrices = "nu", type = "free")

# Indicates that Q10 would improve fit. We can also look at residuals:
residuals(mod_strong)

# Let's try freeing intercept 10:
mod_strong_partial <- mod_strong \%>\% groupfree("nu",10) \%>\% runmodel

# Compare all models:
compare(configural = mod_configural,
        weak = mod_weak,
        strong = mod_strong,
        strong_partial = mod_strong_partial)

# This seems worth it and lead to an acceptable model! It seems that older people find 
# the latest special effects more marvellous!
mod_strong_partial \%>\% getmatrix("nu")

# Now let's investigate strict invariance:
mod_strict <- mod_strong_partial \%>\% groupequal("sigma_epsilon") \%>\% runmodel

# Compare all models:
compare(configural = mod_configural,
        weak = mod_weak,
        strong_partial = mod_strong_partial,
        strict = mod_strict)
# Strict invariance can be accepted!

#  Now we can test for homogeneity!
# Are the latent variances equal?
mod_eqvar <- mod_strict \%>\% groupequal("sigma_zeta") \%>\% runmodel

# Compare:
compare(strict = mod_strict, eqvar = mod_eqvar) 

# This is acceptable. What about the means? (alpha = nu_eta)
mod_eqmeans <- mod_eqvar \%>\% groupequal("nu_eta") \%>\% runmodel

# Compare:
compare(eqvar = mod_eqvar, eqmeans = mod_eqmeans)

# Rejected! We could look at MIs again:
mod_eqmeans \%>\% MIs(matrices = "nu_eta", type = "free")

# Indicates the strongest effect for prequels. Let's see what happens:
eqmeans2 <- mod_eqvar \%>\% 
  groupequal("nu_eta",row = c("Original","Sequels")) \%>\% runmodel

# Compare:
compare(eqvar = mod_eqvar, eqmeans = eqmeans2)
# Questionable, what about the sequels as well?

eqmeans3 <- mod_eqvar \%>\% groupequal("nu_eta", row = "Original") \%>\% runmodel

# Compare:
compare(eqvar = mod_eqvar, eqmeans = eqmeans3)

# Still questionable.. Let's look at the mean differences:
mod_eqvar \%>\% getmatrix("nu_eta")

# Looks like people over 30 like the prequels better and the other two trilogies less!
}
}
