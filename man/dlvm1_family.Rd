\name{dlvm1}
\alias{dlvm1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lag-1 dynamic latent variable family of psychonetrics models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dlvm1(data, vars, lambda, within_latent = c("cov", "chol", "prec", "ggm"), within_residual = c("cov", "chol", "prec", "ggm"), between_latent = c("cov", "chol", "prec", "ggm"), between_residual = c("cov", "chol", "prec", "ggm"), beta = "full", omega_zeta_within = "full", delta_zeta_within = "full", kappa_zeta_within = "full", sigma_zeta_within = "full", lowertri_zeta_within = "full", omega_epsilon_within = "empty", delta_epsilon_within = "empty", kappa_epsilon_within = "empty", sigma_epsilon_within = "empty", lowertri_epsilon_within = "empty", omega_zeta_between = "full", delta_zeta_between = "full", kappa_zeta_between = "full", sigma_zeta_between = "full", lowertri_zeta_between = "full", omega_epsilon_between = "empty", delta_epsilon_between = "empty", kappa_epsilon_between = "empty", sigma_epsilon_between = "empty", lowertri_epsilon_between = "empty", tau, mu_eta, identify = TRUE, identification = c("loadings", "variance"), latents, groups, covs, means, nobs, missing = "listwise", equal = "none", baseline_saturated = TRUE, estimator = "ML", optimizer = "default", storedata = FALSE, sampleStats)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{vars}{
%%     ~~Describe \code{vars} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{within_latent}{
%%     ~~Describe \code{within_latent} here~~
}
  \item{within_residual}{
%%     ~~Describe \code{within_residual} here~~
}
  \item{between_latent}{
%%     ~~Describe \code{between_latent} here~~
}
  \item{between_residual}{
%%     ~~Describe \code{between_residual} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{omega_zeta_within}{
%%     ~~Describe \code{omega_zeta_within} here~~
}
  \item{delta_zeta_within}{
%%     ~~Describe \code{delta_zeta_within} here~~
}
  \item{kappa_zeta_within}{
%%     ~~Describe \code{kappa_zeta_within} here~~
}
  \item{sigma_zeta_within}{
%%     ~~Describe \code{sigma_zeta_within} here~~
}
  \item{lowertri_zeta_within}{
%%     ~~Describe \code{lowertri_zeta_within} here~~
}
  \item{omega_epsilon_within}{
%%     ~~Describe \code{omega_epsilon_within} here~~
}
  \item{delta_epsilon_within}{
%%     ~~Describe \code{delta_epsilon_within} here~~
}
  \item{kappa_epsilon_within}{
%%     ~~Describe \code{kappa_epsilon_within} here~~
}
  \item{sigma_epsilon_within}{
%%     ~~Describe \code{sigma_epsilon_within} here~~
}
  \item{lowertri_epsilon_within}{
%%     ~~Describe \code{lowertri_epsilon_within} here~~
}
  \item{omega_zeta_between}{
%%     ~~Describe \code{omega_zeta_between} here~~
}
  \item{delta_zeta_between}{
%%     ~~Describe \code{delta_zeta_between} here~~
}
  \item{kappa_zeta_between}{
%%     ~~Describe \code{kappa_zeta_between} here~~
}
  \item{sigma_zeta_between}{
%%     ~~Describe \code{sigma_zeta_between} here~~
}
  \item{lowertri_zeta_between}{
%%     ~~Describe \code{lowertri_zeta_between} here~~
}
  \item{omega_epsilon_between}{
%%     ~~Describe \code{omega_epsilon_between} here~~
}
  \item{delta_epsilon_between}{
%%     ~~Describe \code{delta_epsilon_between} here~~
}
  \item{kappa_epsilon_between}{
%%     ~~Describe \code{kappa_epsilon_between} here~~
}
  \item{sigma_epsilon_between}{
%%     ~~Describe \code{sigma_epsilon_between} here~~
}
  \item{lowertri_epsilon_between}{
%%     ~~Describe \code{lowertri_epsilon_between} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
}
  \item{mu_eta}{
%%     ~~Describe \code{mu_eta} here~~
}
  \item{identify}{
%%     ~~Describe \code{identify} here~~
}
  \item{identification}{
%%     ~~Describe \code{identification} here~~
}
  \item{latents}{
%%     ~~Describe \code{latents} here~~
}
  \item{groups}{
%%     ~~Describe \code{groups} here~~
}
  \item{covs}{
%%     ~~Describe \code{covs} here~~
}
  \item{means}{
%%     ~~Describe \code{means} here~~
}
  \item{nobs}{
%%     ~~Describe \code{nobs} here~~
}
  \item{missing}{
%%     ~~Describe \code{missing} here~~
}
  \item{equal}{
%%     ~~Describe \code{equal} here~~
}
  \item{baseline_saturated}{
%%     ~~Describe \code{baseline_saturated} here~~
}
  \item{estimator}{
%%     ~~Describe \code{estimator} here~~
}
  \item{optimizer}{
%%     ~~Describe \code{optimizer} here~~
}
  \item{storedata}{
%%     ~~Describe \code{storedata} here~~
}
  \item{sampleStats}{
%%     ~~Describe \code{sampleStats} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, vars, lambda, within_latent = c("cov", "chol", 
    "prec", "ggm"), within_residual = c("cov", "chol", "prec", 
    "ggm"), between_latent = c("cov", "chol", "prec", "ggm"), 
    between_residual = c("cov", "chol", "prec", "ggm"), beta = "full", 
    omega_zeta_within = "full", delta_zeta_within = "full", kappa_zeta_within = "full", 
    sigma_zeta_within = "full", lowertri_zeta_within = "full", 
    omega_epsilon_within = "empty", delta_epsilon_within = "empty", 
    kappa_epsilon_within = "empty", sigma_epsilon_within = "empty", 
    lowertri_epsilon_within = "empty", omega_zeta_between = "full", 
    delta_zeta_between = "full", kappa_zeta_between = "full", 
    sigma_zeta_between = "full", lowertri_zeta_between = "full", 
    omega_epsilon_between = "empty", delta_epsilon_between = "empty", 
    kappa_epsilon_between = "empty", sigma_epsilon_between = "empty", 
    lowertri_epsilon_between = "empty", tau, mu_eta, identify = TRUE, 
    identification = c("loadings", "variance"), latents, groups, 
    covs, means, nobs, missing = "listwise", equal = "none", 
    baseline_saturated = TRUE, estimator = "ML", optimizer = "default", 
    storedata = FALSE, sampleStats) 
{
    if (missing(lambda)) {
        stop("'lambda' may not be missing")
    }
    if (is.character(lambda)) {
        stop("'lambda' may not be a string")
    }
    within_latent <- match.arg(within_latent)
    between_latent <- match.arg(between_latent)
    within_residual <- match.arg(within_residual)
    between_residual <- match.arg(between_residual)
    identification <- match.arg(identification)
    if (missing(vars)) {
        stop("'vars' argument may not be missing")
    }
    if (!is.matrix(vars)) {
        stop("'vars' must be a design matrix, with rows indicating variables and columns indicating measurements.")
    }
    allVars <- na.omit(as.vector(vars))
    if (missing(sampleStats)) {
        sampleStats <- samplestats(data = data, vars = allVars, 
            groups = groups, covs = covs, means = means, nobs = nobs, 
            missing = ifelse(estimator == "FIML", "pairwise", 
                missing), fimldata = estimator == "FIML", storedata = storedata, 
            weightsmatrix = ifelse(!estimator \%in\% c("WLS", "ULS", 
                "DWLS"), "none", switch(estimator, WLS = "full", 
                ULS = "identity", DWLS = "diag")))
    }
    design <- as(1 * (!is.na(vars)), "sparseMatrix")
    timePerVar <- as.vector(design * row(design))
    timePerVar <- timePerVar[timePerVar != 0]
    nVar <- nrow(vars)
    nTime <- ncol(vars)
    nLat <- ncol(lambda)
    if (is.null(rownames(vars))) {
        rownames(vars) <- paste0("V", seq_len(nrow(vars)))
    }
    varnames <- rownames(vars)
    if (is.null(colnames(vars))) {
        colnames(vars) <- paste0("T", seq_len(ncol(vars)))
    }
    timenames <- colnames(vars)
    if (missing(latents)) {
        latents <- paste0("Eta_", seq_len(nLat))
    }
    if (length(latents) != nLat) {
        stop("Length of 'latents' is not equal to number of latent variables in model.")
    }
    model <- generate_psychonetrics(model = "dlvm1", types = list(within_latent = within_latent, 
        between_latent = between_latent, within_residual = within_residual, 
        between_residual = between_residual), sample = sampleStats, 
        computed = FALSE, equal = equal, optimizer = optimizer, 
        estimator = estimator, distribution = "Gaussian", identification = identification)
    nGroup <- nrow(model@sample@groups)
    nAllVar <- length(allVars)
    model@sample@nobs <- nAllVar * (nAllVar + 1)/2 * nGroup + 
        nAllVar * nGroup
    modMatrices <- list()
    expMeans <- lapply(model@sample@means, function(x) tapply(x, 
        timePerVar, mean, na.rm = TRUE))
    modMatrices$tau <- matrixsetup_mu(tau, nNode = nVar, nGroup = nGroup, 
        labels = varnames, equal = "tau" \%in\% equal, expmeans = expMeans, 
        sampletable = sampleStats, name = "tau")
    modMatrices$mu_eta <- matrixsetup_mu(mu_eta, nNode = nLat, 
        nGroup = nGroup, labels = latents, equal = "mu_eta" \%in\% 
            equal, expmeans = lapply(1:nGroup, function(x) rep(0, 
            nLat)), sampletable = sampleStats, name = "mu_eta")
    firstVars <- apply(vars, 1, function(x) na.omit(x)[1])
    secondVars <- apply(vars, 1, function(x) na.omit(x)[2])
    firstSigma0 <- lapply(sampleStats@covs, function(x) spectralshift(x[firstVars, 
        firstVars]))
    firstSigma1 <- lapply(sampleStats@covs, function(x) spectralshift(x[secondVars, 
        firstVars]))
    modMatrices$lambda <- matrixsetup_lambda(lambda, expcov = firstSigma0, 
        nGroup = nGroup, observednames = sampleStats@variables$label, 
        latentnames = latents, sampletable = sampleStats, name = "lambda")
    prior_bet_cov <- lapply(firstSigma1, function(x) spectralshift(0.5 * 
        (x + t(x))))
    prior_wit_cov <- lapply(seq_along(firstSigma1), function(i) spectralshift(firstSigma0[[i]] - 
        prior_bet_cov[[i]]))
    prior_sig_zeta_within <- lapply(seq_along(firstSigma1), function(i) {
        curLam <- Matrix(as.vector(modMatrices$lambda$start[, 
            , i, drop = FALSE]), nVar, nLat)
        inv <- corpcor::pseudoinverse(as.matrix(kronecker(curLam, 
            curLam)))
        matrix(inv \%*\% as.vector(prior_wit_cov[[i]])/2, nLat, 
            nLat)
    })
    modMatrices <- c(modMatrices, matrixsetup_flexcov(sigma = sigma_zeta_within, 
        lowertri = lowertri_zeta_within, omega = omega_zeta_within, 
        delta = delta_zeta_within, kappa = kappa_zeta_within, 
        type = within_latent, name = "zeta_within", sampleStats = sampleStats, 
        equal = equal, nNode = nLat, expCov = prior_sig_zeta_within, 
        nGroup = nGroup, labels = latents))
    modMatrices$beta <- matrixsetup_beta(beta, name = "beta", 
        nNode = nLat, nGroup = nGroup, labels = latents, equal = "beta" \%in\% 
            equal, sampletable = sampleStats)
    modMatrices <- c(modMatrices, matrixsetup_flexcov(sigma_epsilon_within, 
        lowertri_epsilon_within, omega_epsilon_within, delta_epsilon_within, 
        kappa_epsilon_within, type = within_residual, name = "epsilon_within", 
        sampleStats = sampleStats, equal = equal, nNode = nVar, 
        expCov = lapply(1:nGroup, function(x) diag(0.5, nVar)), 
        nGroup = nGroup, labels = sampleStats@variables$label))
    prior_sig_zeta_between <- lapply(seq_along(firstSigma1), 
        function(i) {
            curLam <- Matrix(as.vector(modMatrices$lambda$start[, 
                , i, drop = FALSE]), nVar, nLat)
            inv <- corpcor::pseudoinverse(as.matrix(kronecker(curLam, 
                curLam)))
            matrix(inv \%*\% as.vector(prior_bet_cov[[i]])/2, nLat, 
                nLat)
        })
    modMatrices <- c(modMatrices, matrixsetup_flexcov(sigma_zeta_between, 
        lowertri_zeta_between, omega_zeta_between, delta_zeta_between, 
        kappa_zeta_between, type = between_latent, name = "zeta_between", 
        sampleStats = sampleStats, equal = equal, nNode = nLat, 
        expCov = prior_sig_zeta_between, nGroup = nGroup, labels = latents))
    modMatrices <- c(modMatrices, matrixsetup_flexcov(sigma_epsilon_between, 
        lowertri_epsilon_between, omega_epsilon_between, delta_epsilon_between, 
        kappa_epsilon_between, type = between_residual, name = "epsilon_between", 
        sampleStats = sampleStats, equal = equal, nNode = nVar, 
        expCov = lapply(1:nGroup, function(x) diag(0.5, nVar)), 
        nGroup = nGroup, labels = sampleStats@variables$label))
    pars <- do.call(generateAllParameterTables, modMatrices)
    model@parameters <- pars$partable
    model@matrices <- pars$mattable
    model@extramatrices <- list(D = psychonetrics::duplicationMatrix(nAllVar), 
        D_y = psychonetrics::duplicationMatrix(nVar), D_eta = psychonetrics::duplicationMatrix(nLat), 
        Dstar_y = psychonetrics::duplicationMatrix(nVar, diag = FALSE), 
        Dstar_eta = psychonetrics::duplicationMatrix(nLat, diag = FALSE), 
        I_y = Diagonal(nVar), I_eta = Diagonal(nLat), A_y = psychonetrics::diagonalizationMatrix(nVar), 
        A_eta = psychonetrics::diagonalizationMatrix(nLat), C_y_y = as(lavaan::lav_matrix_commutation(nVar, 
            nVar), "sparseMatrix"), C_y_eta = as(lavaan::lav_matrix_commutation(nVar, 
            nLat), "sparseMatrix"), C_eta_eta = as(lavaan::lav_matrix_commutation(nLat, 
            nLat), "sparseMatrix"), L_y = psychonetrics::eliminationMatrix(nVar), 
        L_eta = psychonetrics::eliminationMatrix(nLat), design = design)
    muDummy <- matrix(rep(1:nVar, nTime))
    sigDummy <- matrix(0, nVar, nVar)
    sigDummy[lower.tri(sigDummy, diag = TRUE)] <- max(muDummy) + 
        seq_len(nVar * (nVar + 1)/2)
    sigDummy[upper.tri(sigDummy)] <- t(sigDummy)[upper.tri(sigDummy)]
    U <- list(sigDummy)
    for (i in 1:(nTime - 1)) {
        U[[length(U) + 1]] <- matrix(max(unlist(U)) + seq_len(nVar^2), 
            nVar, nVar)
    }
    allSigmas <- blockToeplitz(U)
    totElements <- max(allSigmas)
    subMu <- muDummy[as.vector(design == 1), , drop = FALSE]
    subSigmas <- allSigmas[as.vector(design == 1), as.vector(design == 
        1)]
    inds <- c(as.vector(subMu), subSigmas[lower.tri(subSigmas, 
        diag = TRUE)])
    distVec <- c(as.vector(subMu), subSigmas[lower.tri(subSigmas, 
        diag = TRUE)])
    nTotal <- length(distVec)
    distVecrawts <- seq_along(distVec)[distVec != 0]
    distVec <- distVec[distVec != 0]
    model@extramatrices$P <- sparseMatrix(i = distVecrawts, j = distVec, 
        dims = c(nTotal, totElements))
    model@modelmatrices <- formModelMatrices(model)
    if (is.list(baseline_saturated)) {
        model@baseline_saturated <- baseline_saturated
    }
    else if (isTRUE(baseline_saturated)) {
        model@baseline_saturated$baseline <- varcov(data, type = "chol", 
            lowertri = "empty", vars = allVars, groups = groups, 
            covs = covs, means = means, nobs = nobs, missing = missing, 
            equal = equal, estimator = estimator, baseline_saturated = FALSE, 
            sampleStats = sampleStats)
        model@baseline_saturated$saturated <- varcov(data, type = "chol", 
            lowertri = "full", vars = allVars, groups = groups, 
            covs = covs, means = means, nobs = nobs, missing = missing, 
            equal = equal, estimator = estimator, baseline_saturated = FALSE, 
            sampleStats = sampleStats)
        if (estimator != "FIML") {
            model@baseline_saturated$saturated@computed <- TRUE
            model@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(parVector(model@baseline_saturated$saturated), 
                model@baseline_saturated$saturated)
        }
    }
    if (identify) {
        model <- identify(model)
    }
    return(model)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
