# Default matrices to penalize for PML estimation
# Matches the default matrix selection in prune() for each model family
defaultPenalizeMatrices <- function(x) {
  if (!is(x, "psychonetrics")) stop("input must be a psychonetrics object")

  matrices <- character(0)

  if (x@model == "varcov") {
    if (x@submodel == "ggm") {
      matrices <- "omega"
    } else if (x@submodel == "prec") {
      matrices <- "kappa"
    } else if (x@submodel == "chol") {
      matrices <- "lowertri"
    } else if (x@submodel == "cov") {
      matrices <- "sigma"
    }

  } else if (x@model == "meta_varcov") {
    if (x@types$y == "ggm") {
      matrices <- c(matrices, "omega_y")
    } else if (x@types$y == "prec") {
      matrices <- c(matrices, "kappa_y")
    }
    if (x@types$randomEffects == "ggm") {
      matrices <- c(matrices, "omega_randomEffects")
    } else if (x@types$randomEffects == "prec") {
      matrices <- c(matrices, "kappa_randomEffects")
    }

  } else if (x@model == "lvm") {
    if (x@types$latent == "ggm") {
      matrices <- c(matrices, "omega_zeta")
    } else if (x@types$latent == "prec") {
      matrices <- c(matrices, "kappa_zeta")
    }
    if (x@types$residual == "ggm") {
      matrices <- c(matrices, "omega_epsilon")
    } else if (x@types$residual == "prec") {
      matrices <- c(matrices, "kappa_epsilon")
    }

  } else if (x@model == "lnm") {
    matrices <- "omega_eta"

  } else if (x@model == "rnm") {
    matrices <- "omega_epsilon"

  } else if (x@model == "gvar") {
    matrices <- c("beta", "omega_zeta")

  } else if (x@model == "var1") {
    matrices <- c("beta")
    if (x@types$zeta == "prec") {
      matrices <- c(matrices, "kappa_zeta")
    } else if (x@types$zeta == "ggm") {
      matrices <- c(matrices, "omega_zeta")
    }

  } else if (x@model == "panelvar1") {
    matrices <- c("beta")
    if (x@types$contemporaneous == "prec") {
      matrices <- c(matrices, "kappa_zeta")
    } else if (x@types$contemporaneous == "ggm") {
      matrices <- c(matrices, "omega_zeta")
    }
    if (x@types$between == "prec") {
      matrices <- c(matrices, "kappa_mu")
    } else if (x@types$between == "ggm") {
      matrices <- c(matrices, "omega_mu")
    }

  } else if (x@model %in% c("ml_lvm", "dlvm1")) {
    if (x@model == "dlvm1") {
      matrices <- c("beta")
    }
    if (x@types$within_latent == "prec") {
      matrices <- c(matrices, "kappa_zeta_within")
    } else if (x@types$within_latent == "ggm") {
      matrices <- c(matrices, "omega_zeta_within")
    }
    if (x@types$within_residual == "prec") {
      matrices <- c(matrices, "kappa_epsilon_within")
    } else if (x@types$within_residual == "ggm") {
      matrices <- c(matrices, "omega_epsilon_within")
    }
    if (x@types$between_latent == "prec") {
      matrices <- c(matrices, "kappa_zeta_between")
    } else if (x@types$between_latent == "ggm") {
      matrices <- c(matrices, "omega_zeta_between")
    }
    if (x@types$between_residual == "prec") {
      matrices <- c(matrices, "kappa_epsilon_between")
    } else if (x@types$between_residual == "ggm") {
      matrices <- c(matrices, "omega_epsilon_between")
    }

  } else if (x@model == "tsdlvm1") {
    matrices <- c("beta")
    if (x@types$zeta == "prec") {
      matrices <- c(matrices, "kappa_zeta")
    } else if (x@types$zeta == "ggm") {
      matrices <- c(matrices, "omega_zeta")
    }
    if (x@types$epsilon == "prec") {
      matrices <- c(matrices, "kappa_epsilon")
    } else if (x@types$epsilon == "ggm") {
      matrices <- c(matrices, "omega_epsilon")
    }

  } else if (x@model == "Ising") {
    matrices <- c("omega")

  } else {
    stop("No default penalization matrices defined for model '", x@model, "'.")
  }

  matrices
}
