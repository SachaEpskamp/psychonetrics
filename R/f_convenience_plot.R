# --------------------------------------------------------------------------
# plot() method for psychonetrics models.
# --------------------------------------------------------------------------
# A thin wrapper around qgraph, in the spirit of bootnet's plot method for
# estimateNetwork objects: the same overwritten defaults (cut = 0,
# theme = "colorblind", layout = "spring", parallelEdge = TRUE) and everything
# else passed on through '...'. The model matrices are obtained through
# getmatrix(), so 'group', 'threshold' (and its 'alpha'/'adjust'/'mode'
# companions) and the RI-CLPM 'time' block selector all behave exactly as they
# do there. Unlike getmatrix(), 'matrix' may be a vector: each element becomes
# a panel, in the given order.


# Readable panel titles for the matrix names of the various frameworks:
psychonetrics_paneltitle <- function(mat){
  base <- sub("_(within|between)$", "", mat)
  suffix <- if (grepl("_within$", mat)) "within-person" else if (grepl("_between$", mat)) "between-person" else ""

  main <- if (grepl("^(beta|PDC)_epsilon", base)){
    "Residual temporal"
  } else if (grepl("^(beta|PDC)$", base)){
    "Temporal"
  } else if (grepl("_zeta$", base)){
    "Contemporaneous"
  } else if (grepl("_epsilon$", base)){
    "Residual"
  } else if (grepl("^lambda", base)){
    "Factor loadings"
  } else if (grepl("_y$", base)){
    "Variance-covariance"
  } else ""

  if (main == "" && suffix == "") return(mat)
  if (main == "") return(paste0(toupper(substring(suffix, 1, 1)), substring(suffix, 2)))
  if (suffix == "") return(main)
  # "Contemporaneous (within-person)" is redundant; keep the between/within
  # distinction only where it adds information:
  if (main == "Contemporaneous" && suffix == "within-person") return("Contemporaneous")
  if (main == "Contemporaneous" && suffix == "between-person") return("Between-person")
  paste0(main, " (", suffix, ")")
}


# Which kind of graph does a matrix name imply?
psychonetrics_graphkind <- function(mat, x){
  # The Ising / Blume-Capel 'beta' is the (1 x 1) inverse temperature, not a
  # temporal network:
  if (mat == "beta" && x@model %in% c("Ising", "BlumeCapel")) return("scalar")
  if (grepl("^omega", mat))                     return("omega")
  if (grepl("^kappa", mat))                     return("kappa")
  if (grepl("^sigma", mat))                     return("sigma")
  if (grepl("^rho", mat))                       return("rho")
  if (grepl("^(lowertri|exo_cholesky)", mat))   return("cholesky")
  if (grepl("^(beta|PDC)", mat))                return("directed")
  if (grepl("^lambda", mat))                    return("lambda")
  if (grepl("^(delta|SD)", mat))                return("scaling")
  "other"
}


# Node labels of a model matrix, taken from the parameter table:
psychonetrics_matrixlabels <- function(x, mat){
  nr <- x@matrices$nrow[x@matrices$name == mat]
  nc <- x@matrices$ncol[x@matrices$name == mat]
  if (length(nr) != 1 || length(nc) != 1) return(NULL)

  p  <- x@parameters[x@parameters$matrix == mat, , drop = FALSE]
  rl <- rep(NA_character_, nr)
  cl <- rep(NA_character_, nc)
  ok <- !is.na(p$var1_id) & p$var1_id >= 1 & p$var1_id <= nr
  rl[p$var1_id[ok]] <- p$var1[ok]
  ok <- !is.na(p$var2_id) & p$var2_id >= 1 & p$var2_id <= nc
  cl[p$var2_id[ok]] <- p$var2[ok]

  # Symmetric matrices store only one triangle, so the first row and the last
  # column never appear in the parameter table and their labels are missing.
  # Rows and columns refer to the same variables there, so pool them. This
  # must never happen for a matrix whose rows and columns are different
  # variables (lambda), hence the check that the labels agree wherever both
  # are known:
  sym  <- x@matrices$symmetrical[x@matrices$name == mat]
  both <- !is.na(rl) & !is.na(cl)
  poolable <- nr == nc &&
    (isTRUE(sym) || (any(both) && all(rl[both] == cl[both])))
  if (poolable){
    pooled <- ifelse(is.na(rl), cl, rl)
    rl <- cl <- pooled
  }
  rl[is.na(rl)] <- paste0("V", which(is.na(rl)))
  cl[is.na(cl)] <- paste0("V", which(is.na(cl)))
  list(rows = rl, cols = cl)
}


# Default matrices per model family. Returns the sentinel "semPaths" for
# latent variable models that have no natural network representation.
psychonetrics_defaultmatrices <- function(x){
  cc  <- covchooser
  mod <- x@model

  mats <- if (mod == "varcov"){
    cc(x@submodel)
  } else if (mod == "meta_varcov"){
    cc(x@types$y, "y")
  } else if (mod %in% c("Ising", "BlumeCapel")){
    "omega"
  } else if (mod == "var1"){
    c("beta", cc(x@types$zeta, "zeta"))
  } else if (mod == "tsdlvm1"){
    c("beta", cc(x@types$zeta, "zeta"))
  } else if (mod %in% c("panelvar", "ml_var1")){
    c("beta", cc(x@types$within_latent, "zeta_within"),
      cc(x@types$between_latent, "zeta_between"))
  } else if (mod %in% c("dlvm1", "ml_lvm")){
    c("beta", cc(x@types$within_latent, "zeta_within"),
      cc(x@types$between_latent, "zeta_between"))
  } else if (mod == "ml_varcov"){
    c(cc(x@types$within, "within"), cc(x@types$between, "between"))
  } else if (mod == "lvm"){
    if (identical(x@submodel, "RI_CLPM")){
      rictype <- x@extramatrices$ri_clpm_type
      if (is.null(rictype)) rictype <- x@types$latent
      c("beta", cc(rictype, "zeta_within"), cc(rictype, "zeta_between"))
    } else if (identical(x@submodel, "lnm")){
      cc(x@types$latent, "zeta")
    } else if (identical(x@submodel, "rnm")){
      cc(x@types$residual, "epsilon")
    } else if (identical(x@submodel, "lrnm")){
      c(cc(x@types$latent, "zeta"), cc(x@types$residual, "epsilon"))
    } else {
      # A latent variable model without a network parameterization: the whole
      # path diagram is the natural picture.
      return("semPaths")
    }
  } else {
    stop("No default 'matrix' for a '", mod, "' model; please specify one.")
  }

  # PDC temporal parameterization: the temporal matrix is "PDC", not "beta":
  if ((!is.null(x@types$temporal) && x@types$temporal == "PDC") ||
      (!is.null(x@types$temporal_latent) && x@types$temporal_latent == "PDC")){
    mats[mats == "beta"] <- "PDC"
  }
  if (any(mats == "none")){
    stop("No default 'matrix' could be determined for this model; please specify one.")
  }
  mats
}


plot.psychonetrics <- function(
    x,
    matrix,                 # Matrix/matrices to plot (one panel each)
    ...,                    # Passed to qgraph()
    time,                   # Only for ri_clpm()/ri_clpn(): the wave to plot
    group,                  # Group(s) to plot; defaults to all groups
    threshold = FALSE,      # As in getmatrix()
    alpha = 0.01,
    adjust = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
    mode = c("tested", "all"),
    covariances = FALSE,    # Plot unstandardized (co)variances / precisions
    labels,                 # Node labels (overrides the model's variable names)
    layout = "spring",
    cut = 0,
    theme = "colorblind",
    parallelEdge = TRUE,
    titles = TRUE,          # TRUE / FALSE / a character vector of panel titles
    equalLayout = TRUE,     # Constrain comparable panels to a common layout
    verbose = TRUE
){
  if (!is(x, "psychonetrics")){
    stop("Input must be a 'psychonetrics' object.")
  }
  adjust <- match.arg(adjust)
  mode   <- match.arg(mode)

  if (!x@computed){
    warning("Model has not been computed! Plotting start values.")
  }

  # --- Which matrices? ---
  if (missing(matrix)){
    matrix <- psychonetrics_defaultmatrices(x)
  }
  if (!is.character(matrix)){
    stop("'matrix' must be a character vector of matrix names.")
  }

  # --- Latent variable models without a network: draw the path diagram ---
  if (identical(matrix, "semPaths")){
    if (!requireNamespace("semPlot", quietly = TRUE)){
      stop("Plotting a latent variable model requires the 'semPlot' package. Install it, or specify a 'matrix' to plot as a network (e.g. matrix = \"lambda\").")
    }
    if (verbose){
      message("Drawing the path diagram with semPlot::semPaths(); specify 'matrix' to plot a model matrix as a network instead.")
    }
    lav <- tryCatch(tolavaan(x, verbose = FALSE),
                    error = function(e){
                      stop("Could not convert this model to lavaan for semPaths(): ",
                           conditionMessage(e),
                           ". Specify a 'matrix' to plot as a network instead.", call. = FALSE)
                    })
    sdots <- list(...)
    if (is.null(sdots$whatLabels) && is.null(sdots$what)) sdots$whatLabels <- "est"
    return(invisible(do.call(semPlot::semPaths, c(list(lav), sdots))))
  }

  # --- Groups ---
  allgroups <- x@sample@groups$label
  if (missing(group)){
    group <- allgroups
  } else if (is.numeric(group)){
    group <- x@sample@groups$label[match(group, x@sample@groups$id)]
  }
  if (any(is.na(match(group, allgroups)))){
    stop("'group' does not correspond to a group in the model.")
  }

  # --- Collect the panels ---
  panels <- list()
  for (g in group){
    for (mat in matrix){
      kind <- psychonetrics_graphkind(mat, x)

      if (kind == "cholesky"){
        stop("Matrix '", mat, "' is a Cholesky decomposition and has no network interpretation. Plot the (co)variance or partial correlation matrix instead (e.g. \"", sub("^lowertri", "omega", sub("^exo_cholesky", "omega", mat)), "\").")
      }
      if (kind == "scaling"){
        stop("Matrix '", mat, "' is a diagonal scaling matrix and has no network interpretation.")
      }
      if (kind == "scalar"){
        stop("Matrix 'beta' of an Ising / Blume-Capel model is the (scalar) inverse temperature and has no network interpretation.")
      }

      # Obtain the matrix. For an RI-CLPM the plain matrix names refer to the
      # blocks here (as in CIplot): plotting the full latent-sized matrix,
      # which mixes every wave and the random intercepts, is never what is
      # wanted. getmatrix() keeps those names meaning the full matrix for
      # backward compatibility, so the blocks are requested wave by wave.
      spec <- NULL
      if (is_riclpm(x) && missing(time)){
        spec <- riclpm_resolve(x, mat, implicit = TRUE)
        if (!is.null(spec) && spec$kind %in% c("between", "mean_between")) spec <- NULL
      }
      views <- if (!is.null(spec)){
        stats::setNames(
          lapply(spec$views, function(v)
            getmatrix(x, mat, group = g, threshold = threshold, alpha = alpha,
                      adjust = adjust, mode = mode, time = v$waves[1])),
          vapply(spec$views, "[[", character(1), "label"))
      } else {
        got <- if (missing(time)){
          getmatrix(x, mat, group = g, threshold = threshold, alpha = alpha,
                    adjust = adjust, mode = mode)
        } else {
          getmatrix(x, mat, group = g, threshold = threshold, alpha = alpha,
                    adjust = adjust, mode = mode, time = time)
        }
        if (is.list(got)) got else stats::setNames(list(got), "")
      }

      for (v in seq_along(views)){
        W    <- as.matrix(views[[v]])
        vlab <- names(views)[v]

        if (kind != "lambda" && nrow(W) != ncol(W)){
          stop("Matrix '", mat, "' is not square and cannot be plotted as a network.")
        }
        if (kind == "other"){
          stop("Matrix '", mat, "' has no network interpretation (it is not a (partial) covariance, temporal or loading matrix).")
        }
        if (nrow(W) < 2){
          stop("Matrix '", mat, "' has fewer than two rows and cannot be plotted as a network.")
        }

        # Node labels: the block machinery of ri_clpm() already sets dimnames.
        lab <- if (!is.null(rownames(W)) && !is.null(colnames(W))){
          list(rows = rownames(W), cols = colnames(W))
        } else psychonetrics_matrixlabels(x, mat)
        if (is.null(lab)) lab <- list(rows = paste0("V", seq_len(nrow(W))),
                                      cols = paste0("V", seq_len(ncol(W))))

        panels[[length(panels) + 1]] <- list(
          W = unname(W), kind = kind, matrix = mat, group = g, view = vlab,
          labels = lab
        )
      }
    }
  }

  # --- Standardization and edge type per panel ---
  msg_sigma <- FALSE
  msg_kappa <- FALSE
  for (i in seq_along(panels)){
    p <- panels[[i]]
    W <- p$W
    p$directed     <- FALSE
    p$bidirectional <- FALSE
    p$diag         <- FALSE

    if (p$kind == "sigma"){
      if (covariances){
        p$diag <- TRUE
      } else {
        W <- suppressWarnings(stats::cov2cor(W))
        diag(W) <- 0
        msg_sigma <- TRUE
      }
      p$directed <- TRUE; p$bidirectional <- TRUE
    } else if (p$kind == "rho"){
      diag(W) <- 0
      p$directed <- TRUE; p$bidirectional <- TRUE
    } else if (p$kind == "kappa"){
      if (covariances){
        p$diag <- TRUE
      } else {
        W <- computePCC(W)
        msg_kappa <- TRUE
      }
    } else if (p$kind == "omega"){
      diag(W) <- 0
    } else if (p$kind == "directed"){
      # Temporal networks: the diagonal holds the autoregressive effects, which
      # are drawn as self-loops (as in graphicalVAR / mlVAR).
      p$directed <- TRUE
      p$diag     <- TRUE
    }

    p$W <- W
    panels[[i]] <- p
  }
  if (verbose && msg_sigma){
    message(if (covariances)
      "Plotting the covariance matrix (self-loops show the variances)."
      else "Plotting the covariance matrix as a correlation network; use covariances = TRUE for the unstandardized covariances.")
  }
  if (verbose && msg_kappa){
    message(if (covariances)
      "Plotting the unstandardized precision matrix (self-loops show the diagonal)."
      else "Plotting the precision matrix as a partial correlation network; use covariances = TRUE for the unstandardized precisions.")
  }

  # --- Panel titles ---
  npan <- length(panels)
  if (is.character(titles)){
    if (length(titles) != npan){
      stop("'titles' must have one element per panel (", npan, " panels).")
    }
    for (i in seq_along(panels)) panels[[i]]$title <- titles[i]
  } else if (isTRUE(titles)){
    # Name the matrix only when more than one is shown: with a single matrix
    # across groups or waves, the group / wave alone is the informative label.
    onemat <- length(unique(vapply(panels, "[[", character(1), "matrix"))) == 1
    for (i in seq_along(panels)){
      p  <- panels[[i]]
      tt <- if (npan == 1 || onemat) "" else psychonetrics_paneltitle(p$matrix)
      if (nzchar(p$view)) tt <- paste0(tt, if (nzchar(tt)) " " else "", p$view)
      if (length(group) > 1) tt <- paste0(tt, if (nzchar(tt)) ": " else "", p$group)
      if (!nzchar(tt) && npan > 1) tt <- psychonetrics_paneltitle(p$matrix)
      panels[[i]]$title <- tt
    }
  } else {
    for (i in seq_along(panels)) panels[[i]]$title <- ""
  }

  # --- User-supplied labels ---
  if (!missing(labels)){
    for (i in seq_along(panels)){
      lb <- if (is.list(labels)) labels[[min(i, length(labels))]] else labels
      if (panels[[i]]$kind == "lambda"){
        panels[[i]]$labels$rows <- lb
      } else {
        panels[[i]]$labels <- list(rows = lb, cols = lb)
      }
    }
  }

  # --- A common layout for the panels that share a node set ---
  commonlayout <- NULL
  if (equalLayout && npan > 1 && is.character(layout) && length(layout) == 1){
    net <- vapply(panels, function(p) p$kind != "lambda", logical(1))
    if (sum(net) > 1){
      key <- vapply(panels[net], function(p) paste(p$labels$rows, collapse = "\r"), character(1))
      if (length(unique(key)) == 1){
        commonlayout <- tryCatch(
          do.call(qgraph::averageLayout,
                  c(unname(lapply(panels[net], "[[", "W")), list(layout = layout))),
          error = function(e) NULL)
      }
    }
  }

  # --- Plot ---
  if (npan > 1){
    nrowp <- length(group)
    ncolp <- npan / nrowp
    if (nrowp == 1 && npan > 4){
      ncolp <- min(npan, 3)
      nrowp <- ceiling(npan / ncolp)
    }
    # Only 'mfrow' is touched (and restored); saving the full par() would also
    # restore 'new', which warns on a device that has not been drawn on yet:
    oldmfrow <- graphics::par(mfrow = c(nrowp, ncolp))
    on.exit(graphics::par(oldmfrow), add = TRUE)
  }

  # Anything the user passes through '...' overrides the defaults set here:
  dots     <- list(...)
  dotnames <- names(dots)
  if (is.null(dotnames)) dotnames <- rep("", length(dots))
  named    <- dots[nzchar(dotnames)]
  unnamed  <- dots[!nzchar(dotnames)]

  out <- lapply(panels, function(p){
    if (p$kind == "lambda"){
      # Factor loadings: a reflective measurement diagram. qgraph.loadings()
      # numbers the factor nodes and offers no argument to name them for more
      # than one factor, so the labels are patched into the qgraph object
      # before it is drawn.
      args <- list(model = "reflective", labels = p$labels$rows,
                   cut = cut, theme = theme, title = p$title)
      args <- utils::modifyList(args, named)
      donotplot <- isTRUE(args$DoNotPlot)
      args$DoNotPlot <- TRUE
      q <- do.call(qgraph::qgraph.loadings, c(list(p$W), unnamed, args))
      nman <- nrow(p$W); nfac <- ncol(p$W)
      if (length(q$graphAttributes$Nodes$labels) >= nman + nfac){
        q$graphAttributes$Nodes$labels[seq_len(nman)] <- p$labels$rows
        q$graphAttributes$Nodes$labels[nman + seq_len(nfac)] <- p$labels$cols
      }
      if (!donotplot) plot(q)
      q
    } else {
      lay  <- if (!is.null(commonlayout)) commonlayout else layout
      args <- list(labels = p$labels$rows, directed = p$directed,
                   bidirectional = p$bidirectional, diag = p$diag,
                   parallelEdge = parallelEdge, theme = theme, cut = cut,
                   layout = lay, title = p$title)
      args <- utils::modifyList(args, named)
      do.call(qgraph::qgraph, c(list(p$W), unnamed, args))
    }
  })

  names(out) <- vapply(panels, function(p){
    paste0(p$matrix,
           if (nzchar(p$view)) paste0(" (", p$view, ")") else "",
           if (length(group) > 1) paste0(": ", p$group) else "")
  }, character(1))

  if (length(out) == 1) return(invisible(out[[1]]))
  invisible(out)
}
