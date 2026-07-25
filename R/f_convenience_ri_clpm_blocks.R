# --------------------------------------------------------------------------
# Block extraction for RI-CLPM / RI-CLPN models.
# --------------------------------------------------------------------------
# ri_clpm() is built on the lvm() framework, so all of its model matrices span
# the FULL set of latent variables: the nVar within-person components of each
# of the nTime waves, followed by the nVar random intercepts. The quantities a
# user actually wants are nVar x nVar blocks of those matrices:
#
#   * the temporal network of a wave-to-wave transition   (a block of 'beta')
#   * the contemporaneous network of one wave             (a wave block of the
#                                                          typed '*_zeta' matrix)
#   * the between-person network                          (the random-intercept
#                                                          block of that matrix)
#   * the wave-specific intercepts / random-intercept means
#
# Rather than a new vocabulary, these are addressed with the matrix names that
# the sibling panel frameworks already use -- "beta", "omega_zeta_within",
# "omega_zeta_between", "nu_within", "nu_between" -- plus a 'time' argument
# that selects the wave. getmatrix() and CIplot() call the resolver below to
# translate such a request into (real matrix, row indices, column indices).
#
# The block layout is stored explicitly in @extramatrices$ri_clpm_blocks at
# construction time; riclpm_index() falls back to deriving it from the design
# matrix for models created before psychonetrics 0.17.3.


# Build the block index (also called by ri_clpm() at construction):
riclpm_makeindex <- function(vars){
  nVar  <- nrow(vars)
  nTime <- ncol(vars)
  list(
    nVar        = nVar,
    nTime       = nTime,
    varnames    = rownames(vars),
    timenames   = colnames(vars),
    # Latent variables: wave-1 ... wave-nTime components, then random intercepts
    latent_wave = c(rep(seq_len(nTime), each = nVar), rep(NA_integer_, nVar)),
    latent_var  = c(rep(seq_len(nVar), nTime), seq_len(nVar)),
    # Observed variables: wave-1 ... wave-nTime
    obs_wave    = rep(seq_len(nTime), each = nVar),
    obs_var     = rep(seq_len(nVar), nTime)
  )
}

# Retrieve the block index of a fitted / specified RI-CLPM model, or NULL when
# it cannot be resolved. Bootstrap aggregates carry no extramatrices of their
# own, so the index is taken from one of the aggregated models.
riclpm_index_safe <- function(x){
  if (is(x, "psychonetrics_bootstrap")){
    if (length(x@models) == 0) return(NULL)
    return(riclpm_index_safe(x@models[[1]]))
  }
  if (!is(x, "psychonetrics")) return(NULL)
  idx <- x@extramatrices$ri_clpm_blocks
  if (!is.null(idx)) return(idx)
  # Backward compatibility with models made before 0.17.3:
  vars <- x@extramatrices$vars
  if (is.null(vars) || !is.matrix(vars)) return(NULL)
  if (is.null(rownames(vars))) rownames(vars) <- paste0("V", seq_len(nrow(vars)))
  if (is.null(colnames(vars))) colnames(vars) <- paste0("T", seq_len(ncol(vars)))
  riclpm_makeindex(vars)
}

riclpm_index <- function(x){
  idx <- riclpm_index_safe(x)
  if (is.null(idx)){
    stop("This RI-CLPM model does not carry a design matrix and its blocks cannot be resolved.")
  }
  idx
}

# Is this an RI-CLPM model whose blocks can be resolved?
is_riclpm <- function(x){
  isTRUE(try(identical(x@submodel, "RI_CLPM"), silent = TRUE)) &&
    !is.null(riclpm_index_safe(x))
}

# Latent indices of one wave block, or of the random-intercept block:
riclpm_latents <- function(idx, wave){
  if (is.na(wave)) which(is.na(idx$latent_wave)) else which(!is.na(idx$latent_wave) & idx$latent_wave == wave)
}
# Observed indices of one wave:
riclpm_observed <- function(idx, wave){
  which(idx$obs_wave == wave)
}

# Translate a user-supplied 'time' (integer or wave label) to a wave number:
riclpm_time <- function(time, idx, valid, what){
  if (is.character(time)){
    time <- match(time, idx$timenames)
    if (anyNA(time)) stop("'time' does not match a wave label of the design matrix.")
  }
  time <- as.integer(time)
  if (anyNA(time) || any(!time %in% valid)){
    stop("'time' must be ", what, " (", paste(range(valid), collapse = "-"), ").")
  }
  time
}

# The typed covariance parameterizations that can carry a "_within"/"_between"
# suffix (all of them are block-diagonal over the RI-CLPM latents):
riclpm_zeta_matrices <- paste0(c("sigma", "omega", "kappa", "lowertri", "rho", "SD", "delta"),
                               "_zeta")


# --------------------------------------------------------------------------
# The resolver.
# --------------------------------------------------------------------------
# Returns NULL when 'matrix' is not a block request (so the caller falls back
# to its normal behavior), and otherwise a list with:
#
#   real   : the name of the underlying model matrix ("beta", "omega_zeta", ...)
#   kind    : "temporal" / "PDC" / "within" / "between" / "mean_within" /
#             "mean_between"
#   views   : a list of blocks, each with $label, $rows, $cols and $waves
#   rownames/colnames : the variable labels of a block
#
# Views are collapsed over waves that share exactly the same free parameters
# (i.e. after a stationarity constraint), so a stationary temporal structure
# yields one block instead of nTime - 1 identical copies.
#
# 'implicit' controls whether the plain model-matrix names ("beta", "nu",
# "omega_zeta") count as block requests when no 'time' is given. getmatrix()
# uses implicit = FALSE, so that those names keep returning the full matrix as
# in earlier versions; CIplot() uses implicit = TRUE, because a confidence-
# interval plot of the full latent-sized matrix mixes all waves and the random
# intercepts and is never what is wanted.
riclpm_resolve <- function(x, matrix, time, implicit = FALSE){
  if (!is_riclpm(x)) return(NULL)
  hastime <- !missing(time) && !is.null(time)
  blockable <- hastime || isTRUE(implicit)
  idx     <- riclpm_index(x)
  nVar    <- idx$nVar
  nTime   <- idx$nTime

  # --- Which kind of block is being requested? ---
  real <- NULL; kind <- NULL
  if (matrix %in% paste0(riclpm_zeta_matrices, "_within")){
    real <- sub("_within$", "", matrix); kind <- "within"
  } else if (matrix %in% paste0(riclpm_zeta_matrices, "_between")){
    real <- sub("_between$", "", matrix); kind <- "between"
  } else if (matrix == "nu_within"){
    real <- "nu";     kind <- "mean_within"
  } else if (matrix %in% c("nu_between", "nu_eta_between")){
    real <- "nu_eta"; kind <- "mean_between"
  } else if (matrix == "PDC"){
    real <- "beta";   kind <- "PDC"
  } else if (blockable && matrix == "beta"){
    real <- "beta";   kind <- "temporal"
  } else if (blockable && matrix == "nu"){
    real <- "nu";     kind <- "mean_within"
  } else if (blockable && matrix %in% riclpm_zeta_matrices){
    real <- matrix;   kind <- "within"
  } else {
    # Not a block request: the caller handles 'matrix' as usual.
    return(NULL)
  }

  # 'time' is meaningless for the between-person blocks:
  if (kind %in% c("between", "mean_between") && hastime){
    stop("'time' does not apply to a between-person block: the random intercepts are constant over waves.")
  }

  # --- Which waves? ---
  if (kind %in% c("temporal", "PDC")){
    if (nTime < 2) stop("Fewer than two waves: there are no temporal blocks.")
    waves <- if (hastime) riclpm_time(time, idx, seq_len(nTime - 1),
                                      "a transition index, i.e. the wave the effect starts from") else seq_len(nTime - 1)
  } else if (kind %in% c("within", "mean_within")){
    waves <- if (hastime) riclpm_time(time, idx, seq_len(nTime), "a wave index") else seq_len(nTime)
  } else {
    waves <- NA_integer_
  }

  # --- Row / column indices per wave ---
  blockof <- function(w){
    switch(kind,
           temporal     = list(rows = riclpm_latents(idx, w + 1L), cols = riclpm_latents(idx, w)),
           PDC          = list(rows = riclpm_latents(idx, w + 1L), cols = riclpm_latents(idx, w)),
           within       = { i <- riclpm_latents(idx, w);  list(rows = i, cols = i) },
           between      = { i <- riclpm_latents(idx, NA); list(rows = i, cols = i) },
           mean_within  = list(rows = riclpm_observed(idx, w), cols = 1L),
           mean_between = list(rows = riclpm_latents(idx, NA), cols = 1L))
  }
  labelof <- function(w){
    switch(kind,
           temporal     = paste0(idx$timenames[w], " -> ", idx$timenames[w + 1L]),
           PDC          = paste0(idx$timenames[w], " -> ", idx$timenames[w + 1L]),
           within       = idx$timenames[w],
           mean_within  = idx$timenames[w],
           "between")
  }

  # Which parameters determine a view's contents? For an ordinary block that is
  # simply that block of the underlying matrix; a PDC block is a function of
  # the temporal block AND of the innovation structure of the target wave.
  parblocksof <- function(w){
    if (kind == "PDC"){
      tgt   <- riclpm_latents(idx, w + 1L)
      zmats <- intersect(riclpm_zeta_matrices, unique(x@parameters$matrix))
      c(list(list(matrix = "beta",
                  rows = riclpm_latents(idx, w + 1L), cols = riclpm_latents(idx, w))),
        lapply(zmats, function(zm) list(matrix = zm, rows = tgt, cols = tgt)))
    } else {
      b <- blockof(w)
      list(list(matrix = real, rows = b$rows, cols = b$cols))
    }
  }

  views <- lapply(waves, function(w){
    b <- blockof(w)
    list(label = labelof(w), rows = b$rows, cols = b$cols, waves = w,
         parblocks = if (is.na(w)) NULL else parblocksof(w))
  })

  # --- Collapse waves that share exactly the same free parameters ---
  # (i.e. after ri_clpm_stationary(); avoids returning identical copies.)
  if (length(views) > 1 && !hastime){
    pars <- x@parameters
    # Bootstrap aggregates carry no 'par' column; then no collapsing is done.
    signature <- function(v){
      if (is.null(pars$par) || is.null(v$parblocks)) return(NA_character_)
      paste(vapply(v$parblocks, function(pb){
        sub  <- pars[pars$matrix == pb$matrix, , drop = FALSE]
        ri   <- match(sub$row, pb$rows)
        ci   <- match(sub$col, pb$cols)
        keep <- !is.na(ri) & !is.na(ci)
        key  <- paste0(ri[keep], "_", ci[keep])
        ord  <- order(key)
        paste0(pb$matrix, ":", paste(key[ord], sub$par[keep][ord], collapse = "|"))
      }, character(1)), collapse = "||")
    }
    sigs <- vapply(views, signature, character(1))
    if (!anyNA(sigs)){
      grp    <- match(sigs, unique(sigs))
      merged <- lapply(unique(grp), function(g){
        v  <- views[[which(grp == g)[1]]]
        ws <- unlist(lapply(views[grp == g], "[[", "waves"))
        if (length(ws) > 1){
          v$label <- if (length(ws) == length(views)){
            "stationary"
          } else if (all(diff(ws) == 1)){
            paste0(labelof(ws[1]), " ... ", labelof(ws[length(ws)]))
          } else paste(vapply(ws, labelof, character(1)), collapse = ", ")
          v$waves <- ws
        }
        v
      })
      views <- merged
    }
  }

  list(
    real     = real,
    kind     = kind,
    request  = matrix,
    views    = views,
    rowlabs  = idx$varnames,
    collabs  = if (kind %in% c("mean_within", "mean_between")) NULL else idx$varnames,
    idx      = idx
  )
}


# Compute the PDC (partial directed correlation) block of a transition.
# The RI-CLPM has no PDC model matrix (it is not a parameterization option
# there), so it is derived from the temporal block and the innovation
# precision of the TARGET wave -- the innovation that enters at wave t + 1.
# In a fully stationary RI-CLPM every wave-2+ innovation block is the same,
# so this reproduces exactly the panelvar/panelgvar PDC. Like the PDC of the
# other frameworks, the result has from = row and to = column.
riclpm_PDC <- function(x, view, groupID){
  idx  <- riclpm_index(x)
  w    <- view$waves[1]
  kap  <- riclpm_latents(idx, w + 1L)
  lapply(groupID, function(g){
    mm    <- x@modelmatrices[[g]]
    beta  <- as.matrix(mm$beta)[view$rows, view$cols, drop = FALSE]
    kappa <- as.matrix(mm$kappa_zeta)[kap, kap, drop = FALSE]
    out   <- computePDC(beta, kappa)
    dimnames(out) <- list(idx$varnames, idx$varnames)
    out
  })
}


# Subset a full model matrix to a block and attach variable labels:
riclpm_subset <- function(mat, view, spec){
  out <- mat[view$rows, view$cols, drop = FALSE]
  if (nrow(out) == length(spec$rowlabs)) rownames(out) <- spec$rowlabs
  if (!is.null(spec$collabs) && ncol(out) == length(spec$collabs)){
    colnames(out) <- spec$collabs
  } else {
    colnames(out) <- NULL
  }
  out
}
