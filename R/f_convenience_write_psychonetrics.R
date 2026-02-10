# Write comprehensive psychonetrics output to a text file

write_psychonetrics <- function(x, file = "psychonetrics_output.txt",
                                 matrices = TRUE, MIs = TRUE, logbook = TRUE) {

  # Validation:
  stopifnot(is(x, "psychonetrics"))
  if (!x@computed) {
    warning("Model has not been computed. Output will be incomplete.")
  }

  LINE_WIDTH <- 72

  # Accumulate all output lines:
  lines <- character(0)

  lines <- c(lines, .wp_header(x, LINE_WIDTH))
  lines <- c(lines, .wp_general(x, LINE_WIDTH))
  lines <- c(lines, .wp_sample(x, LINE_WIDTH))
  lines <- c(lines, .wp_model(x, LINE_WIDTH))

  if (x@computed) {
    lines <- c(lines, .wp_parameters(x, LINE_WIDTH))

    if (!is.null(x@fitmeasures)) {
      lines <- c(lines, .wp_fit(x, LINE_WIDTH))
    }

    if (matrices) {
      lines <- c(lines, .wp_matrices(x, LINE_WIDTH))
    }

    if (MIs && any(!is.na(x@parameters$mi))) {
      lines <- c(lines, .wp_MIs(x, LINE_WIDTH))
    }
  }

  if (logbook) {
    lines <- c(lines, .wp_logbook(x, LINE_WIDTH))
  }

  lines <- c(lines, .wp_footer(LINE_WIDTH))

  # Write to file:
  writeLines(lines, con = file)
  message("Output written to: ", file)
  invisible(file)
}


# --- Internal formatting helpers ---

.wp_separator <- function(width = 72, char = "=") {
  paste0(rep(char, width), collapse = "")
}

.wp_section_title <- function(title, width = 72) {
  sep <- .wp_separator(width, "=")
  pad <- max(0, (width - nchar(title)) %/% 2)
  centered <- paste0(paste0(rep(" ", pad), collapse = ""), title)
  c("", sep, centered, sep, "")
}

.wp_indent <- function(lines, n = 2) {
  pad <- paste0(rep(" ", n), collapse = "")
  paste0(pad, lines)
}

# Capture print.data.frame output without line wrapping:
.wp_capture_df <- function(df, indent = 2) {
  old_width <- getOption("width")
  options(width = 10000)
  on.exit(options(width = old_width))
  table_lines <- capture.output(print.data.frame(df, row.names = FALSE))
  paste0(paste0(rep(" ", indent), collapse = ""), table_lines)
}

# Capture print (for matrices) without line wrapping:
.wp_capture_mat <- function(mat, indent = 6) {
  old_width <- getOption("width")
  options(width = 10000)
  on.exit(options(width = old_width))
  mat_lines <- capture.output(print(mat))
  paste0(paste0(rep(" ", indent), collapse = ""), mat_lines)
}

# Wrap a long list of variable names:
.wp_wrap_vars <- function(vars, indent = 24, width = 72) {
  current <- ""
  result <- character(0)
  pad <- paste0(rep(" ", indent), collapse = "")
  for (v in vars) {
    test <- if (nchar(current) == 0) v else paste0(current, " ", v)
    if (nchar(test) + indent > width && nchar(current) > 0) {
      result <- c(result, current)
      current <- v
    } else {
      current <- test
    }
  }
  if (nchar(current) > 0) result <- c(result, current)
  paste(result, collapse = paste0("\n", pad))
}


# --- Section: Header / Logo ---

.wp_header <- function(x, width = 72) {
  version <- read.dcf(file = system.file("DESCRIPTION", package = "psychonetrics"),
                       fields = "Version")
  version_string <- paste0("Version: ", version)

  logo <- c(
    "                       _                      _        _          ",
    "                      | |                    | |      (_)         ",
    "  _ __  ___ _   _  ___| |__   ___  _ __   ___| |_ _ __ _  ___ ___ ",
    " |  _ \\/ __| | | |/ __|  _ \\ / _ \\|  _ \\ / _ \\ __|  __| |/ __/ __|",
    " | |_) \\__ \\ |_| | (__| | | | (_) | | | |  __/ |_| |  | | (__\\__ \\",
    " | .__/|___/\\__, |\\___|_| |_|\\___/|_| |_|\\___|\\__|_|  |_|\\___|___/",
    " | |         __/ |                                                ",
    paste0(" |_|        |___/  ",
           paste0(rep(" ", max(0, 66 - 19 - nchar(version_string))), collapse = ""),
           version_string)
  )

  c(logo, "",
    paste0("  Output generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
}


# --- Section: General info ---

.wp_general <- function(x, width = 72) {
  version <- read.dcf(file = system.file("DESCRIPTION", package = "psychonetrics"),
                       fields = "Version")

  last_time <- if (length(x@log) > 0) {
    as.character(x@log[[length(x@log)]]@time)
  } else {
    "(unknown)"
  }

  lines <- .wp_section_title("GENERAL INFORMATION", width)
  lines <- c(lines,
    paste0("  psychonetrics version : ", version),
    paste0("  Model last edited     : ", last_time),
    paste0("  R version             : ", R.version.string),
    paste0("  Platform              : ", R.version$platform)
  )
  lines
}


# --- Section: Sample info ---

.wp_sample <- function(x, width = 72) {
  groups <- x@sample@groups
  vars <- x@sample@variables

  lines <- .wp_section_title("SAMPLE INFORMATION", width)
  lines <- c(lines,
    paste0("  Number of cases               : ", sum(groups$nobs)),
    paste0("  Number of groups              : ", nrow(groups)),
    paste0("  Number of observed statistics : ", x@sample@nobs),
    paste0("  Correlation input             : ", x@sample@corinput)
  )

  # Groups:
  if (nrow(groups) > 1 || groups$label[1] != "fullsample") {
    lines <- c(lines, "", "  Groups:")
    for (i in seq_len(nrow(groups))) {
      lines <- c(lines,
        paste0("    ", groups$label[i], " (n = ", groups$nobs[i], ")"))
    }
  }

  # Variables:
  ordered_vars <- vars$label[vars$ordered]
  continuous_vars <- vars$label[!vars$ordered]
  lines <- c(lines, "", "  Observed variables:")
  if (length(continuous_vars) > 0) {
    lines <- c(lines,
      paste0("    Continuous (", length(continuous_vars), "): ",
             .wp_wrap_vars(continuous_vars, indent = 24, width = width)))
  }
  if (length(ordered_vars) > 0) {
    lines <- c(lines,
      paste0("    Ordered    (", length(ordered_vars), "): ",
             .wp_wrap_vars(ordered_vars, indent = 24, width = width)))
  }

  lines
}


# --- Section: Model specification ---

.wp_model <- function(x, width = 72) {
  # Model name (same switch as printMethod.R):
  mod <- switch(x@model,
    "gvar" = "Graphical vector-autoregression (GVAR)",
    "varcov" = "Variance-covariance matrix (varcov)",
    "lvm" = "Latent variable model (LVM)",
    "var1" = "Lag-1 vector-autoregression (VAR1)",
    "panelvar1" = "Lag-1 panel vector auto-regression (panelvar1)",
    "dlvm1" = "Lag-1 dynamic latent variable model for panel data (dlvm1)",
    "tsdlvm1" = "Lag-1 dynamic latent variable model for time-series data (tsdlvm1)",
    "meta_varcov" = "Variance-covariance matrix meta analysis",
    "Ising" = "Ising model",
    "ml_lvm" = "Multi-level latent variable model",
    x@model
  )

  submod <- switch(x@submodel,
    "none" = "None",
    "lnm" = "Latent Network Model (LNM)",
    "ggm" = "Gaussian graphical model (GGM)",
    "rnm" = "Residual network model (RNM)",
    "cholesky" = "Cholesky decomposition (cholesky)",
    "sem" = "Structural equation model (SEM)",
    "lrnm" = "Latent & residual network model (LRNM)",
    "gvar" = "Graphical vector-autoregression (GVAR)",
    "var" = "Vector-autoregression (VAR)",
    "ml_lnm" = "Multi-level latent network model",
    "ml_rnm" = "Multi-level residual network model",
    x@submodel
  )
  if (is.null(submod)) submod <- x@submodel

  lines <- .wp_section_title("MODEL SPECIFICATION", width)
  lines <- c(lines,
    paste0("  Model type            : ", mod),
    paste0("  Submodel              : ", submod),
    paste0("  Number of parameters  : ", max(x@parameters$par)),
    paste0("  Distribution          : ", x@distribution),
    paste0("  Identification        : ", x@identification),
    paste0("  Mean structure        : ", ifelse(x@meanstructure, "Modeled", "Not modeled"))
  )

  if (x@computed) {
    estimator <- switch(x@estimator,
      "ML" = "Maximum likelihood estimation (ML)",
      "FIML" = "Full information maximum likelihood (FIML)",
      "ULS" = "Unweighted least squares (ULS)",
      "WLS" = "Weighted least squares (WLS)",
      "DWLS" = "Diagonally weighted least squares (DWLS)",
      x@estimator
    )
    lines <- c(lines,
      paste0("  Estimator             : ", estimator),
      paste0("  Optimizer             : ", x@optim$optimizer),
      paste0("  Optimization message  : ", x@optim$message)
    )
  }

  # Matrix list:
  lines <- c(lines, "", "  Matrices in model:")
  for (i in seq_len(nrow(x@matrices))) {
    mat_info <- x@matrices[i, ]
    props <- c()
    if (mat_info$symmetrical) props <- c(props, "symmetric")
    if (mat_info$diagonal) props <- c(props, "diagonal")
    if (mat_info$sparse) props <- c(props, "sparse")
    prop_str <- if (length(props) > 0) paste0(" (", paste(props, collapse = ", "), ")") else ""
    lines <- c(lines,
      paste0("    - ", mat_info$name, " [", mat_info$nrow, " x ", mat_info$ncol, "]", prop_str))
  }

  lines
}


# --- Section: Parameter estimates ---

.wp_parameters <- function(x, width = 72) {
  lines <- .wp_section_title("PARAMETER ESTIMATES", width)

  parTable <- x@parameters

  # Determine columns (same logic as parameters()):
  has_boots <- !all(is.na(parTable$se_boot))
  if (has_boots) {
    cols <- c("var1", "op", "var2", "est", "se", "p", "se_boot", "p_boot",
              "matrix", "row", "col", "group", "par")
  } else {
    cols <- c("var1", "op", "var2", "est", "se", "p",
              "matrix", "row", "col", "group", "par")
  }

  # Filter non-fixed or non-zero:
  parTable <- parTable %>%
    filter(drop(!.data[["fixed"]] | .data[["est"]] != 0)) %>%
    select(all_of(cols))

  parTable$var2 <- ifelse(is.na(parTable$var2), "", parTable$var2)

  # Format numbers:
  parTable$est <- goodNum2(parTable$est)
  parTable$se <- goodNum(parTable$se)
  parTable$p <- goodNum(parTable$p)
  if (has_boots) {
    parTable$se_boot <- goodNum(parTable$se_boot)
    parTable$p_boot <- goodNum(parTable$p_boot)
  }

  # For each group/matrix:
  for (g in x@sample@groups$label) {
    lines <- c(lines, "", paste0("  Parameters for group: ", g))

    for (mat in unique(parTable$matrix[parTable$group == g])) {
      if (x@matrices$diagonal[x@matrices$name == mat]) {
        brackets <- " (diagonal)"
      } else if (x@matrices$symmetrical[x@matrices$name == mat]) {
        brackets <- " (symmetric)"
      } else {
        brackets <- ""
      }

      lines <- c(lines, "", paste0("    Matrix: ", mat, brackets))

      subTable <- parTable %>%
        filter(.data[["group"]] == g, .data[["matrix"]] == mat) %>%
        select(-.data[["matrix"]], -.data[["group"]])

      lines <- c(lines, .wp_capture_df(subTable, indent = 4))
    }
  }

  lines
}


# --- Section: Fit measures ---

.wp_fit <- function(x, width = 72) {
  lines <- .wp_section_title("FIT MEASURES", width)

  df <- data.frame(
    Measure = names(x@fitmeasures),
    Value = goodNum2(unlist(x@fitmeasures)),
    stringsAsFactors = FALSE
  )

  lines <- c(lines, .wp_capture_df(df, indent = 2))

  lines
}


# --- Section: Model matrices ---

.wp_matrices <- function(x, width = 72) {
  lines <- .wp_section_title("MODEL MATRICES", width)

  varLabels <- x@sample@variables$label
  nObs <- length(varLabels)

  for (g in seq_len(nrow(x@sample@groups))) {
    group_label <- x@sample@groups$label[g]
    lines <- c(lines, "", paste0("  Group: ", group_label),
               paste0("  ", .wp_separator(width - 2, "-")))

    for (mat_name in names(x@modelmatrices[[g]])) {
      mat <- tryCatch(as.matrix(x@modelmatrices[[g]][[mat_name]]),
                       error = function(e) NULL)
      if (is.null(mat)) next

      # Skip very large matrices:
      if (nrow(mat) > 50 || ncol(mat) > 50) {
        lines <- c(lines, "",
          paste0("    ", mat_name, " [", nrow(mat), " x ", ncol(mat),
                 "] -- omitted (matrix too large)"))
        next
      }

      # Infer dimnames from parameter table:
      dimnames_inferred <- .wp_infer_dimnames(x, mat_name, g, nrow(mat), ncol(mat))
      if (!is.null(dimnames_inferred$rows)) rownames(mat) <- dimnames_inferred$rows
      if (!is.null(dimnames_inferred$cols)) colnames(mat) <- dimnames_inferred$cols

      lines <- c(lines, "", paste0("    ", mat_name, ":"))

      # Format: round to 2 digits, replace exact zeros with "."
      formatted <- matrix(
        ifelse(mat == 0, ".", formatC(round(mat, 2), format = "f", digits = 2)),
        nrow = nrow(mat), ncol = ncol(mat),
        dimnames = dimnames(mat)
      )
      lines <- c(lines, .wp_capture_mat(noquote(formatted), indent = 6))
    }
  }

  lines
}

# Infer row/column names from the parameter table:
.wp_infer_dimnames <- function(x, mat_name, group_id, nrow_mat, ncol_mat) {
  pars <- x@parameters[x@parameters$matrix == mat_name &
                         x@parameters$group_id == group_id, ]

  if (nrow(pars) == 0) return(list(rows = NULL, cols = NULL))

  # Row names from var1:
  row_names <- rep(NA_character_, nrow_mat)
  for (i in seq_len(nrow(pars))) {
    r <- pars$row[i]
    if (r >= 1 && r <= nrow_mat && is.na(row_names[r])) {
      row_names[r] <- pars$var1[i]
    }
  }

  # Column names from var2:
  col_names <- rep(NA_character_, ncol_mat)
  for (i in seq_len(nrow(pars))) {
    v2 <- pars$var2[i]
    cc <- pars$col[i]
    if (!is.na(v2) && cc >= 1 && cc <= ncol_mat && is.na(col_names[cc])) {
      col_names[cc] <- v2
    }
  }

  # Fill gaps with indices:
  row_names[is.na(row_names)] <- paste0("[", which(is.na(row_names)), "]")
  col_names[is.na(col_names)] <- paste0("[", which(is.na(col_names)), "]")

  list(rows = row_names, cols = col_names)
}


# --- Section: Modification indices ---

.wp_MIs <- function(x, width = 72, top = 20) {
  lines <- .wp_section_title("MODIFICATION INDICES", width)

  if (all(is.na(x@parameters$mi))) {
    lines <- c(lines, "  No modification indices available. Use addMIs() to compute them.")
    return(lines)
  }

  miTable <- x@parameters %>%
    filter(!is.na(.data[["mi"]])) %>%
    select(.data[["var1"]], .data[["op"]], .data[["var2"]],
           .data[["mi"]], .data[["pmi"]], .data[["epc"]],
           .data[["matrix"]], .data[["group"]]) %>%
    arrange(desc(.data[["mi"]])) %>%
    head(top)

  miTable$var2 <- ifelse(is.na(miTable$var2), "", miTable$var2)
  miTable$mi <- goodNum(miTable$mi)
  miTable$pmi <- goodNum(miTable$pmi)
  miTable$epc <- goodNum2(miTable$epc)

  lines <- c(lines, paste0("  Top ", min(top, nrow(miTable)), " modification indices:"), "")

  lines <- c(lines, .wp_capture_df(miTable, indent = 2))

  lines
}


# --- Section: Logbook ---

.wp_logbook <- function(x, width = 72) {
  lines <- .wp_section_title("LOGBOOK", width)

  log <- x@log

  if (length(log) == 0) {
    lines <- c(lines, "  (no log entries)")
    return(lines)
  }

  for (i in seq_along(log)) {
    entry <- log[[i]]
    time_str <- format(entry@time, "%Y-%m-%d %H:%M:%S")
    lines <- c(lines, paste0("  [", time_str, "] ", entry@event))
  }

  lines
}


# --- Section: Footer ---

.wp_footer <- function(width = 72) {
  c("",
    .wp_separator(width, "="),
    paste0("  End of psychonetrics output - generated ",
           format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    .wp_separator(width, "="),
    "")
}
