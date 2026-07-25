colorblind <- function (n, shift = 0) 
{
  if (n > 7) 
    warning("'colorblind' palette only supports 8 colors.")
  Palette <- rgb(c(230, 86, 0, 240, 204, 213, 0), c(159, 180, 
                                                    158, 228, 121, 94, 114), c(0, 233, 115, 66, 167, 0, 178), 
                 maxColorValue = 255)
  Palette[(((shift + 1:n) - 1)%%8) + 1]
}

# Small inner function:
covchooser <- function(x,name=""){
  matrices <- "none"
  if (x == "ggm"){
    matrices <- "omega"
  } else if (x == "prec"){
    matrices <- "kappa"
  } else if (x == "chol"){
    matrices <- "lowertri"
  } else if (x == "cov"){
    matrices <- "sigma"
  } else if (x == "cor"){
    matrices <- "rho"
  }
  if (name != ""){
    matrices <- paste0(matrices,"_",name)  
  }
  
  return(matrices)
}

CIplot <- function(
    x,
    matrices,
    alpha_ci = 0.05,
    alpha_color = c(0.05,0.01,0.001,0.0001),
    labels,
    labels2,
    labelstart,
    print = TRUE,
    major_break = 0.2,
    minor_break = 0.1,
    split0,
    prop0,
    prop0_cex = 1,
    prop0_alpha = 0.95,
    prop0_minAlpha = 0.25,
    time   # Only for ri_clpm()/ri_clpn() models: the wave (block) to plot
){
  # ggplot2 is only needed for the plotting and is a Suggests dependency:
  if(!requireNamespace("ggplot2")) stop("'ggplot2' package needs to be installed.")

  # Check if the model is a bootstrap aggregate:
  if (is(x,"psychonetrics_bootstrap")){
    boot_agg <- TRUE
  } else {
    stopifnot(is(x,"psychonetrics"))
    # Bootstrap warning:
    if (x@sample@bootstrap){
      boot_warning()
    }
    boot_agg <- FALSE
  }
  
  
  
  if (!boot_agg){
    logs <- sapply(x@log,function(x)x@event)
    if (any(grepl("Pruned",logs)) | any(grepl("step-up",logs))  | any(grepl("modelsearch",logs)) | any(grepl("Post-selection refit",logs))){
      warning("A model search algorithm was used in creating this model. The CIs are likely invalid.")
    }
  }
  
  
  
  
  # 'time' selects a wave (block) and only exists for RI-CLPM models:
  if (!missing(time) && !is_riclpm(x)){
    stop("The 'time' argument is only available for ri_clpm() / ri_clpn() models.")
  }

  if (missing(labels2) && !missing(labels)){
    labels2 <- labels
  }
  
  
  # Stop if model is not psychonetrics:
  # stopifnot(is(x,"psychonetrics"))
  
  # Default for matrices argument:
  if (missing(matrices)){
    matrices <- "none"
    
    # Varcov default:
    if (x@model == "varcov"){
      matrices <- covchooser(x@submodel)
    }  else if (x@model == "meta_varcov"){
      matrices <- covchooser(x@types$y,"y")
    }  else if (x@model == "lvm"){

      if (is_riclpm(x)){

        # RI-CLPM / RI-CLPN: mirror the panelvar default, but on the blocks
        # (temporal, contemporaneous within-person, between-person) rather
        # than on the full latent-sized matrices, which would mix all waves
        # and the random intercepts into one unreadable plot:
        rictype <- if (boot_agg) NULL else x@extramatrices$ri_clpm_type
        if (is.null(rictype)) rictype <- x@types$latent
        matrices <- c("beta",
                      covchooser(rictype,"zeta_within"),
                      covchooser(rictype,"zeta_between"))

      } else if (x@submodel == "rnm"){

        matrices <- covchooser(x@types$residual,"epsilon")

      } else if (x@submodel == "lnm"){
        
        matrices <- covchooser(x@types$latent,"zeta")
        
      } else if (x@submodel == "lrnm"){
        
        matrices <- c(covchooser(x@types$latent,"zeta"), covchooser(x@types$residual,"epsilon"))
        
      } else {
        
        matrices <- covchooser(x@types$latent,"zeta")
        
      }
    } else if (x@model == "var1"){
      
      matrices <- c("beta",covchooser(x@types$zeta,"zeta"))
      
    } else if (x@model %in% c("panelvar","ml_var1")){
      matrices <- c(
        "beta",
        covchooser(x@types$within_latent,"zeta_within"),
        covchooser(x@types$between_latent,"zeta_between")
      )

    } else if (x@model %in% c("ml_lvm","dlvm1")){
      # The between-person block is typed by 'between_latent'; using
      # 'within_residual' here was a bug (invisible whenever the two types
      # coincide, which they do by default):
      matrices <- c(
        "beta",
        covchooser(x@types$within_latent,"zeta_within"),
        covchooser(x@types$between_latent,"zeta_between")
      )
      
    } else if (x@model == "tsdlvm1"){
      matrices <- c(
        "beta",
        covchooser(x@types$zeta,"zeta")
      )

    } else if (x@model == "ml_varcov"){
      matrices <- c(
        covchooser(x@types$within,"within"),
        covchooser(x@types$between,"between")
      )

    }  else if (x@model %in% c("Ising", "BlumeCapel")){
      matrices <- c("omega")
      
    }  else stop("No default argument for 'matrices' for current model.")

    # PDC temporal parameterization: the temporal matrix is "PDC", not "beta":
    if ((!is.null(x@types$temporal) && x@types$temporal == "PDC") ||
        (!is.null(x@types$temporal_latent) && x@types$temporal_latent == "PDC")){
      matrices[matrices == "beta"] <- "PDC"
    }
  }
  if (any(matrices == "none")){
    stop("No default matrix could be selected")
  }
  
  
  # Extract parameters:
  pars <- x@parameters

  # Dimensions of each plotted matrix (used to check 'labels'/'labels2'):
  matdims <- data.frame(name = x@matrices$name, nrow = x@matrices$nrow,
                        ncol = x@matrices$ncol, stringsAsFactors = FALSE)

  # --- RI-CLPM / RI-CLPN: expand block requests ---
  # A matrix name such as "omega_zeta_within" (or a plain "beta") refers to an
  # nVar x nVar block of the full latent-sized matrix. Each requested block
  # becomes its own pseudo-matrix in 'pars', with the wave in its name and the
  # plain variable names as node labels, so that the plotting code below can
  # stay untouched. Waves that share the same free parameters (after a
  # stationarity constraint) are collapsed into a single block.
  if (is_riclpm(x)){
    newmats <- character(0)
    for (mat in matrices){
      spec <- if (missing(time)){
        riclpm_resolve(x, mat, implicit = TRUE)
      } else riclpm_resolve(x, mat, time, implicit = TRUE)

      if (is.null(spec)){
        newmats <- c(newmats, mat)
        next
      }
      if (spec$kind == "PDC"){
        stop("The RI-CLPM has no PDC parameters (PDC is derived from beta and the innovation precision), so its confidence intervals cannot be plotted. Use matrices = 'beta' instead, or getmatrix(x, 'PDC', time = ...) for the point estimates.")
      }
      for (v in spec$views){
        lab <- if (spec$kind %in% c("between","mean_between")){
          spec$request
        } else paste0(spec$request, " (", v$label, ")")

        sub  <- pars[pars$matrix == spec$real, , drop = FALSE]
        ri   <- match(sub$row, v$rows)
        ci   <- match(sub$col, v$cols)
        keep <- !is.na(ri) & !is.na(ci)
        sub  <- sub[keep, , drop = FALSE]
        if (nrow(sub) == 0) next

        sub$matrix  <- lab
        sub$var1_id <- ri[keep]
        sub$var2_id <- ci[keep]
        sub$var1    <- spec$rowlabs[sub$var1_id]
        sub$var2    <- if (is.null(spec$collabs)) "" else spec$collabs[sub$var2_id]

        pars    <- rbind(pars, sub)
        newmats <- c(newmats, lab)
        matdims <- rbind(matdims, data.frame(name = lab, nrow = length(v$rows),
                                             ncol = length(v$cols),
                                             stringsAsFactors = FALSE))
      }
    }
    if (length(newmats) == 0){
      stop("None of the requested matrices contain parameters to plot.")
    }
    matrices <- newmats
  }

  # Make a data frame for each matrix:
  plots <- lapply(matrices,function(mat,labels,labels2,labelstart,split0,prop0){
    
    # Type of edge:
    edge <- "--"
    if (grepl("lambda",mat)){
      edge <- "->"
    } else if (grepl("beta",mat)){
      edge <- "->"
    } else if (grepl("PDC",mat)){
      # PDC matrices encode from = row (var1) and to = column (var2), so the
      # arrow direction follows the default var1 -> var2 node order (no flip
      # as for beta/lambda):
      edge <- "->"
    } else if (grepl("sigma",mat)){
      edge <- "<->"
    } else if (grepl("kappa",mat)){
      edge <- "--"
    } else if (grepl("lowerTri",mat)){
      edge <- "-c-"
    } 
    
    # symmetrical <- x@matrices$symmetrical[x@matrices$name == mat]
    # diagonal  <- x@matrices$diagonal[x@matrices$name == mat]
    # lowertri   <- x@matrices$lowertri[x@matrices$name == mat]
    ncol   <- matdims$ncol[matdims$name == mat]
    nrow   <- matdims$nrow[matdims$name == mat]
    
    # Edgelist:
    if (missing(labels)){
      if (grepl("lambda",mat) |grepl("beta",mat)){
        node1 <- pars$var2[pars$matrix == mat]
      } else {
        node1 <- pars$var1[pars$matrix == mat]
      }      
    } else {
      if (grepl("lambda",mat) |grepl("beta",mat)){
        if (length(labels) != ncol){
          stop("'labels' is not of the correct length")
        }
        node1 <- labels[pars$var2_id[pars$matrix == mat]]
      } else {
        if (length(labels) != nrow){
          stop("'labels' is not of the correct length")
        }
        node1 <- labels[pars$var1_id[pars$matrix == mat]]
      } 
    }
    
    if (missing(labels2)){
      if (grepl("lambda",mat) |grepl("beta",mat)){
        node2 <- pars$var1[pars$matrix == mat]
      } else {
        node2 <- pars$var2[pars$matrix == mat]
      }      
    } else {
      if (grepl("lambda",mat) |grepl("beta",mat)){
        if (length(labels2) != nrow){
          stop("'labels2' is not of the correct length")
        }
        node2 <- labels2[pars$var1_id[pars$matrix == mat]]
      } else {
        if (length(labels2) != ncol){
          stop("'labels2' is not of the correct length")
        }
        node2 <- labels2[pars$var2_id[pars$matrix == mat]]
      } 
    }
    
    # Edges :
    edges <- paste0(node1," ",edge," ",node2)
    
    # FIXME: UGLY REPETATIVE CODE FOR NOW
    
    #  Non bootstrap:
    if (!boot_agg){
      
      
      # Est:
      est <-  pars$est[pars$matrix == mat]
      se <- pars$se[pars$matrix == mat]
      p <- pars$p[pars$matrix == mat]
      
      # Compute CI:
      ci_lower <- est - qnorm(1-alpha_ci/2) * se
      ci_upper <- est + qnorm(1-alpha_ci/2) * se
      
      df_cis <- data.frame(
        edge = edges,
        est = est,
        p = p,
        lower = ci_lower,
        upper = ci_upper
      )
      
      # Fix CI:
      df_cis$lower[is.na(df_cis$lower)] <- df_cis$est[is.na(df_cis$lower)] 
      df_cis$upper[is.na(df_cis$upper)] <- df_cis$est[is.na(df_cis$upper)]
      
      # Order edges by estimate:
      df_cis$edge <- factor(df_cis$edge, levels = df_cis$edge[order(df_cis$est)])
      
      # Bounds:
      bound_upper <- max(df_cis$upper, na.rm=TRUE)
      bound_lower <- min(df_cis$lower, na.rm=TRUE)
      
      # Significance:
      alpha_color <- sort(alpha_color,decreasing = TRUE)
      alpha_color_lab <- format(alpha_color, scientific=FALSE,drop0trailing=TRUE)
      df_cis$sig <- paste0("p > ",alpha_color_lab[1])
      for (i in seq_along(alpha_color)){
        df_cis$sig[df_cis$p < alpha_color[i]] <- paste0("p < ",alpha_color_lab[i])
      }
      
      # Set levels:
      df_cis$sig <- factor(df_cis$sig,levels=c(paste0("p > ",alpha_color_lab[1]),paste0("p < ",alpha_color_lab)))
      
      # Add matrix:
      df_cis$matrix <- mat
      
      # If labelstart is missing, choose it:
      if (missing(labelstart)){
        labelstart <- quantile(df_cis$est,0.95)
      }
      
      colors <- colorblind(7)[c(3,1,6,5,2,4,7)][seq_along(alpha_color)]
      
      df_cis$labstart <- ifelse(df_cis$upper > labelstart, df_cis$lower, df_cis$upper)
      df_cis$hjust <- ifelse(df_cis$upper > labelstart, 1.05, -0.05)
      
      plots <- ggplot2::ggplot(df_cis,
                               ggplot2::aes(x = .data[["edge"]], y = .data[["est"]], colour = .data[["sig"]],
                                                   ymin = .data[["lower"]], ymax = .data[["upper"]])) +
        ggplot2::geom_hline(yintercept = 0, alpha = 0.2) +
        ggplot2::geom_errorbar() + 
        ggplot2::geom_point(cex = 1.5) +
        # geom_point(cex = 1.5, aes(y = edge_pruned)) +
        # geom_point(cex = 0.5*1.5, aes(y = edge_pruned), colour = "white") +
        ggplot2::geom_text(ggplot2::aes(y = .data[["labstart"]], label = .data[["edge"]], hjust = .data[["hjust"]]), colour = rgb(0.2,0.2,0.2), cex = 2.5) +
        ggplot2::coord_flip() +
        # facet_grid( ~  model) +
        ggplot2::scale_y_continuous(breaks = seq(floor(min(df_cis$lower)),ceiling(max(df_cis$upper)),by=major_break),
                                    minor_breaks = seq(floor(min(df_cis$lower)),ceiling(max(df_cis$upper)),by=minor_break)) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::scale_colour_manual("",values = c("black",colors),drop=FALSE) +  ggplot2::theme(legend.position = "top", # axis.text.y = element_text(size = 7),
                                                                                                 axis.text.y=ggplot2::element_blank(),
                                                                                                 axis.ticks.y=ggplot2::element_blank(),
                                                                                                 panel.grid.major.y = ggplot2::element_blank()) + 
        ggplot2::xlab("") + 
        ggplot2::ylab("") 
      
      # Return plot:
      return(plots)
      
    } else {
      
      # Check for model selection for split0 and plot0:
      if (missing(split0)){
        logs <- sapply(x@models[[1]]@log,function(x)x@event)
        if (any(grepl("Pruned",logs)) | any(grepl("step-up",logs))  | any(grepl("modelsearch",logs)) | any(grepl("Post-selection refit",logs))){
          split0 <- TRUE
        } else {
          split0 <- FALSE
        }
      }
      if (missing(prop0)){
        prop0 <- split0
      }
      
      
      # Bootstraps:  
      if (alpha_ci != 0.05){
        stop("Only alpha = 0.05 supported for bootstrap samples.")
      }
      
      prop0_boots <- 1 - pars$prop_non0[pars$matrix == mat]
      sample_est <- pars$est_sample[pars$matrix == mat]
      
      if (split0){
        avg <- pars$avg_non0[pars$matrix == mat]
        ci_lower <- pars$q2.5_non0[pars$matrix == mat]
        ci_upper <- pars$q97.5_non0[pars$matrix == mat]
        
        # Alpha:
        alpha <- prop0_minAlpha + (1-prop0_minAlpha) * (1-prop0_boots)
        alpha[!is.finite(avg)] <- 0
        
        avg[!is.finite(avg)] <- 0
        ci_lower[!is.finite(ci_lower)] <- 0
        ci_upper[!is.finite(ci_upper)] <- 0
        
        
        
      } else {
        avg <- pars$avg[pars$matrix == mat]
        ci_lower <- pars$q2.5[pars$matrix == mat]
        ci_upper <- pars$q97.5[pars$matrix == mat]
        alpha <- 1
      }
      
      # bootstrap DF:
      df_cis_boot <- data.frame(
        edge = edges,
        est = avg,
        lower = ci_lower,
        upper = ci_upper,
        prop0=prop0_boots,
        alpha=alpha,
        type = "boot"
      )
      
      # sample DF:
      df_cis_sample <- data.frame(
        edge = edges,
        est = sample_est,
        alpha=1,
        type = "sample"
      )
      
      # Combine:
      df_cis <- bind_rows(df_cis_boot,df_cis_sample)
      
      # Fix CI:
      df_cis$lower[is.na(df_cis$lower)] <- df_cis$est[is.na(df_cis$lower)] 
      df_cis$upper[is.na(df_cis$upper)] <- df_cis$est[is.na(df_cis$upper)]
      
      # Order edges by estimate:
      df_cis$edge <- factor(df_cis$edge, levels = df_cis_sample$edge[order(df_cis_sample$est)])
      
      # Bounds:
      bound_upper <- max(df_cis$upper, na.rm=TRUE)
      bound_lower <- min(df_cis$lower, na.rm=TRUE)
      
      
      # Add matrix:
      df_cis$matrix <- mat
      
      # If labelstart is missing, choose it:
      if (missing(labelstart)){
        labelstart <- quantile(df_cis$est,0.95,na.rm=TRUE)
      }
      
      
      df_cis$labstart <- ifelse(df_cis$upper > labelstart, df_cis$lower, df_cis$upper)
      # Nudge zero labels:
      if (split0){
        df_cis$labstart <- ifelse(abs(df_cis$labstart) < max(abs(df_cis$labstart))/50, max(abs(df_cis$labstart))/50, df_cis$labstart)  
      }
      
      # Types:
      df_cis$type <- factor(df_cis$type, levels = c("boot","sample"),labels = c("Bootstrap mean","Sample"))
      
      df_cis$hjust <- ifelse(df_cis$upper > labelstart, 1.05, -0.05)
      
      plots <- ggplot2::ggplot(df_cis %>% filter(.data[["type"]] == "Bootstrap mean"),
                               ggplot2::aes(x = .data[["edge"]], y = .data[["est"]],
                                                   ymin = .data[["lower"]], ymax = .data[["upper"]])) +
        ggplot2::geom_hline(yintercept = 0, alpha = 0.2) +
        ggplot2::geom_errorbar(ggplot2::aes(alpha = .data[["alpha"]])) +
        ggplot2::geom_point(ggplot2::aes(alpha = .data[["alpha"]], colour = .data[["type"]]),cex = 1.5, data = df_cis) +
        # geom_point(cex = 1.5, aes(y = edge_pruned)) +
        # geom_point(cex = 0.5*1.5, aes(y = edge_pruned), colour = "white") +
        ggplot2::geom_text(ggplot2::aes(y = .data[["labstart"]], label = .data[["edge"]], hjust = .data[["hjust"]]), colour = rgb(0.2,0.2,0.2), cex = 2.5) +
        ggplot2::coord_flip() + 
        # facet_grid( ~  model) +
        ggplot2::scale_y_continuous(breaks = seq(floor(min(df_cis$lower)),ceiling(max(df_cis$upper)),by=major_break),
                                    minor_breaks = seq(floor(min(df_cis$lower)),ceiling(max(df_cis$upper)),by=minor_break)) + 
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::scale_color_manual("",values = c("black","darkred"), labels = c("Bootstrap mean","Sample")) +
        ggplot2::theme(legend.position = "top", # axis.text.y = element_text(size = 7),
                          axis.text.y=ggplot2::element_blank(),
                          axis.ticks.y=ggplot2::element_blank(),
                          panel.grid.major.y = ggplot2::element_blank()) + 
        ggplot2::xlab("") + 
        ggplot2::ylab("")  + 
        ggplot2::scale_alpha_continuous(guide="none", range = c(0,1)) + 
        ggplot2::theme(legend.position = "top") 
      
      
      
      if (prop0){
        
        plots <- plots + ggplot2::geom_label(ggplot2::aes(y=0,label=format(round(prop0, 2), nsmall = 2)), cex = prop0_cex * 2, 
                                             label.padding = ggplot2::unit(0.1, "lines"),
                                             label.size = 0.1, alpha = prop0_alpha,
                                             colour = "black")
      }
      # Return plot:
      return(plots)
      
    }    
    
    
  }, labels=labels,labels2=labels2,labelstart=labelstart,split0=split0,prop0=prop0)
  
  
  # If only one plot, plot it:
  if (length(matrices)==1){
    if (print) print(plots[[1]])
    invisible(plots[[1]])
  } else {
    for (i in seq_along(plots)){
      plots[[i]] <- plots[[i]] + ggplot2::ggtitle(matrices[i])
    }
    return(plots)
  }
  
}