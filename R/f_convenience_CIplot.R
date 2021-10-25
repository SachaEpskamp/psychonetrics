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
  minor_break = 0.1
){

  logs <- sapply(x@log,function(x)x@event)
  if (any(grepl("Pruned",logs)) | any(grepl("step-up",logs))){
    warning("A model search algorithm was used in creating this model. The CIs are likely invalid.")
  }
  
  if (missing(labels2) && !missing(labels)){
    labels2 <- labels
  }

  
  # Stop if model is not psychonetrics:
  stopifnot(is(x,"psychonetrics"))
  
  # Default for matrices argument:
  if (missing(matrices)){
    matrices <- "none"
    
    # Varcov default:
    if (x@model == "varcov"){
      matrices <- covchooser(x@submodel)
    }  else if (x@model == "meta_varcov"){
      matrices <- covchooser(x@types$y,"y")
    }  else if (x@model == "lvm"){
      
      if (x@submodel == "rnm"){
        
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
      
    } else if (x@model %in% c("ml_lvm","dlvm1")){
      matrices <- c(
        "beta",
        covchooser(x@types$within_latent,"zeta_within"),
        covchooser(x@types$within_residual,"zeta_between")
      )
      
    } else if (x@model == "tsdlvm1"){
      matrices <- c(
        "beta",
        covchooser(x@types$zeta,"zeta")
      )
      
      
    }  else if (x@model == "Ising"){
      matrices <- c("omega")
      
    }  else stop("No default argument for 'matrices' for current model.")
  }
  if (any(matrices == "none")){
    stop("No default matrix could be selected")
  }
  

  # Extract parameters:
  pars <- x@parameters
  
  # Make a data frame for each matrix:
  plots <- lapply(matrices,function(mat,labels,labels2,labelstart){
    
    # Type of edge:
    edge <- "--"
    if (grepl("lambda",mat)){
      edge <- "->"
    } else if (grepl("beta",mat)){
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
    ncol   <- x@matrices$ncol[x@matrices$name == mat]
    nrow   <- x@matrices$nrow[x@matrices$name == mat]
    
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
        node2 <- labels[pars$var1_id[pars$matrix == mat]]
      } else {
        if (length(labels2) != ncol){
          stop("'labels2' is not of the correct length")
        }
        node2 <- labels2[pars$var2_id[pars$matrix == mat]]
      } 
    }

    # Edges :
    edges <- paste0(node1," ",edge," ",node2)
    
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
                             ggplot2::aes_string(x = "edge", y = "est", colour = "sig",
                        ymin = "lower", ymax = "upper")) +
      ggplot2::geom_hline(yintercept = 0, alpha = 0.2) +
      ggplot2::geom_errorbar() + 
      ggplot2::geom_point(cex = 1.5) +
      # geom_point(cex = 1.5, aes(y = edge_pruned)) +
      # geom_point(cex = 0.5*1.5, aes(y = edge_pruned), colour = "white") +
      ggplot2::geom_text(ggplot2::aes_string(y = "labstart", label = "edge", hjust = "hjust"), colour = rgb(0.2,0.2,0.2), cex = 2.5) +
      ggplot2::coord_flip() + 
      # facet_grid( ~  model) +
      ggplot2::scale_y_continuous(breaks = seq(floor(min(df_cis$lower)),ceiling(max(df_cis$lower)),by=major_break),
                         minor_breaks = seq(floor(min(df_cis$lower)),ceiling(max(df_cis$lower)),by=minor_break)) + 
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::scale_colour_manual("",values = c("black",colors),drop=FALSE) +  ggplot2::theme(legend.position = "top", # axis.text.y = element_text(size = 7),
                                                                                              axis.text.y=ggplot2::element_blank(),
                                                                                              axis.ticks.y=ggplot2::element_blank(),
                                                                                              panel.grid.major.y = ggplot2::element_blank()) + 
      ggplot2::xlab("") + 
      ggplot2::ylab("") 
    
    # Return plot:
    return(plots)
    
  }, labels=labels,labels2=labels2,labelstart=labelstart)
  
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