# Inner function:
# von Oertzen, T., Schmiedek, F., & Voelkle, M. C. (2020). Ergodic Subspace Analysis. Journal of Intelligence, 8(1), 3.
esa_manual <- function(sigma_wp, sigma_bp, cutoff = 0.1){
  # Step 1: Combine both matrices:
  Sigma <- 0.5 * (sigma_wp + sigma_bp) # (8)
  
  # Step 2: whitening. First an eigenvalue decomposition:
  eig <- eigen(Sigma)
  Q1 <- t(eig$vectors)
  D1 <- diag(eig$values)
  
  # Whiten:
  Qwhite <- diag(1/sqrt(eig$values)) %*% Q1 # (10)
  
  # Step 3: eigen decomposition of between covmatrix:
  eig2 <- eigen((Qwhite %*% (0.5 * sigma_bp) %*% t(Qwhite)))
  Q2 <- t(eig2$vectors)
  D_BP <- diag(eig2$values)
  
  # Total transformation of ESA transformation:
  Qesa <- Q2 %*% Qwhite
  
  # So the BP lambdas are:
  lambda_BP <- diag(D_BP)
  # And the WP lambdas:
  lambda_WP <- 1 - diag(D_BP)
  
  # Wich makes the ergodicity:
  erg <- lambda_BP - lambda_WP
  
  # Step 4: select an arbitrary c:
  c <- cutoff
  
  # Now obtain the subspaces:
  V_BP <- Qesa[erg > c, , drop=FALSE]
  V_ergodic <- Qesa[erg >= -c & erg <= c, , drop=FALSE]
  V_WP <- Qesa[erg < -c, , drop=FALSE]
  
  # Now return a full list:
  res <- list(
    ergodicity = erg,
    Q_esa = Qesa,
    V_bp = V_BP,
    V_ergodic = V_ergodic,
    V_wp = V_WP,
    cutoff = c
  )
  class(res) <- c("esa_manual","list")
  return(res)
}

# Ergodic subspace analysis
# von Oertzen, T., Schmiedek, F., & Voelkle, M. C. (2020). Ergodic Subspace Analysis. Journal of Intelligence, 8(1), 3.
esa <- function(x, cutoff = 0.1, between = c("crosssection","between")){
  between <- match.arg(between)
  if (x@model != "dlvm1"){
    stop("Only implemented for 'dlvm1' model family")
  }
  
  # Number of groups:
  nGroup <- nrow(x@sample@groups)
  
  
  ### Observed variables:
  obs <- lapply(1:nGroup,function(g){
    # Within matrix:
    sigma_wp <- x@modelmatrices[[g]]$sigma_within
    
    # Between matrix:
    if (between == "crosssection"){
      sigma_bp <- x@modelmatrices[[g]]$sigma_crosssection
    } else {
      sigma_bp <- x@modelmatrices[[g]]$sigma_between
    }
    
    
    # Result:
    esa_manual(sigma_wp, sigma_bp, cutoff = cutoff)
  })
  names(obs) <- x@sample@groups$label
  
  ### Latent variables:
  lat <- lapply(1:nGroup,function(g){
    # Within matrix:
    sigma_wp <- x@modelmatrices[[g]]$sigma_eta_within
    
    # Between matrix:
    sigma_bp <- x@modelmatrices[[g]]$sigma_zeta_between
    
    # Result:
    esa_manual(sigma_wp, sigma_bp, cutoff = cutoff)
  })
  names(lat) <- x@sample@groups$label
  
  # Final result:
  res <- list(
    observed = obs,
    latent = lat
  )
  class(res) <- c("esa", "list")
  return(res)
}

# Print function:
print.esa_manual <- function(x, printref = TRUE, ...){
  cat("Ergodicity coefficients:\n")
  cat(round(x$ergodicity,2))
  cat("\n\n")
  cat("Between persons subspace:\n")
  print(round(x$V_bp,2))
  cat("\n")
  cat("Ergodic subspace:\n")
  print(round(x$V_ergodic,2))
  cat("\n")
  cat("Within persons subspace:\n")
  print(round(x$V_wp,2))
  cat("\n")
  if (printref){
    cat("More information: von Oertzen, T., Schmiedek, F., & Voelkle, M. C. (2020). Ergodic Subspace Analysis. Journal of Intelligence, 8(1), 3.")
  }
}

print.esa <- function(x, printref = TRUE, ...){
 
  nGroups <- length(x$observed)
  
  cat("ESA results for observed variables:\n\n")
  
  for (g in  1:nGroups){
    # if (nGroups > 1){
      cat("\t- Group:",names(x$observed)[g],"\n\n")      
    # }
    
    print(x$observed[[g]], printref = FALSE)
    
  }
  
  
  cat("\n\nESA results for latent variables:\n\n")
  
  for (g in  1:nGroups){
    # if (nGroups > 1){
    cat("\t- Group:",names(x$latent)[g],"\n\n")      
    # }
    
    print(x$latent[[g]], printref = FALSE)
    
  }
  
  if (printref){
    cat("More information: von Oertzen, T., Schmiedek, F., & Voelkle, M. C. (2020). Ergodic Subspace Analysis. Journal of Intelligence, 8(1), 3.")
  }

}

# Plot method:
plot.esa_manual <- function(x,...){
  # Scree plots:
  df <- data.frame(
    ev = seq_along(x$ergodicity),
    ergodicity = x$ergodicity
  )
  
  
  g <- ggplot2::ggplot(df, ggplot2::aes_string(x = "ev", y = "ergodicity")) + 
    ggplot2::geom_line(lwd = 1.5) + ggplot2::geom_point(cex = 3) + ggplot2::ylim(-1,1) + 
    ggplot2::geom_text(aes(x=mean(seq_along(x$ergodicity)), y = 1, label = "Dominantly between-subject"), colour = "black") + 
    ggplot2::geom_text(aes(x=mean(seq_along(x$ergodicity)), y = -1, label = "Dominantly Within-subject"), colour = "black") + 
    ggplot2::geom_text(aes(x=1.5, y = -0, label = "Ergodic"), colour = "black") + 
    ggplot2::theme_bw() + 
    ggplot2::ylab("") + ggplot2::xlab("Component") +
    ggplot2::geom_hline(yintercept = -0.1, lwd = 0.5) + 
    ggplot2::geom_hline(yintercept = 0.1, lwd = 0.5) +
    ggplot2::theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
    ) + 
    ggplot2::scale_color_discrete("")
  
  return(g) 
}

plot.esa <- function(x, plot = c("observed","latent"),...){
  plot <- match.arg(plot)
  nGroups <- length(x$observed)
  
  # Create data frame:
  df <- do.call(rbind, lapply(1:nGroups,function(g){
    data.frame(
      ev = seq_along(x[[plot]][[g]]$ergodicity),
      ergodicity = x[[plot]][[g]]$ergodicity,
      group = names(x[[plot]])[g]
    ) }))
  
  erg <- df$ergodicity
    
    # Create base plot:
    if (nGroups == 1){
      g <- ggplot2::ggplot(df, aes_string(x = "ev", y = "ergodicity"))
    } else {
      g <- ggplot2::ggplot(df, aes_string(x = "ev", y = "ergodicity", colour = "factor(group)"))
    }
    
    g <- g + 
      ggplot2::geom_line(lwd = 1.5) + ggplot2::geom_point(cex = 3) + ggplot2::ylim(-1,1) + 
      ggplot2::geom_text(aes(x=mean(seq_along(erg)), y = 1, label = "Dominantly between-subject"), colour = "black") + 
      ggplot2::geom_text(aes(x=mean(seq_along(erg)), y = -1, label = "Dominantly Within-subject"), colour = "black") + 
      ggplot2::geom_text(aes(x=1.5, y = -0, label = "Ergodic"), colour = "black") + 
      ggplot2::theme_bw() + 
      ggplot2::ylab("") + ggplot2::xlab("Component") +
      ggplot2::geom_hline(yintercept = -0.1, lwd = 0.5) + 
      ggplot2::geom_hline(yintercept = 0.1, lwd = 0.5) +
      ggplot2::theme(
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
      ) + 
      ggplot2::scale_color_discrete("")
    
    return(g) 
}
