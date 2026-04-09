partial_resid_plot <- function(modobj, term, int.term = NULL, group.term = NULL, add.intercept = FALSE, 
                       xlabel = NULL, ylabel = NULL, plot.title = NULL,show.data=FALSE,newdata = NULL) {
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(mgcv)
  
  terms=c(term,int.term)
  # Identify model type and extract data
  if ("gam" %in% class(modobj)) {
    model_type <- "gam"
    df <- if (!is.null(newdata)) newdata else modobj$model
  } else if ("gamm" %in% class(modobj) | "list" %in% class(modobj)) {
    model_type <- "gamm"
    df <- if (!is.null(newdata)) newdata else modobj$gam$model  # Extract model from the 'gam' component
    modobj <- modobj$gam  # Set the modobj to the 'gam' component for plotting
  } else if ("lm" %in% class(modobj)) {
    model_type <- "lm"
    df <- if (!is.null(newdata)) newdata else model.frame(modobj)
  } else if ("lmerMod" %in% class(modobj)) {
    model_type <- "lmer"
    df <- if (!is.null(newdata)) newdata else modobj@frame  # Use @frame for lmerMod objects
  } else {
    stop("Unsupported model type. Use gam, gamm, lm, or lmer.")
  }
  
  # Extract intercept if present
  mod.intercept <- if (add.intercept) {
    if (!"(Intercept)" %in% names(fixef(modobj))) {
      warning("Intercept not found in model coefficients; setting to 0.")
      0
    } else {
      fixef(modobj)["(Intercept)"]
    }
  } else {
    0
  }
  
  # Handling gam or gamm models
  if (model_type == "gam" || model_type == "gamm") {
    pterms <- predict(modobj, type = "terms", se.fit = TRUE,unconditional=TRUE)

    # Handle interaction terms if specified
    if (!is.null(int.term)) {
      intcol <- df %>% select(int.term)
      int_type <- ifelse(any(class(intcol[[1]]) == "numeric"), "numeric", "factor")
      if (int_type=="numeric") {
      
      }
    }
    
    # Adjust for intercept if needed
    if (add.intercept) {
      pterms.fit <- pterms$fit + mod.intercept
    } else {
      pterms.fit <- pterms$fit
    }
    
    pterms.sefit <- pterms$se.fit
    
    # Clean term names (e.g., removing "s(" and ")")
    colnames(pterms.fit) <- gsub(x = colnames(pterms.fit), pattern = "s\\(", replacement = "") %>%
      gsub(pattern = "\\)", replacement = "")
    colnames(pterms.sefit) <- gsub(x = colnames(pterms.sefit), pattern = "s\\(", replacement = "") %>%
      gsub(pattern = "\\)", replacement = "")
    
    # Create the pterms dataframe
    pterms.df <- data.frame(pterms.fit) %>%
      select(any_of(terms)) %>%
      mutate(fit = rowSums(across(where(is.numeric)))) %>%
      select(fit) %>%
      cbind(data.frame(pterms.sefit) %>%
              select(any_of(term)) %>%
              mutate(se.fit = rowSums(across(where(is.numeric))))) %>%
      mutate(upr = fit + 1.96 * se.fit,
             lwr = fit - 1.96 * se.fit) %>%
      mutate(partial.residuals = fit + resid(modobj)) %>%
      select(fit, se.fit, lwr, upr, partial.residuals)
    
    partial.residuals.df <- cbind(pterms.df, df %>%
                                    select(matches(term)) %>%
                                    rename_all(.funs = function(x) paste0(x, ".raw")))
    
    # Extract residuals and fits
    pterms <- partial.residuals.df
  } 
  
  # For lm models, extract terms and adjust for intercept if needed, and calculate residuals
  else if (model_type == "lm") {
    terms_pred <- predict(modobj, type = "terms", se.fit = TRUE)
    pterms.fit <- as.data.frame(terms_pred$fit)
    pterms.sefit <- as.data.frame(terms_pred$se.fit)
    
    if (add.intercept) {
      pterms.fit <- pterms + mod.intercept
    }
    
    # Compute residuals
    residuals <- resid(modobj)
    pterms.df <- data.frame(pterms.fit) %>%
      select(matches(terms)) %>%
      mutate(fit = rowSums(across(where(is.numeric)))) %>%
      select(fit) %>%
      cbind(data.frame(pterms.sefit) %>%
              select(matches(terms)) %>%
              mutate(se.fit = rowSums(across(where(is.numeric))))) %>%
      mutate(upr = fit + 1.96 * se.fit,
             lwr = fit - 1.96 * se.fit) %>%
      mutate(partial.residuals = fit + residuals) %>%
      select(fit, se.fit, lwr, upr, partial.residuals)
    pterms <- pterms.df
    
  } 
  
  # For lmer models, calculate residuals and combine with fixed-effects predictions for partial residuals
  else if (model_type == "lmer") {
    # Calculate fixed-effect contributions
    fixed_effects <- fixef(modobj)
    design_matrix <- model.matrix(terms(modobj), df)
    pterms <- as.data.frame(design_matrix %*% fixed_effects)  # Ensure pterms is a data frame
    
    # Add intercept if requested
    pterms <- pterms + if (add.intercept) mod.intercept else 0
    names(pterms) <- "fit"
    
    # Standard errors are not computed for lmer fixed effects in this implementation
    se_terms <- data.frame(rep(NA, nrow(df)))
    
    # Compute residuals
    residuals <- residuals(modobj)
    
    # Calculate partial residuals: fitted values + residuals
    pterms$partial.residuals <- pterms$fit + residuals
  } 
  
  else {
    stop("Unsupported model type.")
  }
  
  # Prepare plot data
  plot_data <- df %>%
    mutate(partial_residuals = pterms$partial.residuals,
           fit = pterms$fit, 
           lwr = if (model_type == "lmer") NA else pterms$lwr,  # Set NA for lwr and upr if model is lmer
           upr = if (model_type == "lmer") NA else pterms$upr) %>%
    select(all_of(c(term, int.term, group.term)), partial_residuals, fit, lwr, upr)
  
  # Start plotting
  gg <- ggplot(plot_data, aes(x = .data[[term]], y = partial_residuals)) +
    theme(
      axis.text.x = element_blank(),   # remove x-axis labels
      axis.ticks.x = element_blank(),  # remove x-axis ticks
      axis.title.x = element_blank()   # remove x-axis title
    )
  
  # Add points if desired
  if (show.data) {
    # Add points, colored by int.term if specified
    if (!is.null(int.term)) {
      gg <- gg + geom_point(aes(color = .data[[int.term]], fill = .data[[int.term]]), 
                            alpha = 0.75, shape = 21)
      # Plot connected lines (if group.term is specified)
      if (!is.null(group.term)) {
        gg <- gg + geom_line(aes(group = .data[[group.term]], color = .data[[int.term]]), 
                             alpha = 0.6, size = 1)  # Connect points by group.term with color reflecting int.term
      }
    } else {
      gg <- gg + geom_point(alpha = 0.75, color = "black", shape = 21,fill="black")  # Default black points if no int.term
      # Plot connected lines (if group.term is specified)
      if (!is.null(group.term)) {
        gg <- gg + geom_line(aes(group = .data[[group.term]]), 
                             alpha = 0.6, size = 1,color="black")  # Connect points by group.term with color reflecting int.term
      }
    }
    
  }
  
  
  # Skip the ribbon if model is 'lmer' because standard errors are not available
  if (model_type != "lmer") {
    # Add confidence ribbons, reflecting int.term color if specified
    if (!is.null(int.term)) {
      gg <- gg + geom_ribbon(aes(y = fit, ymin = lwr, ymax = upr, fill = .data[[int.term]]), 
                             alpha = 0.25, color = NA)  # Ribbon color reflects int.term
    } else {
      gg <- gg + geom_ribbon(aes(y = fit, ymin = lwr, ymax = upr), 
                             fill = "gray", alpha = 0.4, color = NA)
    }
  }
  
  # Add the fit lines, colored by int.term if specified
  if (!is.null(int.term)) {
    gg <- gg + geom_line(aes(y = fit, color = .data[[int.term]]), 
                         size = 1.5, alpha = 0.7)
  } else {
    gg <- gg + geom_line(aes(y = fit), color = "black", size = 1.5, alpha = 1)
    
  }
  
  # Customize labels
  if (!is.null(xlabel)) gg <- gg + xlab(xlabel)
  if (!is.null(ylabel)) gg <- gg + ylab(ylabel)
  if (!is.null(plot.title)) gg <- gg + ggtitle(plot.title)
  
  ## ============ PART 2: Derivative plot ============
  derv <- derivatives(modobj, term = sprintf('s(%s)', term), partial_match = TRUE)
  derv <- derv %>%
    mutate(sig = !(0 > .lower_ci & 0 < .upper_ci)) %>%
    mutate(sig_deriv = ifelse(sig == TRUE, .derivative, 0))
  dplot <- ggplot(data = derv) +
    geom_tile(aes(x = !!sym(term), y = 0.5),
              fill = NA,                  # transparent fill
              color = "black",            # black border
              size = 2,                 # border thickness
              height = 1) +               # height
    geom_tile(aes(x = !!sym(term), y = .5, fill = sig_deriv),
              show.legend = FALSE) +
    theme(panel.spacing = unit(-.01, "cm")) +
    scale_fill_gradient2(low = "#8c96c6", midpoint = 0, mid = "white",
                         high = "#440154FF",
                         limits = c(min(derv$sig_deriv), max(derv$sig_deriv))) +
    labs(x = xlabel, fill = sprintf("\u0394%s", term)) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 14),
          axis.line = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_line(),
          strip.background = element_blank(),
          strip.text.y = element_blank(),
          legend.position = "bottom",
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  ## ============ Combine both ============
  combined_plot <- gg / dplot + plot_layout(heights = c(3, 0.5))
  
  # ============ SIGNIFICANT RANGE REPORT ============
  # Handle factor-level smooths if present
  if (".smooth" %in% names(derv)) {
    derv <- derv %>%
      mutate(fac_levels = ifelse(str_detect(.smooth, ":"),
                                 sub("^.*:", "", .smooth),
                                 "all"))
  } else {
    derv$fac_levels <- "all"
  }
  
  smooth_var <- term
  factor_levels <- unique(derv$fac_levels)
  
  for (this_level in factor_levels) {
    this_derv <- derv %>% filter(fac_levels == this_level)
    sig_df <- this_derv %>% filter(sig == TRUE)
    cat(sprintf("\nSig change in %s:\n", this_level))
    
    if (nrow(sig_df) == 0) {
      cat("No change\n")
      next
    }
    
    sig_df <- sig_df %>%
      mutate(is_pos = .derivative > 0,
             blocks = data.table::rleid(is_pos))
    
    sig_ranges <- sig_df %>%
      group_by(blocks) %>%
      summarize(
        min = min(.data[[smooth_var]], na.rm = TRUE),
        max = max(.data[[smooth_var]], na.rm = TRUE),
        positive = first(is_pos),
        .groups = "drop"
      ) %>%
      arrange(min)
    
    for (r in seq_len(nrow(sig_ranges))) {
      if (isTRUE(sig_ranges$positive[r])) {
        cat(sprintf("Increasing from %1.3f to %1.3f\n",
                    sig_ranges$min[r], sig_ranges$max[r]))
      } else {
        cat(sprintf("Decreasing from %1.3f to %1.3f\n",
                    sig_ranges$min[r], sig_ranges$max[r]))
      }
    }
  }
  
  return(combined_plot)
  
}


