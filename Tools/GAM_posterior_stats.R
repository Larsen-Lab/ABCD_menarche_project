#POSTERIOR DISTRIBUTION SMOOTH FUNCTIONS
##Function to simulate the posterior distribution from a fitted GAM, calculate smooths for individual posterior draws, and return smooth max and min values + 95% credible intervals

gam.posterior.smooths2 <- function(gam.model,
                                   smooth_var,
                                   factor_var = NULL,
                                   set_fx = FALSE,
                                   draws = 1000,
                                   increments = 200,
                                   return_draws = TRUE,
                                   UNCONDITIONAL = FALSE,
                                   make_plots = FALSE,
                                   newdata = NULL,
                                   seed = 42) {
  
  # --- package checks
  required_pkgs <- c("gratia", "dplyr", "tidyr", "ggplot2")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop("Please install these packages first: ", paste(missing_pkgs, collapse = ", "))
  }
  
  # Make sure we have a gam or gamm
  if (any(class(gam.model) == "gam")) {
    gam_obj <- gam.model
  } else if (!is.null(gam.model$gam) && any(class(gam.model$gam) == "gam")) {
    gam_obj <- gam.model$gam
  } else {
    stop("Can't find a 'gam' object. Pass an mgcv::gam or an object containing $gam")
  }
  
  # reproducible random draws
  set.seed(seed)
  
  # extract model frame & term labels
  df <- gam_obj$model
  all_terms <- attr(terms(df), "term.labels")
  
  # Make sure smooth_var exists in model frame
  if (!(smooth_var %in% names(df))) {
    stop("smooth_var '", smooth_var, "' not found in the model frame. Check the variable name.")
  }
  
  # generate a prediction dataframe if newdata not provided
  if (is.null(newdata)) {
    # safer range extraction
    x_range <- range(df[[smooth_var]], na.rm = TRUE)
    if (!all(is.finite(x_range))) {
      stop("smooth_var '", smooth_var, "' has no finite values in the model data (all NA / non-finite).")
    }
    
    # handle constant variable by adding a tiny jitter range
    if (diff(x_range) == 0) {
      EPS <- ifelse(x_range[1] == 0, 1e-6, abs(x_range[1]) * 1e-6)
      x_min <- x_range[1] - EPS
      x_max <- x_range[2] + EPS
      warning("smooth_var '", smooth_var, "' is constant in the model data — creating a tiny range for prediction.")
    } else {
      x_min <- x_range[1]
      x_max <- x_range[2]
    }
    
    xseq <- seq(from = x_min, to = x_max, length.out = increments)
    
    if (is.null(factor_var)) {
      newdata <- data.frame(tmp_x = xseq, stringsAsFactors = FALSE)
      names(newdata) <- smooth_var
    } else {
      if (!(factor_var %in% names(df))) stop("factor_var '", factor_var, "' not found in model frame.")
      flevs <- unique(df[[factor_var]])
      newdata <- expand.grid(tmp_x = xseq, f = flevs, stringsAsFactors = FALSE)
      names(newdata) <- c(smooth_var, factor_var)
    }
    
    # set other covariates to typical values (mean or modal) to avoid missing columns required for predict
    covariates <- setdiff(all_terms, c(smooth_var, factor_var))
    if (length(covariates) > 0) {
      cov_vals <- lapply(df[covariates], function(col) {
        if (is.numeric(col)) {
          mean(col, na.rm = TRUE)
        } else {
          # pick the most frequent level (works for factor/character)
          as.character(names(sort(table(col), decreasing = TRUE))[1])
        }
      })
      for (cov in names(cov_vals)) newdata[[cov]] <- cov_vals[[cov]]
    }
  } else {
    # if user supplied newdata, verify necessary columns exist
    if (!(smooth_var %in% names(newdata))) stop("Provided newdata must contain column '", smooth_var, "'.")
    if (!is.null(factor_var) && !(factor_var %in% names(newdata))) stop("Provided newdata must contain column '", factor_var, "'.")
    # ensure order by factor then smooth var so groupwise operations are contiguous
    if (is.null(factor_var)) {
      newdata <- newdata[order(newdata[[smooth_var]]), , drop = FALSE]
    } else {
      newdata <- newdata[order(newdata[[factor_var]], newdata[[smooth_var]]), , drop = FALSE]
    }
    rownames(newdata) <- NULL
  }
  
  # --- identify relevant smooth label via gratia::smooths()
  all_smooths <- gratia::smooths(gam_obj)
  # try to match the smooth that contains the smooth_var name
  matched <- all_smooths[grepl(smooth_var, all_smooths, fixed = TRUE)]
  if (length(matched) == 0) {
    # fallback: try partial grep without fixed
    matched <- all_smooths[grepl(smooth_var, all_smooths)]
  }
  if (length(matched) == 0) {
    stop("Couldn't find a smooth matching '", smooth_var, "' in the model's smooth labels. Available smooths: ",
         paste(all_smooths, collapse = ", "))
  }
  # if multiple matches, pick first but warn
  if (length(matched) > 1) {
    warning("Multiple smooths matched the string '", smooth_var, "'. Using the first match: '", matched[1], "'.")
  }
  sm_label <- matched[1]
  
  # --- Posterior smooth draws via gratia::smooth_samples
  smooth_samples <- tryCatch({
    gratia::smooth_samples(model = gam_obj,
                           select = sm_label,
                           data = newdata,
                           n = draws,
                           unconditional = UNCONDITIONAL)
  }, error = function(e) {
    stop("gratia::smooth_samples failed: ", conditionMessage(e))
  })
  
  # standardize names: gratia output typically contains 'draw', smooth_var, and 'value'
  smooth_draws <- dplyr::rename(smooth_samples, predicted = .value)
  
  # --- Posterior derivatives via gratia::derivatives
  deriv_samples <- tryCatch({
    # derivatives(..., interval="none") returns draws if n specified
    gratia::derivatives(object = gam_obj,
                        select = sm_label,
                        data = newdata,
                        n = draws,
                        unconditional = UNCONDITIONAL)
  }, error = function(e) {
    stop("gratia::derivatives failed: ", conditionMessage(e))
  })
  
  # ensure a consistent naming convention for derivative column
  if ("derivative" %in% names(deriv_samples)) {
    deriv_samples <- dplyr::rename(deriv_samples, deriv = derivative)
  } else if ("estimate" %in% names(deriv_samples)) {
    deriv_samples <- dplyr::rename(deriv_samples, deriv = estimate)
  }
  
  # --- compute extrema (x value at which predicted y is max/min per draw)
  # smooth_draws: columns draw, <smooth_var>, predicted (and possibly factor_var)
  if (!is.null(factor_var)) {
    extrema <- smooth_draws %>%
      dplyr::group_by(.draw, dplyr::across(dplyr::all_of(factor_var))) %>%
      dplyr::summarise(
        x_at_max = .data[[smooth_var]][which.max(predicted)],
        x_at_min = .data[[smooth_var]][which.min(predicted)],
        .groups = "drop"
      )
  } else {
    extrema <- smooth_draws %>%
      dplyr::group_by(.draw) %>%
      dplyr::summarise(
        x_at_max = .data[[smooth_var]][which.max(predicted)],
        x_at_min = .data[[smooth_var]][which.min(predicted)],
        .groups = "drop"
      )
  }
  
  # summarise extrema across draws to median + 95% CI
  if (!is.null(factor_var)) {
    extrema_summary <- extrema %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(factor_var))) %>%
      dplyr::summarise(
        max.median = median(x_at_max, na.rm = TRUE),
        max.lower = stats::quantile(x_at_max, 0.025, na.rm = TRUE),
        max.upper = stats::quantile(x_at_max, 0.975, na.rm = TRUE),
        min.median = median(x_at_min, na.rm = TRUE),
        min.lower = stats::quantile(x_at_min, 0.025, na.rm = TRUE),
        min.upper = stats::quantile(x_at_min, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    extrema_summary <- data.frame(
      max.median = median(extrema$x_at_max, na.rm = TRUE),
      max.lower = as.numeric(stats::quantile(extrema$x_at_max, 0.025, na.rm = TRUE)),
      max.upper = as.numeric(stats::quantile(extrema$x_at_max, 0.975, na.rm = TRUE)),
      min.median = median(extrema$x_at_min, na.rm = TRUE),
      min.lower = as.numeric(stats::quantile(extrema$x_at_min, 0.025, na.rm = TRUE)),
      min.upper = as.numeric(stats::quantile(extrema$x_at_min, 0.975, na.rm = TRUE))
    )
  }
  
  # --- optional plotting: histograms of extrema per factor level
  if (isTRUE(make_plots)) {
    if (!is.null(factor_var)) {
      print(ggplot2::ggplot(extrema, ggplot2::aes_string(x = "x_at_max", fill = factor_var)) +
              ggplot2::geom_histogram(position = "identity", alpha = 0.5) +
              ggplot2::ggtitle("Posterior maxima of smooth"))
    } else {
      print(ggplot2::ggplot(extrema, ggplot2::aes(x = x_at_max)) +
              ggplot2::geom_histogram() + ggplot2::ggtitle("Posterior maxima of smooth"))
    }
  }
  
  # --- assemble output
  out <- list(
    smooth_draws = smooth_draws,
    deriv_draws = deriv_samples,
    extrema = extrema,
    extrema_summary = extrema_summary,
    newdata = newdata,
    smooth_label = sm_label
  )
  
  if (isTRUE(return_draws)) return(out)
  return(extrema_summary)
}


# Function to compute KDE mode to be used for derivative
get_mode_kde <- function(x, ...) {
  d <- density(x, ...)
  d$x[which.max(d$y)]
}

derivative_stats <- function(gam.model, smooth_var, factor_var = NULL, set_fx = FALSE, EPS = .01, draws = 1000, increments = 200, return_draws = FALSE, UNCONDITIONAL = FALSE, make_plots = FALSE, newdata = NULL) {
  libararies <- c("gratia", "dplyr", "ggplot2", "tidyr", "purrr", "rlang","HDInterval")
  # Load the libraries
  invisible(lapply(libararies, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "needed but not installed. Please install it."))
    }
    library(pkg, character.only = TRUE)
  }))
  # Check the model
  if (any(class(gam.model) == "gam")) {
    gam.model <- gam.model
  } else if (class(gam.model$gam) == "gam") {
    gam.model <- gam.model$gam
  } else {
    stop("Can't find a gam object")
  }
  
  
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  all_terms <- attr(terms(df),"term.labels")
  
  if (is.null(newdata)) {
    if (is.null(factor_var)) {
      newdata <- with(df,expand.grid(x=seq(from=quantile(df[,smooth_var],.05),to=quantile(df[,smooth_var],.95),length.out=increments)))
      colnames(newdata)=smooth_var
    } else {
      newdata <- with(df,expand.grid(x=seq(from=min(df[,smooth_var]),to=max(df[,smooth_var]),length.out=increments),f=unique(df[,factor_var])))
      colnames(newdata)=c(smooth_var,factor_var)
    }
    covariates <- all_terms[!(all_terms%in%c(smooth_var,factor_var))]
    if (!is_empty(covariates)) {
      covariate_values <- df %>% dplyr::summarize(across(.cols = all_of(covariates),.fns = function(x){ifelse(is.numeric(x),yes = mean(x,na.rm=T),no = levels(x)[2])}))
      for (cov in covariates) {
        newdata[,cov]=covariate_values[,cov]
      }
    }
  }
  
  fd_samp <- derivative_samples(gam.model,data = newdata, type = "central",focal = smooth_var, eps = EPS, seed = 12, n_sim = draws)
  sd_samp <- derivative_samples(gam.model,data = newdata, order = 2, type = "central",focal = smooth_var, eps = EPS, seed = 12, n_sim = draws)
  
  # --- Calculate Maximum Curvature (k) ---
  # We need to join fd and sd samples by draw and x location
  # Note: We rely on the fact that both samples used the same 'data' grid and seed

  # Join and calculate k
  k_samp <- left_join(
    fd_samp %>% select(.draw, .row, .derivative) %>% rename(fd = .derivative),
    sd_samp %>% select(.draw, .row, .derivative) %>% rename(sd = .derivative),
    by = c(".draw", ".row")
  ) %>%
    # Bring back the x-values (smooth_var) from the original data using .row
    mutate(!!sym(smooth_var) := newdata[[smooth_var]][.row]) %>%
    mutate(k = abs(sd) / (1 + fd^2)^1.5)

  # Get the value of smooth var at the Max Curvature
  # We only care about MAX curvature (sharpest turn)
  k.extrema <- k_samp %>%
    group_by(.draw) %>%
    dplyr::summarize(
      max.ind = which.max(k),
      max.val = .data[[smooth_var]][max.ind],
      .groups = "drop"
    )

  # Summarize K
  summary.k.max <- tibble(
    mode_max = get_mode_kde(k.extrema$max.val),
    lower_ci = hdi(k.extrema$max.val)[1],
    upper_ci = hdi(k.extrema$max.val)[2],
    median_max = median(k.extrema$max.val)
  )
  
  # Calculate True (Point Estimate) Curvature for reference
  # We can reuse the true.fd and true.sd calculated later
  # (Might need to move the true.fd/true.sd calculation up higher in the script to use here, 
  #  or just calculate it temporarily here)
  
  # Get the value of smooth var at the fd (first derivative extremes)
  fd.extrema <- fd_samp %>%
    group_by(.draw) %>%
    dplyr::summarize(
      max.ind = which.max(.derivative),
      max.val = .data[[smooth_var]][max.ind],
      min.ind = which.min(.derivative),
      min.val = .data[[smooth_var]][min.ind],
      .groups = "drop"
    )
  
  # dplyr::summarize the values across the draws
  # FD summary
  summary.fd.max <- tibble(
    mode_max = get_mode_kde(fd.extrema$max.val),
    lower_ci = hdi(fd.extrema$max.val)[1],
    upper_ci = hdi(fd.extrema$max.val)[2],
    median_max = median(fd.extrema$max.val)
  )
  
  summary.fd.min <- tibble(
    mode_min = get_mode_kde(fd.extrema$min.val),
    lower_ci = hdi(fd.extrema$min.val)[1],
    upper_ci = hdi(fd.extrema$min.val)[2],
    median_min = median(fd.extrema$min.val)
  )
  
  # OLD VERSION USING QUANTILES -- we now use the HDInterval package for more robust credible intervals
  # Get the median and 95% credible interval of the smooth var at the peak of the fd
  # summary.fd.peak <- fd.max %>%
  #   dplyr::summarize(
  #     median_peak = median(max.val),
  #     lower_ci = quantile(max.val, 0.025),
  #     upper_ci = quantile(max.val, 0.975)
  #   )
  

  # Get the value of smooth var at the extremes of the sd (second derivative)
  sd.extrema <- sd_samp %>%
    group_by(.draw) %>%
    dplyr::summarize(
      max.ind = which.max(.derivative),
      max.val = .data[[smooth_var]][max.ind],
      min.ind = which.min(.derivative),
      min.val = .data[[smooth_var]][min.ind],
      .groups = "drop"
    )
  
  # dplyr::summarize the values across the draws
  # SD summary
  summary.sd.max <- tibble(
    mode_max = get_mode_kde(sd.extrema$max.val),
    lower_ci = hdi(sd.extrema$max.val)[1],
    upper_ci = hdi(sd.extrema$max.val)[2],
    median_max = median(sd.extrema$max.val)
  )
  
  summary.sd.min <- tibble(
    mode_min = get_mode_kde(sd.extrema$min.val),
    lower_ci = hdi(sd.extrema$min.val)[1],
    upper_ci = hdi(sd.extrema$min.val)[2],
    median_min = median(sd.extrema$min.val)
  )
  
  
  # # OLD VERSION USING QUANTILES -- we now use the HDInterval package for more robust credible intervals
  # # Get the median and 95% credible interval of the smooth var at the peak of the sd
  # summary.sd.peak <- sd.max %>%
  #   dplyr::summarize(
  #     median_peak = median(max.val),
  #     lower_ci = quantile(max.val, 0.025),
  #     upper_ci = quantile(max.val, 0.975)
  #   )
  
  # First derivative: get both max and min (non-simulated)
  true.fd <- derivatives(gam.model, data = newdata, select = smooth_var, partial_match = T, type = "central", eps = EPS)
  
  true.fd.extrema <- true.fd %>%
    filter(str_detect(.smooth,smooth_var))%>%
    dplyr::summarize(
      max.ind = which.max(.derivative),
      max.val = .data[[smooth_var]][max.ind],
      min.ind = which.min(.derivative),
      min.val = .data[[smooth_var]][min.ind],
      .groups = "drop"
    )
  
  # Second derivative: get both max and min (non-simulated)
  true.sd <- derivatives(gam.model, data = newdata, order = 2, type = "central", select = smooth_var, partial_match = T, eps = EPS)
  true.sd.extrema <- true.sd %>%
    filter(str_detect(.smooth,smooth_var))%>%
    dplyr::summarize(
      max.ind = which.max(.derivative),
      max.val = .data[[smooth_var]][max.ind],
      min.ind = which.min(.derivative),
      min.val = .data[[smooth_var]][min.ind],
      .groups = "drop"
    )
  
  # Calculate True Curvature (Point Estimate)
  
  # We need to join the true FD and true SD data frames.
  true.k <- inner_join(
    true.fd %>% 
      filter(stringr::str_detect(.smooth, smooth_var)) %>%
      select(all_of(smooth_var), fd = .derivative),
    true.sd %>% 
      filter(stringr::str_detect(.smooth, smooth_var)) %>%
      select(all_of(smooth_var), sd = .derivative),
    by = smooth_var
  ) %>%
    mutate(k = abs(sd) / (1 + fd^2)^1.5)
  
  # Find the single x-value where true curvature is maximized
  true.k.extrema <- true.k %>%
    dplyr::slice_max(k, n = 1, with_ties = FALSE) %>%
    dplyr::rename(max.val = !!sym(smooth_var))
    
  # Output summary table
  deriv_summary <- tibble(
    fd_max_value     = true.fd.extrema$max.val,
    fd_max_median    = summary.fd.max$median_max,
    fd_max_mode      = summary.fd.max$mode_max,
    fd_max_lower_ci  = summary.fd.max$lower_ci,
    fd_max_upper_ci  = summary.fd.max$upper_ci,
    
    fd_min_value     = true.fd.extrema$min.val,
    fd_min_median    = summary.fd.min$median_min,
    fd_min_mode      = summary.fd.min$mode_min,
    fd_min_lower_ci  = summary.fd.min$lower_ci,
    fd_min_upper_ci  = summary.fd.min$upper_ci,
    
    sd_max_value     = true.sd.extrema$max.val,
    sd_max_median    = summary.sd.max$median_max,
    sd_max_mode      = summary.sd.max$mode_max,
    sd_max_lower_ci  = summary.sd.max$lower_ci,
    sd_max_upper_ci  = summary.sd.max$upper_ci,
    
    sd_min_value     = true.sd.extrema$min.val,
    sd_min_median    = summary.sd.min$median_min,
    sd_min_mode      = summary.sd.min$mode_min,
    sd_min_lower_ci  = summary.sd.min$lower_ci,
    sd_min_upper_ci  = summary.sd.min$upper_ci,
    
    # Curvature
    k_max_value      = true.k.extrema$max.val,
    k_max_median     = summary.k.max$median_max,
    k_max_mode       = summary.k.max$mode_max,
    k_max_lower_ci   = summary.k.max$lower_ci,
    k_max_upper_ci   = summary.k.max$upper_ci
  )
  
  # cli::cli_h1("Posterior Derivative Summary")
  # 
  # cli::cli_text("First derivative extrema (change rate peaks/troughs):")
  # cli::cli_text("  ▸ Max slope at {true.fd.extrema$max.val} [{round(summary.fd.max$lower_ci, 2)} – {round(summary.fd.max$upper_ci, 2)}]")
  # cli::cli_text("  ▸ Min slope at {true.fd.extrema$min.val} [{round(summary.fd.min$lower_ci, 2)} – {round(summary.fd.min$upper_ci, 2)}]")
  # 
  # cli::cli_text("")
  # cli::cli_text("Second derivative extrema (curvature):")
  # cli::cli_text("  ▸ Max curvature at {true.sd.extrema$max.val} [{round(summary.sd.max$lower_ci, 2)} – {round(summary.sd.max$upper_ci, 2)}]")
  # cli::cli_text("  ▸ Min curvature at {true.sd.extrema$min.val} [{round(summary.sd.min$lower_ci, 2)} – {round(summary.sd.min$upper_ci, 2)}]")
  
  # Posterior Derivative Summary (mimicking cli::cli_h1)
  cat("================================\n")
  cat("Posterior Derivative Summary\n")
  cat("================================\n\n")
  
  # First derivative extrema (change rate peaks/troughs):
  cat("First derivative extrema (change rate peaks/troughs):\n")
  
  # Format for the CI range: [XX.XX - YY.YY]
  ci_format <- "[%.2f – %.2f]"
  
  # Max slope output
  max_slope_ci <- sprintf(ci_format, 
                          summary.fd.max$lower_ci, 
                          summary.fd.max$upper_ci)
  cat(sprintf("  > Max slope at %s %s\n", 
              true.fd.extrema$max.val, 
              max_slope_ci))
  
  # Min slope output
  min_slope_ci <- sprintf(ci_format, 
                          summary.fd.min$lower_ci, 
                          summary.fd.min$upper_ci)
  cat(sprintf("  > Min slope at %s %s\n", 
              true.fd.extrema$min.val, 
              min_slope_ci))
  
  cat("\n")
  
  # Second derivative extrema (curvature):
  cat("Second derivative extrema (curvature):\n")
  
  # Max curvature output
  max_curvature_ci <- sprintf(ci_format, 
                              summary.sd.max$lower_ci, 
                              summary.sd.max$upper_ci)
  cat(sprintf("  > Max curvature at %s %s\n", 
              true.sd.extrema$max.val, 
              max_curvature_ci))
  
  # Min curvature output
  min_curvature_ci <- sprintf(ci_format, 
                              summary.sd.min$lower_ci, 
                              summary.sd.min$upper_ci)
  cat(sprintf("  > Min curvature at %s %s\n", 
              true.sd.extrema$min.val, 
              min_curvature_ci))
  

  if (make_plots) {
    library(ggridges)
    library(patchwork) 
    library(ggplot2)
    library(dplyr)
    
    # --- 1. PREPARE POINT ESTIMATES FOR PLOTTING ON CURVE ---
    # Create a dataframe of the specific X-locations (modes) for all metrics
    if (!rlang::is_empty(covariates)) {
      peak_pts_full <- newdata[rep(1, 5), , drop = FALSE] 
      peak_pts_full[[smooth_var]] <- c(summary.fd.max$mode_max, summary.fd.min$mode_min,
                                       summary.sd.max$mode_max, summary.sd.min$mode_min,
                                       summary.k.max$mode_max) 
    } else {
      peak_pts_full <- setNames(
        data.frame(c(summary.fd.max$mode_max, summary.fd.min$mode_min,
                     summary.sd.max$mode_max, summary.sd.min$mode_min,
                     summary.k.max$mode_max)), 
        smooth_var
      )
    }
    
    # Assign types
    peak_pts_full$type <- c("fd max", "fd min", "sd max", "sd min", "k max")
    
    # Get the predicted Y values (fitted smooth) at these peak X locations
    # We use partial_match = TRUE to handle complex model formulas
    peak_vals <- smooth_estimates(gam.model, data = peak_pts_full, select = smooth_var, partial_match = TRUE) %>%
      mutate(type = peak_pts_full$type)
    
    
    # --- 2. PREPARE POSTERIOR DRAWS FOR RIDGE PLOTS ---
    posterior_locs <- bind_rows(
      tibble(type = "fd max", value = fd.extrema$max.val),
      tibble(type = "fd min", value = fd.extrema$min.val),
      tibble(type = "sd max", value = sd.extrema$max.val),
      tibble(type = "sd min", value = sd.extrema$min.val),
      tibble(type = "k max",  value = k.extrema$max.val)
    )
    
    # Set Factor Order for plotting (Curvature at top, then SD, then FD)
    posterior_locs$type <- factor(posterior_locs$type, 
                                  levels = c("k max", "sd min", "sd max", "fd min", "fd max"))
    
    # Define Color Palette
    my_colors <- c(
      "fd max" = "#FF6666",  # Light Red
      "fd min" = "#6699FF",  # Light Blue
      "sd max" = "#CC0000",  # Dark Red
      "sd min" = "#003399",  # Dark Blue
      "k max"  = "#800080"   # Purple
    )
    
    # Get background smooth fit for plotting
    smooth.fit <- smooth_estimates(gam.model, data = newdata, select = smooth_var, partial_match = TRUE)
    
    # Define shared x-axis range to align plots
    xrange <- range(c(smooth.fit[[smooth_var]], posterior_locs$value), na.rm = TRUE)
    
    
    # --- 3. PLOT: SMOOTH CURVE WITH PEAK MARKERS ---
    smooth_plot <- draw(smooth.fit) +
      geom_point(data = peak_vals, aes(x = !!sym(smooth_var), y = .estimate, color = type), size = 3) +
      scale_color_manual(values = my_colors) +
      labs(title = "Smooth Fit with Transition Points", y = "Response") +
      coord_cartesian(xlim = xrange, expand = FALSE) +
      theme_classic() +
      theme(legend.position = "none") 
    
    
    # --- 4. PLOT: POSTERIOR RIDGE DISTRIBUTIONS ---
    posterior_ridges <- ggplot(posterior_locs, aes(x = value, y = type, fill = type)) +
      geom_density_ridges(alpha = 0.7, scale = 1.2, quantile_lines = TRUE, quantiles = 2) + 
      scale_fill_manual(values = my_colors) +
      labs(title = "Posterior Uncertainty of Transition Points", x = smooth_var, y = NULL) +
      coord_cartesian(xlim = xrange, expand = FALSE) +
      theme_classic() +
      theme(legend.position = "none")
    
    # Combine and Print Main Plot
    final_plot <- (posterior_ridges / smooth_plot) + 
      patchwork::plot_layout(heights = c(1, 1))
    print(final_plot)
    
    
    # --- 5. PLOT: SIGNIFICANT INTERVALS (FD) ---
    # Filter for significance
    true.fd <- true.fd %>%
      filter(stringr::str_detect(.smooth, smooth_var)) %>%
      mutate(significant = ifelse(.lower_ci > 0 | .upper_ci < 0, "yes", "no")) %>%
      mutate(significant.deriv = ifelse(significant=="yes", .derivative, NA))
    
    fd_ci_plot <- ggplot(true.fd, aes(x = !!sym(smooth_var), y = 1, fill = significant.deriv)) +
      geom_tile(height = 0.2) +
      scale_fill_gradient2(low = "#6699FF", mid = "white", high = "#FF6666", na.value = "grey95", name = "Slope") +
      # Add Curvature Point (Purple Diamond) for reference
      geom_point(data = peak_vals %>% filter(type == "k max"), 
                 aes(x = !!sym(smooth_var), y = 1), color = "#800080", size = 3, shape = 18, inherit.aes = FALSE) +
      # Error bar for FD Max
      geom_errorbarh(data = deriv_summary, 
                     aes(xmin = fd_max_lower_ci, xmax = fd_max_upper_ci, x = fd_max_value, y = 1), 
                     color = "black", height = 0.1, size = 1, inherit.aes = FALSE) +
      theme_classic() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      labs(title = "Significant Slope Regions (Purple Diamond = Max Curvature)", x = smooth_var, y = NULL)
    
    print(fd_ci_plot)
    
    
    # --- 6. PLOT: SIGNIFICANT INTERVALS (SD) ---
    # Filter for significance
    true.sd <- true.sd %>%
      filter(stringr::str_detect(.smooth, smooth_var)) %>%
      mutate(significant = ifelse(.lower_ci > 0 | .upper_ci < 0, "yes", "no")) %>%
      mutate(significant.deriv = ifelse(significant=="yes", .derivative, NA))
    
    sd_ci_plot <- ggplot(true.sd, aes(x = !!sym(smooth_var), y = 1, fill = significant.deriv)) +
      geom_tile(height = 0.2) +
      scale_fill_gradient2(low = "#003399", mid = "white", high = "#CC0000", na.value = "grey95", name = "Accel") +
      # Add Curvature Point (Purple Diamond) for reference
      geom_point(data = peak_vals %>% filter(type == "k max"), 
                 aes(x = !!sym(smooth_var), y = 1), color = "#800080", size = 3, shape = 18, inherit.aes = FALSE) +
      # Error bar for SD Max
      geom_errorbarh(data = deriv_summary, 
                     aes(xmin = sd_max_lower_ci, xmax = sd_max_upper_ci, x = sd_max_value, y = 1), 
                     color = "black", height = 0.1, size = 1, inherit.aes = FALSE) +
      theme_classic() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      labs(title = "Significant Acceleration (Purple Diamond = Max Curvature)", x = smooth_var, y = NULL)
    
    print(sd_ci_plot)
  }
  
  if (return_draws) {
    return(list(summary = deriv_summary, fd_draws = fd_samp, sd_draws = sd_samp,sd_extrema = sd.extrema, fd_extrema = fd.extrema))
  } else {
    return(deriv_summary)
  }
}