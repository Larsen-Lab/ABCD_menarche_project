library(gratia)
library(tidyverse)

# Function to extract significant derivative periods
get_derivative_table <- function(modobj, smooth_var) {
  
  # Handle both gam and gamm objects
  if (any(class(modobj) == "gam")) {
    model <- modobj
  } else if (class(modobj$gam) == "gam") {
    model <- modobj$gam
  } else {
    stop("Can't find a gam object to process")
  }
  
  # Get derivatives
  derv <- derivatives(model, term = sprintf('s(%s)', smooth_var), partial_match = TRUE)
  
  # Identify significant regions
  derv <- derv %>%
    mutate(
      sig = !(0 > .lower_ci & 0 < .upper_ci),
      direction = case_when(
        !sig ~ "non-significant",
        .derivative > 0 ~ "increase",
        .derivative < 0 ~ "decrease",
        TRUE ~ "non-significant"
      )
    )
  
  # Find continuous significant regions
  derv <- derv %>%
    mutate(
      region_id = cumsum(sig != lag(sig, default = FALSE) | 
                           direction != lag(direction, default = ""))
    )
  
  # Summarize significant regions
  sig_regions <- derv %>%
    filter(sig == TRUE) %>%
    group_by(region_id, direction) %>%
    summarise(
      start = min(!!sym(smooth_var)),
      end = max(!!sym(smooth_var)),
      mean_derivative = mean(.derivative),
      max_derivative = max(abs(.derivative)),
      n_points = n(),
      .groups = "drop"
    ) %>%
    select(-region_id) %>%
    arrange(start)
  
  return(sig_regions)
}

# Example usage with multiple models:
# derivative_stats <- brain_models_onset_output %>%
#   mutate(sig_changes = map(model, ~get_derivative_table(., smooth_var = "Time_to_from_Menarche_y"))) %>%
#   unnest(sig_changes) %>%
#   select(-model)