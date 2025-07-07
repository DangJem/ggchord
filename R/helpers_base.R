#' Missing value handling operator
#'
#' Used to safely handle NULL values: returns y if x is NULL, otherwise returns x
#'
#' @param x Any R object (may be NULL)
#' @param y Default value to return when x is NULL
#' @return x if x is not NULL, otherwise y
#' @keywords internal
if_null_else <- function (x, y)
{
  if (is.null(x)) y else x
}


#' Process panel margin parameters
#'
#' Standardizes input margin parameters into a list containing t (top), r (right), b (bottom), l (left). Supports single-value or list input.
#'
#' @param arg_list Numeric (single value) or list (named/unnamed), margin parameters
#' @return List containing four elements: t, r, b, l (numeric, margin sizes)
#' @keywords internal
process_panel_margin <- function(arg_list) {
  # Initialize result list with default values of 0
  result <- list(t = 0, r = 0, b = 0, l = 0)

  # Check if input is a list
  if (!is.list(arg_list)) {
    # Handle single-value input
    if (is.numeric(arg_list) && length(arg_list) == 1) {
      value <- arg_list
      result <- list(t = value, r = value, b = value, l = value)
      return(result)
    } else {
      warning("Input is not a valid list or single numeric value; default values will be used")
      return(result)
    }
  }

  # Handle empty list
  if (length(arg_list) == 0) {
    return(result)
  }

  # Handle named list
  if (!is.null(names(arg_list)) && all(names(arg_list) != "")) {
    valid_names <- c("t", "r", "b", "l")

    # Iterate over each element of the input list
    for (name in names(arg_list)) {
      if (name %in% valid_names) {
        # Check if value is numeric
        if (is.numeric(arg_list[[name]]) && length(arg_list[[name]]) == 1) {
          result[[name]] <- arg_list[[name]]
        } else {
          warning(paste("Parameter", name, "is not a single numeric value; default 0 will be used"))
        }
      } else {
        warning(paste("Unknown parameter", name, "will be ignored"))
      }
    }
  }
  # Handle unnamed list
  else {
    param_order <- c("t", "r", "b", "l")
    num_args <- length(arg_list)

    # Handle single-element unnamed list
    if (num_args == 1 && is.numeric(arg_list[[1]])) {
      value <- arg_list[[1]]
      result <- list(t = value, r = value, b = value, l = value)
      return(result)
    }

    # Assign values in order
    for (i in 1:min(num_args, length(param_order))) {
      if (is.numeric(arg_list[[i]]) && length(arg_list[[i]]) == 1) {
        result[[param_order[i]]] <- arg_list[[i]]
      } else {
        warning(paste("Parameter at position", i, "is not a single numeric value; default 0 will be used"))
      }
    }
  }

  return(result)
}


#' Calculate plot extremes
#'
#' Extracts x/y coordinate extremes from all plot elements (sequence arcs, ribbons, gene arrows, etc.) for adjusting the plot range
#'
#' @param allRibbon data.frame, ribbon data (with x, y columns), default NULL
#' @param seqArcs List, sequence arc data (each element is a data frame with x, y, seq_id), default NULL
#' @param axisLines data.frame, axis line data (with x, y, seq_id columns), default NULL
#' @param axisTicks data.frame, tick mark data (with x0, y0, x1, y1, label_x, label_y columns), default NULL
#' @param gene_arrows data.frame, gene label data (with text_x, text_y columns), default NULL
#' @param gene_polys data.frame, gene arrow polygon data (with x, y columns), default NULL
#' @param show_axis Logical, whether to include extreme value calculation for axis-related elements, default FALSE
#' @return List containing x_min (minimum x), x_max (maximum x), y_min (minimum y), y_max (maximum y)
#' @keywords internal
get_plot_extremes <- function(allRibbon=NULL, seqArcs=NULL, axisLines=NULL, axisTicks=NULL, gene_arrows=NULL, gene_polys = NULL, show_axis = F) {
  # Initialize vectors to store x and y coordinates
  x_coords <- numeric(0)
  y_coords <- numeric(0)

  # 1. Process ribbons (allRibbon)
  if (!is.null(allRibbon) && nrow(allRibbon) > 0) {
    x_coords <- c(x_coords, allRibbon$x)
    y_coords <- c(y_coords, allRibbon$y)
  }
  # 2. Process sequence arcs (seqArcs, convert list to data frame)
  if (!is.null(seqArcs) && length(seqArcs) > 0) {
    seq_df <- do.call(rbind, seqArcs)
    if (nrow(seq_df) > 0) {
      x_coords <- c(x_coords, seq_df$x)
      y_coords <- c(y_coords, seq_df$y)
    }
  }
  # 3. Process gene labels (gene_arrows)
  if (!is.null(gene_arrows) && nrow(gene_arrows) > 0) {
    x_coords <- c(x_coords, gene_arrows$text_x)
    y_coords <- c(y_coords, gene_arrows$text_y)
  }
  # 4. Process axis lines (axisLines)
  if (show_axis && !is.null(axisLines) && nrow(axisLines) > 0) {
    x_coords <- c(x_coords, axisLines$x)
    y_coords <- c(y_coords, axisLines$y)
  }
  # 5. Process tick marks (axisTicks)
  if (show_axis && !is.null(axisTicks) && nrow(axisTicks) > 0) {
    # x coordinates: x0, x1, label_x
    x_coords <- c(x_coords, axisTicks$x0, axisTicks$x1, axisTicks$label_x)
    # y coordinates: y0, y1, label_y
    y_coords <- c(y_coords, axisTicks$y0, axisTicks$y1, axisTicks$label_y)
  }
  # 6. Process gene arrow polygons (gene_polys): contains x, y coordinates (polygon vertices)
  if (!is.null(gene_polys) && nrow(gene_polys) > 0) {
    x_coords <- c(x_coords, gene_polys$x)
    y_coords <- c(y_coords, gene_polys$y)
  }

  # Filter missing values (NA)
  x_coords <- x_coords[!is.na(x_coords)]
  y_coords <- y_coords[!is.na(y_coords)]

  # Calculate extremes (return list containing min/max for x and y)
  list(
    x_min = min(x_coords), # Left extreme
    x_max = max(x_coords), # Right extreme
    y_min = min(y_coords), # Bottom extreme
    y_max = max(y_coords) # Top extreme
  )
}
