#' Process sequence-related parameters
#'
#' Standardizes sequence parameters (e.g., radius, gap) from various input formats (single value/vector/named vector) into a vector named by sequence IDs
#'
#' @param param Input parameter (can be NULL, single value, vector, named vector)
#' @param seqs Character vector, list of sequence IDs
#' @param param_name Character, name of the parameter (used in error messages)
#' @param default_value Default value when param is NULL, default NULL
#' @param allow_null Logical, whether to allow param to be NULL, default FALSE
#' @return Named vector (names are seq_ids), standardized parameter values
#' @keywords internal
process_sequence_param <- function(param, seqs, param_name, default_value = NULL, allow_null = FALSE) {
  n <- length(seqs)

  # Handle case where parameter is NULL
  if (is.null(param)) {
    if (allow_null) return(NULL)

    # Ensure default value is a single value or vector matching sequence length
    if (!is.null(default_value)) {
      if (length(default_value) == 1) {
        return(setNames(rep(default_value, n), seqs))
      } else if (length(default_value) == n) {
        # Ensure default value vector has correct naming
        if (is.null(names(default_value))) {
          return(setNames(default_value, seqs))
        } else {
          # Check if names match sequence IDs
          if (all(names(default_value) %in% seqs)) {
            return(default_value[seqs])
          } else {
            warning(paste0("Names of default ", param_name, " do not fully match sequence IDs; using positional matching"))
            return(setNames(default_value, seqs))
          }
        }
      } else {
        stop(paste0("Length of default ", param_name, " must be 1 or equal to the number of sequences (", n, ")"))
      }
    }

    stop(paste0(param_name, " cannot be NULL and no default value is specified"))
  }

  # Handle case where parameter is a single value
  if (length(param) == 1) {
    return(setNames(rep(param, n), seqs))
  }

  # Handle case where parameter is an unnamed vector with length matching number of sequences
  if (length(param) == n && is.null(names(param))) {
    return(setNames(param, seqs))
  }

  # Handle case where parameter is a named vector
  if (!is.null(names(param))) {
    # Check if all names are in sequence IDs
    if (!all(names(param) %in% seqs)) {
      stop(paste0(param_name, " contains unknown sequence IDs: ",
                  paste(setdiff(names(param), seqs), collapse = ", ")))
    }

    # Return parameter values in sequence ID order
    return(param[seqs])
  }

  stop(paste0("Unprocessable format for ", param_name, ". Please provide a single value, named vector, or unnamed vector with length equal to the number of sequences"))
}


#' Process gene-related parameters
#'
#' Standardizes gene parameters (e.g., offset, width) from various input formats (single value/vector/list) into a list seperated by sequence and strand
#'
#' @param param Input parameter (can be NULL, single value, vector, list)
#' @param seqs Character vector, list of sequence IDs
#' @param param_name Character, name of the parameter (used in error messages)
#' @param default_value Default value when param is NULL
#' @param is_logical Logical, whether the parameter is logical (TRUE/FALSE), default FALSE
#' @return List (named by seq_id), where each element is a vector with "+"/"-" (parameter values for positive/negative strands)
#' @keywords internal
process_gene_param <- function(param, seqs, param_name, default_value, is_logical = FALSE) {
  n <- length(seqs)
  # Initialize result list: each sequence corresponds to a vector containing "+" and "-"
  result <- setNames(lapply(seqs, function(id) {
    if (is_logical) {
      c("+" = default_value, "-" = default_value)  # Logical parameter
    } else {
      c("+" = default_value, "-" = default_value)  # Numeric parameter
    }
  }), seqs)

  if (is.null(param)) {
    return(result)
  }

  # 1. Single value (same value for all sequences and strands)
  if (length(param) == 1 && !is.list(param)) {
    # Validate parameter type
    if (is_logical && !is.logical(param)) {
      stop(paste(param_name, "is a logical parameter and must be TRUE/FALSE"))
    }
    if (!is_logical && !is.numeric(param)) {
      stop(paste(param_name, "is a numeric parameter and must be a number"))
    }
    val <- param
    # Key modification: remove names with unname(val)
    return(setNames(lapply(seqs, function(id) c("+" = unname(val), "-" = unname(val))), seqs))
  }

  # 2. Vector input (not a list)
  if (is.vector(param) && !is.list(param)) {
    # 2.1 Named vector: validate sequence names
    if (!is.null(names(param))) {
      unknown <- setdiff(names(param), seqs)
      if (length(unknown) > 0) {
        warning(paste(param_name, "contains non-existent sequence IDs:", paste(unknown, collapse = ", ")))
      }
      valid_names <- intersect(names(param), seqs)
      for (id in valid_names) {
        val <- param[id]
        # Validate type
        if (is_logical && !is.logical(val)) {
          stop(paste("Elements of named vector for", param_name, "must be TRUE/FALSE"))
        }
        if (!is_logical && !is.numeric(val)) {
          stop(paste("Elements of named vector for", param_name, "must be numeric"))
        }
        # Key modification: remove names with unname(val)
        result[[id]] <- c("+" = unname(val), "-" = unname(val))  # Vector does not distinguish strands; same value for both strands
      }
      return(result)
    }

    # 2.2 Unnamed vector: validate length match
    if (length(param) != n) {
      stop(paste("Length of unnamed vector for", param_name, "must match the number of sequences (current number of sequences:", n, ")"))
    }
    for (i in seq_along(seqs)) {
      val <- param[i]
      # Validate type
      if (is_logical && !is.logical(val)) {
        stop(paste("Elements of unnamed vector for", param_name, "must be TRUE/FALSE"))
      }
      if (!is_logical && !is.numeric(val)) {
        stop(paste("Elements of unnamed vector for", param_name, "must be numeric"))
      }
      # Key modification: remove names with unname(val)
      result[[seqs[i]]] <- c("+" = unname(val), "-" = unname(val))  # Vector does not distinguish strands
    }
    return(result)
  }

  # 3. List input (unchanged as issues do not involve lists)
  if (is.list(param)) {
    # 3.1 Single-value list (same values for positive/negative strands across all sequences)
    if (length(param) == 1 && is.null(names(param))) {
      elem <- param[[1]]
      # Check if list element contains "+" and "-"
      if (!all(c("+", "-") %in% names(elem))) {
        stop(paste("Single-value list element for", param_name, "must be a named vector containing '+' and '-'"))
      }
      # Validate type
      if (is_logical && (!is.logical(elem["+"]) || !is.logical(elem["-"]))) {
        stop(paste("List elements for", param_name, "must be logical values (TRUE/FALSE)"))
      }
      if (!is_logical && (!is.numeric(elem["+"]) || !is.numeric(elem["-"]))) {
        stop(paste("List elements for", param_name, "must be numeric"))
      }
      return(setNames(lapply(seqs, function(id) elem[c("+", "-")]), seqs))
    }

    # 3.2 Named list (specify values for positive/negative strands of specific sequences)
    if (!is.null(names(param))) {
      unknown <- setdiff(names(param), seqs)
      if (length(unknown) > 0) {
        warning(paste("List for", param_name, "contains non-existent sequence IDs:", paste(unknown, collapse = ", ")))
      }
      valid_names <- intersect(names(param), seqs)
      for (id in valid_names) {
        elem <- param[[id]]
        if (!all(c("+", "-") %in% names(elem))) {
          stop(paste("List element", id, "for", param_name, "must be a named vector containing '+' and '-'"))
        }
        # Validate type
        if (is_logical && (!is.logical(elem["+"]) || !is.logical(elem["-"]))) {
          stop(paste("List element", id, "for", param_name, "must be logical values (TRUE/FALSE)"))
        }
        if (!is_logical && (!is.numeric(elem["+"]) || !is.numeric(elem["-"]))) {
          stop(paste("List element", id, "for", param_name, "must be numeric"))
        }
        result[[id]] <- elem[c("+", "-")]
      }
      return(result)
    }

    # 3.3 Unnamed list (specify values for positive/negative strands in sequence order)
    if (length(param) != n) {
      stop(paste("Length of unnamed list for", param_name, "must match the number of sequences (current number of sequences:", n, ")"))
    }
    for (i in seq_along(seqs)) {
      elem <- param[[i]]
      if (!all(c("+", "-") %in% names(elem))) {
        stop(paste("Unnamed list element", i, "for", param_name, "must be a named vector containing '+' and '-'"))
      }
      # Validate type
      if (is_logical && (!is.logical(elem["+"]) || !is.logical(elem["-"]))) {
        stop(paste("Unnamed list element", i, "for", param_name, "must be logical values (TRUE/FALSE)"))
      }
      if (!is_logical && (!is.numeric(elem["+"]) || !is.numeric(elem["-"]))) {
        stop(paste("Unnamed list element", i, "for", param_name, "must be numeric"))
      }
      result[[seqs[i]]] <- elem[c("+", "-")]
    }
    return(result)
  }

  # Invalid format
  stop(paste("Invalid format for", param_name, ". Supported input methods:\n",
             "1. Single value (shared by all sequences/strands)\n",
             "2. Named vector (names are sequence IDs, shared by all strands)\n",
             "3. Unnamed vector (length matches number of sequences, shared by all strands)\n",
             "4. Single-value list (element is a named vector with '+'/'-', shared by all sequences)\n",
             "5. Named list (names are sequence IDs, elements are named vectors with '+'/'-')\n",
             "6. Unnamed list (length matches number of sequences, elements are named vectors with '+'/'-')"))
}


#' Process axis label orientation parameters
#'
#' Standardizes axis label orientation parameters in various formats (character/numeric/vector) into a named vector (mapped by sequence ID)
#'
#' @param param Character ("horizontal"), numeric (angle), vector (length matches number of sequences), or named vector, label orientation parameter
#' @param seqs Character vector, list of sequence IDs
#' @return Named vector (names are seq_id), values are "horizontal" or numeric angles
#' @keywords internal
process_axis_orientation <- function(param, seqs) {
  n <- length(seqs)

  # Handle single value "horizontal"
  if (is.character(param) && length(param) == 1 && tolower(param) == "horizontal") {
    return(setNames(rep("horizontal", n), seqs))
  }

  # Handle single numeric value
  if (is.numeric(param) && length(param) == 1) {
    return(setNames(rep(param, n), seqs))
  }

  # Handle vector input (including mixed types)
  if (is.vector(param) && length(param) == n) {
    result <- character(n)
    names(result) <- seqs

    for (i in seq_along(param)) {
      val <- param[i]
      seq_id <- seqs[i]

      if (is.character(val) && tolower(val) == "horizontal") {
        result[seq_id] <- "horizontal"
      } else if (is.numeric(val)) {
        result[seq_id] <- as.character(val)
      } else {
        stop(paste("Element", i, "of axis_label_orientation has incorrect format; must be numeric or 'horizontal'"))
      }
    }

    return(result)
  }

  # Handle named vector
  if (!is.null(names(param)) && all(names(param) %in% seqs)) {
    result <- setNames(rep("0", n), seqs)
    for (id in names(param)) {
      val <- param[id]
      if (is.character(val) && tolower(val) == "horizontal") {
        result[id] <- "horizontal"
      } else if (is.numeric(val)) {
        result[id] <- as.character(val)
      } else {
        stop(paste("Format of", id, "in axis_label_orientation is incorrect; must be numeric or 'horizontal'"))
      }
    }
    return(result)
  }

  stop("Incorrect format for axis_label_orientation parameter. Please provide:\n",
       "- 'horizontal' (default)\n",
       "- A single numeric value\n",
       "- A vector with length matching the number of sequences (can mix numeric values and 'horizontal')\n",
       "- A named vector (names correspond to sequence IDs)")
}


#' Process gene color parameters in strand mode
#'
#' Standardizes gene color parameters in strand mode (color by strand direction) into a named vector with "+"/"-"
#'
#' @param gene_colors Color vector (can be NULL, single value, vector of length 2, named vector with "+"/"-")
#' @return Named vector (names are "+"/"-"), standardized color values (default "+" is red, "-" is blue)
#' @keywords internal
process_strand_colors <- function(gene_colors) {
  # Default values
  default <- c("+" = "#E41A1C", "-" = "#377EB8")
  if (is.null(gene_colors)) {
    return(default)
  }
  # Handle named vector
  if (!is.null(names(gene_colors))) {
    if (!all(names(gene_colors) %in% c("+", "-"))) {
      stop("In 'strand' mode, named vectors for gene_colors can only contain '+' and '-'")
    }
    res <- default
    res[names(gene_colors)] <- gene_colors
    return(res)
  } else {
    # Handle unnamed vector
    len <- length(gene_colors)
    if (len == 1) {
      return(c("+" = gene_colors[1], "-" = gene_colors[1]))
    } else if (len == 2) {
      return(c("+" = gene_colors[1], "-" = gene_colors[2]))
    } else {
      stop("In 'strand' mode, unnamed gene_colors must have length 1 or 2")
    }
  }
}

#' Process gene color parameters in manual mode
#'
#' Standardizes gene color parameters in manual mode (color by gene annotation) into a vector named by gene annotation
#'
#' @param gene_colors Color vector (can be NULL, single value, vector, named vector with gene annotations)
#' @param unique_anno Character vector, unique gene annotation names
#' @param gene_order Character vector, display order of genes in the legend, default NULL (order of appearance)
#' @return Named vector (names are gene annotations), standardized color values (default uses RColorBrewer Set1)
#' @keywords internal
process_manual_colors <- function(gene_colors, unique_anno, gene_order) {
  # Determine final gene order
  if (!is.null(gene_order)) {
    # Validate all elements in gene_order exist in unique_anno
    unknown <- setdiff(gene_order, unique_anno)
    if (length(unknown) > 0) {
      stop("'gene_order' contains unknown gene annotations: ", paste(unknown, collapse = ","))
    }
    # Ensure all unique_anno are included, ordered by gene_order, with unordered ones at the end
    final_order <- c(gene_order, setdiff(unique_anno, gene_order))
  } else {
    final_order <- unique_anno
  }
  n_anno <- length(final_order)

  # Generate default colors (for supplementation)
  if (n_anno <= 9) {
    default_pal <- brewer.pal(n_anno, "Set1")
  } else {
    default_pal <- colorRampPalette(brewer.pal(9, "Set1"))(n_anno)
  }
  names(default_pal) <- final_order

  if (is.null(gene_colors)) {
    return(default_pal)
  }

  # Handle named vector (matching annotations)
  if (!is.null(names(gene_colors))) {
    unknown <- setdiff(names(gene_colors), unique_anno)
    if (length(unknown) > 0) {
      warning("In 'manual' mode, gene_colors contains unknown annotations: ", paste(unknown, collapse = ","))
    }
    res <- default_pal
    # Only update colors present in final_order
    common_names <- intersect(names(gene_colors), final_order)
    res[common_names] <- gene_colors[common_names]
    return(res)
  } else {
    # Handle unnamed vector (use user-provided colors first, supplement with defaults if insufficient)
    len <- length(gene_colors)
    res <- character(n_anno)

    # Fill with user-provided colors
    if (len >= 1) {
      res[1:min(len, n_anno)] <- gene_colors[1:min(len, n_anno)]
    }

    # Supplement remaining with default colors
    if (len < n_anno) {
      res[(len + 1):n_anno] <- default_pal[(len + 1):n_anno]
    }

    names(res) <- final_order
    return(res)
  }
}
