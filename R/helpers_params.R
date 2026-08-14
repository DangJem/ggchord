#' Process sequence-related parameters
#'
#' Standardizes sequence parameters (e.g., radius, gap) from flexible input
#' formats into a vector named by sequence IDs. Supported formats:
#'
#' 1. A single value, recycled to every sequence.
#' 2. A named vector by sequence ID, e.g. `c("MT108731.1" = 3, "OR222515.1" = 1)`.
#' 3. An unnamed vector with length equal to the number of sequences (matched
#'    by sequence order).
#' 4. A list named by sequence ID.
#' 5. A list named by sequence order ("1", "2", ...).
#' 6. An unnamed list matched by sequence order (length 1 or equal to the
#'    number of sequences); a length-one list recycles.
#'
#' @param param Input parameter (can be NULL, single value, vector, named vector, or list)
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
    if (!is.null(default_value)) {
      if (length(default_value) == 1) {
        return(setNames(rep(default_value, n), seqs))
      } else if (length(default_value) == n) {
        if (is.null(names(default_value))) {
          return(setNames(default_value, seqs))
        } else {
          if (all(names(default_value) %in% seqs)) {
            return(default_value[seqs])
          } else {
            warning(paste0("Names of default ", param_name, " do not fully match sequence IDs; using positional matching"))
            return(setNames(default_value, seqs))
          }
        }
      } else {
        ggchord_stop(paste0("Length of default ", param_name, " must be 1 or equal to the number of sequences (", n, ")"))
      }
    }
    ggchord_stop(paste0(param_name, " cannot be NULL and no default value is specified"))
  }

  # Build the base vector used for unspecified sequences
  default_base <- function() {
    if (is.null(default_value)) return(rep(NA, n))
    if (length(default_value) == 1) return(rep(default_value, n))
    if (length(default_value) == n) return(default_value)
    ggchord_stop(paste0("Length of default ", param_name, " must be 1 or equal to the number of sequences (", n, ")"))
  }

  # Handle a single unnamed value: recycled to every sequence
  if (length(param) == 1 && is.null(names(param)) && !is.list(param)) {
    return(setNames(rep(param, n), seqs))
  }

  # Handle a plain unnamed vector with length matching the number of sequences
  if (length(param) == n && is.null(names(param)) && !is.list(param)) {
    return(setNames(param, seqs))
  }

  # Handle a named vector (names are sequence IDs)
  if (!is.null(names(param)) && !is.list(param)) {
    if (!all(names(param) %in% seqs)) {
      ggchord_stop(paste0(param_name, " contains unknown sequence IDs: ",
                  paste(setdiff(names(param), seqs), collapse = ", ")))
    }
    # Missing sequences keep the default value (or NA when no default is set)
    out <- setNames(default_base(), seqs)
    out[names(param)] <- param
    return(out)
  }

  # Handle list input (named by sequence ID, named by order index, or unnamed)
  if (is.list(param)) {
    nm <- names(param)
    idx_names <- as.character(seq_along(seqs))

    # length-one list: recycle the single element
    if (length(param) == 1 && (is.null(nm) || all(nm == ""))) {
      val <- param[[1]]
      if (length(val) != 1) {
        ggchord_stop(paste0("Elements of the list for ", param_name, " must be single values"))
      }
      return(setNames(rep(val, n), seqs))
    }

    # named by sequence ID
    if (!is.null(nm) && length(nm) > 0 && all(nm %in% seqs)) {
      out <- setNames(default_base(), seqs)
      for (id in nm) {
        val <- param[[id]]
        if (length(val) != 1) {
          ggchord_stop(paste0("Elements of the list for ", param_name, " must be single values"))
        }
        out[id] <- val
      }
      return(out)
    }

    # named by sequence order ("1", "2", ...)
    if (!is.null(nm) && length(nm) > 0 && all(nm %in% idx_names)) {
      out <- setNames(default_base(), seqs)
      for (i in seq_along(nm)) {
        val <- param[[i]]
        if (length(val) != 1) {
          ggchord_stop(paste0("Elements of the list for ", param_name, " must be single values"))
        }
        out[seqs[as.integer(nm[i])]] <- val
      }
      return(out)
    }

    # unnamed list matched by sequence order
    if (is.null(nm) || all(nm == "")) {
      if (length(param) != n) {
        ggchord_stop(paste0("Length of unnamed list for ", param_name,
                    " must be 1 or match the number of sequences (", n, ")"))
      }
      vals <- unlist(param)
      if (length(vals) != n) {
        ggchord_stop(paste0("Elements of the list for ", param_name, " must be single values"))
      }
      return(setNames(vals, seqs))
    }

    ggchord_stop(paste0("Names of the list for ", param_name,
                " must all be sequence IDs, all be sequence order indices ('1', '2', ...), or be omitted"))
  }

  ggchord_stop(paste0("Unprocessable format for ", param_name, ". Please provide a single value, named vector, unnamed vector with length equal to the number of sequences, or a list (named by sequence ID, named by order, or unnamed)"))
}

#' Process gene-related parameters
#'
#' Standardizes gene parameters (e.g., offset, width, label rotation) from
#' flexible input formats into a list separated by sequence and strand. All of
#' the following formats are supported:
#'
#' 1. A single value: applied to every sequence and strand, e.g. `20`.
#' 2. A named vector by strand, applied to every sequence:
#'    `c("+" = -15, "-" = -45)` (a missing strand keeps the default).
#' 3. A named vector by sequence ID (same value on both strands), e.g.
#'    `c("MT118296.1" = 20, "OR222515.1" = 30)`.
#' 4. An unnamed vector with length equal to the number of sequences (same
#'    value on both strands, matched by sequence order).
#' 5. A list named by sequence ID, each element a `+`/`-` named vector.
#' 6. A list named by sequence order ("1", "2", ...), each element a `+`/`-`
#'    named vector.
#' 7. An unnamed list with length equal to the number of sequences, each
#'    element a `+`/`-` named vector (matched by sequence order).
#' 8. A length-one list that recycles: `list(20)` applies 20 to everything,
#'    `list(c("+" = -15, "-" = -45))` applies the per-strand values to every
#'    sequence.
#'
#' @param param Input parameter (can be NULL, single value, vector, or list)
#' @param seqs Character vector, list of sequence IDs
#' @param param_name Character, name of the parameter (used in error messages)
#' @param default_value Default value when param is NULL
#' @param is_logical Logical, whether the parameter is logical (TRUE/FALSE), default FALSE
#' @return List (named by seq_id), where each element is a vector with "+"/"-"
#' @keywords internal
process_gene_param <- function(param, seqs, param_name, default_value, is_logical = FALSE) {
  n <- length(seqs)
  result <- setNames(lapply(seqs, function(id) {
    c("+" = default_value, "-" = default_value)
  }), seqs)

  if (is.null(param)) return(result)

  # Validate and strip a single value
  check_value <- function(val, where) {
    if (length(val) != 1 || is.na(val)) {
      ggchord_stop(param_name, " must be a single non-missing value (", where, ")")
    }
    if (is_logical && !is.logical(val)) {
      ggchord_stop(param_name, " is a logical parameter and must be TRUE/FALSE (", where, ")")
    }
    if (!is_logical && (!is.numeric(val) || !is.finite(val))) {
      ggchord_stop(param_name, " is a numeric parameter and must be finite (", where, ")")
    }
    unname(val)
  }

  # Check that an element is a per-strand specification (named "+"/"-" or an
  # unnamed length-2 vector interpreted as "+", "-")
  strand_spec <- function(elem, where) {
    if (is.null(elem) || length(elem) == 0) {
      ggchord_stop("List element for ", param_name, " (", where, ") cannot be empty")
    }
    if (is.null(names(elem))) {
      if (length(elem) == 2 && (is.numeric(elem) || is.logical(elem))) {
        elem <- c("+" = elem[1], "-" = elem[2])
      } else {
        ggchord_stop("List element for ", param_name, " (", where,
             ") must be a named vector containing '+' and '-' (or a length-2 vector)")
      }
    } else if (!all(names(elem) %in% c("+", "-"))) {
      ggchord_stop("List element for ", param_name, " (", where,
           ") can only be named '+' and/or '-'")
    }
    out <- c("+" = default_value, "-" = default_value)
    for (st in names(elem)) {
      out[st] <- check_value(elem[[st]], paste0(where, ", strand '", st, "'"))
    }
    out
  }

  # ---- 1. single value: everything gets the same value ----
  if (length(param) == 1 && !is.list(param) && is.null(names(param))) {
    val <- check_value(param, "single value")
    return(setNames(lapply(seqs, function(id) c("+" = val, "-" = val)), seqs))
  }

  # ---- 2. plain vector (not a list) ----
  if (is.vector(param) && !is.list(param)) {
    nm <- names(param)
    # 2.0 named by strand: global per-strand values
    if (!is.null(nm) && all(nm %in% c("+", "-"))) {
      vals <- strand_spec(param, "named vector")
      return(setNames(lapply(seqs, function(id) vals), seqs))
    }
    # 2.1 named by sequence ID (same value on both strands)
    if (!is.null(nm)) {
      unknown <- setdiff(nm, seqs)
      if (length(unknown) > 0) {
        ggchord_stop(param_name, " contains unknown sequence IDs: ",
             paste(unknown, collapse = ", "))
      }
      for (id in intersect(nm, seqs)) {
        result[[id]] <- c("+" = check_value(param[[id]], id),
                          "-" = check_value(param[[id]], id))
      }
      return(result)
    }
    # 2.2 unnamed vector matched by sequence order
    if (length(param) != n) {
      ggchord_stop("Length of unnamed vector for ", param_name,
           " must match the number of sequences (current number of sequences: ", n, ")")
    }
    for (i in seq_along(seqs)) {
      result[[seqs[i]]] <- c("+" = check_value(param[i], seqs[i]),
                             "-" = check_value(param[i], seqs[i]))
    }
    return(result)
  }

  # ---- 3. list ----
  if (is.list(param)) {
    nm <- names(param)
    idx_names <- as.character(seq_along(seqs))

    # 3.0 named by sequence ID
    if (!is.null(nm) && length(nm) > 0 && all(nm %in% seqs)) {
      for (id in nm) {
        result[[id]] <- strand_spec(param[[id]], id)
      }
      return(result)
    }
    # 3.1 named by sequence order ("1", "2", ...)
    if (!is.null(nm) && length(nm) > 0 && all(nm %in% idx_names)) {
      for (i in seq_along(nm)) {
        result[[seqs[as.integer(nm[i])]]] <- strand_spec(param[[i]], nm[i])
      }
      return(result)
    }
    # 3.2 unnamed list
    if (is.null(nm) || all(nm == "")) {
      # length-one list: recycle
      if (length(param) == 1) {
        elem <- param[[1]]
        if (length(elem) == 1 && is.null(names(elem))) {
          # list(20): scalar recycled everywhere
          val <- check_value(elem, "single-value list")
          return(setNames(lapply(seqs, function(id) c("+" = val, "-" = val)), seqs))
        }
        # list(c("+"=.., "-"=..)) or list(c(.., ..)): global per-strand
        vals <- strand_spec(elem, "single-value list")
        return(setNames(lapply(seqs, function(id) vals), seqs))
      }
      if (length(param) != n) {
        ggchord_stop("Length of unnamed list for ", param_name,
             " must be 1 or match the number of sequences (current number of sequences: ", n, ")")
      }
      for (i in seq_along(seqs)) {
        result[[seqs[i]]] <- strand_spec(param[[i]], as.character(i))
      }
      return(result)
    }
    # 3.3 mixed / unknown names
    ggchord_stop("Names of the list for ", param_name,
         " must all be sequence IDs, all be sequence order indices ('1', '2', ...), or be omitted")
  }

  ggchord_stop("Invalid format for ", param_name, ". Supported input methods:\n",
       "1. Single value (shared by all sequences/strands)\n",
       "2. Named vector by strand (c('+' = .., '-' = ..), shared by all sequences)\n",
       "3. Named vector by sequence ID (shared by both strands)\n",
       "4. Unnamed vector (length matches the number of sequences, shared by both strands)\n",
       "5. List named by sequence ID (elements are '+'/'-' named vectors)\n",
       "6. List named by sequence order '1', '2', ... (elements are '+'/'-' named vectors)\n",
       "7. Unnamed list (length matches the number of sequences, elements are '+'/'-' named vectors)\n",
       "8. Length-one list (scalar or '+'/'-' vector, recycled to all sequences)")
}

#' Process axis label orientation parameters
#'
#' Standardizes axis label orientation parameters in various formats (character/numeric/vector) into a named vector (mapped by sequence ID)
#'
#' @param param Character ("horizontal", "parallel" or "perpendicular"),
#'   numeric (angle), vector (length matches number of sequences), or named
#'   vector, label orientation parameter
#' @param seqs Character vector, list of sequence IDs
#' @return Named vector (names are seq_id), values are "horizontal",
#'   "parallel", "perpendicular" or numeric angles
#' @keywords internal
process_axis_orientation <- function(param, seqs) {
  n <- length(seqs)

  is_keyword <- function(x) {
    is.character(x) && tolower(x) %in% c("horizontal", "parallel", "perpendicular")
  }

  # Handle a single keyword ("horizontal", "parallel" or "perpendicular")
  if (is.character(param) && length(param) == 1 && is_keyword(param)) {
    return(setNames(rep(tolower(param), n), seqs))
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

      num_val <- suppressWarnings(as.numeric(val))
      if (is_keyword(val)) {
        result[seq_id] <- tolower(val)
      } else if (is.numeric(val) || (!is.na(num_val) && num_val == val)) {
        result[seq_id] <- as.character(num_val)
      } else {
        ggchord_stop(paste("Element", i, "of axis_label_orientation has incorrect format; must be numeric, 'horizontal', 'parallel' or 'perpendicular'"))
      }
    }

    return(result)
  }

  # Handle named vector
  if (!is.null(names(param)) && all(names(param) %in% seqs)) {
    result <- setNames(rep("0", n), seqs)
    for (id in names(param)) {
      val <- param[id]
      num_val <- suppressWarnings(as.numeric(val))
      if (is_keyword(val)) {
        result[id] <- tolower(val)
      } else if (is.numeric(val) || (!is.na(num_val) && num_val == val)) {
        result[id] <- as.character(num_val)
      } else {
        ggchord_stop(paste("Format of", id, "in axis_label_orientation is incorrect; must be numeric, 'horizontal', 'parallel' or 'perpendicular'"))
      }
    }
    return(result)
  }

  ggchord_stop("Incorrect format for axis_label_orientation parameter. Please provide:\n",
       "- A keyword: 'horizontal', 'parallel' (default, parallel to the axis) or 'perpendicular'\n",
       "- A single numeric value\n",
       "- A vector with length matching the number of sequences (can mix keywords and numeric values)\n",
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
      ggchord_stop("In 'strand' mode, named vectors for gene_colors can only contain '+' and '-'")
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
      ggchord_stop("In 'strand' mode, unnamed gene_colors must have length 1 or 2")
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
#' @return Named vector (names are gene annotations), standardized color values (default uses the built-in Set1 palette)
#' @keywords internal
process_manual_colors <- function(gene_colors, unique_anno, gene_order) {
  # Determine final gene order
  if (!is.null(gene_order)) {
    # Validate all elements in gene_order exist in unique_anno
    unknown <- setdiff(gene_order, unique_anno)
    if (length(unknown) > 0) {
      ggchord_stop("'gene_order' contains unknown gene annotations: ", paste(unknown, collapse = ","))
    }
    # Ensure all unique_anno are included, ordered by gene_order, with unordered ones at the end
    final_order <- c(gene_order, setdiff(unique_anno, gene_order))
  } else {
    final_order <- unique_anno
  }
  n_anno <- length(final_order)

  # Generate default colors (for supplementation)
  default_pal <- chord_default_palette(n_anno)
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

# ---------------------------------------------------------------------------
# Sequence grouping helper (v0.8.0)
# ---------------------------------------------------------------------------

#' Resolve a seq_group specification into a named character vector
#'
#' @param seq_data data.frame containing at least `seq_id` and, optionally,
#'   a `seq_group` column.
#' @param seqs Character vector of sequence IDs in drawing order.
#' @param seq_group NULL, a column name in `seq_data`, or a parameter accepted
#'   by [process_sequence_param()] (single value, named vector, unnamed vector
#'   matching the sequences, or a list).
#' @return A named character vector with one element per sequence (names are
#'   sequence IDs), or NULL when grouping is disabled.
#' @keywords internal
resolve_ggchord_seq_group <- function(seq_data, seqs, seq_group) {
  if (is.null(seq_group)) {
    if (is.data.frame(seq_data) && "seq_group" %in% colnames(seq_data)) {
      seq_group <- stats::setNames(as.character(seq_data$seq_group),
                                   as.character(seq_data$seq_id))
      seq_group <- seq_group[seqs]
      return(seq_group)
    }
    return(NULL)
  }

  # A single unnamed string is treated as a column name when that column
  # exists in seq_data; otherwise it is a group label recycled to every
  # sequence (which is rarely useful but is a valid process_sequence_param()
  # input).
  if (length(seq_group) == 1 && is.character(seq_group) &&
      is.null(names(seq_group)) && !is.list(seq_group) &&
      is.data.frame(seq_data) && seq_group %in% colnames(seq_data)) {
    seq_group <- stats::setNames(as.character(seq_data[[seq_group]]),
                                 as.character(seq_data$seq_id))
    seq_group <- seq_group[seqs]
    return(seq_group)
  }

  out <- process_sequence_param(seq_group, seqs, "seq_group", allow_null = TRUE)
  if (is.null(out)) return(NULL)
  as.character(out)
}

#' Normalise a named vector of group colours
#'
#' @param colors Named vector, or NULL. When NULL a default grey is used later.
#' @param groups Character vector of group names in display order.
#' @return Named character vector with names equal to `groups`; unnamed colour
#'   vectors are recycled positionally.
#' @keywords internal
resolve_ggchord_group_colors <- function(colors, groups) {
  if (is.null(colors)) return(NULL)
  if (is.null(names(colors))) {
    if (length(colors) == 1) colors <- setNames(rep(colors, length(groups)), groups)
    else if (length(colors) == length(groups)) colors <- setNames(as.character(colors), groups)
    else ggchord_stop("seq_group_colors must be length 1 or match the number of groups")
  } else {
    unknown <- setdiff(names(colors), groups)
    if (length(unknown) > 0) {
      ggchord_stop("seq_group_colors contains unknown group name(s): ",
                   paste(unknown, collapse = ", "))
    }
    out <- rep("grey20", length(groups))
    names(out) <- groups
    out[names(colors)] <- as.character(colors)
    colors <- out
  }
  colors
}
