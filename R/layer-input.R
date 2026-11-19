# R partially matches named arguments before the function sees its dots.
# Inspect the original call without evaluating mappings or data callbacks.
ggchord_check_argument_names <- function(call, fn, forwarded = NULL) {
  # A dots-only signature expands a caller's ... without matching or forcing it.
  supplied <- names(match.call(function(...) NULL, call, envir = parent.frame(2L))[-1L])
  formal <- setdiff(c(names(formals(fn)), if (!is.null(forwarded)) names(formals(forwarded))), "...")
  unknown <- setdiff(supplied[nzchar(supplied)], formal)
  for (key in unknown) {
    matches <- formal[startsWith(formal, key)]
    if (length(matches)) ggchord_stop("Spell out argument '", key,
      "'. Matching names: ", paste(matches, collapse = ", "))
  }
}

# Validate the names before a discarded parameter can become a silent no-op.
ggchord_check_param_names <- function(dots, allowed, caller) {
  if (!length(dots)) return(invisible(NULL))
  keys <- names(dots)
  if (is.null(keys) || any(!nzchar(keys)))
    ggchord_stop(caller, " needs named parameters in `...`")
  canonical <- ggplot2::standardise_aes_names(keys)
  if (anyDuplicated(canonical))
    ggchord_stop(caller, " received duplicate parameter '", keys[anyDuplicated(canonical)], "'")
  bad <- which(!canonical %in% allowed)
  if (length(bad)) {
    hints <- vapply(keys[bad], function(key) {
      distance <- utils::adist(key, allowed)
      closest <- which(distance == min(distance))
      if (length(closest) == 1L && distance[closest] <= 2L)
        paste0("'", key, "' (did you mean '", allowed[closest], "'?)")
      else paste0("'", key, "'")
    }, character(1))
    ggchord_stop(caller, " received unknown parameter(s): ", paste(hints, collapse = ", "))
  }
  invisible(NULL)
}

ggchord_check_geom_params <- function(dots, geom, caller, position = "identity",
    stat = ggplot2::StatIdentity, extra = character()) {
  if (!length(dots)) return(invisible(NULL))
  if (identical(position, "identity")) position <- ggplot2::PositionIdentity
  if (!inherits(geom, "Geom") || !inherits(stat, "Stat") || !inherits(position, "Position")) {
    # Let ggplot2 resolve named extensions before asking what they consume.
    receiver <- ggplot2::layer(geom = geom, stat = stat, position = position)
    geom <- receiver$geom
    stat <- receiver$stat
    position <- receiver$position
  }
  geom_names <- union(geom$aesthetics(), geom$parameters(TRUE))
  stat_names <- stat$parameters(TRUE)
  allowed <- unique(c(geom_names, stat_names, position$aesthetics(),
    if (isTRUE(geom$rename_size)) "size", "layout", "key_glyph", "name", extra))
  ggchord_check_param_names(dots, allowed, caller)
  shared <- intersect(ggplot2::standardise_aes_names(names(dots)), intersect(geom_names, stat_names))
  if (length(shared) && !identical(stat, StatGgchordCoverage) &&
      !identical(stat, get0(class(stat)[1L], asNamespace("ggplot2"), inherits = FALSE))) {
    ggchord_stop(caller, ": ambiguous parameter '", shared[1L], "' is accepted by both ",
      class(stat)[1L], " (Stat) and ", class(geom)[1L], " (Geom). ",
      "Give the Stat parameter a distinct name.")
  }
  invisible(geom$aesthetics())
}

#' Capture user data and mappings separately from the computed placeholder
#' @noRd
ggchord_capture_layer_input <- function(lyr, data, mapping, roles) {
  if (ggchord_local_position(lyr$position) && is.null(lyr$ggchord_params$feature_position))
    ggchord_stop("Layer '", lyr$ggchord_type, "' cannot use a feature position")
  mapping <- ggchord_normalize_mapping(mapping)
  if (rlang::is_formula(data)) data <- rlang::as_function(data)
  lyr$ggchord_input_data <- ggchord_normalize_accver(data)
  lyr$ggchord_input_mapping <- mapping
  lyr$ggchord_role_aes <- roles
  lyr
}

# Capture the first source result on the render-owned layer. Input projection
# and display-instance expansion still run in their respective contexts.
ggchord_render_data_once <- function(callback) {
  force(callback)
  resolved <- FALSE
  value <- NULL
  reader <- function(input) {
    if (!resolved) {
      value <<- callback(input)
      resolved <<- TRUE
    }
    value
  }
  attr(reader, "ggchord_render_source") <- TRUE
  reader
}

ggchord_native_facet_input <- function(layer) {
  copy <- ggplot2::ggproto(NULL, layer)
  # Materialize the subclass method on this render's clone, including any
  # transformation it applies after reading data or running a data callback.
  read <- ggchord_render_data_once(copy$layer_data)
  copy$layer_data <- function(plot_data) read(plot_data)
  copy$ggchord_facet_input_cached <- TRUE
  copy
}

ggchord_render_inputs <- function(plot) {
  for (i in seq_along(plot$layers)) {
    layer <- plot$layers[[i]]
    # A native layer may introduce or rename samples. Facet selection and
    # ggplot2 drawing must see the same result from its layer_data() method.
    if (is.null(layer$ggchord_type) && !inherits(plot$facet, "FacetNull") &&
        !isTRUE(layer$ggchord_facet_input_cached)) {
      plot$layers[[i]] <- ggchord_native_facet_input(layer)
      next
    }
    reader <- layer$ggchord_input_data
    if (!is.function(reader) || isTRUE(attr(reader, "ggchord_render_source"))) next
    copy <- reconstruct_layer(layer, layer$data, mapping = layer$mapping)
    copy$ggchord_input_data <- ggchord_render_data_once(reader)
    plot$layers[[i]] <- copy
  }
  plot
}

#' Normalise public US spelling aliases before geometry or scale inference
#' @noRd
ggchord_normalize_mapping <- function(mapping) {
  if (is.null(mapping)) return(mapping)
  if ("seq_id" %in% names(mapping)) {
    ggchord_stop("`seq_id` mapping was removed in v0.13.1; use `accver`")
  }
  if (anyDuplicated(names(mapping))) {
    duplicated_names <- unique(names(mapping)[duplicated(names(mapping))])
    ggchord_stop(
      "Duplicated aesthetic after colour/color normalisation: ",
      paste(duplicated_names, collapse = ", ")
    )
  }
  aliases <- c(seq_color = "seq_colour", ribbon_color = "ribbon_colour", link_color = "link_colour")
  for (alias in names(aliases)) {
    canonical <- aliases[[alias]]
    if (alias %in% names(mapping) && canonical %in% names(mapping)) {
      ggchord_stop("Use only one of `", canonical, "` and `", alias, "`")
    }
    if (alias %in% names(mapping)) names(mapping)[names(mapping) == alias] <- canonical
  }
  mapping
}

#' Resolve a fixed colour/color alias without ambiguous precedence
#' @noRd
ggchord_colour_alias <- function(colour, dots, caller, call = sys.call(-1L)) {
  raw <- names(as.list(call)[-1L]) %||% character()
  has_color <- "color" %in% names(dots)
  has_colour <- "colour" %in% raw || "colour" %in% names(dots)
  if (has_color && has_colour) {
    ggchord_stop(caller, ": use only one of colour and color")
  }
  if (has_color) {
    colour <- dots$color
    dots$color <- NULL
  }
  if ("color" %in% names(dots)) names(dots)[names(dots) == "color"] <- "colour"
  list(colour = colour, dots = dots)
}

#' Normalise colour/color inside variadic fixed aesthetics
#' @noRd
ggchord_colour_dots <- function(dots, caller) {
  if ("color" %in% names(dots) && "colour" %in% names(dots))
    ggchord_stop(caller, ": use only one of colour and color")
  if ("color" %in% names(dots)) names(dots)[names(dots) == "color"] <- "colour"
  dots
}

ggchord_layer_input_data <- function(layer, fallback) {
  input <- layer$ggchord_input_data %||% fallback
  if (!is.function(input)) return(input)
  tryCatch(input(fallback), error = function(e) rlang::abort(
    paste0("Data callback for layer '", layer$ggchord_layer_id, "' failed."),
    parent = e, call = NULL))
}

ggchord_layer_fallback <- function(layer, data) {
  p <- layer$ggchord_params
  if (p$type %in% c("ribbon", "alignment")) return(data$ribbon_data)
  if (p$type %in% c("gene", "gene_label", "gene_label_repel"))
    return(ggchord_annotation_fallback(data))
  if (identical(p$native_stat, "coverage")) return(data$feature_data)
  if (p$type == "seq_region") return(if (is.function(p$regions)) data$seq_data else p$regions)
  if (p$type %in% c("link", "restriction_site")) return(NULL)
  data$seq_data
}

ggchord_layer_input_context <- function(layer) {
  paste0("Layer '", layer$ggchord_layer_id, "', ", layer$ggchord_params$input_description)
}

ggchord_eval_input_mapping <- function(role, mapping, data, caller) {
  tryCatch(rlang::eval_tidy(mapping, data = data),
    error = function(e) rlang::abort(paste0(caller, ": cannot evaluate `", role, " = ",
      rlang::as_label(mapping), "`."), parent = e, call = NULL))
}

# Row numbers refer to the supplied table, before panels and display instances.
ggchord_input_rows <- function(bad) {
  rows <- which(bad)
  if (!length(rows)) return("")
  paste0(" (input rows: ", paste(utils::head(rows, 8L), collapse = ", "),
    if (length(rows) > 8L) paste0(", +", length(rows) - 8L, " more"), ")")
}

ggchord_check_reference_ids <- function(data, fields, ids, caller) {
  ggchord_require_columns(data, fields, caller)
  for (field in fields) {
    values <- as.character(data[[field]])
    bad <- is.na(values) | !values %in% ids
    if (any(bad)) ggchord_stop(caller, ": unknown reference accession in ", field, " ",
      paste(unique(values[bad]), collapse = ", "), ggchord_input_rows(bad))
  }
}

#' Read and adapt a layer's input before evaluating mappings
#' @noRd
ggchord_layer_source <- function(lyr, fallback) {
  data <- ggchord_layer_input_data(lyr, fallback)
  if (is.null(data)) return(NULL)
  type <- lyr$ggchord_params$type %||% ""
  if (inherits(data, "GRanges") || inherits(data, "Seqinfo")) {
    data <- if (type %in% c("seq", "seq_label", "seq_break", "seq_center_label")) {
      as_ggchord_seq(data)
    } else if (inherits(data, "GRanges")) as_ggchord_feature(data) else {
      ggchord_stop("Seqinfo supplies sequence lengths; use GRanges for feature layers")
    }
  } else if (inherits(data, "GInteractions")) {
    data <- as_ggchord_link(data)
    if (identical(type, "link")) {
      data$qpos <- (data$qstart + data$qend) / 2
      data$spos <- (data$sstart + data$send) / 2
      data$.point_projection <- "anchor_midpoint"
    }
  }
  if (!is.data.frame(data)) {
    ggchord_stop("ggchord layer data must be a data.frame, a supported interval object, ",
      "or a function returning supported data")
  }
  data <- ggchord_normalize_accver(as.data.frame(data, stringsAsFactors = FALSE))
  capture <- lyr$ggchord_params$stage_capture
  if (!is.null(capture) && "source" %in% capture$requested && is.null(capture$source))
    capture$source <- data
  data
}

# Inspect calls rather than printed labels, which may contain a stage name in
# a column, a string, or an abbreviated expression. Match positional stage()
# arguments to the same formals as ggplot2.
ggchord_has_stage <- function(x, stages = c("after_stat", "after_scale")) {
  if (missing(x)) return(FALSE)
  x <- rlang::get_expr(x)
  # ggplot2 still recognises these older spellings of computed aesthetics.
  if ("after_stat" %in% stages) {
    if (rlang::is_symbol(x)) return(grepl("^\\.\\.[A-Za-z._]+\\.\\.$", rlang::as_string(x)))
    if (rlang::is_call(x, "stat")) return(TRUE)
  }
  if (!is.call(x)) return(FALSE)
  if (rlang::is_call(x, "stage")) {
    args <- as.list(match.call(ggplot2::stage, x))[-1L]
    if (any(!vapply(args[intersect(names(args), stages)], is.null, logical(1)))) return(TRUE)
  }
  rlang::is_call(x, stages) || any(vapply(as.list(x)[-1L], ggchord_has_stage, logical(1), stages = stages))
}

# Starting expression available for grouping and scale training.
ggchord_start_mapping <- function(q) {
  expr <- rlang::get_expr(q)
  if (rlang::is_call(expr, "stage")) {
    first <- match.call(ggplot2::stage, expr)$start
    if (is.null(first)) return(NULL)
    return(rlang::new_quosure(first, rlang::get_env(q)))
  }
  if (ggchord_has_stage(expr)) NULL else q
}

# A render owns one source table per layer. The row map follows segment
# expansion, instance selection and faceting independently of geometry roles.
ggchord_visual_rows <- function(visual, rows) {
  if (is.null(visual$data)) return(visual)
  local <- new.env(parent = emptyenv())
  local$data <- visual$data
  local$values <- visual$values
  local$rows <- visual$rows[rows]
  local
}

# Classify data references without running the mapping. A computed subscript
# needs the combined mask because its target column is known only at runtime.
ggchord_visual_fields <- function(expr, env) {
  if (rlang::is_symbol(expr)) return(rlang::as_string(expr))
  if (!is.call(expr)) return(character())
  if (rlang::is_call(expr, c("$", "[["))) {
    if (rlang::is_symbol(expr[[2L]], ".env")) return(character())
    if (rlang::is_symbol(expr[[2L]], ".data")) {
      key <- expr[[3L]]
      if (rlang::is_call(expr, "$")) return(as.character(key))
      if (rlang::is_symbol(key)) key <- get0(as.character(key), env, inherits = TRUE)
      return(if (rlang::is_string(key)) key else NA_character_)
    }
  }
  unique(unlist(lapply(as.list(expr)[-1L], ggchord_visual_fields, env = env), use.names = FALSE))
}

# A source expression is evaluated once per layer, before replicating its
# values across geometry. Mixed expressions need the owning source row and
# the generated vertex. An expression using only generated columns stays
# with ggplot2. A list distinguishes an evaluated NULL from deferred work.
ggchord_visual_value <- function(visual, key, start, geometry, rows, projected = TRUE) {
  if (is.null(visual$data) || is.null(start) || !rlang::is_quosure(start)) return(NULL)
  fields <- ggchord_visual_fields(rlang::get_expr(start), rlang::get_env(start))
  source <- visual$data
  external <- setdiff(fields, names(source))
  if (!anyNA(external) && all(!external %in% names(geometry)) &&
      all(vapply(external, exists, logical(1), envir = rlang::get_env(start), inherits = TRUE))) {
    if (!exists(key, visual$values, inherits = FALSE)) {
      value <- rlang::eval_tidy(start, source)
      if (!is.null(value) && !length(value) %in% c(1L, nrow(source)))
        ggchord_stop("Mapped ", key, " must be scalar or supply one value per input row")
      if (length(value) == 1L) value <- rep(value, nrow(source))
      visual$values[[key]] <- list(value)
    }
    return(list(visual$values[[key]][[1L]][visual$rows[rows]]))
  }
  columns <- if (anyNA(fields)) names(source) else intersect(fields, names(source))
  if (!projected || !length(columns)) return(NULL)
  for (column in columns) geometry[[column]] <- source[[column]][visual$rows[rows]]
  list(rlang::eval_tidy(start, geometry))
}

ggchord_draw_layer <- function(layer, geometry, input = NULL) {
  mapping <- ggchord_effective_mapping(layer)
  visual <- layer$ggchord_params$visual_source
  prepared <- reconstruct_layer(layer, geometry, mapping)
  if (!nrow(geometry) || is.null(visual$data)) return(prepared)
  rows <- geometry$source_row
  if (is.null(rows)) {
    key <- intersect(c("display_instance_id", "accver"), intersect(names(geometry), names(input)))[1L]
    if (is.na(key)) return(prepared)
    rows <- match(geometry[[key]], input[[key]])
  }
  keys <- intersect(setdiff(names(layer$ggchord_input_mapping), layer$ggchord_role_aes), names(mapping))
  columns <- character()
  for (key in keys) {
    value <- ggchord_visual_value(visual, key, ggchord_start_mapping(mapping[[key]]), geometry, rows)
    if (is.null(value)) next
    if (is.null(value[[1L]])) {
      columns[key] <- NA_character_
      next
    }
    column <- tail(make.unique(c(names(geometry), paste0(".ggchord_start_", key))), 1L)
    geometry[[column]] <- value[[1L]]
    columns[key] <- column
  }
  if (!length(columns)) return(prepared)
  prepared$data <- geometry
  prepared$ggchord_params$visual_columns <- columns
  # Use source values for the initial evaluation only. Keeping the declared
  # mapping on the layer preserves guide titles and later stat/scale stages.
  ggplot2::ggproto(NULL, prepared, compute_aesthetics = function(self, data, plot) {
    start <- self$computed_mapping
    for (key in names(columns)) {
      column <- columns[[key]]
      expr <- if (is.na(column)) NULL else rlang::sym(column)
      start[[key]] <- rlang::new_quosure(expr, emptyenv())
    }
    evaluation <- ggplot2::ggproto(NULL, self, computed_mapping = start)
    ggplot2::ggproto_parent(prepared, evaluation)$compute_aesthetics(data, plot)
  })
}

#' Resolve a layer's input table and role mappings
#' @noRd
ggchord_resolve_layer_input <- function(lyr, fallback = NULL,
                                        data = ggchord_layer_source(lyr, fallback)) {
  if (!is.null(lyr$ggchord_params$panel_input)) return(lyr$ggchord_params$panel_input)
  if (!is.null(lyr$ggchord_params$backbone_input)) return(lyr$ggchord_params$backbone_input)
  if (is.null(data)) return(NULL)
  type <- lyr$ggchord_params$type %||% ""
  source <- lyr$ggchord_params$input_source %||%
    if (is.null(lyr$ggchord_input_data)) attr(fallback, "ggchord_source") else NULL
  # A data function may subset while retaining arbitrary data-frame attributes.
  # Only a resolver-owned row map may be used as constructor provenance.
  attr(data, "ggchord_source") <- NULL
  if (is.null(lyr$ggchord_input_data) && isTRUE(attr(fallback, "ggchord_annotation_input")) &&
      type %in% c("gene", "gene_label", "gene_label_repel") &&
      !isTRUE(lyr$ggchord_params$is_feature) && !isTRUE(lyr$ggchord_params$is_feature_label)) {
    selected <- ggchord_annotation_gene_rows(data)
    data <- data[selected, , drop = FALSE]
    source$rows <- source$rows[selected]
  }
  mapping <- lyr$ggchord_input_mapping
  roles <- intersect(lyr$ggchord_role_aes %||% character(0), names(mapping))
  input <- data
  for (role in roles) {
    expr <- rlang::get_expr(mapping[[role]])
    if (ggchord_has_stage(expr)) {
      ggchord_stop(
        ggchord_layer_input_context(lyr), ": role aesthetic `", role,
        "` controls geometry and cannot use after_stat() or after_scale()"
      )
    }
    value <- ggchord_eval_input_mapping(role, mapping[[role]], input, ggchord_layer_input_context(lyr))
    if (length(value) == 1 && nrow(data) != 1) value <- rep(value, nrow(data))
    if (length(value) != nrow(data)) {
      ggchord_stop(ggchord_layer_input_context(lyr), ": mapped `", role, " = ",
        rlang::as_label(mapping[[role]]), "` must return one value per input row (got ",
        length(value), " for ", nrow(data), " rows)")
    }
    data[[role]] <- value
  }
  reference_ids <- lyr$ggchord_params$reference_ids
  if (!is.null(reference_ids)) {
    fields <- if (type %in% c("link", "ribbon")) c("qaccver", "saccver") else "accver"
    ggchord_check_reference_ids(data, fields, reference_ids, ggchord_layer_input_context(lyr))
  }
  for (nm in setdiff(names(mapping), c(roles, names(lyr$aes_params)))) {
    expr <- rlang::get_expr(mapping[[nm]])
    if (!ggchord_has_stage(expr)) {
      columns <- intersect(all.vars(expr), c("pident", "length", "anno", "type", "label"))
      ggchord_require_columns(data, columns, paste0(ggchord_layer_input_context(lyr),
        ": mapping `", nm, " = ", rlang::as_label(mapping[[nm]]), "`"))
    }
  }
  capture <- lyr$ggchord_params$stage_capture
  if (!is.null(capture) && "mapped" %in% capture$requested &&
      is.null(capture$mapped) && is.null(lyr$ggchord_params$native_geom))
    capture$mapped <- data
  if (type %in% c("gene", "gene_label", "gene_label_repel") &&
      !isTRUE(lyr$ggchord_params$is_feature) && !isTRUE(lyr$ggchord_params$is_feature_label) &&
      (is.null(data$strand) || anyNA(data$strand) || any(!data$strand %in% c("+", "-")))) ggchord_stop(
    "geom_feature(preset = 'gene'): directional genes need '+'/'-' strand; use geom_feature() for nondirectional intervals"
  )
  transform <- lyr$ggchord_input_transform
  visual <- lyr$ggchord_params$visual_source
  capture_visual <- !is.null(visual) && is.null(visual$data) &&
    length(setdiff(names(mapping), c(roles, names(lyr$aes_params)))) > 0L
  if (capture_visual) {
    visual$data <- input
    visual$values <- new.env(parent = emptyenv())
    visual$rows <- seq_len(nrow(input))
    if (!is.null(transform)) {
      row_key <- tail(make.unique(c(names(data), ".ggchord_visual_row")), 1L)
      data[[row_key]] <- visual$rows
    }
  }
  if (!is.null(transform)) {
    data <- transform(data)
    if (capture_visual) {
      visual$rows <- data[[row_key]]
      data[[row_key]] <- NULL
    }
  }
  if (!is.null(source)) {
    if (is.null(lyr$ggchord_params$input_source) && isTRUE(lyr$ggchord_params$is_feature) && "segments" %in% names(data) &&
        ".biological_source_row" %in% names(data))
      source$rows <- source$rows[data$.biological_source_row]
    attr(data, "ggchord_source") <- source
  }

  data
}

#' Resolve omitted fixed styles against local and inherited user mappings
#' @noRd
ggchord_resolve_style_defaults <- function(layer, mapping) {
  defaults <- layer$ggchord_params$style_defaults
  if (is.null(defaults)) return(layer)
  user <- layer$ggchord_input_mapping
  if (isTRUE(layer$inherit.aes)) {
    inherited <- ggchord_normalize_mapping(mapping)
    for (nm in setdiff(names(inherited), names(user))) user[[nm]] <- inherited[[nm]]
  }
  layer <- reconstruct_layer(layer, layer$data)
  # Ribbon outlines and line types use role names only when that role is mapped.
  if (identical(layer$ggchord_type, "ribbon")) {
    for (nm in c("colour", "linetype")) {
      role <- paste0("ribbon_", nm)
      if (role %in% names(user) && nm %in% names(layer$aes_params)) {
        layer$aes_params[role] <- layer$aes_params[nm]
        layer$aes_params[nm] <- NULL
        names(defaults)[names(defaults) == nm] <- role
      }
    }
    old_key <- layer$geom$draw_key
    layer$geom <- make_ribbon_geom(
      outline = "ribbon_colour" %in% names(user),
      linetype = "ribbon_linetype" %in% names(user))
    layer$geom$draw_key <- old_key
  }
  mapped_defaults <- intersect(names(defaults)[defaults], names(user))
  layer$aes_params[mapped_defaults] <- NULL
  layer
}

#' Combine computed and user visual mappings
#' @noRd
ggchord_effective_mapping <- function(lyr) {
  out <- lyr$mapping
  user <- lyr$ggchord_input_mapping
  input <- lyr$ggchord_resolved_input
  if (identical(lyr$ggchord_type, "ribbon") && !is.null(input) &&
      !"pident" %in% names(input) && !"ribbon_fill" %in% names(user)) {
    out$ribbon_fill <- NULL
  }
  if (lyr$ggchord_type %in% c("gene_poly") && !is.null(input) &&
      (!"anno" %in% names(input) || all(is.na(input$anno)))) {
    # Gene fill defaults to strand, which remains valid without optional anno.
    # Only the annotation-based feature default needs to disappear here.
    if (!"feature_fill" %in% names(user)) out$feature_fill <- NULL
  }
  if (!is.null(user)) {
    visual <- setdiff(names(user), lyr$ggchord_role_aes %||% character(0))
    for (nm in visual) out[[nm]] <- user[[nm]]
  }
  # Fixed aesthetics bypass scales, including automatically supplied mappings.
  out[intersect(names(out), names(lyr$aes_params))] <- NULL
  out
}

#' Add original input columns to expanded geometry
#' @noRd
ggchord_attach_input_columns <- function(geometry, input) {
  if (is.null(input) || !is.data.frame(input) || nrow(geometry) == 0) {
    return(geometry)
  }
  protected <- c(
    "source_row",
    "x", "y", "x0", "y0", "x1", "y1", "xend", "yend", "group",
    "text_x", "text_y", "text", "label_x", "label_y", "label",
    "text_angle", "label_angle", "hjust", "vjust", "size", "alpha",
    "colour", "fill", "zfill", "zcolour", "zregionfill", "zoutline",
    "zlinetype", "outline_col", "linetype_val", "seq_colour",
    "ribbon_fill", "ribbon_alpha", "ribbon_colour",
    "ribbon_linetype", "gene_fill", "feature_fill", "feature_shape",
    "scale_length", "bp_length", "sequence_unit", ".ggchord_adjacent",
    "feature_label_colour", "feature_label_orientation", "region_fill", ".component",
    "track_id", "label_text_original", "label_size_requested", "label_size_final",
    "label_size_min", "label_fit_reason", ".ggchord_resolved_size"
  )
  if ("source_row" %in% names(geometry)) {
    idx <- geometry$source_row
  } else if ("accver" %in% names(geometry) && "accver" %in% names(input)) {
    idx <- match(as.character(geometry$accver), as.character(input$accver))
    geometry <- geometry[!is.na(idx), , drop = FALSE]
    idx <- idx[!is.na(idx)]
  } else {
    return(geometry)
  }
  valid <- !is.na(idx) & idx >= 1 & idx <= nrow(input)
  # Protect computed geometry, not names that the geometry does not contain.
  # A source column such as label may still supply an ordinary visual mapping.
  cols <- setdiff(names(input), intersect(protected, names(geometry)))
  for (nm in cols) {
    source <- input[[nm]]
    if (is.list(source) && !is.data.frame(source)) {
      value <- rep(list(NULL), nrow(geometry))
      value[valid] <- source[idx[valid]]
    } else {
      value <- rep(ggchord_typed_na(source), nrow(geometry))
      value[valid] <- source[idx[valid]]
    }
    geometry[[nm]] <- value
  }
  geometry
}

# Bind an explicitly used ggplot2 aesthetic to its standard scale and guide.
# Domain-specific aesthetics still supply independently trained default styles.
ggchord_bind_visual_mapping <- function(layer, plot) {
  if (is.null(layer$ggchord_type)) return(layer)
  if (identical(layer$ggchord_type, "track")) {
    keys <- layer$geom_params$track_keys
    user <- layer$ggchord_input_mapping
    use <- vapply(names(keys), function(name) {
      key <- keys[[name]]
      key %in% names(user) && plot$scales$has_scale(name) && !plot$scales$has_scale(key)
    }, logical(1))
    if (any(use)) {
      target <- keys[use]
      to_standard <- stats::setNames(names(target), unname(target))
      names(user)[match(target, names(user))] <- names(target)
      layer$ggchord_input_mapping <- ggchord_track_stage_mapping(user,
        layer$ggchord_params$track_id, replacements = to_standard)
      defaults <- layer$geom$default_aes
      names(defaults)[match(target, names(defaults))] <- names(target)
      layer$geom <- ggchord_geom_defaults(layer$geom, defaults)
      layer$geom_params$track_keys[use] <- names(target)
    }
    # ggplot2 matches a named legend flag to the aesthetic trained by its scale.
    # Resolve public names after choosing between plot and track-specific scales.
    if (!is.null(names(layer$show.legend))) {
      flags <- ggplot2::standardise_aes_names(names(layer$show.legend))
      for (name in names(keys)) {
        aliases <- c(name, paste0("track_", name), keys[[name]])
        flags[flags %in% aliases] <- layer$geom_params$track_keys[[name]]
      }
      names(layer$show.legend) <- flags
    }
    return(layer)
  }
  line_roles <- c(colour = "link_colour", alpha = "link_alpha",
    linewidth = "link_linewidth", linetype = "link_linetype")
  roles <- switch(layer$ggchord_type,
    alignment = line_roles, link = line_roles,
    gene_poly = c(fill = if (isTRUE(layer$ggchord_params$is_feature)) "feature_fill" else "gene_fill"),
    seq = c(colour = "seq_colour"), ribbon = c(fill = "ribbon_fill", alpha = "ribbon_alpha"),
    seq_region = c(fill = "region_fill"), NULL)
  user <- layer$ggchord_input_mapping
  if (isTRUE(layer$inherit.aes)) for (nm in setdiff(names(plot$mapping), names(user)))
    user[[nm]] <- plot$mapping[[nm]]
  for (nm in names(roles)) {
    role <- roles[[nm]]
    if (role %in% names(user)) next
    mapped <- nm %in% names(user)
    fixed <- nm %in% names(layer$aes_params)
    if (!mapped && !fixed) next
    if (mapped && !is.null(layer$aes_params[[role]])) {
      if (isTRUE(unname(layer$ggchord_params$style_defaults[role]))) layer$aes_params[[role]] <- NULL else {
        if (nm %in% names(layer$ggchord_input_mapping))
          ggchord_stop("Layer '", layer$ggchord_layer_id, "' both maps and fixes ", nm)
        layer$ggchord_input_mapping[[nm]] <- NULL
        next
      }
    }
    if (mapped && fixed && nm %in% names(layer$ggchord_input_mapping))
      ggchord_stop("Layer '", layer$ggchord_layer_id, "' both maps and fixes ", nm)
    defaults <- layer$geom$default_aes
    names(defaults)[names(defaults) == role] <- nm
    layer$geom <- ggchord_geom_defaults(layer$geom, defaults)
    if (!fixed) layer$mapping[[nm]] <- user[[nm]]
    layer$mapping[[role]] <- NULL
    if (role %in% names(layer$show.legend)) names(layer$show.legend)[names(layer$show.legend) == role] <- nm
    if (identical(role, "feature_fill")) {
      layer$ggchord_params$user_feature_fill <- TRUE
      layer$ggchord_params$gene_colors <- NULL
      layer$ggchord_params$feature_fill_mapping <- layer$mapping[[nm]]
      layer$ggchord_params$feature_fill_fixed <- layer$aes_params[[nm]]
    }
  }
  layer
}

# ggproto retains the parent expression. Force it before replacing layer$geom.
ggchord_geom_defaults <- function(parent, defaults) {
  force(parent)
  ggplot2::ggproto(NULL, parent, default_aes = defaults)
}
