#' Capture user data and mappings separately from the computed placeholder
#' @noRd
ggchord_capture_layer_input <- function(lyr, data, mapping, roles) {
  mapping <- ggchord_normalize_mapping(mapping)
  if (rlang::is_formula(data)) data <- rlang::as_function(data)
  lyr$ggchord_input_data <- ggchord_normalize_accver(data)
  lyr$ggchord_input_mapping <- mapping
  lyr$ggchord_role_aes <- roles
  lyr
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

#' Resolve a layer's input table and role mappings
#' @noRd
ggchord_resolve_layer_input <- function(lyr, fallback = NULL) {
  data <- lyr$ggchord_input_data %||% fallback
  if (is.null(data)) return(NULL)
  if (is.function(data)) {
    data_function <- data
    data <- tryCatch(
      data_function(fallback),
      error = function(e) ggchord_stop(
        "ggchord layer data function failed: ", conditionMessage(e)
      )
    )
  }
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
    contains_stage <- function(x) {
      if (!is.call(x)) return(FALSE)
      head <- x[[1L]]
      name <- if (is.symbol(head)) as.character(head) else if (
        is.call(head) && as.character(head[[1L]]) %in% c("::", ":::")) {
        as.character(head[[3L]])
      } else ""
      if (name %in% c("after_stat", "after_scale")) return(TRUE)
      if (identical(name, "stage")) {
        args <- as.list(x)[-1L]
        late <- args[intersect(names(args), c("after_stat", "after_scale"))]
        if (any(!vapply(late, is.null, logical(1)))) return(TRUE)
      }
      any(vapply(as.list(x)[-1L], contains_stage, logical(1)))
    }
    if (contains_stage(expr)) {
      ggchord_stop(
        "Role aesthetic `", role,
        "` controls geometry and cannot use after_stat() or after_scale()"
      )
    }
    value <- tryCatch(
      rlang::eval_tidy(mapping[[role]], data = input),
      error = function(e) ggchord_stop(
        "Cannot evaluate `", role, "` in layer mapping: ", conditionMessage(e)
      )
    )
    if (length(value) == 1 && nrow(data) != 1) value <- rep(value, nrow(data))
    if (length(value) != nrow(data)) {
      ggchord_stop("Mapped `", role, "` must return one value per input row")
    }
    data[[role]] <- value
  }
  for (nm in setdiff(names(mapping), c(roles, names(lyr$aes_params)))) {
    expr <- rlang::get_expr(mapping[[nm]])
    text <- rlang::as_label(expr)
    if (!grepl("after_stat|after_scale", text)) {
      columns <- intersect(all.vars(expr), c("pident", "length", "anno", "type", "label"))
      ggchord_require_columns(data, columns, paste0("Mapping `", nm, "`"))
    }
  }
  transform <- lyr$ggchord_input_transform
  if (!is.null(transform)) data <- transform(data)
  if (!is.null(source)) {
    if (is.null(lyr$ggchord_params$input_source) && isTRUE(lyr$ggchord_params$is_feature) && "segments" %in% names(data) &&
        ".biological_source_row" %in% names(data))
      source$rows <- source$rows[data$.biological_source_row]
    attr(data, "ggchord_source") <- source
  }
  if (type %in% c("gene", "gene_label", "gene_label_repel") &&
      !isTRUE(lyr$ggchord_params$is_feature) && !isTRUE(lyr$ggchord_params$is_feature_label) &&
      (is.null(data$strand) || anyNA(data$strand) || any(!data$strand %in% c("+", "-")))) ggchord_stop(
    "geom_gene(): directional genes need '+'/'-' strand; use geom_feature() for nondirectional intervals"
  )
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
