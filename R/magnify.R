ggchord_magnify_spec <- function(from, to, region, highlight = TRUE,
    fill = "#187B77", alpha = 0.12, colour = "#187B77", linewidth = 0.4, id = NULL,
    connector = "interval", linetype = NULL) {
  ggchord_local_console_errors()
  caller <- "Fusion Magnify"
  connector <- ggchord_match_choice(connector, c("interval", "border"), "connector")
  if (!is.null(linetype)) {
    valid_codes <- c("blank", "solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
    if (length(linetype) != 1L || anyNA(linetype) ||
        !(is.numeric(linetype) && linetype %in% 0:6 ||
          is.character(linetype) && (linetype %in% valid_codes ||
            grepl("^([1-9A-Fa-f]{2}){1,4}$", linetype))))
      ggchord_stop(caller, ": use a single grid line type or NULL")
  }
  from_key <- ggchord_fusion_check_path(from, caller)
  to_key <- ggchord_fusion_check_path(to, caller)
  if (identical(from_key, to_key)) ggchord_stop(caller, ": provide distinct panels")
  ggchord_require_columns(region, c("accver", "start", "end"), caller)
  if (!is.character(region$accver) || anyNA(region$accver) || any(!nzchar(region$accver)))
    ggchord_stop(caller, ": accver must contain nonempty source IDs")
  ggchord_project_integer(region$start, "start", caller)
  ggchord_project_integer(region$end, "end", caller)
  if (!nrow(region) || any(region$start < 1 | region$start > region$end)) ggchord_stop(caller, ": intervals must have ordered endpoints")
  for (name in intersect(c("from_instance", "to_instance"), names(region))) {
    if (!is.character(region[[name]]) || any(!is.na(region[[name]]) & !nzchar(region[[name]])))
      ggchord_stop(caller, ": instance selectors must be strings")
  }
  for (name in intersect(c("from_facet", "to_facet"), names(region))) {
    if (!is.list(region[[name]])) ggchord_stop(caller, ": facet selectors must be a list column")
    for (value in region[[name]]) ggchord_fusion_facet(value, caller)
  }
  if (anyDuplicated(region[intersect(c("accver", "start", "end", "from_instance",
      "to_instance", "from_facet", "to_facet"), names(region))]))
    ggchord_stop(caller, ": duplicate region declarations require distinct instance selectors")
  if (!is.logical(highlight) || length(highlight) != 1L || is.na(highlight))
    ggchord_stop(caller, ": highlight must be TRUE or FALSE")
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1 ||
      !is.numeric(linewidth) || length(linewidth) != 1L || !is.finite(linewidth) || linewidth < 0)
    ggchord_stop(caller, ": alpha must be in [0,1] and linewidth must be nonnegative")
  for (value in list(fill, colour)) {
    if (length(value) != 1L) ggchord_stop(caller, ": colours must be scalar")
    tryCatch(col2rgb(value), error = function(e) ggchord_stop(caller, ": invalid colour"))
  }
  structure(list(from = from, to = to, from_key = from_key, to_key = to_key,
    region = as.data.frame(region), highlight = highlight, fill = fill,
    alpha = alpha, colour = colour, linewidth = linewidth, id = id,
    connector = connector, linetype = linetype),
    class = "ggchord_fusion_zoom")
}

#' Compose an overview and a source-aligned detail
#'
#' @param overview An independent source plot. Leave NULL when adding to an
#'   existing Fusion with `from` and `to`.
#' @param detail For Auto, an optional editor for the derived plot. Set
#'   `level = "residue"` to edit a registered character plot.
#'   For standalone Manual, an independent detail plot.
#' @param region Nonempty data frame of root source `accver`, `start` and `end`
#'   (one-based closed intervals). Use `from_instance` and `to_instance` to
#'   select display instances and named list columns `from_facet` and
#'   `to_facet` to select actual facet values.
#' @param detail_coord Sequence-path Coord for Feature detail.
#' @param level `"feature"` derives a sequence window, `"residue"` displays
#'   registered real characters, and `"auto"` uses a residue view when the
#'   actual panel fits its glyphs. An explicit editor or detail Coord selects
#'   the Feature contract in Auto mode.
#' @param expand Nonnegative integer padding in source coordinate units,
#'   limited by the source view. The highlighted interval stays unchanged.
#' @param layout Optional layout for a standalone composition.
#' @param guides Guide policy for a standalone composition.
#' @param ... Named fixed decoration parameters: `highlight`, `fill`,
#'   `colour`, `alpha`, `linewidth`, `linetype`, `connector`, and `id`.
#'   Set `connector = "border"` to lead to the facing detail border, or use
#'   the default `"interval"` to join source and detail interval fronts.
#' @param from,to Paths when adding to an existing Fusion. Auto creates a new
#'   direct child named `to`; Manual connects two existing leaf views.
#' @return A Fusion for standalone calls, or a Magnify declaration to add to
#'   a Fusion.
#' @details Auto derives display windows over the original reference and source
#'   observations. Feature detail has the original visual and Track domains.
#'   Residue detail uses the reference `sequence` column or a registered
#'   character plot. Set `alphabet` in the reference when its unit and
#'   characters cannot identify DNA, RNA or amino acids. Observed alignment
#'   cells retain their original columns and gaps, and reference windows are
#'   shown on independent rows. Split selections at a circular origin.
#'   Use Manual for complex faceted or third-party feature geometry.
#'   Manual accepts a registered character detail and requires the same source
#'   interval in both views. A standalone composition gives overview and detail
#'   two columns, in a 2:1 ratio.
#' @name fusion_magnify
#' @examples
#' seqs <- data.frame(accver = "A", length = 1000, sequence_unit = "bp")
#' overview <- ggchord(seqs) + geom_seq() + coord_collinear()
#' region <- data.frame(accver = "A", start = 100, end = 240)
#' fusion_magnify_auto(overview, region, expand = 40)
#'
#' detail <- overview + layout_sequences(windows = data.frame(
#'   instance_id = "A_window", accver = "A", start = 60, end = 280))
#' fusion_magnify_manual(overview, detail, region, connector = "border")
#'
#' cells <- flagship_fusion_example$base_cells
#' characters <- ggplot2::ggplot(cells,
#'   ggplot2::aes(alignment_column, accver, label = symbol)) +
#'   geom_residue(ggplot2::aes(fill = symbol)) +
#'   scale_fill_residue("DNA")
#' characters <- register_ggchord_alignment(characters, cells)
#' first <- min(cells$source_position, na.rm = TRUE)
#' bases <- data.frame(accver = cells$source_accver[1], start = first, end = first + 9)
#' fusion_magnify_auto(characters, bases, level = "residue")
NULL

#' @rdname fusion_magnify
#' @export
fusion_magnify_auto <- function(overview = NULL, region,
    detail_coord = coord_collinear(), expand = 0, layout = NULL,
    guides = "collect", ..., detail = NULL, from = NULL, to = NULL,
    level = "auto") {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  if (!is.null(detail) && !is.function(detail))
    ggchord_stop("fusion_magnify_auto(): detail must be a function of the derived plot")
  if (!is.numeric(expand) || length(expand) != 1L || !is.finite(expand) ||
      expand < 0 || expand != trunc(expand))
    ggchord_stop("fusion_magnify_auto(): expand must be a nonnegative integer")
  level <- ggchord_match_choice(level, c("auto", "feature", "residue"), "level")
  coord_explicit <- !missing(detail_coord)
  if (level == "residue" && coord_explicit)
    ggchord_stop("fusion_magnify_auto(): detail_coord applies to Feature detail")
  ggchord_magnify("auto", overview, detail, region, from, to, layout, guides,
    !missing(layout) || !missing(guides), list(...), detail_coord, expand, level, coord_explicit)
}

#' @rdname fusion_magnify
#' @export
fusion_magnify_manual <- function(overview = NULL, detail = NULL, region,
    layout = NULL, guides = "keep", ..., from = NULL, to = NULL) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  ggchord_magnify("manual", overview, detail, region, from, to, layout, guides,
    !missing(layout) || !missing(guides), list(...))
}

ggchord_magnify <- function(mode, overview, detail, region, from, to, layout,
    guides, root_options, styles, detail_coord = NULL, expand = NULL,
    level = NULL, coord_explicit = FALSE) {
  standalone <- !is.null(overview)
  if (standalone && (!is.null(from) || !is.null(to)))
    ggchord_stop("Fusion Magnify: use plots or from/to paths")
  if (!standalone && (is.null(from) || is.null(to) || root_options ||
      (mode == "manual" && !is.null(detail))))
    ggchord_stop("Fusion Magnify: supply from and to in + mode; set layout and guides on the Fusion")
  if (standalone && mode == "manual" && is.null(detail))
    ggchord_stop("fusion_magnify_manual(): supply a detail plot")
  if (standalone) {
    from <- "overview"
    to <- "detail"
  }
  allowed <- c("highlight", "fill", "alpha", "colour", "linewidth", "id", "connector", "linetype")
  if (length(styles) && (is.null(names(styles)) || anyDuplicated(names(styles)) ||
      any(!names(styles) %in% allowed)))
    ggchord_stop("Fusion Magnify: use named decoration parameters in ...")
  spec <- do.call(ggchord_magnify_spec, c(list(from = from, to = to, region = region), styles))
  spec$mode <- mode
  spec$detail_coord <- detail_coord
  spec$expand <- expand
  spec$level <- level
  spec$coord_explicit <- coord_explicit
  if (mode == "auto") spec$edit <- detail
  if (!standalone) return(spec)
  panels <- list(overview = overview)
  if (mode == "manual") panels$detail <- detail
  layout <- layout %||% layout_fusion_grid(ncol = 2, widths = c(2, 1))
  ggchord_fusion(panels, layout = layout, guides = guides) + spec
}

ggchord_fusion_leaf <- function(plot, path, caller) {
  plot <- ggchord_fusion_get(plot, path, caller)
  if (inherits(plot, "ggchord_fusion"))
    ggchord_stop(caller, ": select a leaf view")
  plot
}

ggchord_magnify_register <- function(plot, spec) {
  caller <- paste0("fusion_magnify_", spec$mode, "()")
  source <- ggchord_fusion_leaf(plot, spec$from, caller)
  if (spec$mode == "manual") {
    ggchord_fusion_leaf(plot, spec$to, caller)
    return(list(plot = plot, spec = spec))
  }
  if (length(spec$to) != 1L || spec$to_key %in% names(plot$panels))
    ggchord_stop(caller, ": to must be a new direct child ID")
  derived <- ggchord_magnify_derive(source, spec, caller)
  plot$panels[[spec$to_key]] <- list(plot = derived$plot, derived_from = spec$from,
    residue_recipe = derived$residue_recipe)
  list(plot = plot, spec = derived$spec)
}

ggchord_magnify_window_ids <- function(region, caller) {
  target <- region$to_instance %||% rep(NA_character_, nrow(region))
  for (i in which(is.na(target))) {
    id <- paste0("window-", i)
    while (id %in% target) id <- paste0(id, "-detail")
    target[i] <- id
  }
  if (anyDuplicated(target))
    ggchord_stop(caller, ": each detail window needs a distinct to_instance")
  target
}

ggchord_magnify_register_windows <- function(spec, windows) {
  spec$requested_region <- spec$region
  spec$region$from_instance <- windows$pred
  spec$region$to_instance <- windows$instance_id
  spec$detail_window <- windows
  spec
}

ggchord_magnify_residue_detail <- function(spec, recipe, caller, reason) {
  if (is.null(recipe$alphabet))
    ggchord_stop(caller, ": residue level requires real, source-mapped characters (",
      recipe$reason, "). Register a character alignment or supply the reference sequence and alphabet.")
  view <- ggchord_magnify_residue_plot(recipe)
  registered <- view$ggchord_alignment
  if (!is.null(spec$edit)) view <- spec$edit(view)
  if (!inherits(view, "ggplot") || inherits(view, "ggchord_fusion") ||
      !identical(view$ggchord_alignment, registered))
    ggchord_stop(caller, ": a residue detail editor must preserve the source cell registration")
  spec$edit <- NULL
  spec$effective_level <- "residue"
  spec$level_reason <- reason
  list(plot = view, spec = ggchord_magnify_register_windows(spec, recipe$windows))
}

ggchord_magnify_derive <- function(source, spec, caller) {
  if (!is.null(source$ggchord_alignment)) {
    if (spec$level == "feature" || spec$coord_explicit)
      ggchord_stop(caller, ": a registered character source supports residue detail; use Manual for Feature")
    recipe <- ggchord_magnify_residue_alignment(source, spec, caller)
    return(ggchord_magnify_residue_detail(spec, recipe, caller, "registered_source"))
  }
  if (!inherits(source, "ggchord") || !inherits(source$facet, "FacetNull") ||
      isTRUE(source$coordinates$ggchord_dotplot) || length(spec$from) > 1L)
    ggchord_stop(caller, ": Auto requires a direct, unfaceted sequence-path ggchord or a registered character view; use Manual")
  source <- ggchord_render_plot_copy(source)
  configured <- ggchord_configure_layout(ggchord_resolve_reference(source))
  reference <- ggchord_query_reference(configured, NULL)
  region <- spec$region
  target <- ggchord_magnify_window_ids(region, caller)
  if (length(Filter(Negate(is.null), c(region$from_facet, region$to_facet))))
    ggchord_stop(caller, ": facet selection is available in Manual")
  windows <- lapply(seq_len(nrow(region)), function(i) {
    r <- region[i, , drop = FALSE]
    local_start <- r$start - reference$source_origin + 1
    local_end <- r$end - reference$source_origin + 1
    selected <- ggchord_fusion_selector(r, "from_instance")
    hit <- which(reference$source_accver == r$accver &
      local_start >= reference$local_start & local_end <= reference$local_end &
      (is.null(selected) | reference$instance_id %in% selected))
    if (length(hit) != 1L)
      ggchord_stop(caller, ": region row ", i, " must select one fully visible source window; candidates: ",
        paste(reference$instance_id[reference$source_accver == r$accver], collapse = ", "))
    origin <- reference$source_origin[hit]
    lo <- max(reference$local_start[hit], local_start[hit] - spec$expand)
    hi <- min(reference$local_end[hit], local_end[hit] + spec$expand)
    parent <- reference$instance_id[hit]
    orientation <- configured$ggchord$sequence_layout$orientation
    orientation_key <- if (isTRUE(source$coordinates$ggchord_sankey)) sub("^(left|right):", "", parent) else parent
    data.frame(instance_id = target[i], accver = reference$accver[hit],
      start = lo, end = hi, pred = parent, source_origin = origin,
      source_accver = r$accver, sequence_unit = reference$sequence_unit[hit],
      orientation = unname(orientation[orientation_key]),
      expand_clipped = lo > local_start[hit] - spec$expand || hi < local_end[hit] + spec$expand)
  })
  windows <- do.call(rbind, windows)
  seqs <- configured$ggchord$data$seq_data
  claimed_unit <- region[["sequence_unit"]]
  if (!is.null(claimed_unit) && any(!is.na(claimed_unit) & claimed_unit != windows$sequence_unit))
    ggchord_stop(caller, ": region unit differs from the source reference")
  recipe <- if (spec$level != "feature" && !spec$coord_explicit &&
      (is.null(spec$edit) || spec$level == "residue"))
    ggchord_magnify_residue_reference(ggchord_magnify_reference_strings(source, seqs), windows)
  if (spec$level == "residue")
    return(ggchord_magnify_residue_detail(spec, recipe, caller, "explicit_residue"))
  spec$effective_level <- "feature"
  spec$level_reason <- if (spec$level == "feature") "explicit_feature" else
    recipe$reason %||% if (is.null(recipe)) "feature_editor_or_coord" else "unmeasured"
  if (is.null(recipe$alphabet)) recipe <- NULL
  coord <- spec$detail_coord
  if (!inherits(coord, "Coord") || isTRUE(coord$ggchord_dotplot) || isTRUE(coord$ggchord_sankey) ||
      !isTRUE(coord$ggchord_coord))
    ggchord_stop(caller, ": choose a Collinear, Chord, Circular or Radial Feature detail Coord, or use Manual")
  if (isTRUE(coord$ggchord_circular) && (nrow(windows) != 1L || windows$start != 1 ||
      windows$end != seqs$length[match(windows$accver, seqs$accver)] ||
      isTRUE(ggchord_focus_source(seqs)$origin[match(windows$accver, seqs$accver)] != 1)))
    ggchord_stop(caller, ": Circular detail requires a complete reference; use a linear detail for windows")
  if (nrow(windows) > 1L && !isTRUE(coord$ggchord_collinear) && anyDuplicated(windows$accver))
    ggchord_stop(caller, ": use Collinear for repeated windows of one reference")
  original <- source$ggchord$sequence_layout
  arrangement <- layout_sequences(order = windows$instance_id,
    windows = windows[c("instance_id", "accver", "start", "end", "orientation")],
    length_mode = original$length_mode %||% "proportional")
  parent_layout <- configured$ggchord$sequence_layout
  match_parent <- match(sub("^(left|right):", "", windows$pred), parent_layout$order)
  same_coord <- identical(class(source$coordinates)[1L], class(coord)[1L])
  compatible <- if (!same_coord) character() else if (isTRUE(coord$ggchord_circular)) {
    c("angular_gap", "radius")
  } else if (isTRUE(coord$ggchord_radial)) {
    c("inner_radius", "curvature")
  } else if (isTRUE(coord$ggchord_collinear)) {
    c("row_gap", "align", "curvature")
  } else c("angular_gap", "radius", "curvature")
  for (field in intersect(compatible, names(Filter(Negate(is.null), original)))) {
    values <- switch(field, angular_gap = parent_layout$gap * 360,
      row_gap = original$row_gap, parent_layout[[field]])
    arrangement[[field]] <- if (field %in% c("angular_gap", "curvature", "radius"))
      unname(values[match_parent]) else if (field == "row_gap")
        unname(parent_layout$row_gap[match_parent]) else values
  }
  spec$omitted_layout <- setdiff(names(Filter(Negate(is.null), original)),
    c("order", "orientation", "windows", "length_mode", compatible))
  detail <- ggchord_render_plot_copy(source)
  detail$coordinates <- ggproto(NULL, coord)
  detail$ggchord$sequence_layout <- arrangement
  detail$ggchord$layout_configured <- NULL
  detail$ggchord$dotplot_layout <- NULL
  detail$ggchord$position_scales_prepared <- NULL
  detail$ggchord$magnify_instances <- windows
  detail$ggchord$magnify_reference <- reference
  axis <- detail$scales$get_scales("seq_position")
  axis <- if (is.null(axis)) scale_seq_continuous(coordinates = "source") else axis$clone()
  axis$ggchord_position_coordinates <- "source"
  selected <- axis$ggchord_axis_sequences
  if (length(selected)) {
    selected <- windows$instance_id[windows$pred %in% selected | windows$accver %in% selected]
    axis$ggchord_axis_sequences <- axis$ggchord_axis_guide$sequences <- selected
  }
  detail$scales$add(axis)
  omitted <- list()
  for (i in seq_along(detail$layers)) {
    layer <- detail$layers[[i]]
    p <- layer$ggchord_params
    if (is.null(layer$ggchord_type))
      ggchord_stop(caller, ": unsupported layer '", layer$ggchord_layer_id, "'; use Manual")
    if (identical(p$label_kind, "sequence") && identical(p$anchor, "center") &&
        !isTRUE(coord$ggchord_circular)) {
      omitted[[layer$ggchord_layer_id]] <- "circular_center_decoration"
      detail$layers[i] <- list(NULL)
    }
  }
  dropped <- names(omitted)
  detail$layers <- Filter(Negate(is.null), detail$layers)
  detail$ggchord$scoped_scales <- Filter(function(scope) length(scope$layers) || length(scope$tracks),
    lapply(detail$ggchord$scoped_scales, function(scope) {
      if (!is.null(scope$layers)) scope$layers <- setdiff(scope$layers, dropped)
      scope
    }))
  before <- detail
  if (!is.null(spec$edit)) detail <- spec$edit(detail)
  if (!inherits(detail, "ggchord") || !inherits(detail$facet, "FacetNull") ||
      !identical(detail$coordinates, before$coordinates) ||
      !identical(detail$ggchord$data, before$ggchord$data) ||
      !identical(detail$ggchord$sequence_layout, arrangement))
    ggchord_stop(caller, ": detail editor must preserve the derived reference, Coord and windows")
  spec <- ggchord_magnify_register_windows(spec, windows)
  spec$omissions <- omitted
  spec$edit <- NULL
  scale_channels <- unique(unlist(lapply(c(before$scales$scales, detail$scales$scales),
    `[[`, "aesthetics")))
  scale_changed <- scale_channels[vapply(scale_channels, function(key) {
    old <- before$scales$get_scales(key)
    new <- detail$scales$get_scales(key)
    if (is.null(old) || is.null(new)) return(!identical(old, new))
    fields <- function(sc) {
      value <- as.list(sc)
      value[c("super", "range")] <- NULL
      value
    }
    !identical(fields(old), fields(new))
  }, logical(1))]
  unchanged <- Filter(function(layer) any(vapply(detail$layers,
    identical, logical(1), y = layer)), before$layers)
  unchanged_ids <- vapply(unchanged, function(layer) layer$ggchord_layer_id, character(1))
  detail$ggchord$magnify_parent <- list(unchanged_layers = unchanged_ids, scale_changed = scale_changed,
    scopes = before$ggchord$scoped_scales, position_scales = before$ggchord$position_scales,
    preset = before$ggchord$palette_preset, mapping = before$mapping)
  list(plot = detail, spec = spec, residue_recipe = recipe)
}

ggchord_magnify_source <- function(layer) {
  region <- layer$region
  region$component_type <- "decoration"
  region$magnify_mode <- layer$mode
  region$layer_id <- layer$id
  region$composition_id <- ""
  region$source_row <- seq_len(nrow(region))
  region$from_panel_path <- rep(list(layer$from), nrow(region))
  region$to_panel_path <- rep(list(layer$to), nrow(region))
  if (identical(layer$mode, "auto")) {
    windows <- layer$detail_window
    region$detail_start <- windows$start + windows$source_origin - 1
    region$detail_end <- windows$end + windows$source_origin - 1
    region$detail_local_start <- windows$start
    region$detail_local_end <- windows$end
    region$source_origin <- windows$source_origin
    region$expand <- layer$expand
    region$requested_level <- layer$level
    region$effective_level <- layer$effective_level
    region$level_reason <- layer$level_reason
    region$window_status <- ifelse(windows$expand_clipped, "expand_clipped", "complete")
    region$omitted_layers <- rep(list(layer$omissions), nrow(region))
    region$omitted_layout <- rep(list(layer$omitted_layout), nrow(region))
    region$omitted_observations <- rep(list(layer$omitted_observations), nrow(region))
    region$scale_inheritance <- rep(list(layer$scale_inheritance), nrow(region))
    region$statistic_inheritance <- if (identical(layer$effective_level, "residue"))
      "not_applicable" else "source_domain"
  }
  region
}

ggchord_render_zoom <- function(layer, panels) {
  connector <- layer$connector %||% "interval"
  line <- layer$linetype %||% if (connector == "border") "dashed" else "solid"
  outline <- gpar(fill = alpha(layer$fill, layer$alpha),
    col = layer$colour, lwd = layer$linewidth * ggplot2::.pt, lty = line)
  resolve <- function(region, panel, side, row) {
    instance <- ggchord_fusion_selector(region, paste0(side, "_instance"))
    facet <- ggchord_fusion_selector(region, paste0(side, "_facet"))
    ggchord_resolve_panel_interval(panels, region$accver, region$start, region$end,
      panel, instance, paste0("Fusion Magnify '", layer$id, "' row ", row, " ", side, ": "), facet)
  }
  source_overlay <- function(front) {
    half <- front$cell_half_height_mm %||% rep(0.7, nrow(front))
    x <- front$x + front$normal_x * half
    y <- front$y + front$normal_y * half
    xx <- front$x - front$normal_x * half
    yy <- front$y - front$normal_y * half
    polygonGrob(c(x, rev(xx)), c(y, rev(yy)), default.units = "mm",
      gp = gpar(fill = alpha(layer$fill, layer$alpha),
        col = layer$colour, lwd = layer$linewidth * ggplot2::.pt))
  }
  regions <- ggchord_magnify_source(layer)
  parts <- lapply(seq_len(nrow(regions)), function(i) {
    region <- regions[i, , drop = FALSE]
    q <- resolve(region, layer$from_key, "from", i)
    s <- resolve(region, layer$to_key, "to", i)
    ggchord_check_character_origin(q, s, "Fusion Magnify")
    if (q$unit != s$unit) ggchord_stop("Fusion Magnify: source and detail must use the same unit")
    ggchord_check_interval_unit(region, q$unit, "Fusion Magnify")
    ggchord_fusion_reference_match(q, s, "Fusion Magnify")
    region$from_facet <- list(q$facet_key)
    region$to_facet <- list(s$facet_key)
    qf <- ggchord_fusion_front(panels[[q$panel]], q)
    sf <- ggchord_fusion_front(panels[[s$panel]], s)
    region$from_instance <- q$instance
    region$to_instance <- s$instance
    region$from_panel_path <- list(layer$from)
    region$to_panel_path <- list(layer$to)
    region$sequence_unit <- q$unit
    region$layer_id <- layer$id
    region$source_row <- i
    region$connector <- connector
    region$component_type <- "decoration"
    region$composition_id <- ""
    path <- grob <- NULL
    if (connector == "interval") {
      xy <- ggchord_fusion_path(qf, sf, curvature = 0)
      path <- data.frame(x = xy[, 1L], y = xy[, 2L], layer_id = layer$id,
        source_row = i, component_type = "decoration", composition_id = "", units = "mm")
      grob <- polygonGrob(path$x, path$y, default.units = "mm", gp = outline)
    }
    list(region = region, path = path, grob = grob, front = qf,
      source_facet = q$facet %||% 1L, target_facet = s$facet %||% 1L,
      overlay = if (layer$highlight) setNames(list(source_overlay(qf)), q$panel) else list())
  })
  if (connector == "interval" || !length(parts)) return(parts)

  # Draw one envelope for all source intervals feeding the same two facets.
  # Preserve each source highlight and interval in the exported region table.
  groups <- split(seq_along(parts), vapply(parts, function(part)
    paste(part$source_facet, part$target_facet, sep = ":"), character(1)))
  bounds <- function(panel, facet) {
    context <- panel$context
    result <- context$panel_bounds_mm
    leaf <- if (length(context$facet_contexts) >= facet) context$facet_contexts[[facet]]
    if (!is.null(leaf))
      result <- leaf$panel_bounds_mm - (context$device_panel_bounds_mm - result)
    matrix(result, nrow = 2L)
  }
  for (where in groups) {
    first <- parts[[where[1L]]]
    from <- bounds(panels[[layer$from_key]], first$source_facet)
    to <- bounds(panels[[layer$to_key]], first$target_facet)
    direction <- colMeans(from) - colMeans(to)
    separating <- pmax(from[1L, ] - to[2L, ], to[1L, ] - from[2L, ], 0)
    normal <- if (any(separating > 0)) which.max(separating) else which.max(abs(direction))
    side <- if (direction[normal] > 0) 2L else 1L
    along <- 3L - normal
    tips <- do.call(rbind, lapply(where, function(i) {
      front <- parts[[i]]$front
      data.frame(x = front$x, y = front$y, source_row = i)
    }))
    pick <- c(which.min(tips[[along]]), which.max(tips[[along]]))
    if (pick[1L] == pick[2L]) pick <- c(which.min(tips[[normal]]), which.max(tips[[normal]]))
    roots <- tips[pick, , drop = FALSE]
    edge <- matrix(rep(to[side, ], each = 2L), nrow = 2L)
    edge[, along] <- to[, along]
    xy <- rbind(as.matrix(roots[c("x", "y")]), edge)
    path <- data.frame(x = xy[, 1L], y = xy[, 2L], layer_id = layer$id,
      source_row = rep(roots$source_row, 2L), path_part = rep(1:2, 2L),
      component_type = "decoration", composition_id = "", units = "mm")
    parts[[where[1L]]]$path <- path
    parts[[where[1L]]]$grob <- segmentsGrob(roots$x, roots$y, edge[, 1L], edge[, 2L],
      default.units = "mm", gp = outline)
    for (i in where[-1L]) parts[[i]]$grob <- nullGrob()
  }
  parts
}
