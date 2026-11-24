# A recipe carries existing strings or registered observations. No alignment
# columns are inferred from interval links.
ggchord_magnify_reference_strings <- function(plot, seqs) {
  for (layer in plot$layers) {
    if (!identical(layer$ggchord_type, "seq")) next
    claim <- layer$ggchord_input_data
    if (!is.data.frame(claim)) next
    at <- match(claim$accver, seqs$accver)
    for (key in intersect(c("sequence", "alphabet", "reference_namespace"), names(claim))) {
      supplied <- as.character(claim[[key]])
      current <- as.character(seqs[[key]] %||% rep(NA_character_, nrow(seqs)))
      known <- !is.na(supplied)
      conflict <- known & !is.na(current[at]) & current[at] != supplied
      if (any(conflict))
        ggchord_stop("Auto Magnify: conflicting ", key, " in registered reference declarations")
      current[at[known]] <- supplied[known]
      seqs[[key]] <- current
    }
  }
  seqs
}

ggchord_magnify_residue_reference <- function(seqs, windows) {
  at <- match(windows$accver, seqs$accver)
  strings <- seqs[["sequence"]][at]
  if (is.null(strings) || anyNA(strings) || any(!nzchar(strings)))
    return(list(reason = "missing_reference_characters"))
  if (!is.character(strings) || any(nchar(strings) != seqs$length[at]) ||
      any(grepl("[-.[:space:]]", strings)))
    return(list(reason = "invalid_reference_characters"))
  declared <- seqs[["alphabet"]][at] %||% rep(NA_character_, length(at))
  unit <- windows$sequence_unit
  has_u <- grepl("U", strings, ignore.case = TRUE)
  has_t <- grepl("T", strings, ignore.case = TRUE)
  inferred <- ifelse(unit == "aa", "AA",
    ifelse(unit == "bp" & xor(has_u, has_t), ifelse(has_u, "RNA", "DNA"), NA_character_))
  alphabet <- ifelse(is.na(declared), inferred, toupper(declared))
  if (anyNA(alphabet)) return(list(reason = "reference_alphabet_required"))
  if (any(!alphabet %in% c("DNA", "RNA", "AA")) ||
      any(unit != ifelse(alphabet == "AA", "aa", "bp")))
    return(list(reason = "reference_alphabet_unit_conflict"))
  if (length(unique(alphabet)) != 1L) return(list(reason = "mixed_residue_alphabets"))
  pattern <- switch(alphabet[1L], DNA = "[^ACGTRYSWKMBDHVN]",
    RNA = "[^ACGURYSWKMBDHVN]", AA = "[^ACDEFGHIKLMNPQRSTVWYBXZJUO*]")
  if (any(grepl(pattern, strings, ignore.case = TRUE)))
    return(list(reason = "invalid_reference_alphabet_characters"))
  root <- ggchord_focus_source(seqs)
  list(alphabet = alphabet[1L], strings = strings, windows = windows, width = max(windows$end - windows$start + 1),
    reference_length = root$length[at], namespace = seqs[["reference_namespace"]][at],
    circular = seqs$topology[at] %in% "circular", row = at,
    provenance = attr(seqs, "ggchord_import") %||% attr(seqs, "ggchord_conversion"))
}

ggchord_magnify_residue_alignment <- function(source, spec, caller) {
  layout <- ggchord_alignment_layout(ggplot_build(source), source$ggchord_alignment)
  data <- layout$character_data
  ids <- ggchord_magnify_window_ids(spec$region, caller)
  windows <- vector("list", nrow(spec$region))
  cells <- lapply(seq_len(nrow(spec$region)), function(i) {
    region <- spec$region[i, , drop = FALSE]
    instance <- ggchord_fusion_selector(region, "from_instance")
    matches <- ggchord_alignment_interval(list(layout = layout), spec$from_key,
      region$accver, region$start, region$end, instance)
    matches <- Filter(function(match) !isTRUE(match$fragment), matches)
    if (length(matches) != 1L)
      ggchord_stop(caller, ": residue interval must match one source-backed character row; candidates: ",
        paste(unique(data$display_row[data$source_accver == region$accver]), collapse = ", "))
    endpoint <- matches[[1L]]
    ggchord_check_interval_unit(region, endpoint$unit, caller)
    row <- data[data$display_row == endpoint$instance, , drop = FALSE]
    observed <- sort(row$source_position[!row$is_gap])
    block <- split(observed, cumsum(c(TRUE, diff(observed) != 1)))
    chosen <- Filter(function(x) min(x) <= region$start && max(x) >= region$end, block)
    lo <- max(region$start - spec$expand, min(chosen[[1L]]))
    hi <- min(region$end + spec$expand, max(chosen[[1L]]))
    boundaries <- range(row$alignment_column[!row$is_gap & row$source_position %in% c(lo, hi)])
    row <- row[row$alignment_column >= boundaries[1L] & row$alignment_column <= boundaries[2L], , drop = FALSE]
    row$parent_alignment_row <- as.character(row$accver)
    row$accver <- ids[i]
    windows[[i]] <<- data.frame(instance_id = ids[i], accver = row$parent_alignment_row[1L],
      start = lo, end = hi, pred = endpoint$instance, source_origin = 1,
      source_accver = region$accver, sequence_unit = endpoint$unit,
      orientation = if (row$strand[1L] == "+") 1 else -1,
      expand_clipped = lo > region$start - spec$expand || hi < region$end + spec$expand)
    row[c("display_column", "display_row")] <- NULL
    row
  })
  cells <- ggchord_rbind_fill(cells)
  alphabet <- unique(cells$alphabet)
  if (length(alphabet) != 1L) ggchord_stop(caller, ": split mixed alphabets across residue views")
  attr(cells, "ggchord_conversion") <- attr(source$ggchord_alignment, "ggchord_conversion")
  list(cells = cells, windows = do.call(rbind, windows), alphabet = alphabet,
    width = diff(range(cells$alignment_column)) + 1)
}

ggchord_magnify_residue_plot <- function(recipe) {
  cells <- recipe$cells
  if (is.null(cells)) {
    windows <- recipe$windows
    cells <- lapply(seq_len(nrow(windows)), function(i) {
      window <- windows[i, , drop = FALSE]
      observed <- substring(recipe$strings[i], window$start, window$end)
      reversed <- isTRUE(window$orientation == -1)
      if (reversed) {
        observed <- paste(rev(strsplit(observed, "", fixed = TRUE)[[1L]]), collapse = "")
        observed <- ggchord_alignment_complement(observed, recipe$alphabet)
      }
      source_start <- window$source_origin + (if (reversed) window$end else window$start) - 1
      source <- data.frame(accver = window$instance_id, source_accver = window$source_accver,
        start = source_start, strand = if (reversed) "-" else "+",
        source_length = recipe$reference_length[i])
      if (window$source_origin == 1 && isTRUE(nchar(recipe$strings[i]) == recipe$reference_length[i]))
        source$sequence <- recipe$strings[i]
      out <- as_ggchord_alignment(setNames(observed, window$instance_id),
        alphabet = recipe$alphabet, source_data = source, circular = recipe$circular[i] &&
          is.finite(recipe$reference_length[i]))
      out$.source_row <- recipe$row[i]
      if (length(recipe$namespace)) out$reference_namespace <- recipe$namespace[i]
      out
    })
    cells <- ggchord_rbind_fill(cells)
    attr(cells, "ggchord_conversion") <- list(format = "reference windows",
      source_windows = windows, reference_input = recipe$provenance)
  }
  cells$accver <- factor(cells$accver, levels = rev(recipe$windows$instance_id))
  labels <- with(recipe$windows, setNames(paste0(source_accver, " ",
    start + source_origin - 1, "-", end + source_origin - 1), instance_id))
  positions <- if (nrow(recipe$windows) == 1L)
    cells$source_position[match(unique(cells$alignment_column), cells$alignment_column)]
  x_scale <- if (is.null(positions)) ggplot2::scale_x_continuous(breaks = NULL,
    expand = ggplot2::expansion(add = .5)) else {
    columns <- unique(cells$alignment_column)
    ticks <- function(limits) {
      wanted <- round(scales::breaks_extended(n = 5)(range(columns)))
      sort(unique(c(min(columns), wanted[wanted %in% columns], max(columns))))
    }
    ggplot2::scale_x_continuous(breaks = ticks,
      labels = function(x) ifelse(is.na(positions[match(x, columns)]), "", positions[match(x, columns)]),
      expand = ggplot2::expansion(add = .5))
  }
  plot <- ggplot(cells, aes(x = alignment_column, y = accver, label = symbol)) +
    geom_residue(aes(fill = toupper(symbol)), id = "residues", cell_colour = "grey60") +
    scale_fill_residue(recipe$alphabet) + x_scale +
    ggplot2::scale_y_discrete(labels = labels, expand = ggplot2::expansion(add = .55)) +
    ggplot2::labs(x = if (!is.null(positions)) "Source position") +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::theme_minimal() + theme(axis.title.y = element_blank(), panel.grid = element_blank())
  register_ggchord_alignment(plot, cells)
}

ggchord_magnify_residue_fit <- function(recipe, dimensions) {
  # Rule out impossible widths before materialising a long reference window.
  defaults <- GeomResidue$default_aes
  symbols <- strsplit(switch(recipe$alphabet, DNA = "ACGTRYSWKMBDHVN",
    RNA = "ACGURYSWKMBDHVN", AA = "ACDEFGHIKLMNPQRSTVWYBXZJUO*"), "")[[1L]]
  glyphs <- data.frame(label = symbols, size = as.numeric(get_expr(defaults$size)),
    family = as.character(get_expr(defaults$family)), fontface = as.numeric(get_expr(defaults$fontface)),
    lineheight = as.numeric(get_expr(defaults$lineheight)))
  metrics <- ggchord_residue_glyph_metrics(glyphs) * 1.15
  if (recipe$width * min(metrics[, "width"]) > dimensions[1L] ||
      nrow(recipe$windows) * min(metrics[, "height"]) > dimensions[2L])
    return(NULL)
  plot <- ggchord_magnify_residue_plot(recipe)
  built <- ggplot_build(plot)
  table <- ggplot_gtable(built)
  context <- ggchord_panel_context(table, built)
  if (!identical(context$status, "measured")) return(NULL)
  data <- built$data[[1L]]
  coord <- built$layout$coord
  coord_params <- built$layout$panel_params[[1L]]
  lower <- coord$transform(data.frame(x = data$xmin, y = data$ymin), coord_params)
  upper <- coord$transform(data.frame(x = data$xmax, y = data$ymax), coord_params)
  fits <- ggchord_residue_glyph_fit(data, abs(upper$x - lower$x) * context$panel_mm[1L],
    abs(upper$y - lower$y) * context$panel_mm[2L])
  if (all(fits)) plot
}

ggchord_fusion_residue_levels <- function(plot, rectangles, area, shift) {
  changed <- FALSE
  for (id in names(plot$zooms)) {
    spec <- plot$zooms[[id]]
    child <- plot$panels[[spec$to_key]]
    recipe <- child$residue_recipe
    if (is.null(recipe)) next
    key <- spec$to_key
    rect <- rectangles[[key]]
    size <- c(diff(rect[1:2]), diff(rect[3:4])) * area
    center <- c(mean(rect[1:2]), mean(rect[3:4])) * area + shift
    alternative <- if (!isTRUE(child$plot$ggchord_fusion_level_lock)) local({
      vp <- grid::viewport(x = center[1L], y = center[2L],
        width = size[1L], height = size[2L], default.units = "mm")
      pushViewport(vp)
      on.exit(popViewport())
      ggchord_magnify_residue_fit(recipe, size)
    })
    plot$panels[[key]]$residue_recipe <- NULL
    if (is.null(alternative)) {
      spec$level_reason <- if (isTRUE(child$plot$ggchord_fusion_level_lock))
        "explicit_visual_scope" else "glyphs_do_not_fit"
    } else {
      plot$panels[[key]]$plot <- alternative
      spec$effective_level <- "residue"
      spec$level_reason <- "glyphs_fit"
      changed <- TRUE
    }
    plot$zooms[[id]] <- spec
  }
  list(plot = plot, changed = changed)
}
