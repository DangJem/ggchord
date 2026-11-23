# A source entity can own several display rows. The legacy engines consume
# unique private keys, while public data and scales retain biological IDs.
ggchord_instance_param <- function(value, instances) {
  if (is.null(value) || is.unit(value)) return(value)
  if (is.list(value) && all(names(value) %in% c("source", "instance"))) {
    source <- value$source
    local <- value$instance
    if (anyDuplicated(names(source)) || anyDuplicated(names(local)) ||
        any(!names(source) %in% instances$accver) ||
        any(!names(local) %in% instances$instance_id))
      ggchord_stop("Unknown or repeated source/instance parameter keys")
    out <- vector("list", nrow(instances))
    for (i in seq_len(nrow(instances))) {
      out[[i]] <- list(local[[instances$instance_id[i]]] %||% source[[instances$accver[i]]])
    }
    out <- lapply(out, `[[`, 1L)
    names(out) <- instances$key
    out <- Filter(Negate(is.null), out)
    if (all(lengths(out) == 1L)) return(unlist(out, use.names = TRUE))
    return(out)
  }
  if (is.null(names(value))) return(value)
  keys <- unique(c(instances$accver, instances$instance_id))
  if (!all(names(value) %in% keys)) return(value)
  out <- list()
  for (i in seq_len(nrow(instances))) {
    available <- c(instances$instance_id[i], instances$accver[i])
    key <- available[available %in% names(value)][1L]
    if (!is.na(key)) out[[instances$key[i]]] <- value[[key]]
  }
  if (is.list(value)) out else unlist(out, use.names = TRUE)
}

ggchord_instance_position <- function(position, instances) {
  if (!inherits(position, "Position")) return(position)
  changes <- list()
  for (field in c("offset", "width", "gap", "spacing")) if (!is.null(position[[field]]))
    changes[[field]] <- ggchord_instance_param(position[[field]], instances)
  if (!is.null(position$base_position)) changes$base_position <-
    ggchord_instance_position(position$base_position, instances)
  do.call(ggproto, c(list(NULL, position), changes))
}

ggchord_instance_expand <- function(input, instances, reference = FALSE, location = NULL) {
  if (is.null(input)) return(list(data = NULL, rows = integer()))
  ggchord_require_columns(input, "accver", "display instances")
  point <- isTRUE(location %in% c("site", "cut"))
  interval <- !point && all(c("start", "end") %in% names(input))
  rows <- integer(); targets <- integer()
  for (i in seq_len(nrow(input))) {
    at <- which(instances$accver == as.character(input$accver[i]))
    selected <- input$instance_id[i] %||% input$display_instance_id[i]
    if (!is.null(selected) && !is.na(selected)) {
      at <- at[instances$instance_id[at] == selected]
      if (!length(at)) ggchord_stop("No display instance matches input row ", i)
    }
    if ("window_start" %in% names(instances)) {
      low <- if (interval) input$start[i] else input$pos[i] %||% input$x[i]
      high <- if (interval) input$end[i] else low
      if (!is.null(low)) at <- at[is.finite(low) & is.finite(high) &
        high >= instances$window_start[at] & low <= instances$window_end[at]]
    }
    rows <- c(rows, rep(i, length(at))); targets <- c(targets, at)
  }
  out <- input[rows, , drop = FALSE]
  out$accver <- instances$key[targets]
  out$display_instance_id <- instances$instance_id[targets]
  if ("window_start" %in% names(instances)) {
    from <- instances$window_start[targets]
    to <- instances$window_end[targets]
    fields <- c("pos", "x", "xmin", "xmax", "cut_top", "cut_bottom",
      "recognition_start", "recognition_end")
    if (!point) fields <- c("start", "end", fields)
    for (field in intersect(fields, names(out))) {
      boundary <- field %in% c("xmin", "xmax", "cut_top", "cut_bottom")
      out[[paste0(".window_source_", field)]] <- out[[field]]
      if (field %in% c("cut_top", "cut_bottom")) {
        observed <- out[[field]] >= from - 1 & out[[field]] <= to
        out[[field]] <- ifelse(observed, out[[field]] - from + 1, NA_real_)
      } else {
        out[[field]] <- pmax(from - boundary, pmin(to, out[[field]])) - from + 1
      }
    }
    if (identical(location, "cut") && "cut_status" %in% names(out)) {
      out$.window_source_cut_status <- out$cut_status
      out$cut_status <- ggchord_annotation_cut_status(out)
    }
    if (reference && "length" %in% names(out)) out$length <- as.numeric(to - from + 1L)
    nested <- c("segments", ".feature_runs", ".feature_boundaries", "cleavage_arrows")
    if (!point && nrow(out) && any(nested %in% names(out))) {
      parts <- lapply(seq_len(nrow(out)), function(i) {
        ggchord_focus_feature_parts(out[i, , drop = FALSE], from[i], to[i],
          origin = 1, caller = "layout_sequences(windows)")
      })
      keep <- !vapply(parts, is.null, logical(1))
      out <- if (any(keep)) ggchord_rbind_fill(parts[keep]) else out[FALSE, , drop = FALSE]
      rows <- rows[keep]
    }
  }
  list(data = out, rows = rows)
}

ggchord_route_instances <- function(input, instances, between, layer_id, mirrored = FALSE, source_rows = seq_len(nrow(input))) {
  if (is.null(input)) return(list(data = NULL, report = data.frame()))
  ggchord_require_columns(input, c("qaccver", "saccver"), "instance routing")
  if (!is.null(between)) {
    if (length(between) != 2L || anyNA(between)) ggchord_stop("between needs exactly two display instances")
    if (is.numeric(between)) {
      if (any(between != floor(between)) || any(!between %in% instances$declared_order)) ggchord_stop("Unknown between index")
      between <- instances$instance_id[match(between, instances$declared_order)]
    }
    if (any(!between %in% instances$instance_id) || anyDuplicated(between)) ggchord_stop("between names unknown or repeated display instances")
  }
  report <- list(); out <- input
  omitted <- logical(nrow(input))
  out$qaccver <- as.character(out$qaccver); out$saccver <- as.character(out$saccver)
  for (i in seq_len(nrow(input))) {
    q <- which(instances$accver == as.character(input$qaccver[i]))
    s <- which(instances$accver == as.character(input$saccver[i]))
    if (mirrored) { q <- q[instances$side[q] == "left"]; s <- s[instances$side[s] == "right"] }
    if ("window_start" %in% names(instances)) {
      endpoint <- function(at, side) {
        low <- (input[[paste0(side, "start")]] %||% input[[paste0(side, "pos")]])[i]
        high <- (input[[paste0(side, "end")]] %||% input[[paste0(side, "pos")]])[i]
        at[pmax(low, high) >= instances$window_start[at] &
          pmin(low, high) <= instances$window_end[at]]
      }
      q <- endpoint(q, "q")
      s <- endpoint(s, "s")
    }
    if (!is.null(input$qinstance) && !is.na(input$qinstance[i]))
      q <- q[instances$instance_id[q] == input$qinstance[i]]
    if (!is.null(input$sinstance) && !is.na(input$sinstance[i]))
      s <- s[instances$instance_id[s] == input$sinstance[i]]
    candidates <- expand.grid(q = q, s = s)
    if (!is.null(between)) candidates <- candidates[
      instances$instance_id[candidates$q] %in% between & instances$instance_id[candidates$s] %in% between, , drop = FALSE]
    if (!is.null(between) && input$qaccver[i] == input$saccver[i])
      candidates <- candidates[instances$instance_id[candidates$q] == between[1], , drop=FALSE]
    if (!nrow(candidates) && "window_start" %in% names(instances) &&
        is.null(between) && is.null(input$qinstance) && is.null(input$sinstance)) {
      omitted[i] <- TRUE
      next
    }
    if (nrow(candidates) != 1L) ggchord_stop("HSP row ", i,
      " has an absent or ambiguous instance route; specify between = c(instance1, instance2)")
    qi <- candidates$q; si <- candidates$s
    out$qaccver[i] <- instances$key[qi]; out$saccver[i] <- instances$key[si]
    if ("window_start" %in% names(instances)) {
      for (name in c("q", "s")) {
        window <- if (name == "q") qi else si
        from <- instances$window_start[window]
        to <- instances$window_end[window]
        for (field in intersect(paste0(name, c("start", "end", "pos")), names(input)))
          out[[field]][i] <- pmax(from, pmin(to, input[[field]][i])) - from + 1
      }
    }
    stable <- "source_hsp_id" %in% names(input)
    source_id <- if (stable) as.character(input$source_hsp_id[i]) else paste0(layer_id, ":source:", i)
    if (is.na(source_id) || !nzchar(source_id)) ggchord_stop("source_hsp_id must be nonempty when supplied")
    report[[i]] <- data.frame(layer_id = layer_id, input_row = source_rows[i],
      source_hsp_id = source_id, source_scope = if (stable) "declared" else "layer",
      displayed_relation_id = paste0(layer_id, ":display:", i),
      qaccver = as.character(input$qaccver[i]), saccver = as.character(input$saccver[i]),
      q_instance = instances$instance_id[qi], s_instance = instances$instance_id[si],
      repeated_display = FALSE, stringsAsFactors = FALSE)
  }
  list(data = out[!omitted, , drop = FALSE], rows = which(!omitted),
    report = ggchord_rbind_fill(report) %||% data.frame())
}

ggchord_reference_table <- function(data, caller) {
  ggchord_require_columns(data, c("accver", "length"), caller)
  if (!nrow(data)) ggchord_stop(caller, ": the reference table is empty")
  out <- as_ggchord_seq(data)
  ggchord_require_columns(out, c("accver", "length"), caller)
  if (anyDuplicated(out$accver) || anyNA(out$accver) || any(!nzchar(trimws(out$accver))))
    ggchord_stop(caller, ": reference IDs must be unique and nonempty")
  if (!is.numeric(out$length) || any(!is.finite(out$length) | out$length < 1 |
      out$length != trunc(out$length)))
    ggchord_stop(caller, ": length must contain complete positive integer sequence lengths")
  if (is.null(out$sequence_unit)) out$sequence_unit <- rep("bp", nrow(out))
  if (is.null(out$topology)) out$topology <- rep(NA_character_, nrow(out))
  out
}

ggchord_check_reference_claim <- function(reference, claim, facet_fields, subset, caller) {
  at <- match(claim$accver, reference$accver)
  if (anyNA(at) || (!subset && length(at) != nrow(reference)))
    ggchord_stop(caller, ": conflicting reference accession sets")
  original <- ggchord_focus_source(reference)
  proposed <- ggchord_focus_source(claim)
  for (field in c("accver", "origin", "length")) {
    left <- original[[field]]
    right <- proposed[[field]]
    if (!identical(unname(left[at]), unname(right)))
      ggchord_stop(caller, ": conflicting reference source ", field)
  }
  fields <- unique(c("length", "sequence_unit", "topology", facet_fields))
  for (name in fields) {
    left <- reference[[name]]
    right <- claim[[name]]
    # A display subset may omit ancillary facet columns and inherit them only
    # for reference membership. It cannot redefine any value that it supplies.
    if (subset && is.null(right) && name %in% facet_fields) next
    if (!isTRUE(all.equal(left[at], right, check.attributes = FALSE)))
      ggchord_stop(caller, ": conflicting reference ", name)
  }
  invisible(at)
}

ggchord_resolve_reference <- function(plot) {
  if (isTRUE(plot$ggchord$reference_resolved)) return(plot)
  reference <- if (!is_waiver(plot$data)) plot$data
  global <- !is.null(reference)
  facet <- plot$facet$params
  fields <- names(c(facet$facets, facet$rows, facet$cols))
  for (layer in plot$layers) {
    if (!identical(layer$ggchord_params$type, "seq")) next
    input <- layer$ggchord_input_data
    if (!is.data.frame(input)) next
    claim <- ggchord_reference_table(input, "geom_seq()")
    if (is.null(reference)) reference <- claim else
      ggchord_check_reference_claim(reference, claim, fields, global, "geom_seq()")
  }
  if (is.null(reference)) ggchord_stop(
    "ggchord(): supply a complete reference in ggchord(data = ...) or geom_seq(data = ...)")
  reference <- reference[order(reference$accver, method = "radix"), , drop = FALSE]
  rownames(reference) <- NULL
  plot$ggchord$data$seq_data <- reference
  plot$ggchord$complete_reference_ids <- reference$accver
  plot$ggchord$reference_resolved <- TRUE
  plot
}

ggchord_sequence_reference_data <- function(plot) {
  ggchord_resolve_reference(plot)$ggchord$data$seq_data
}

ggchord_prepare_instances <- function(plot, seqs) {
  if (isTRUE(plot$ggchord$instances_prepared)) return(plot)
  arrangement <- plot$ggchord$sequence_layout
  windows <- arrangement$windows
  order <- arrangement$order %||% seqs$accver
  configured <- any(vapply(plot$layers, function(l) !is.null(l$ggchord_params$between), logical(1)))
  mirrored <- isTRUE(plot$coordinates$ggchord_sankey)
  if (is.null(windows) && !mirrored && !configured) return(plot)
  refs <- windows$accver %||% order
  ids <- windows$instance_id %||% order
  if (mirrored) {
    ids <- c(paste0("left:", ids), paste0("right:", ids))
    refs <- rep(refs, 2L)
  }
  side <- if (mirrored) rep(c("left", "right"), each = length(ids) / 2L) else rep(NA_character_, length(ids))
  instances <- data.frame(instance_id = ids, accver = refs,
    declared_order = seq_along(ids), key = ids, side = side, stringsAsFactors = FALSE)
  original <- plot$ggchord$data
  virtual <- seqs[match(instances$accver, seqs$accver), , drop = FALSE]
  if (!is.null(windows)) {
    at <- rep(seq_len(nrow(windows)), if (mirrored) 2L else 1L)
    instances$window_start <- windows$start[at]
    instances$window_end <- windows$end[at]
    if (!is.null(windows$ring_id)) instances$ring_id <- windows$ring_id[at]
    provenance <- ggchord_focus_source(virtual)
    virtual$.source_accver <- provenance$accver
    virtual$.source_origin <- provenance$origin + instances$window_start - 1L
    virtual$.source_length <- provenance$length
    virtual$length <- as.numeric(instances$window_end - instances$window_start + 1L)
    for (field in setdiff(names(windows), c("instance_id", "accver", "start", "end", "orientation", "ring_id", "row_id")))
      virtual[[field]] <- windows[[field]][at]
  }
  plot$ggchord$facet_reference <- virtual
  virtual$accver <- instances$key
  routes <- omissions <- list()
  facet <- plot$facet$params
  facet_keys <- names(c(facet$facets, facet$rows, facet$cols))
  numeric_fields <- c("seq_radius", "seq_curvature", "seq_labels", "seq_thickness_mm",
    "seq_label_radius", "seq_label_rotation", "seq_label_size", "seq_label_hjust", "seq_label_vjust",
    "gene_width", "gene_label_radial_offset", "gene_label_circum_offset", "gene_label_circum_limit",
    "gene_label_rotation", "ribbon_gap", "front_gap_mm", "link_gap", "label_size", "label_offset")
  for (i in seq_along(plot$layers)) {
    layer <- plot$layers[[i]]; p <- layer$ggchord_params
    if (is.null(p)) next
    type <- p$type
    fallback <- if (type %in% c("gene", "gene_label", "gene_label_repel")) ggchord_annotation_fallback(original) else
      if (type == "ribbon") original$ribbon_data else if (type == "seq_region") p$regions else original$seq_data
    input <- if (type == "link") ggchord_resolve_layer_input(layer) else ggchord_resolve_layer_input(layer, fallback)
    source <- attr(input, "ggchord_source")
    # Routing rewrites accessions and can omit or replicate rows. Facets
    # belong to the source observations, including their factor levels.
    keys <- intersect(facet_keys, names(input))
    if (length(keys)) p$instance_facet_source <- input[keys]
    if (type %in% c("ribbon", "link")) {
      input_rows <- seq_len(nrow(input))
      retain <- function(keep,reason) {
        if (any(!keep)) omissions[[length(omissions)+1L]] <<- data.frame(
          layer_id=layer$ggchord_layer_id,input_row=input_rows[!keep],reason=reason,
          qaccver=as.character(input$qaccver[!keep]),saccver=as.character(input$saccver[!keep]))
        input <<- input[keep,,drop=FALSE]; input_rows <<- input_rows[keep]
      }
      if (!is.null(input) && nrow(input)) {
        if (any(!input$qaccver %in% seqs$accver | !input$saccver %in% seqs$accver))
          ggchord_stop("Instance routing: unknown source sequence ID")
        retain(input$qaccver %in% instances$accver & input$saccver %in% instances$accver,"sequence_not_displayed")
        if (mirrored && !isTRUE(plot$coordinates$sankey_self_links))
          retain(as.character(input$qaccver) != as.character(input$saccver),"self_link_filtered")
        if (!is.null(p$pairs)) {
          key <- function(q,s) paste(pmin(q,s), pmax(q,s), sep="\r")
          retain(key(input$qaccver,input$saccver) %in% key(p$pairs$qaccver,p$pairs$saccver),"pair_not_selected")
        }
      }
      routed <- ggchord_route_instances(input, instances, p$between, layer$ggchord_layer_id, mirrored,input_rows)
      input <- routed$data
      if (!is.null(routed$rows)) input_rows <- input_rows[routed$rows]
      p$visual_source <- ggchord_visual_rows(p$visual_source, input_rows)
      field <- paste0(type, "_ctrl_point")
      p[[field]] <- ggchord_link_control_rows(p[[field]], input_rows)
      if (length(keys)) p$instance_facet_rows <- input_rows
      if (nrow(routed$report)) {
        routed$report$allow_repeated <- isTRUE(p$allow_repeated)
        routes[[length(routes)+1L]] <- routed$report
      }
      if (!is.null(p$between)) {
        if (isTRUE(p$adjacent_only)) {
          delta <- abs(match(input$qaccver,instances$key)-match(input$saccver,instances$key))
          adjacent <- delta == 1L | ((!isTRUE(plot$coordinates$ggchord_collinear) || isTRUE(plot$coordinates$ggchord_radial)) & delta == nrow(instances)-1L)
          if (any(!adjacent)) ggchord_stop("between conflicts with explicit adjacent_only=TRUE")
        } else p$adjacent_only <- FALSE
      }
      p$pairs <- NULL
      if (mirrored) p$adjacent_only <- FALSE
    } else if (!is.null(input) && "accver" %in% names(input)) {
      location <- p$annotation_location %||% p$anchor_location
      if (is.null(location) && type == "domain_label" && NROW(input) &&
          ggchord_label_kind(p, input) == "site") location <- "site"
      expanded <- ggchord_instance_expand(input, instances,
        reference = ggchord_reference_layer(type, p$label_kind), location = location)
      # A separate keyed label table need not contain coordinates. Decide
      # which windows contain its target using the target's original rows.
      if (type == "domain_label" && !isTRUE(p$anchor_default_data) &&
          NROW(p$anchor_data) && NROW(expanded$data) && "window_start" %in% names(instances)) {
        matching <- ggchord_label_match(input, p$anchor_data, p$label_kind,
          isTRUE(p$anchor_same_source), p$anchor_location)
        members <- attr(matching, "entity_rows") %||% as.list(matching)
        targets <- ggchord_instance_expand(p$anchor_data, instances, location = p$anchor_location)
        visible <- split(targets$rows, targets$data$accver)
        keep <- vapply(seq_along(expanded$rows), function(row) {
          any(members[[expanded$rows[row]]] %in% visible[[expanded$data$accver[row]]])
        }, logical(1))
        expanded$data <- expanded$data[keep, , drop = FALSE]
        expanded$rows <- expanded$rows[keep]
      }
      input <- expanded$data
      p$visual_source <- ggchord_visual_rows(p$visual_source, expanded$rows)
      if (length(keys)) p$instance_facet_rows <- expanded$rows
      if (!is.null(source)) {
        source$rows <- source$rows[expanded$rows]
        p$input_source <- source
      }
    }
    if (type %in% c("gene", "gene_label", "gene_label_repel"))
      layer$ggchord_annotation_origin <- layer$ggchord_input_data
    layer$ggchord_input_data <- input
    p$reference_ids <- NULL  # validated before replacing accessions with instance keys
    p$reference <- NULL
    layer$ggchord_input_transform <- NULL
    p$source_stat <- NULL
    for (role in intersect(layer$ggchord_role_aes, names(layer$ggchord_input_mapping)))
      layer$ggchord_input_mapping[[role]] <- new_quosure(sym(role), emptyenv())
    for (field in intersect(names(p), numeric_fields)) p[[field]] <- ggchord_instance_param(p[[field]], instances)
    if (!is.null(p$feature_position))
      p$feature_position <- ggchord_instance_position(p$feature_position, instances)
    if (type == "seq") {
      p$backbone_input <- NULL
      p$seq_labels <- p$seq_labels %||% setNames(instances$accver, instances$key)
      entities <- unique(instances$accver)
      colours <- setNames(chord_default_palette(length(entities)), entities)
      p$seq_colors <- setNames(unname(colours[instances$accver]), instances$key)
    }
    layer$ggchord_params <- p
    plot$layers[[i]] <- layer
  }
  report <- ggchord_rbind_fill(routes) %||% data.frame()
  coord <- plot$coordinates
  if (!is.null(coord$collinear_anchor) || !is.null(coord$collinear_segments)) {
    source_coord <- coord
    coord <- ggproto(NULL, source_coord)
    if (!is.null(coord$collinear_anchor) && any(!coord$collinear_anchor$accver %in% instances$key))
      coord$collinear_anchor <- ggchord_instance_expand(coord$collinear_anchor, instances)$data
    plot$coordinates <- coord
  }
  plot$ggchord$complete_display_instances <- instances
  arrangement$order <- instances$key
  if (mirrored) for (field in c("orientation", "gap", "row_gap", "radius", "curvature")) {
    value <- arrangement[[field]]
    if (!is.null(value)) arrangement[[field]] <- setNames(rep(value, 2L), instances$key)
  }
  plot$ggchord$sequence_layout <- arrangement
  plot$ggchord$data$seq_data <- virtual
  plot$ggchord$display_instances <- instances
  plot$ggchord$instance_routing <- report
  plot$ggchord$instance_omissions <- ggchord_rbind_fill(omissions) %||% data.frame()
  plot$ggchord$instances_prepared <- TRUE
  plot
}

# Validate actual displayed source rows after pair/adjacency/direction filtering.
# This also runs when ordinary Chord/Collinear need no instance expansion.
ggchord_check_displayed_hsps <- function(layout, plot) {
  identities <- list()
  routes <- layout$instance_routing
  if (is.data.frame(routes) && nrow(routes)) {
    routes$retained <- FALSE
    routes$repeated_display <- FALSE
  }
  for (layer in plot$layers) {
    component <- layer$ggchord_type
    if (is.null(component) || !component %in% c("ribbon", "link")) next
    id <- layer$ggchord_layer_id
    input <- layout$layer_inputs[[id]][[component]]
    geometry <- layout$layer_geometry[[id]][[component]]
    visible <- unique(geometry$source_row)
    if (identical(layer$ggchord_params$ribbon_stat$type, "bundle")) {
      members <- layout$layer_layouts[[id]]$ribbon_stat_report$source_rows[visible]
      visible <- as.integer(unlist(strsplit(members, ",", fixed = TRUE)))
      input <- ggchord_resolve_layer_input(layer, plot$ggchord$data$ribbon_data)
    }
    if (is.data.frame(routes) && nrow(routes)) {
      at <- which(routes$layer_id == id)[visible]
      routes$retained[at] <- TRUE
    }
    if (!"source_hsp_id" %in% names(input)) next
    stable <- as.character(input$source_hsp_id)
    if (anyNA(stable) || any(!nzchar(stable)))
      ggchord_stop("source_hsp_id must be nonempty when supplied")
    if (!length(visible)) next
    identities[[length(identities) + 1L]] <- data.frame(
      source_hsp_id = stable[visible], allow_repeated = isTRUE(layer$ggchord_params$allow_repeated))
  }
  identity <- ggchord_rbind_fill(identities)
  if (is.data.frame(identity) && nrow(identity)) {
    repeated <- duplicated(identity$source_hsp_id) | duplicated(identity$source_hsp_id, fromLast = TRUE)
    if (any(repeated & !identity$allow_repeated)) ggchord_stop(
      "A declared source_hsp_id is displayed more than once; use disjoint layer data or allow_repeated=TRUE explicitly")
    if (is.data.frame(routes) && nrow(routes)) routes$repeated_display <-
      routes$retained & routes$source_scope == "declared" &
      routes$source_hsp_id %in% identity$source_hsp_id[repeated]
  }
  layout$instance_routing <- routes
  layout
}

ggchord_restore_instance_data <- function(data, instances, fields = names(data)) {
  if (!is.data.frame(data) || !nrow(data) || is.null(instances)) return(data)
  for (pair in list(c("accver", "display_instance_id"), c("qaccver", "q_instance"), c("saccver", "s_instance"))) {
    if (!pair[1] %in% fields || !pair[1] %in% names(data)) next
    at <- match(as.character(data[[pair[1]]]), instances$key)
    known <- !is.na(at)
    if (!any(known)) next
    if (!pair[2] %in% names(data)) data[[pair[2]]] <- rep(NA_character_, nrow(data))
    data[[pair[2]]][known] <- instances$instance_id[at[known]]
    data[[pair[1]]][known] <- instances$accver[at[known]]
  }
  data
}
