# A source entity can own several display rows. The legacy engines consume
# unique private keys, while public data and scales retain biological IDs.
ggchord_instance_param <- function(value, instances) {
  if (is.null(value) || grid::is.unit(value) || is.null(names(value))) return(value)
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
  for (field in c("offset", "width", "gap")) if (!is.null(position[[field]]))
    changes[[field]] <- ggchord_instance_param(position[[field]], instances)
  if (!is.null(position$base_position)) changes$base_position <-
    ggchord_instance_position(position$base_position, instances)
  do.call(ggplot2::ggproto, c(list(NULL, position), changes))
}

ggchord_instance_expand <- function(input, instances) {
  if (is.null(input)) return(list(data = NULL, rows = integer()))
  ggchord_require_columns(input, "accver", "display instances")
  rows <- integer(); targets <- integer()
  for (i in seq_len(nrow(input))) {
    at <- which(instances$accver == as.character(input$accver[i]))
    if ("display_instance_id" %in% names(input) && !is.na(input$display_instance_id[i]))
      at <- at[instances$instance_id[at] == as.character(input$display_instance_id[i])]
    if (!length(at) && "display_instance_id" %in% names(input) && !is.na(input$display_instance_id[i]))
      ggchord_stop("No display instance matches input row ", i, " (", input$accver[i], ")")
    rows <- c(rows, rep(i, length(at))); targets <- c(targets, at)
  }
  out <- input[rows, , drop = FALSE]
  out$accver <- instances$key[targets]
  out$display_instance_id <- instances$instance_id[targets]
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
  out$qaccver <- as.character(out$qaccver); out$saccver <- as.character(out$saccver)
  for (i in seq_len(nrow(input))) {
    q <- which(instances$accver == as.character(input$qaccver[i]))
    s <- which(instances$accver == as.character(input$saccver[i]))
    if (mirrored) { q <- q[instances$side[q] == "left"]; s <- s[instances$side[s] == "right"] }
    candidates <- expand.grid(q = q, s = s)
    candidates <- candidates[candidates$q != candidates$s, , drop = FALSE]
    if (!is.null(between)) candidates <- candidates[
      instances$instance_id[candidates$q] %in% between & instances$instance_id[candidates$s] %in% between, , drop = FALSE]
    if (!is.null(between) && input$qaccver[i] == input$saccver[i])
      candidates <- candidates[instances$instance_id[candidates$q] == between[1], , drop=FALSE]
    if (nrow(candidates) != 1L) ggchord_stop("HSP row ", i,
      " has an absent or ambiguous instance route; specify between = c(instance1, instance2)")
    qi <- candidates$q; si <- candidates$s
    out$qaccver[i] <- instances$key[qi]; out$saccver[i] <- instances$key[si]
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
  list(data = out, report = ggchord_rbind_fill(report) %||% data.frame())
}

ggchord_prepare_instances <- function(plot) {
  if (isTRUE(plot$ggchord$instances_prepared)) return(plot)
  skeleton <- which(vapply(plot$layers, function(l) identical(l$ggchord_params$type, "seq"), logical(1)))
  order <- if (length(skeleton)) plot$layers[[tail(skeleton, 1)]]$ggchord_params$seq_order else NULL
  configured <- any(vapply(plot$layers, function(l) !is.null(l$ggchord_params$between), logical(1)))
  mirrored <- isTRUE(plot$coordinates$ggchord_sankey)
  activate <- mirrored || configured || (!is.null(order) && (anyDuplicated(order) ||
    (!is.null(names(order)) && any(names(order) != unname(order)))))
  if (!activate) return(plot)
  if (isTRUE(plot$coordinates$ggchord_circular)) ggchord_stop("coord_circular() does not support repeated display instances")
  seqs <- if (length(skeleton)) ggchord_resolve_layer_input(plot$layers[[tail(skeleton,1)]],
    plot$ggchord$data$seq_data) else plot$ggchord$data$seq_data
  order <- order %||% seqs$accver
  if (!is.character(order) || !length(order) || anyNA(order) || any(!order %in% seqs$accver))
    ggchord_stop("seq_order must contain known source sequence IDs")
  if (mirrored) {
    if (anyDuplicated(order)) ggchord_stop("coord_sankey() owns left/right copies; use each source once in seq_order")
    if (configured) ggchord_stop("coord_sankey() routes query left and subject right; between is for repeated Chord/Collinear views")
    labels <- names(order) %||% unname(order)
    prefix <- ""
    while (any(c(paste0(prefix,"left:",labels),paste0(prefix,"right:",labels)) %in% seqs$accver)) prefix <- paste0(".",prefix)
    order <- stats::setNames(rep(unname(order),2),c(paste0(prefix,"left:",labels),paste0(prefix,"right:",labels)))
  }
  ids <- names(order)
  if (is.null(ids)) {
    if (!anyDuplicated(order)) ids <- unname(order) else {
      prefix <- "instance-"
      while (any(paste0(prefix, seq_along(order)) %in% seqs$accver)) prefix <- paste0(".", prefix)
      ids <- paste0(prefix, seq_along(order))
    }
  }
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) ggchord_stop("Display instance names must be nonempty and unique")
  for (i in seq_along(ids)) if (ids[i] %in% seqs$accver &&
      (ids[i] != order[i] || sum(order == ids[i]) > 1L)) ggchord_stop(
    "Instance name '", ids[i], "' conflicts with a source ID; give repeated instances distinct names")
  instances <- data.frame(instance_id = ids, accver = unname(order),
    declared_order = seq_along(order), key = ids, stringsAsFactors = FALSE)
  instances$side <- if (mirrored) rep(c("left","right"),each=nrow(instances)/2) else NA_character_
  original <- plot$ggchord$data
  virtual <- seqs[match(instances$accver, seqs$accver), , drop = FALSE]
  virtual$accver <- instances$key
  routes <- omissions <- list()
  numeric_fields <- c("seq_orientation", "seq_gap", "seq_radius", "seq_curvature", "seq_labels",
    "seq_label_radius", "seq_label_rotation", "seq_label_size", "seq_label_hjust", "seq_label_vjust",
    "gene_width", "gene_label_radial_offset", "gene_label_circum_offset", "gene_label_circum_limit",
    "gene_label_rotation", "ribbon_gap", "link_gap", "label_size", "label_offset")
  for (i in seq_along(plot$layers)) {
    layer <- plot$layers[[i]]; p <- layer$ggchord_params
    if (is.null(p)) next
    type <- p$type
    fallback <- if (type %in% c("gene", "gene_label", "gene_label_repel")) ggchord_annotation_fallback(original) else
      if (type == "ribbon") original$ribbon_data else if (type == "seq_region") p$regions else original$seq_data
    input <- if (type == "link") ggchord_resolve_layer_input(layer) else ggchord_resolve_layer_input(layer, fallback)
    source <- attr(input, "ggchord_source")
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
      expanded <- ggchord_instance_expand(input, instances)
      input <- expanded$data
      if (!is.null(source)) {
        source$rows <- source$rows[expanded$rows]
        p$input_source <- source
      }
    }
    layer$ggchord_input_data <- input
    layer$ggchord_input_transform <- NULL
    for (role in intersect(layer$ggchord_role_aes, names(layer$ggchord_input_mapping)))
      layer$ggchord_input_mapping[[role]] <- rlang::new_quosure(rlang::sym(role), emptyenv())
    for (field in intersect(names(p), numeric_fields)) p[[field]] <- ggchord_instance_param(p[[field]], instances)
    for (field in c("feature_position", "track_position")) if (!is.null(p[[field]]))
      p[[field]] <- ggchord_instance_position(p[[field]], instances)
    if (type == "seq") {
      p$seq_order <- instances$key
      p$seq_labels <- p$seq_labels %||% stats::setNames(instances$accver, instances$key)
      entities <- unique(instances$accver)
      colours <- stats::setNames(chord_default_palette(length(entities)), entities)
      p$seq_colors <- stats::setNames(unname(colours[instances$accver]), instances$key)
    }
    layer$ggchord_params <- p
    plot$layers[[i]] <- layer
  }
  report <- ggchord_rbind_fill(routes) %||% data.frame()
  if (nrow(report)) {
    repeated <- duplicated(report$source_hsp_id) | duplicated(report$source_hsp_id, fromLast = TRUE)
    declared <- report$source_scope == "declared"
    if (any(repeated & declared & !report$allow_repeated)) ggchord_stop(
      "A declared source_hsp_id is displayed more than once; use disjoint layer data or allow_repeated=TRUE explicitly")
    report$repeated_display <- repeated & declared
  }
  coord <- plot$coordinates
  if (!is.null(coord$collinear_anchor) || !is.null(coord$collinear_segments)) {
    source_coord <- coord
    coord <- ggplot2::ggproto(NULL, source_coord)
    if (is.numeric(coord$collinear_anchor)) coord$collinear_anchor <- ggchord_instance_param(coord$collinear_anchor, instances) else
      if (is.data.frame(coord$collinear_anchor)) coord$collinear_anchor <- ggchord_instance_expand(coord$collinear_anchor, instances)$data
    if (!is.null(coord$collinear_segments)) coord$collinear_segments <- ggchord_instance_expand(coord$collinear_segments, instances)$data
    plot$coordinates <- coord
  }
  axis <- plot$scales$get_scales("seq_position")
  if (!is.null(axis$ggchord_axis_sequences)) {
    plot$scales <- plot$scales$clone()
    axis <- plot$scales$get_scales("seq_position")
    chosen <- axis$ggchord_axis_sequences
    axis$ggchord_axis_sequences <- instances$key[instances$accver %in% chosen | instances$instance_id %in% chosen]
  }
  plot$ggchord$data$seq_data <- virtual
  plot$ggchord$display_instances <- instances
  plot$ggchord$instance_routing <- report
  plot$ggchord$instance_omissions <- ggchord_rbind_fill(omissions) %||% data.frame()
  plot$ggchord$instances_prepared <- TRUE
  plot
}

ggchord_restore_instance_data <- function(data, instances) {
  if (!is.data.frame(data) || !nrow(data) || is.null(instances)) return(data)
  for (pair in list(c("accver", "display_instance_id"), c("qaccver", "q_instance"), c("saccver", "s_instance"))) {
    if (!pair[1] %in% names(data)) next
    at <- match(as.character(data[[pair[1]]]), instances$key)
    known <- !is.na(at)
    if (!any(known)) next
    if (!pair[2] %in% names(data)) data[[pair[2]]] <- rep(NA_character_, nrow(data))
    data[[pair[2]]][known] <- instances$instance_id[at[known]]
    data[[pair[1]]][known] <- instances$accver[at[known]]
  }
  data
}
