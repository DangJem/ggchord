# Composition owns device-local panels; source plots remain independent.

#' Compose independent genomic panels
#'
#' A Fusion is a static composition, not a coordinate system. Add up to four
#' [fusion_panel()] objects and independent [fusion_ribbon()] layers. Ordinary
#' [ggplot2::ggsave()] rebuilds panels and interval fronts at the target size.
#' Source plots and their scales remain independent. Composition scales and
#' themes apply only to cross-panel links and their separate legend.
#'
#' @param center A ggchord plot using Circular, Chord or horizontal/vertical Collinear coordinates.
#' @param center_size Optional width and height as fractions of the composition
#'   area, excluding its cross-panel legend. NULL allocates available space.
#' @param layout `"fixed"` uses fixed panel allocations; `"adaptive"` also
#'   adjusts the central allocation to the physical aspect ratio. Both rebuild
#'   labels and guides at actual panel sizes. Neither rescales explicit fonts.
#' @return A Fusion object, printable and accepted by [ggplot2::ggsave()].
#' @export
#' @examples
#' library(dplyr)
#' data(seq_data_example)
#' seqs <- seq_data_example |> select(accver, length) |> slice_head(n = 2)
#' center <- ggchord(seqs) + geom_seq() + coord_chord()
#' detail <- ggchord(seqs) + geom_seq() + coord_collinear(seq_axis = "vertical")
#' ggchord_fusion(center) + fusion_panel(detail, id = "detail", side = "right")
ggchord_fusion <- function(center, center_size = NULL, layout = c("fixed", "adaptive")) {
  ggchord_fusion_check_plot(center, center = TRUE)
  ggchord_fusion_size(center_size)
  # Retain the ggplot carrier so ggplot2 4's S7 Ops bridge can dispatch
  # additions of native scales/themes without conflicting S3 '+' methods.
  x <- ggplot2::ggplot()
  x$panels <- list(center = list(plot = center, side = "center", size = center_size))
  x$ribbons <- list(); x$styles <- list(); x$fusion_layout <- match.arg(layout)
  class(x) <- c("ggchord_fusion", class(x))
  ggplot2::set_last_plot(x)
  x
}

#' Add a genomic panel to a Fusion
#' @param plot A horizontal or vertical Collinear ggchord plot. Sankey,
#'   Radial, Dotplot, facets and external ggplots are not supported here.
#' @param id Unique nonempty panel name; `"center"` is reserved.
#' @param side One of right, top, left or bottom; one panel per side.
#' @param size Optional fractional width and height; overlapping or out-of-frame
#'   allocations fail with a diagnostic instead of moving explicit panels.
#' @return A component for [ggchord_fusion()].
#' @export
#' @examples
#' data(seq_data_example)
#' p <- ggchord(seq_data_example) + geom_seq()
#' ggchord_fusion(p) + fusion_panel(p + coord_collinear(), "right")
fusion_panel <- function(plot, id, side = c("right", "top", "left", "bottom"), size = NULL) {
  ggchord_fusion_check_plot(plot)
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id) || id == "center")
    ggchord_stop("fusion_panel(): id must be one nonempty name other than 'center'")
  ggchord_fusion_size(size)
  structure(list(plot = plot, id = id, side = match.arg(side), size = size), class = "ggchord_fusion_panel")
}

#' Connect real intervals across Fusion panels
#'
#' Uses ordinary qaccver/saccver and qstart/qend/sstart/send columns. Optional
#' qpanel/spanel and q_instance/s_instance select exact panels and display
#' instances. Missing selectors are inferred only when a unique candidate
#' exists. Query/subject roles never depend on the placement of a panel.
#' Intervals are one-based inclusive, with ordered fronts spanning `[min-1,max]`;
#' a one-position interval has nonzero width. Fully visible intervals are
#' required; split/focus the source data explicitly before composition.
#'
#' @param data A relationship data frame, preserving original endpoints.
#' @param mapping Standard ggplot2 aesthetic mapping: fill, colour, alpha,
#'   linewidth and linetype. With pident and no explicit fill, Identity is
#'   mapped by default with 0/50/80/90/95/100 breaks. Composition uses ordinary
#'   scale_fill_* etc.; panel role scales remain local to each panel.
#' @param qpanel,spanel Optional fixed panel IDs. Do not also supply their columns.
#' @param gradient Opt in to query-to-subject sequence-colour interpolation.
#'   Cannot be combined with an explicit fill or fill mapping.
#' @param tip `"none"` or a decorative `"subject"` tip. Source fronts remain
#'   exported in full even when the tip is drawn.
#' @param allow_repeated Explicitly allow repeated source_hsp_id display.
#' @param curvature Control distance relative to endpoint separation, in `[0,1]`.
#' @param fill,colour,alpha,linewidth,linetype Fixed styles; explicit values
#'   take priority over mappings. The default fill without pident is steelblue.
#' @param show.legend Standard ggplot2 legend selection: NA follows mappings;
#'   TRUE/FALSE or a named logical vector controls guide participation.
#' @return A component for [ggchord_fusion()].
#' @export
#' @examples
#' data(seq_data_example)
#' data(ribbon_data_example)
#' p <- ggchord(seq_data_example) + geom_seq()
#' ggchord_fusion(p) +
#'   fusion_panel(p + coord_collinear(), "detail") +
#'   fusion_ribbon(ribbon_data_example, qpanel = "center", spanel = "detail")
fusion_ribbon <- function(data, mapping = NULL, qpanel = NULL, spanel = NULL,
    gradient = FALSE, tip = c("none", "subject"), allow_repeated = FALSE,
    curvature = .3, fill = NULL, colour = NA, alpha = .35, linewidth = .25,
    linetype = 1, show.legend = NA) {
  ggchord_require_columns(data, c("qaccver", "saccver", "qstart", "qend", "sstart", "send"), "fusion_ribbon()")
  if (!is.null(mapping) && !inherits(mapping, "ggplot2::mapping") && !inherits(mapping, "uneval"))
    ggchord_stop("fusion_ribbon(): mapping must be aes()")
  if (length(setdiff(names(mapping), c("fill", "colour", "alpha", "linewidth", "linetype"))))
    ggchord_stop("fusion_ribbon(): only fill, colour, alpha, linewidth and linetype can be mapped")
  for (value in list(gradient, allow_repeated)) if (!is.logical(value) || length(value) != 1 || is.na(value))
    ggchord_stop("fusion_ribbon(): flags must be TRUE or FALSE")
  if (!is.numeric(curvature) || length(curvature) != 1 || !is.finite(curvature) || curvature < 0 || curvature > 1)
    ggchord_stop("fusion_ribbon(): curvature must be in [0,1]")
  if (gradient && (!is.null(fill) || !is.null(mapping$fill)))
    ggchord_stop("fusion_ribbon(): gradient conflicts with explicit fill; choose one")
  selectors <- list(qpanel = qpanel, spanel = spanel)
  for (name in names(selectors)) if (!is.null(selectors[[name]])) {
    if (name %in% names(data)) ggchord_stop("fusion_ribbon(): use either the ", name, " column or argument")
    if (!is.character(selectors[[name]]) || length(selectors[[name]]) != 1 || is.na(selectors[[name]]) || !nzchar(selectors[[name]]))
      ggchord_stop("fusion_ribbon(): panel selectors must be nonempty names")
  }
  fixed <- list()
  for (name in c("fill", "colour", "alpha", "linewidth", "linetype")) {
    if (!switch(name, fill = missing(fill), colour = missing(colour),
        alpha = missing(alpha), linewidth = missing(linewidth), linetype = missing(linetype))) {
      value <- get(name)
      if (length(value) != 1L) ggchord_stop("fusion_ribbon(): fixed styles must be scalar; use aes() for data")
      fixed[[name]] <- value
    }
  }
  structure(list(data = data, mapping = mapping, selectors = selectors, fixed = fixed,
    gradient = gradient, tip = match.arg(tip), allow_repeated = allow_repeated,
    curvature = curvature, show.legend = show.legend), class = "ggchord_fusion_ribbon")
}

#' @export
`+.ggchord_fusion` <- function(e1, e2) {
  if (inherits(e2, "ggchord_fusion_panel")) {
    if (e2$id %in% names(e1$panels)) ggchord_stop("Fusion: duplicate panel id '", e2$id, "'")
    if (e2$side %in% vapply(e1$panels, `[[`, "", "side")) ggchord_stop("Fusion: only one panel per side is supported")
    e1$panels[[e2$id]] <- unclass(e2)
  } else if (inherits(e2, "ggchord_fusion_ribbon")) {
    e1$ribbons[[length(e1$ribbons) + 1L]] <- unclass(e2)
  } else if (inherits(e2, "Scale") || inherits(e2, "theme") || inherits(e2, "Guides") || inherits(e2, "ggplot2::labels")) {
    e1$styles[[length(e1$styles) + 1L]] <- e2
  } else ggchord_stop("Fusion accepts fusion_panel(), fusion_ribbon(), scales, guides, labels and themes")
  ggplot2::set_last_plot(e1)
  e1
}

ggchord_fusion_size <- function(size) {
  if (!is.null(size) && (!is.numeric(size) || length(size) != 2 || any(!is.finite(size)) || any(size <= 0 | size >= 1)))
    ggchord_stop("Fusion: size must contain width and height fractions strictly between 0 and 1")
}

ggchord_fusion_check_plot <- function(plot, center = FALSE) {
  if (!inherits(plot, "ggchord")) ggchord_stop("Fusion panels must be ggchord plots")
  coord <- plot$coordinates
  if (center) {
    if (!isTRUE(coord$ggchord_coord) || isTRUE(coord$ggchord_sankey) || isTRUE(coord$ggchord_radial) || isTRUE(coord$ggchord_dotplot))
      ggchord_stop("Fusion center requires Circular, Chord or Collinear coordinates")
  } else if (!isTRUE(coord$ggchord_collinear) || isTRUE(coord$ggchord_sankey) || isTRUE(coord$ggchord_radial))
    ggchord_stop("Fusion peripheral panels require horizontal or vertical Collinear coordinates")
  if (!inherits(plot$facet, "FacetNull")) ggchord_stop("Fusion panels do not support facets")
}

# Train cross-panel styles through ggplot2; no second colour mapping engine.
ggchord_fusion_styles <- function(x) {
  p <- ggplot2::ggplot() + ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(legend.position = "right", plot.margin = ggplot2::margin(0,0,0,0))
  default_identity <- FALSE
  for (layer in x$ribbons) {
    input <- layer$data
    d <- input[rep(seq_len(nrow(input)), each = 3L), , drop = FALSE]
    d$.fusion_x <- rep(c(0,1,0), nrow(input)); d$.fusion_y <- rep(c(0,0,1), nrow(input))
    d$.fusion_row <- rep(seq_len(nrow(input)), each = 3L)
    mapping <- layer$mapping %||% ggplot2::aes()
    if (is.null(mapping$fill) && is.null(layer$fixed$fill) && !layer$gradient && "pident" %in% names(input)) {
      mapping$fill <- ggplot2::aes(fill = pident)$fill
      default_identity <- TRUE
    }
    mapping$x <- ggplot2::aes(x = .fusion_x)$x; mapping$y <- ggplot2::aes(y = .fusion_y)$y
    mapping$group <- ggplot2::aes(group = .fusion_row)$group
    defaults <- list(fill = "steelblue", colour = NA, alpha = .35, linewidth = .25, linetype = 1)
    fixed <- layer$fixed
    for (name in names(defaults)) if (is.null(mapping[[name]]) && is.null(fixed[[name]])) fixed[[name]] <- defaults[[name]]
    p <- p + do.call(ggplot2::geom_polygon, c(list(data = d, mapping = mapping, inherit.aes = FALSE,
      show.legend = if (layer$gradient) FALSE else layer$show.legend), fixed))
  }
  if (default_identity) p <- p + ggplot2::scale_fill_stepsn(name = "Identity (%)",
    colours = c("#440154", "#482878", "#3E4A89", "#31688E", "#26828E", "#1F9E89", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725"),
    limits = c(0,100), breaks = c(0,50,80,90,95,100),
    guide = ggchord_fusion_identity_guide(p,x$styles))
  for (style in x$styles) p <- suppressMessages(p + style)
  built <- ggplot2::ggplot_build(p)
  table <- ggplot2::ggplot_gtable(built)
  at <- which(grepl("^guide-box", table$layout$name) & !vapply(table$grobs, inherits, logical(1), "zeroGrob"))
  if (length(at) > 1L) ggchord_stop("Fusion: use one common legend position for cross-panel guides")
  guide <- if (length(at)) table$grobs[[at]] else grid::nullGrob()
  side <- if (length(at)) sub("guide-box-", "", table$layout$name[at]) else "none"
  if (!side %in% c("none", "left", "right", "top", "bottom")) ggchord_stop("Fusion: cross-panel legends must be outside the panels")
  styles <- lapply(built$data, function(d) d[seq.int(1L, length.out = nrow(d) / 3L, by = 3L), , drop = FALSE])
  list(data = styles, guide = guide, side = side)
}

ggchord_fusion_rectangles <- function(x, canvas) {
  sides <- vapply(x$panels, `[[`, "", "side")
  available <- c(if ("left" %in% sides) .3 else .02,
    if ("right" %in% sides) .7 else .98,
    if ("bottom" %in% sides) .3 else .02,
    if ("top" %in% sides) .7 else .98)
  center <- c(diff(available[1:2]),diff(available[3:4]))
  center_xy <- c(mean(available[1:2]),mean(available[3:4]))
  if (identical(x$fusion_layout, "adaptive")) {
    physical <- center * canvas
    center <- pmin(center, min(physical) / canvas)
  }
  rects <- lapply(x$panels, function(panel) {
    side <- panel$side
    size <- panel$size %||% switch(side, center = center,
      right = c(.26,.7), left = c(.26,.7), top = c(.4,.26), bottom = c(.4,.26))
    xy <- switch(side, center = center_xy, right = c(.99-size[1]/2,.5),
      left = c(.01+size[1]/2,.5), top = c(.5,.99-size[2]/2), bottom = c(.5,.01+size[2]/2))
    c(xmin=xy[1]-size[1]/2,xmax=xy[1]+size[1]/2,ymin=xy[2]-size[2]/2,ymax=xy[2]+size[2]/2)
  })
  for (i in seq_along(rects)) {
    a <- rects[[i]]
    if (any(a < 0 | a > 1)) ggchord_stop("Fusion: panel allocation lies outside canvas")
    for (j in seq_len(i-1L)) {
      b <- rects[[j]]
      if (min(a[2],b[2]) > max(a[1],b[1]) && min(a[4],b[4]) > max(a[3],b[3]))
        ggchord_stop("Fusion: explicit panel allocations overlap: ", names(rects)[i], " and ", names(rects)[j])
    }
  }
  rects
}

# Partition backgrounds without deleting explicit fills or repainting foreground.
ggchord_fusion_partition <- function(g, background = FALSE, inherited = FALSE) {
  is_background <- inherited || grepl("^(plot[.]background|panel[.]background|background)([.-]|$)", g$name %||% "")
  if (!is.null(g$grobs)) {
    g$grobs <- lapply(g$grobs, ggchord_fusion_partition, background = background, inherited = is_background)
  } else if (!is.null(g$children) && length(g$children)) {
    children <- lapply(g$children, ggchord_fusion_partition, background = background, inherited = is_background)
    g$children <- do.call(grid::gList, children)
  } else if (background != is_background) return(grid::nullGrob(name = g$name))
  g
}

# Resolve entity + instance + fully visible source window before projection.
ggchord_fusion_endpoint <- function(panels, row, prefix, layer_id, source_row) {
  entity <- as.character(row[[paste0(prefix,"accver")]])
  selected_panel <- row[[paste0(prefix,"panel")]]
  selected_instance <- row[[paste0(prefix,"_instance")]]
  start <- row[[paste0(prefix,"start")]]; end <- row[[paste0(prefix,"end")]]
  where <- paste0("Fusion ", layer_id, " row ", source_row, " ", prefix, ": ")
  ggchord_resolve_panel_interval(panels, entity, start, end, selected_panel, selected_instance, where)
}

# Shared locus resolver: decorative focus regions are not HSP relationships.
ggchord_resolve_panel_interval <- function(panels, entity, start, end,
    selected_panel, selected_instance, where) {
  if (length(entity)!=1 || is.na(entity) || !nzchar(entity) ||
      !is.numeric(start) || !is.numeric(end) || any(!is.finite(c(start,end))) || min(start,end) < 1 ||
      any(c(start,end) != floor(c(start,end)))) ggchord_stop(where, "invalid accession or one-based interval")
  if (!is.null(selected_panel) && (is.na(selected_panel) || !selected_panel %in% names(panels)))
    ggchord_stop(where, "unknown panel '", selected_panel, "'")
  candidates <- list()
  for (id in selected_panel %||% names(panels)) {
    layout <- panels[[id]]$layout
    keys <- layout$seqs
    instances <- layout$display_instances
    entities <- if (is.null(instances)) keys else instances$accver[match(keys,instances$key)]
    instance_ids <- if (is.null(instances)) keys else instances$instance_id[match(keys,instances$key)]
    seq_data <- panels[[id]]$seq_data
    original_entities <- entities
    if (".source_accver" %in% names(seq_data))
      original_entities <- as.character(seq_data$.source_accver[match(entities,seq_data$accver)])
    for (at in which(original_entities == entity)) {
      if (!is.null(selected_instance) && (is.na(selected_instance) || instance_ids[at] != selected_instance)) next
      key <- keys[at]
      origin <- 1
      if (".source_origin" %in% names(seq_data))
        origin <- seq_data$.source_origin[match(entities[at],seq_data$accver)]
      if (is.data.frame(layout$sequence_layout) && "source_origin" %in% names(layout$sequence_layout))
        origin <- layout$sequence_layout$source_origin[match(key, layout$sequence_layout$accver)]
      length <- unname(layout$sequence_reference$lens[key])
      if (min(start,end) < origin || max(start,end) > origin+length-1) next
      candidates[[length(candidates)+1L]] <- list(panel=id, key=key, instance=instance_ids[at],
        entity=entity, origin=origin, start=start, end=end, unit=layout$sequence_unit %||% "bp")
    }
  }
  if (length(candidates) != 1L) ggchord_stop(where,
    if (!length(candidates)) "no fully visible matching interval/instance" else "ambiguous panel or instance; specify panel and instance")
  candidates[[1L]]
}

ggchord_fusion_front <- function(panel, endpoint) {
  layout <- panel$layout; ref <- layout$sequence_reference; id <- endpoint$key
  boundary <- c(min(endpoint$start,endpoint$end)-1,max(endpoint$start,endpoint$end)) - endpoint$origin + 1
  if (endpoint$end < endpoint$start) boundary <- rev(boundary)
  pos <- seq(boundary[1], boundary[2], length.out = 40)
  if (isTRUE(layout$collinear)) {
    xy <- cbind(x=ggchord_collinear_x(pos,id,ref), y=rep(unname(ref$y[id]),length(pos)))
    normal <- cbind(x=rep(0,length(pos)),y=rep(1,length(pos)))
    if (identical(layout$seq_axis,"vertical")) {
      xy <- cbind(x=-xy[,2],y=xy[,1]); normal <- cbind(x=-normal[,2],y=normal[,1])
    }
  } else {
    # Track's continuous frame preserves even a single-base front. The
    # legacy nearest-reference Ribbon mapper would quantize short intervals.
    frame <- ggchord_track_frame(pos,id,layout)
    xy <- frame$base; normal <- frame$normal
  }
  coord <- panel$built$layout$coord
  project <- function(points) {
    # ggplot transforms position scales before Coord; reflecting a circular
    # panel must also reflect its interval anchors, even with symmetric limits.
    scaled <- data.frame(
      x=panel$built$layout$panel_scales_x[[1L]]$transform(points[,1]),
      y=panel$built$layout$panel_scales_y[[1L]]$transform(points[,2]))
    coord$transform(scaled,panel$built$layout$panel_params[[1L]])
  }
  projected <- project(xy)
  if (any(!is.finite(c(projected$x,projected$y))) || any(projected$x < -1e-7 | projected$x > 1+1e-7 | projected$y < -1e-7 | projected$y > 1+1e-7))
    ggchord_stop("Fusion: source front is outside the visible coordinate limits in panel '", endpoint$panel, "'")
  bounds <- panel$context$panel_bounds_mm
  normal_tip <- project(xy+normal)
  normal <- cbind((normal_tip$x-projected$x)*diff(bounds[1:2]),
    (normal_tip$y-projected$y)*diff(bounds[3:4]))
  normal <- normal / sqrt(rowSums(normal^2))
  data.frame(position_boundary = pos+endpoint$origin-1, local_x=xy[,1],local_y=xy[,2],
    x=bounds[1]+projected$x*diff(bounds[1:2]),y=bounds[3]+projected$y*diff(bounds[3:4]),
    normal_x=normal[,1],normal_y=normal[,2],tangent_x=-normal[,2],tangent_y=normal[,1])
}

ggchord_fusion_path <- function(q, s, curvature) {
  qc <- colMeans(q[c("x","y")]); sc <- colMeans(s[c("x","y")])
  normal <- function(front, toward) {
    middle <- floor(nrow(front)/2)
    tangent <- unlist(front[middle+1L,c("x","y")])-unlist(front[middle,c("x","y")])
    if (sum(tangent^2) == 0) ggchord_stop("Fusion: degenerate interval front")
    n <- c(-tangent[2],tangent[1])/sqrt(sum(tangent^2))
    if (sum(n*toward)<0) n <- -n
    n
  }
  nq <- normal(q,sc-qc); ns <- normal(s,qc-sc)
  distance <- sqrt(sum((sc-qc)^2))*curvature
  edge <- function(a,b) bezier_pts(a,b,a+nq*distance,b+ns*distance,n=60)
  a <- as.matrix(q[c("x","y")]); b <- as.matrix(s[c("x","y")])
  edge1 <- edge(a[1,],b[1,]); edge2 <- edge(a[nrow(a),],b[nrow(b),])
  path <- rbind(a,edge2,b[nrow(b):1,,drop=FALSE],edge1[nrow(edge1):1,,drop=FALSE])
  attr(path,"endpoint_normals") <- rbind(q=nq,s=ns)
  path
}

#' Export a Fusion at a specified physical size
#'
#' Builds at the requested size on a temporary measurement device, restoring
#' the previous device. Geometry reflects that device's font metrics. A later
#' save rebuilds on its actual output device, so minor metric differences are
#' possible. Source coordinates and identities remain unchanged.
#' @param plot A Fusion object.
#' @param width,height Positive canvas dimensions in inches.
#' @return A list with panels, instances, connections, fronts, paths,
#'   diagnostics and per-panel source layouts. Compositions made with
#'   [magnify_ggchord()] also populate focus_region_data, focus_regions, focus_fronts and decorative
#'   leaders; ordinary Fusion objects return empty tables for these fields.
#'   All global geometry uses mm.
#' @export
#' @examples
#' data(seq_data_example)
#' p <- ggchord(seq_data_example) + geom_seq()
#' f <- ggchord_fusion(p) + fusion_panel(p + coord_collinear(), "detail")
#' export_ggchord_fusion(f, width = 10, height = 8)$panels
export_ggchord_fusion <- function(plot, width = 10, height = 8) {
  if (!inherits(plot,"ggchord_fusion")) ggchord_stop("export_ggchord_fusion(): expected a Fusion")
  if (!is.numeric(width) || !is.numeric(height) || length(width)!=1 || length(height)!=1 ||
      any(!is.finite(c(width,height))) || min(width,height)<=0) ggchord_stop("Fusion: width/height must be positive inches")
  close <- ggchord_measurement_device(width=width,height=height)
  on.exit(close())
  grid::grid.newpage()
  ggchord_build_fusion(plot)$export
}

#' @export
print.ggchord_fusion <- function(x, ...) {
  ggplot2::set_last_plot(x)
  grid::grid.newpage()
  grid::grid.draw(x)
  invisible(x)
}

#' @export
#' @importFrom grid grid.draw
grid.draw.ggchord_fusion <- function(x, recording = TRUE) {
  grid::grid.draw(ggchord_build_fusion(x)$grob, recording = recording)
  invisible(x)
}

#' @export
ggplot_build.ggchord_fusion <- function(plot, ...) {
  ggchord_stop("Fusion is a composition: use grid.draw(), ggsave() or export_ggchord_fusion(), not ggplot_build()")
}

# The default cross-panel ruler uses its final orientation and real text metrics.
ggchord_fusion_identity_guide <- function(plot, styles) {
  theme <- plot$theme
  explicit <- list()
  for (style in styles) if (inherits(style,"theme")) {
    theme <- theme + style
    for (field in c("legend.key.width","legend.key.height","legend.key.size"))
      if (!is.null(style[[field]])) explicit[[field]] <- style[[field]]
  }
  horizontal <- (theme$legend.position %||% "right") %in% c("top","bottom")
  text <- ggplot2::calc_element("legend.text",ggplot2::theme_get()+theme)
  values <- c(0,50,80,90,95,100)
  dimensions <- if(inherits(text,"element_blank")) rep(0,length(values)) else
    vapply(as.character(values),function(label) {
      grob <- grid::textGrob(label,gp=grid::gpar(fontsize=text$size %||% 9,
        fontfamily=text$family %||% "",fontface=text$face %||% "plain"))
      if(horizontal) grid::convertWidth(grid::grobWidth(grob),"mm",TRUE) else
        grid::convertHeight(grid::grobHeight(grob),"mm",TRUE)
    },numeric(1))
  required <- max(((head(dimensions,-1)+tail(dimensions,-1))/2+1.2)/(diff(values)/100))
  field <- if(horizontal) "legend.key.width" else "legend.key.height"
  thickness <- if(horizontal) "legend.key.height" else "legend.key.width"
  length <- explicit[[field]] %||% explicit$legend.key.size
  if(is.null(length)) length <- grid::unit(max(required,if(horizontal) 100 else 85),"mm") else
    if(grid::is.unit(length) && grid::convertUnit(length,"mm",valueOnly=TRUE)<required)
      warning("Fusion Identity labels exceed the explicitly fixed guide length; increase its length or reduce its explicit text size",call.=FALSE)
  width <- explicit[[thickness]] %||% explicit$legend.key.size %||% grid::unit(4,"mm")
  guide_theme <- do.call(ggplot2::theme,stats::setNames(list(length,width),c(field,thickness)))
  ggplot2::guide_coloursteps(even.steps=FALSE,show.limits=TRUE,
    direction=if(horizontal) "horizontal" else "vertical",theme=guide_theme)
}
