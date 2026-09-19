# Coord-owned annotation resources for one circular sequence.

ggchord_empty_circular_annotation_registry <- function() {
  data.frame(
    registry_id = character(), side = character(), kind = character(),
    owner = character(), component = character(), source_row = integer(),
    lane = integer(), track = integer(), radius = numeric(),
    radial_width = numeric(), radius_inner = numeric(),
    radius_outer = numeric(), region = character(), sector = integer(),
    band = integer(), slot = integer(), bbox_xmin = numeric(),
    bbox_xmax = numeric(), bbox_ymin = numeric(), bbox_ymax = numeric(),
    leader_corridor = character(), leader_track_start = integer(),
    leader_track_end = integer(), stringsAsFactors = FALSE
  )
}

ggchord_circular_region <- function(x, y) {
  ifelse(abs(x) >= abs(y), ifelse(x < 0, "left", "right"),
    ifelse(y < 0, "bottom", "top"))
}

ggchord_circular_sector <- function(x, y, count = 24L) {
  angle <- (atan2(y, x) + 2 * pi) %% (2 * pi)
  as.integer(floor(angle / (2 * pi / count))) + 1L
}

ggchord_registry_row <- function(
    id, side, kind, owner, component, source_row = NA_integer_,
    lane = NA_integer_, track = NA_integer_, radius = NA_real_,
    radial_width = NA_real_, radius_inner = NA_real_,
    radius_outer = NA_real_, region = NA_character_, sector = NA_integer_,
    band = NA_integer_, slot = NA_integer_, bbox = NULL,
    leader_corridor = NA_character_, leader_track_start = NA_integer_,
    leader_track_end = NA_integer_) {
  if (is.null(bbox) || !nrow(bbox)) {
    bbox <- data.frame(xmin = NA_real_, xmax = NA_real_,
      ymin = NA_real_, ymax = NA_real_)
  }
  data.frame(
    registry_id = id, side = side, kind = kind, owner = owner,
    component = component, source_row = as.integer(source_row),
    lane = as.integer(lane), track = as.integer(track), radius = radius,
    radial_width = radial_width, radius_inner = radius_inner,
    radius_outer = radius_outer, region = region, sector = as.integer(sector),
    band = as.integer(band), slot = as.integer(slot),
    bbox_xmin = bbox$xmin[1L], bbox_xmax = bbox$xmax[1L],
    bbox_ymin = bbox$ymin[1L], bbox_ymax = bbox$ymax[1L],
    leader_corridor = leader_corridor,
    leader_track_start = as.integer(leader_track_start),
    leader_track_end = as.integer(leader_track_end),
    stringsAsFactors = FALSE
  )
}

ggchord_annotation_text_boxes <- function(x, rows, units_per_inch) {
  if (!length(rows)) return(data.frame())
  labels <- x[rows, , drop = FALSE]
  if (!"text" %in% names(labels)) labels$text <- labels$label
  if (!"text_x" %in% names(labels)) labels$text_x <- labels$x
  if (!"text_y" %in% names(labels)) labels$text_y <- labels$y
  if (!"text_angle" %in% names(labels)) {
    labels$text_angle <- labels$angle %||% 0
  }
  ggchord_text_boxes(labels, units_per_inch = units_per_inch,
    box_padding = .015)
}

# Registry joins cross geom-owned tables whose rows may legitimately be absent
# on a particular device (for example, a hidden short-feature label has no
# label-only resource).  Named-vector indexing returns length zero when the key
# itself is missing, which must never leak into an `if` condition.  Keep every
# lookup scalar and use an explicit fallback supplied by the originating label
# geometry when one is available.
ggchord_registry_track <- function(index, key, fallback = NA_integer_) {
  value <- if (length(key) == 1L && !is.na(key) && nzchar(key)) {
    unname(index[as.character(key)])
  } else integer()
  value <- value[is.finite(value)]
  if (!length(value)) {
    fallback <- as.integer(fallback)
    fallback <- fallback[is.finite(fallback)]
    value <- fallback
  }
  if (length(value)) as.integer(value[1L]) else NA_integer_
}

# Build physical inner bands from the final geometry. Position lanes remain
# local hints; only this coord pass assigns globally meaningful track IDs.
ggchord_resolve_inner_tracks <- function(layer_geometry, layout) {
  units_per_inch <- layout$text_units_per_inch %||% .30
  candidates <- list()
  feature_lookup <- list()

  seq_radii <- numeric()
  for (owner in names(layer_geometry)) {
    geometry <- layer_geometry[[owner]]$seq
    if (is.data.frame(geometry) && nrow(geometry) &&
        all(c("x", "y") %in% names(geometry))) {
      seq_radii <- c(seq_radii, sqrt(geometry$x^2 + geometry$y^2))
    }
  }
  seq_radii <- seq_radii[is.finite(seq_radii)]
  if (!length(seq_radii) && length(layout$seq_arcs)) {
    arc <- layout$seq_arcs[[1L]]
    seq_radii <- sqrt(arc$x^2 + arc$y^2)
  }
  backbone <- if (is.data.frame(layout$backbone_bounds) &&
      nrow(layout$backbone_bounds)) {
    c(min(layout$backbone_bounds$inner), max(layout$backbone_bounds$outer))
  } else if (length(seq_radii)) range(seq_radii) else c(1, 1)
  candidates[[length(candidates) + 1L]] <- data.frame(
    key = "backbone", kind = "backbone", lane = NA_integer_,
    inner = backbone[1L], outer = backbone[2L], mid = mean(backbone),
    stringsAsFactors = FALSE
  )

  for (owner in names(layer_geometry)) {
    geometry <- layer_geometry[[owner]]$gene_poly
    if (!is.data.frame(geometry) || !nrow(geometry) ||
        !all(c("x", "y", "source_row") %in% names(geometry))) next
    polygon <- if (".component" %in% names(geometry)) {
      geometry$.component %in% c("polygon", "point")
    } else rep(TRUE, nrow(geometry))
    rows <- which(polygon & is.finite(geometry$x) & is.finite(geometry$y))
    if (!length(rows)) next
    split_rows <- split(rows, geometry$source_row[rows])
    for (source in names(split_rows)) {
      index <- split_rows[[source]]
      radius <- sqrt(geometry$x[index]^2 + geometry$y[index]^2)
      lane <- if ("lane" %in% names(geometry)) {
        as.integer(geometry$lane[index[1L]])
      } else NA_integer_
      bounds <- range(radius[is.finite(radius)])
      # Arrowheads and compact fallbacks may extend by slightly different
      # amounts, but objects sharing a Position lane still consume one
      # physical band. The coord takes the union of those envelopes below.
      key <- paste("feature", lane, sprintf("%.5f", mean(bounds)),
        sep = ":")
      candidates[[length(candidates) + 1L]] <- data.frame(
        key = key, kind = "feature_band", lane = lane,
        inner = bounds[1L], outer = bounds[2L], mid = mean(bounds),
        stringsAsFactors = FALSE
      )
      feature_lookup[[paste(owner, source, sep = "\r")]] <- list(
        key = key, lane = lane, inner = bounds[1L], outer = bounds[2L]
      )
    }
  }

  label_lookup <- list()
  for (owner in names(layer_geometry)) {
    geometry <- layer_geometry[[owner]]$gene_label_repel
    if (!is.data.frame(geometry) || !nrow(geometry) ||
        !"feature_label_mode" %in% names(geometry)) next
    rows <- which(geometry$.component %in% "text" &
      geometry$feature_label_mode %in% c("inside", "adjacent") &
      !is.na(geometry$label) & nzchar(geometry$label))
    boxes <- ggchord_annotation_text_boxes(geometry, rows, units_per_inch)
    for (j in seq_along(rows)) {
      row <- rows[j]
      source <- geometry$source_row[row]
      radius <- sqrt(geometry$x[row]^2 + geometry$y[row]^2)
      feature_candidates <- Filter(function(value) {
        identical(as.integer(value$lane), as.integer(geometry$lane[row]))
      }, feature_lookup)
      feature <- if (length(feature_candidates)) {
        feature_candidates[[which.min(vapply(feature_candidates, function(x) {
          abs(mean(c(x$inner, x$outer)) - radius)
        }, numeric(1L)))]]
      } else NULL
      if (identical(geometry$feature_label_mode[row], "inside") &&
          !is.null(feature)) {
        key <- feature$key
      } else {
        half_height <- max(.004, boxes$h[j] / 2)
        bounds <- c(radius - half_height, radius + half_height)
        key <- paste("label", format(radius, digits = 10), sep = ":")
        candidates[[length(candidates) + 1L]] <- data.frame(
          key = key, kind = "label_track", lane = NA_integer_,
          inner = bounds[1L], outer = bounds[2L], mid = radius,
          stringsAsFactors = FALSE
        )
      }
      label_lookup[[paste(owner, source, sep = "\r")]] <- list(
        key = key, feature_key = feature$key %||% NA_character_,
        bbox = boxes[j, , drop = FALSE], radius = radius
      )
    }
  }

  axis <- layout$axis_ticks %||% data.frame()
  if (nrow(axis) && all(c("x0", "y0", "x1", "y1") %in% names(axis))) {
    radius <- c(sqrt(axis$x0^2 + axis$y0^2),
      sqrt(axis$x1^2 + axis$y1^2))
    radius <- radius[is.finite(radius)]
    if (length(radius)) candidates[[length(candidates) + 1L]] <- data.frame(
      key = "axis-reserved", kind = "axis_reserved", lane = NA_integer_,
      inner = min(radius), outer = max(radius), mid = mean(range(radius)),
      stringsAsFactors = FALSE
    )
  }

  resources <- ggchord_rbind_fill(candidates)
  resources <- do.call(rbind, lapply(split(resources, resources$key), function(x) {
    x$inner[1L] <- min(x$inner)
    x$outer[1L] <- max(x$outer)
    x$mid[1L] <- mean(range(c(x$inner, x$outer)))
    x[1L, , drop = FALSE]
  }))
  rownames(resources) <- NULL
  backbone_row <- resources$key == "backbone"
  resources$track <- NA_integer_
  resources$track[backbone_row] <- 0L
  inner <- which(!backbone_row)
  resources$track[inner] <- seq_along(inner)[order(order(
    -resources$outer[inner], -resources$mid[inner], resources$key[inner]
  ))]
  track_by_key <- stats::setNames(resources$track, resources$key)

  registry <- list()
  for (i in seq_len(nrow(resources))) {
    item <- resources[i, ]
    registry[[length(registry) + 1L]] <- ggchord_registry_row(
      paste0("inner-resource-", sprintf("%04d", i)), "inner", item$kind,
      ".coord", "resource", lane = item$lane, track = item$track,
      radius = item$mid, radial_width = item$outer - item$inner,
      radius_inner = item$inner, radius_outer = item$outer
    )
  }

  for (owner in names(layer_geometry)) {
    geometry <- layer_geometry[[owner]]$gene_poly
    if (is.data.frame(geometry) && nrow(geometry) &&
        "source_row" %in% names(geometry)) {
      geometry$annotation_side <- "inner"
      geometry$annotation_kind <- "feature_band"
      geometry$physical_track <- NA_integer_
      geometry$track_radius <- NA_real_
      geometry$track_inner <- NA_real_
      geometry$track_outer <- NA_real_
      for (source in unique(geometry$source_row)) {
        item <- feature_lookup[[paste(owner, source, sep = "\r")]]
        if (is.null(item)) next
        resource <- resources[match(item$key, resources$key), , drop = FALSE]
        rows <- geometry$source_row == source
        geometry$physical_track[rows] <- unname(track_by_key[item$key])
        geometry$track_radius[rows] <- resource$mid
        geometry$track_inner[rows] <- resource$inner
        geometry$track_outer[rows] <- resource$outer
      }
      layer_geometry[[owner]]$gene_poly <- geometry
    }

    geometry <- layer_geometry[[owner]]$gene_label_repel
    if (!is.data.frame(geometry) || !nrow(geometry) ||
        !"feature_label_mode" %in% names(geometry)) next
    geometry$annotation_side <- ifelse(
      geometry$feature_label_mode %in% "external", "outer", "inner")
    geometry$annotation_kind <- ifelse(geometry$.component %in% "segment",
      "leader", ifelse(geometry$.component %in% "arc_text", "label_glyph",
        "label"))
    if (!"leader_corridor" %in% names(geometry)) {
      geometry$leader_corridor <- NA_character_
    }
    if (!"leader_track_start" %in% names(geometry)) {
      geometry$leader_track_start <- NA_integer_
    }
    if (!"leader_track_end" %in% names(geometry)) {
      geometry$leader_track_end <- NA_integer_
    }
    for (source in unique(geometry$source_row)) {
      label <- label_lookup[[paste(owner, source, sep = "\r")]]
      if (is.null(label)) next
      rows <- geometry$source_row == source & geometry$annotation_side == "inner"
      if (!any(rows)) next
      source_label_track <- unique(geometry$label_track[rows])
      source_label_track <- source_label_track[is.finite(source_label_track)]
      source_feature_track <- unique(geometry$feature_track[rows])
      source_feature_track <- source_feature_track[
        is.finite(source_feature_track)
      ]
      label_track <- ggchord_registry_track(
        track_by_key, label$key, source_label_track
      )
      feature_track <- ggchord_registry_track(
        track_by_key, label$feature_key, source_feature_track
      )
      geometry$label_track[rows] <- label_track
      geometry$feature_track[rows] <- feature_track
      segment_rows <- rows & geometry$.component %in% "segment"
      corridor_valid <- isTRUE(any(segment_rows)) &&
        length(feature_track) == 1L && is.finite(feature_track) &&
        length(label_track) == 1L && is.finite(label_track) &&
        label_track > feature_track + 1L
      if (corridor_valid) {
        geometry$leader_track_start[segment_rows] <- feature_track + 1L
        geometry$leader_track_end[segment_rows] <- label_track - 1L
        geometry$leader_corridor[segment_rows] <- paste0(
          "inner:", feature_track + 1L, "-", label_track - 1L
        )
      }
      text_rows <- rows & geometry$.component %in% "text"
      if (any(text_rows)) {
        registry[[length(registry) + 1L]] <- ggchord_registry_row(
          paste0("inner-label-", owner, "-", source), "inner", "label",
          owner, "gene_label_repel", source_row = source,
          lane = geometry$lane[which(text_rows)[1L]], track = label_track,
          radius = label$radius,
          bbox = label$bbox
        )
      }
      if (corridor_valid) {
        registry[[length(registry) + 1L]] <- ggchord_registry_row(
          paste0("inner-leader-", owner, "-", source), "inner", "leader",
          owner, "gene_label_repel", source_row = source,
          track = label_track,
          leader_corridor = geometry$leader_corridor[which(segment_rows)[1L]],
          leader_track_start = feature_track + 1L,
          leader_track_end = label_track - 1L
        )
      }
    }
    layer_geometry[[owner]]$gene_label_repel <- geometry
  }

  # The centre title is a real measured obstacle, not a hard-coded radius.
  # Feature candidate generation already keeps a conservative centre reserve;
  # registering the final box here makes the actual occupied area available to
  # fitting, diagnostics and the next allocator stage without inventing a
  # concentric track for a Cartesian title block.
  for (owner in names(layer_geometry)) {
    centre <- layer_geometry[[owner]]$seq_center_label
    if (!is.data.frame(centre) || !nrow(centre) ||
        !all(c("x", "y", "label") %in% names(centre))) next
    labels <- centre
    labels$text <- labels$label
    labels$text_x <- labels$x
    labels$text_y <- labels$y
    labels$text_angle <- labels$angle %||% 0
    boxes <- ggchord_text_boxes(
      labels, units_per_inch = units_per_inch, box_padding = .025
    )
    combined <- data.frame(
      xmin = min(boxes$xmin), xmax = max(boxes$xmax),
      ymin = min(boxes$ymin), ymax = max(boxes$ymax)
    )
    registry[[length(registry) + 1L]] <- ggchord_registry_row(
      paste0("inner-centre-", owner), "inner", "center_reserved",
      owner, "seq_center_label", bbox = combined
    )
  }

  list(layer_geometry = layer_geometry,
    registry = ggchord_rbind_fill(registry))
}

ggchord_resolve_outer_annotations <- function(layer_geometry, layout) {
  units_per_inch <- layout$text_units_per_inch %||% .30
  registry <- list()

  for (owner in names(layer_geometry)) {
    for (component in c("restriction_site", "gene_label_repel")) {
      geometry <- layer_geometry[[owner]][[component]]
      if (!is.data.frame(geometry) || !nrow(geometry)) next
      rows <- if (identical(component, "restriction_site")) {
        which(geometry$.component %in% "label" & !is.na(geometry$label))
      } else if ("feature_label_mode" %in% names(geometry)) {
        which(geometry$.component %in% "text" &
          geometry$feature_label_mode %in% "external" &
          !is.na(geometry$label) & nzchar(geometry$label))
      } else integer()
      if (!length(rows)) next
      boxes <- ggchord_annotation_text_boxes(geometry, rows, units_per_inch)
      region <- ggchord_circular_region(boxes$cx, boxes$cy)
      sector <- ggchord_circular_sector(
        geometry$x[rows], geometry$y[rows]
      )
      band <- if ("outer_track" %in% names(geometry)) {
        as.integer(geometry$outer_track[rows])
      } else rep(1L, length(rows))
      band[is.na(band)] <- 1L
      keys <- paste(region, sector, band, sep = "\r")
      slot <- integer(length(rows))
      for (key in unique(keys)) {
        members <- which(keys == key)
        order_value <- order(
          if (region[members[1L]] %in% c("left", "right")) {
            boxes$cy[members]
          } else boxes$cx[members],
          geometry$anchor_position[rows[members]] %||% rows[members]
        )
        slot[members[order_value]] <- seq_along(members)
      }
      character_columns <- c(
        "annotation_side", "annotation_kind", "annotation_region",
        "leader_corridor"
      )
      integer_columns <- c(
        "annotation_sector", "annotation_band", "annotation_slot"
      )
      numeric_columns <- c(
        "bbox_xmin", "bbox_xmax", "bbox_ymin", "bbox_ymax"
      )
      for (column in setdiff(character_columns, names(geometry))) {
        geometry[[column]] <- NA_character_
      }
      for (column in setdiff(integer_columns, names(geometry))) {
        geometry[[column]] <- NA_integer_
      }
      for (column in setdiff(numeric_columns, names(geometry))) {
        geometry[[column]] <- NA_real_
      }
      geometry$annotation_side[rows] <- "outer"
      geometry$annotation_kind[rows] <- if (
        identical(component, "restriction_site")) "restriction_label" else
        "feature_callout"
      geometry$annotation_region[rows] <- region
      geometry$annotation_sector[rows] <- sector
      geometry$annotation_band[rows] <- band
      geometry$annotation_slot[rows] <- slot
      geometry$bbox_xmin[rows] <- boxes$xmin
      geometry$bbox_xmax[rows] <- boxes$xmax
      geometry$bbox_ymin[rows] <- boxes$ymin
      geometry$bbox_ymax[rows] <- boxes$ymax
      corridor <- if (identical(component, "restriction_site")) {
        as.character(geometry$cluster_id[rows])
      } else paste0("feature-sector-", sprintf("%02d", sector))
      geometry$leader_corridor[rows] <- corridor

      for (j in seq_along(rows)) {
        source <- as.integer(geometry$source_row[rows[j]])
        # Generated label fragments can lack an attached input row.  They still
        # need a stable build-local identity so their leader and registry row
        # can be joined without NA comparisons poisoning the logical mask.
        generated_source <- length(source) != 1L || is.na(source)
        if (generated_source) {
          source <- rows[j]
          geometry$source_row[rows[j]] <- source
        }
        id <- paste("outer", owner, component, source, sep = "-")
        registry[[length(registry) + 1L]] <- ggchord_registry_row(
          id, "outer", geometry$annotation_kind[rows[j]], owner, component,
          source_row = source, region = region[j], sector = sector[j],
          band = band[j], slot = slot[j], bbox = boxes[j, , drop = FALSE],
          leader_corridor = corridor[j]
        )
        leader_rows <- if (generated_source && "group" %in% names(geometry)) {
          geometry$group == geometry$group[rows[j]] &
            (geometry$.component %in% c("segment", "path"))
        } else {
          !is.na(geometry$source_row) & geometry$source_row == source &
            (geometry$.component %in% c("segment", "path"))
        }
        leader_rows[is.na(leader_rows)] <- FALSE
        if (generated_source && any(leader_rows)) {
          geometry$source_row[leader_rows] <- source
        }
        if (identical(component, "restriction_site") &&
            "restriction_component" %in% names(geometry)) {
          leader_rows <- leader_rows & geometry$restriction_component == "leader"
        }
        if (isTRUE(any(leader_rows))) {
          geometry$annotation_side[leader_rows] <- "outer"
          geometry$annotation_kind[leader_rows] <- "leader"
          geometry$annotation_region[leader_rows] <- region[j]
          geometry$annotation_sector[leader_rows] <- sector[j]
          geometry$annotation_band[leader_rows] <- band[j]
          geometry$annotation_slot[leader_rows] <- slot[j]
          geometry$leader_corridor[leader_rows] <- corridor[j]
        }
      }
      layer_geometry[[owner]][[component]] <- geometry
    }
  }
  if (!length(registry)) return(list(
    layer_geometry = layer_geometry,
    registry = ggchord_empty_circular_annotation_registry()
  ))

  registry <- ggchord_rbind_fill(registry)
  slot_key <- paste(registry$region, registry$sector, registry$band, sep = "\r")
  for (key in unique(slot_key)) {
    rows <- which(slot_key == key)
    horizontal <- registry$region[rows[1L]] %in% c("top", "bottom")
    coordinate <- if (horizontal) {
      (registry$bbox_xmin[rows] + registry$bbox_xmax[rows]) / 2
    } else {
      (registry$bbox_ymin[rows] + registry$bbox_ymax[rows]) / 2
    }
    registry$slot[rows[order(coordinate, registry$owner[rows],
      registry$source_row[rows])]] <- seq_along(rows)
  }
  # Copy the globally assigned slot back to both the label and every segment
  # in its leader corridor. This is the point at which slots become shared
  # coord resources instead of layer-local ordinal values.
  for (i in seq_len(nrow(registry))) {
    owner <- registry$owner[i]
    component <- registry$component[i]
    geometry <- layer_geometry[[owner]][[component]]
    if (!is.data.frame(geometry) || !nrow(geometry)) next
    rows <- !is.na(geometry$source_row) &
      geometry$source_row == registry$source_row[i] &
      geometry$annotation_side %in% "outer"
    geometry$annotation_slot[rows] <- registry$slot[i]
    layer_geometry[[owner]][[component]] <- geometry
  }
  list(layer_geometry = layer_geometry, registry = registry)
}

ggchord_resolve_circular_annotation_registry <- function(
    layer_geometry, layout, coord = NULL) {
  # First preserve the two geom-specific candidate systems while sharing their
  # final exterior obstacles. The following passes only assign coord resources.
  layer_geometry <- ggchord_share_external_annotations(layer_geometry, layout)
  inner <- ggchord_resolve_inner_tracks(layer_geometry, layout)
  outer <- ggchord_resolve_outer_annotations(inner$layer_geometry, layout)
  backbone <- inner$registry[
    inner$registry$kind == "backbone" &
      inner$registry$component == "resource", , drop = FALSE
  ]
  collar <- ggchord_empty_circular_annotation_registry()
  if (nrow(backbone)) {
    collar_inner <- backbone$radius_outer[1L]
    tick_radius <- numeric()
    for (owner in names(outer$layer_geometry)) {
      sites <- outer$layer_geometry[[owner]]$restriction_site
      if (!is.data.frame(sites) || !nrow(sites) ||
          !"restriction_component" %in% names(sites)) next
      tick <- sites$restriction_component %in% "tick"
      tick_radius <- c(tick_radius,
        sqrt(sites$x[tick]^2 + sites$y[tick]^2))
    }
    tick_radius <- tick_radius[is.finite(tick_radius)]
    collar_outer <- max(c(collar_inner, tick_radius))
    collar <- ggchord_registry_row(
      "outer-resource-anchor-collar", "outer", "anchor_collar",
      ".coord", "resource", radius = mean(c(collar_inner, collar_outer)),
      radial_width = collar_outer - collar_inner,
      radius_inner = collar_inner, radius_outer = collar_outer,
      region = "collar", band = 0L
    )
  }
  annotation_registry <- ggchord_rbind_fill(list(
    inner$registry, collar, outer$registry
  ))
  if (!nrow(annotation_registry)) {
    annotation_registry <- ggchord_empty_circular_annotation_registry()
  }
  list(layer_geometry = outer$layer_geometry,
    annotation_registry = annotation_registry)
}
