# Dense circular maps share measured tangent space between annotation classes.
# A crowded upper arc may use several short rows; sparse maps keep their
# nearby polar layout.

ggchord_perimeter_side <- function(x, y) {
  radius <- sqrt(x^2 + y^2)
  if (!is.finite(radius) || radius <= 1e-10) return("top")
  x <- x / radius
  y <- y / radius
  if (y < -.75 && abs(x) < .7) return("bottom")
  if (x > .215) return("right")
  if (x < -.05) return("left")
  if (y >= 0) "top" else "bottom"
}

ggchord_perimeter_tasks <- function(registry, units_per_inch) {
  tasks <- list()
  for (id in names(registry)) {
    sites <- registry[[id]]$restriction_site
    if (is.data.frame(sites) && nrow(sites)) {
      rows <- which(sites$restriction_component %in% "label" &
        !is.na(sites$label) & nzchar(sites$label))
      for (row in rows) {
        leaders <- which(sites$restriction_component %in% "leader" &
          sites$source_row == sites$source_row[row])
        if (!length(leaders)) next
        root <- leaders[1L]
        tasks[[length(tasks) + 1L]] <- data.frame(
          id = id, kind = "restriction", row = row,
          source = sites$source_row[row],
          root_x = sites$x[root], root_y = sites$y[root],
          text = sites$label[row], text_x = sites$x[row],
          text_y = sites$y[row], size = sites$size[row],
          angle = 0, hjust = .5, vjust = .5,
          family = if ("family" %in% names(sites))
            as.character(sites$family[row]) else "",
          fontface = "plain", lineheight = 1.2,
          enzyme_label = if ("enzyme_label" %in% names(sites))
            as.character(sites$enzyme_label[row]) else "",
          coordinate_label = if ("coordinate_label" %in% names(sites))
            as.character(sites$coordinate_label[row]) else "",
          enzyme_fontface = if ("enzyme_fontface" %in% names(sites))
            as.character(sites$enzyme_fontface[row]) else "plain",
          stringsAsFactors = FALSE
        )
      }
    }
    geometry <- registry[[id]]$gene_label_repel
    if (!is.data.frame(geometry) || !nrow(geometry)) next
    external <- if ("feature_label_mode" %in% names(geometry)) {
      geometry$feature_label_mode %in% "external"
    } else rep(FALSE, nrow(geometry))
    if ("annotation_class" %in% names(geometry)) {
      external <- external | geometry$annotation_class %in% "primer"
    }
    rows <- which(geometry$.component %in% "text" & external &
      !is.na(geometry$label) & nzchar(geometry$label))
    for (row in rows) {
      leaders <- which(geometry$.component %in% "segment" &
        geometry$source_row == geometry$source_row[row])
      if (!length(leaders)) next
      root <- leaders[which.max(
        (geometry$x[leaders] - geometry$text_x[row])^2 +
          (geometry$y[leaders] - geometry$text_y[row])^2
      )]
      label <- geometry$label[row]
      tasks[[length(tasks) + 1L]] <- data.frame(
        id = id, kind = if ("annotation_class" %in% names(geometry) &&
            geometry$annotation_class[row] %in% "primer") "primer" else
            "feature",
        row = row, source = geometry$source_row[row],
        root_x = geometry$x[root], root_y = geometry$y[root],
        text = label, text_x = geometry$text_x[row],
        text_y = geometry$text_y[row], size = geometry$size[row],
        angle = 0, hjust = .5, vjust = .5,
        family = if ("family" %in% names(geometry))
          as.character(geometry$family[row]) else "",
        fontface = if ("fontface" %in% names(geometry)) {
          face <- geometry$fontface[row]
          if (length(face) != 1L || is.na(face) ||
              (is.character(face) && !nzchar(face))) {
            "plain"
          } else if (is.numeric(face)) {
            c("plain", "bold", "italic", "bold.italic")[
              max(1L, min(4L, as.integer(face)))
            ]
          } else as.character(face)
        } else "plain",
        lineheight = if ("lineheight" %in% names(geometry))
          as.numeric(geometry$lineheight[row]) else 1.2,
        enzyme_label = "", coordinate_label = "",
        enzyme_fontface = "plain",
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(tasks)) return(data.frame())
  tasks <- do.call(rbind, tasks)
  tasks$family[is.na(tasks$family)] <- ""
  tasks$enzyme_fontface[is.na(tasks$enzyme_fontface) |
    !nzchar(tasks$enzyme_fontface)] <- "plain"
  tasks$lineheight[!is.finite(tasks$lineheight) |
    tasks$lineheight <= 0] <- 1.2
  boxes <- ggchord_text_boxes(
    tasks, units_per_inch = units_per_inch, box_padding = .035
  )
  tasks$w <- boxes$bw
  tasks$h <- boxes$bh
  restriction <- which(tasks$kind == "restriction")
  if (length(restriction)) {
    enzyme <- ggchord_restriction_text_metrics(
      tasks$enzyme_label[restriction], tasks$size[restriction],
      units_per_inch, parse = FALSE,
      fontface = tasks$enzyme_fontface[restriction]
    )
    coordinate <- ggchord_restriction_text_metrics(
      tasks$coordinate_label[restriction], tasks$size[restriction],
      units_per_inch, parse = FALSE
    )
    both <- nzchar(tasks$enzyme_label[restriction]) &
      nzchar(tasks$coordinate_label[restriction])
    tasks$w[restriction] <- enzyme$width + coordinate$width +
      ifelse(both, .65 / 25.4 * units_per_inch, 0) +
      .07 * units_per_inch
    tasks$h[restriction] <- pmax(enzyme$height, coordinate$height) +
      .07 * units_per_inch
  }
  callout <- which(tasks$kind == "feature")
  if (length(callout)) {
    line_height <- tasks$size[callout] / 25.4 *
      tasks$lineheight[callout] * units_per_inch
    tasks$w[callout] <- tasks$w[callout] + .24 * line_height
    tasks$h[callout] <- tasks$h[callout] + .24 * line_height
  }
  tasks$pack_w <- tasks$w * ifelse(
    tasks$kind == "restriction", 1.16, 1.08
  ) + .012
  tasks$side <- mapply(ggchord_perimeter_side,
    tasks$root_x, tasks$root_y, USE.NAMES = FALSE)
  root_radius <- sqrt(tasks$root_x^2 + tasks$root_y^2)
  root_radius[!is.finite(root_radius) | root_radius < 1e-8] <- 1
  tasks$direction_x <- tasks$root_x / root_radius
  tasks$direction_y <- tasks$root_y / root_radius
  tasks
}

ggchord_perimeter_rail_crossings <- function(tasks, rows, edge_x, edge_y) {
  if (length(rows) < 2L) return(0L)
  crossings <- 0L
  for (a in seq_len(length(rows) - 1L)) {
    i <- rows[a]
    first <- data.frame(x1 = tasks$root_x[i], y1 = tasks$root_y[i],
      x2 = edge_x, y2 = edge_y[a])
    for (b in seq.int(a + 1L, length(rows))) {
      j <- rows[b]
      if (abs(tasks$root_x[i] - tasks$root_x[j]) < 1e-8 &&
          abs(tasks$root_y[i] - tasks$root_y[j]) < 1e-8) next
      second <- data.frame(x1 = tasks$root_x[j], y1 = tasks$root_y[j],
        x2 = edge_x, y2 = edge_y[b])
      crossings <- crossings + as.integer(
        ggchord_external_leader_crosses(first, second))
    }
  }
  crossings
}

ggchord_perimeter_line_hits_box <- function(line, box) {
  t_min <- 1e-5
  t_max <- 1 - 1e-5
  for (axis in c("x", "y")) {
    start <- line[[paste0("root_", axis)]]
    delta <- line[[paste0("edge_", axis)]] - start
    low <- box[[paste0("target_", axis)]] -
      box[[if (axis == "x") "w" else "h"]] / 2
    high <- box[[paste0("target_", axis)]] +
      box[[if (axis == "x") "w" else "h"]] / 2
    if (abs(delta) <= 1e-9) {
      if (start <= low || start >= high) return(FALSE)
    } else {
      hits <- sort(c((low - start) / delta,
        (high - start) / delta))
      t_min <- max(t_min, hits[1L])
      t_max <- min(t_max, hits[2L])
    }
  }
  t_min + 1e-7 < t_max
}

ggchord_perimeter_pair_blocked <- function(a, b, gap) {
  boxes_meet <- abs(a$target_x - b$target_x) <
    (a$w + b$w) / 2 + gap &&
    abs(a$target_y - b$target_y) < (a$h + b$h) / 2 + gap
  if (boxes_meet) return(TRUE)
  if (ggchord_perimeter_line_hits_box(a, b) ||
      ggchord_perimeter_line_hits_box(b, a)) return(TRUE)
  rx <- a$edge_x - a$root_x
  ry <- a$edge_y - a$root_y
  sx <- b$edge_x - b$root_x
  sy <- b$edge_y - b$root_y
  denominator <- rx * sy - ry * sx
  if (!is.finite(denominator) || abs(denominator) <= 1e-9) {
    return(FALSE)
  }
  dx <- b$root_x - a$root_x
  dy <- b$root_y - a$root_y
  t <- (dx * sy - dy * sx) / denominator
  u <- (dx * ry - dy * rx) / denominator
  t > 1e-9 && t < 1 - 1e-9 && u > 1e-9 && u < 1 - 1e-9
}

# Compress a crowded upper arc into several nearby rows. The root order is
# retained along the tangent; a bounded search assigns height only where a
# real text box or leader needs the extra room. A single horizontal row would
# need the sum of all text widths and turn a local cluster into long rays.
ggchord_perimeter_top_fan <- function(tasks, rows, side_rows,
                                     backbone_outer, gap) {
  if (length(rows) < 8L) return(NULL)
  order <- order(tasks$target_x[rows], tasks$root_x[rows])
  rows <- rows[order]
  span <- diff(range(tasks$target_x[rows]))
  side_limit <- if (length(side_rows)) {
    min(abs(tasks$edge_x[side_rows]), na.rm = TRUE) - .20
  } else Inf
  half_width <- min(.72 * backbone_outer, side_limit)
  if (!is.finite(half_width) || half_width <= .4 ||
      span <= 2.3 * half_width) return(NULL)
  n <- length(rows)
  x <- seq(-half_width, half_width, length.out = n)
  base <- max(backbone_outer + .14,
    max(tasks$root_y[rows], na.rm = TRUE) + .08)
  step <- max(.055, .18 * max(tasks$h[rows], na.rm = TRUE))
  levels <- seq_len(17L) - 1L
  fan_gap <- min(gap, .004)
  best <- NULL
  for (base_offset in c(0, step, 2 * step)) {
    candidates <- lapply(seq_len(n), function(k) {
      lapply(levels, function(level) {
        candidate <- tasks[rows[k], , drop = FALSE]
        candidate$target_x <- x[k]
        candidate$edge_x <- x[k]
        candidate$edge_y <- base + base_offset + level * step
        candidate$target_y <- candidate$edge_y + candidate$h / 2
        candidate
      })
    })
    allowed <- lapply(candidates, function(options) {
      vapply(options, function(candidate) {
        !any(vapply(side_rows, function(j) {
          ggchord_perimeter_pair_blocked(candidate,
            tasks[j, , drop = FALSE], fan_gap)
        }, logical(1L)))
      }, logical(1L))
    })
    blocked <- vector("list", n)
    for (k in seq_len(n)) {
      blocked[[k]] <- vector("list", k - 1L)
      if (k == 1L) next
      for (j in seq_len(k - 1L)) {
        blocked[[k]][[j]] <- outer(seq_along(levels),
          seq_along(levels), Vectorize(function(a, b) {
            ggchord_perimeter_pair_blocked(candidates[[k]][[a]],
              candidates[[j]][[b]], fan_gap)
          }))
      }
    }
    beam <- list(list(levels = integer(), max_y = -Inf,
      total_y = 0))
    for (k in seq_len(n)) {
      next_beam <- list()
      for (state in beam) {
        for (option in which(allowed[[k]])) {
          if (k > 1L && any(vapply(seq_len(k - 1L), function(j) {
              blocked[[k]][[j]][option, state$levels[j]]
            }, logical(1L)))) next
          y <- candidates[[k]][[option]]$target_y
          next_beam[[length(next_beam) + 1L]] <- list(
            levels = c(state$levels, option),
            max_y = max(state$max_y, y),
            total_y = state$total_y + y
          )
        }
      }
      if (!length(next_beam)) {
        beam <- list()
        break
      }
      rank <- order(vapply(next_beam, function(state) {
        10 * state$max_y + state$total_y / k
      }, numeric(1L)))
      beam <- next_beam[head(rank, 180L)]
    }
    if (!length(beam)) next
    chosen <- beam[[1L]]
    score <- 10 * chosen$max_y + chosen$total_y / n
    if (is.null(best) || score < best$score) {
      best <- list(score = score, rows = rows,
        selected = lapply(seq_len(n), function(k) {
          candidates[[k]][[chosen$levels[k]]]
        }))
    }
  }
  best
}

ggchord_place_dense_perimeter <- function(registry, layout,
                                          tasks = NULL,
                                          allow_top_fan = TRUE) {
  units_per_inch <- layout$text_units_per_inch %||% .30
  if (is.null(tasks)) {
    tasks <- ggchord_perimeter_tasks(registry, units_per_inch)
  }
  if (!nrow(tasks)) return(registry)
  backbone_outer <- if (is.data.frame(layout$backbone_bounds) &&
      nrow(layout$backbone_bounds)) {
    max(layout$backbone_bounds$outer, na.rm = TRUE)
  } else 1
  if (!is.finite(backbone_outer)) backbone_outer <- 1

  # A tangent is the axis along which physical label order is monotone.
  # Packing every annotation class together prevents primers from being
  # squeezed through an already fixed enzyme fan.
  tasks$target_x <- NA_real_
  tasks$target_y <- NA_real_
  tasks$edge_x <- NA_real_
  tasks$edge_y <- NA_real_
  gap <- max(.014, .045 * units_per_inch)
  for (side in c("right", "left")) {
    rows <- which(tasks$side == side)
    if (!length(rows)) next
    # Restriction ticks and primer tips can start at slightly different
    # radii. A small physical-y term breaks near-equal angular ties without
    # reversing the broader order around the right side of the circle.
    axis <- ggchord_pack_label_axis(
      tasks$root_y[rows], if (side == "left")
        tasks$root_y[rows] else
          tasks$direction_y[rows] + .15 *
            (tasks$root_y[rows] - tasks$direction_y[rows]),
      tasks$h[rows] / 2, tasks$h[rows] / 2, gap = gap
    )
    boundary <- max(backbone_outer + .12,
      max(abs(tasks$root_x[rows]), na.rm = TRUE) + .07)
    boundary <- if (side == "right") boundary else -boundary
    if (ggchord_perimeter_rail_crossings(
        tasks, rows, boundary, axis) > 0L) {
      step <- max(.008, .04 * units_per_inch)
      shifts <- as.vector(rbind(-seq_len(60L) * step,
        seq_len(60L) * step))
      for (shift in shifts) {
        if (ggchord_perimeter_rail_crossings(
            tasks, rows, boundary, axis + shift) == 0L) {
          axis <- axis + shift
          break
        }
      }
    }
    tasks$target_x[rows] <- boundary +
      if (side == "right") tasks$w[rows] / 2 else -tasks$w[rows] / 2
    tasks$target_y[rows] <- axis
    tasks$edge_x[rows] <- boundary
    tasks$edge_y[rows] <- axis
  }
  side_rows <- which(tasks$side %in% c("right", "left"))
  for (side in c("top", "bottom")) {
    rows <- which(tasks$side == side)
    if (!length(rows)) next
    axis <- ggchord_pack_label_axis(
      tasks$root_x[rows], tasks$direction_x[rows],
      tasks$pack_w[rows] / 2, tasks$pack_w[rows] / 2,
      gap = gap
    )
    boundary <- max(backbone_outer + .14,
      max(abs(tasks$root_y[rows]), na.rm = TRUE) + .08)
    # A small, dense upper arc uses height only where measured text or leader
    # conflicts require it. The ordinary one-row route remains the fallback.
    if (side == "top") {
      proposed <- tasks
      proposed$target_x[rows] <- axis
      proposed$edge_x[rows] <- axis
      fan <- if (allow_top_fan) {
        ggchord_perimeter_top_fan(proposed, rows, side_rows,
          backbone_outer, gap)
      } else NULL
      if (!is.null(fan)) {
        for (k in seq_along(fan$rows)) {
          i <- fan$rows[k]
          candidate <- fan$selected[[k]]
          tasks$target_x[i] <- candidate$target_x
          tasks$target_y[i] <- candidate$target_y
          tasks$edge_x[i] <- candidate$edge_x
          tasks$edge_y[i] <- candidate$edge_y
        }
        next
      }
      # Very long composite enzyme names can leave no measured short-fan
      # placement. Reassign the upper-right roots to the adjacent tangent
      # queue and solve the whole perimeter again, before writing any geom.
      if (allow_top_fan && length(rows) >= 8L &&
          diff(range(axis)) > 2.3 * .72 * backbone_outer) {
        retry <- tasks
        move <- retry$side == "top" & retry$direction_x > .13
        if (any(move)) {
          retry$side[move] <- "right"
          return(ggchord_place_dense_perimeter(registry, layout,
            tasks = retry, allow_top_fan = FALSE))
        }
      }
      if (length(side_rows)) {
        boundary <- max(boundary,
          max(tasks$edge_y[side_rows] + tasks$h[side_rows] / 2,
            na.rm = TRUE) + gap)
      }
    }
    local_boundary <- rep(boundary, length(rows))
    if (length(side_rows)) {
      for (k in seq_along(rows)) {
        i <- rows[k]
        x <- axis[k]
        intersects_x <- abs(x - tasks$target_x[side_rows]) <
          (tasks$w[i] + tasks$w[side_rows]) / 2 + gap
        if (any(intersects_x)) {
          occupied <- side_rows[intersects_x]
          required <- if (side == "top") {
            max(tasks$target_y[occupied] + tasks$h[occupied] / 2 + gap)
          } else {
            max(-tasks$target_y[occupied] + tasks$h[occupied] / 2 + gap)
          }
          local_boundary[k] <- max(local_boundary[k], required)
        }
      }
    }
    # Every leader ends on the circle-facing edge of this one row. Stacking
    # rows would send an upper-row leader through lower-row glyphs, so width
    # pressure is resolved along the tangent while retaining root order. A
    # shallow outer arch restores a circular reading rhythm: middle labels
    # sit farther from the backbone, while the outer labels terminate early.
    spread <- diff(range(axis))
    arch_height <- if (length(rows) >= 4L && spread > 1e-8) {
      min(.23, .025 * (length(rows) - 1L))
    } else 0
    centre <- mean(range(axis))
    arch <- if (arch_height > 0) {
      arch_height * (1 - pmin(1, abs(axis - centre) /
        (spread / 2))^2)
    } else rep(0, length(rows))
    boundary <- if (side == "top") local_boundary else -local_boundary
    tasks$target_x[rows] <- axis
    tasks$target_y[rows] <- boundary +
      if (side == "top") arch + tasks$h[rows] / 2 else
        -arch - tasks$h[rows] / 2
    tasks$edge_x[rows] <- axis
    tasks$edge_y[rows] <- boundary +
      if (side == "top") arch else -arch
  }

  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i, , drop = FALSE]
    id <- task$id
    row <- task$row
    if (task$kind == "restriction") {
      sites <- registry[[id]]$restriction_site
      sites$x[row] <- task$target_x
      sites$y[row] <- task$target_y
      sites$label_attachment_x[row] <- task$edge_x
      sites$label_attachment_y[row] <- task$edge_y
      sites$label_contour_radius[row] <- sqrt(
        task$edge_x^2 + task$edge_y^2
      )
      sites$label_layout[row] <- "local_tangent"
      sites$label_connection_side[row] <- switch(task$side,
        right = "left", left = "right", top = "bottom", bottom = "top")
      sites$outer_track[row] <- 1L
      leaders <- which(sites$restriction_component %in% "leader" &
        sites$source_row == task$source)
      if (length(leaders) >= 2L) {
        sites$x[leaders[-1L]] <- task$edge_x
        sites$y[leaders[-1L]] <- task$edge_y
      }
      registry[[id]]$restriction_site <- sites
    } else {
      geometry <- registry[[id]]$gene_label_repel
      for (nm in c("outer_annotation_region", "outer_spill_reason")) {
        if (!nm %in% names(geometry)) geometry[[nm]] <- NA_character_
      }
      for (nm in c("outer_track", "outer_slot", "dominant_outer_band",
          "leader_crossing_count")) {
        if (!nm %in% names(geometry)) geometry[[nm]] <- NA_integer_
      }
      route <- data.frame(
        x1 = task$root_x, y1 = task$root_y,
        x2 = task$edge_x, y2 = task$edge_y,
        move_endpoint = FALSE
      )
      root_radius <- sqrt(task$root_x^2 + task$root_y^2)
      if (is.finite(root_radius) && root_radius > 1e-8 &&
          root_radius < backbone_outer - .05) {
        # A callout originating on an inner feature must pass the backbone
        # before it joins the exterior fan. A short radial shoulder preserves
        # the source angle and keeps it clear of nearby restriction ticks.
        bend <- c(task$root_x, task$root_y) *
          (backbone_outer + .055) / root_radius
        route <- data.frame(
          x1 = c(task$root_x, bend[1L]),
          y1 = c(task$root_y, bend[2L]),
          x2 = c(bend[1L], task$edge_x),
          y2 = c(bend[2L], task$edge_y),
          move_endpoint = c(FALSE, FALSE)
        )
      }
      candidate <- list(
        row = geometry[row, , drop = FALSE], band = 1L, slot = 0L,
        collision = FALSE, box_collision = FALSE,
        leader_collision = FALSE,
        leader_segments = route
      )
      candidate$row$text_x <- task$target_x
      candidate$row$text_y <- task$target_y
      geometry <- ggchord_shift_feature_callout(geometry, row, candidate)
      # The shared queue measures a centred rectangle. Retaining the
      # natural callout's directional hjust/vjust would move the rendered
      # glyph box after the solve, creating false clearance and collisions.
      geometry$hjust[row] <- .5
      geometry$vjust[row] <- .5
      registry[[id]]$gene_label_repel <- geometry
    }
  }
  ggchord_measure_external_leader_crossings(registry)
}
