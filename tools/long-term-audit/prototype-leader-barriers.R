# Isolated route research on recorded physical geometry, not package behaviour.
# No source feature, text position, label or scientific endpoint is changed.
pkgload::load_all(quiet = TRUE)
shortest <- function(start, end, boxes, barriers, canvas, clearance = .25) {
  padded <- ggchord_route_inflate(boxes, clearance)
  corners <- do.call(rbind, lapply(seq_len(nrow(padded)), function(i)
    as.matrix(expand.grid(x = c(padded$xmin[i], padded$xmax[i]),
      y = c(padded$ymin[i], padded$ymax[i])))))
  endpoints <- rbind(as.matrix(barriers[c('x0', 'y0')]), as.matrix(barriers[c('x1', 'y1')]))
  near <- do.call(rbind, lapply(seq_len(nrow(endpoints)), function(i)
    sweep(as.matrix(expand.grid(x = c(-1, 1)*clearance,
      y = c(-1, 1)*clearance)), 2, endpoints[i, ], '+')))
  nodes <- unique(rbind(start, end, corners, near))
  nodes <- nodes[nodes[, 1] >= 0 & nodes[, 1] <= canvas[1] &
    nodes[, 2] >= 0 & nodes[, 2] <= canvas[2], , drop = FALSE]
  if (!all(nodes[1, ] == start) || !all(nodes[2, ] == end)) return(NULL)
  blocked <- function(a, b) any(ggchord_route_segment_hits(a, b, padded)) ||
    any(ggchord_segments_cross(a[1], a[2], b[1], b[2],
      barriers$x0, barriers$y0, barriers$x1, barriers$y1))
  distance <- rep(Inf, nrow(nodes)); distance[1] <- 0
  previous <- rep(NA_integer_, nrow(nodes)); visited <- rep(FALSE, nrow(nodes))
  for (step in seq_len(nrow(nodes))) {
    available <- which(!visited); at <- available[which.min(distance[available])]
    if (!is.finite(distance[at])) return(NULL)
    if (at == 2L) break
    visited[at] <- TRUE
    for (to in which(!visited)) {
      cost <- sqrt(sum((nodes[to, ]-nodes[at, ])^2))
      if (distance[at]+cost >= distance[to] || blocked(nodes[at, ], nodes[to, ])) next
      distance[to] <- distance[at]+cost; previous[to] <- at
    }
  }
  if (!is.finite(distance[2L])) return(NULL)
  route <- 2L
  while (route[1L] != 1L) route <- c(previous[route[1L]], route)
  nodes[route, , drop = FALSE]
}
count <- function(leaders) {
  ids <- unique(leaders$leader_id)
  if (length(ids) < 2) return(0L)
  pairs <- combn(ids, 2)
  sum(vapply(seq_len(ncol(pairs)), function(i) {
    a <- leaders[leaders$leader_id == pairs[1, i], ]; b <- leaders[leaders$leader_id == pairs[2, i], ]
    any(vapply(seq_len(nrow(a)), function(j) any(ggchord_segments_cross(
      a$x0[j], a$y0[j], a$x1[j], a$y1[j], b$x0, b$y0, b$x1, b$y1)), logical(1)))
  }, logical(1)))
}
input <- jsonlite::fromJSON('tools/long-term-audit/v0.27.0-label-size-results.json')
output <- list()
for (name in c('pTRE_Tight_BI/compact', 'pETDuet_1/compact')) {
  case <- input$observations[[name]]
  leaders <- case$quality$leaders
  text <- case$quality$text
  visibility <- case$visibility
  before <- leaders; attempts <- list()
  candidates <- visibility[!is.na(visibility$side) & visibility$side == 'feature', ]
  for (i in seq_len(nrow(candidates))) {
    source <- candidates[i, ]
    own <- leaders$layer_id == source$layer_id & leaders$source_row == source$source_row
    path <- leaders[own, ]; if (!nrow(path)) next
    baseline_count <- count(leaders)
    others <- leaders[!own, ]
    own_text <- text$layer_id == source$layer_id & !is.na(text$source_row) & text$source_row == source$source_row
    boxes <- text[!own_text, c('xmin_mm', 'xmax_mm', 'ymin_mm', 'ymax_mm')]
    names(boxes) <- c('xmin', 'xmax', 'ymin', 'ymax')
    route <- shortest(c(path$x0[1], path$y0[1]), c(tail(path$x1, 1), tail(path$y1, 1)),
      boxes, others, c(case$width_in, case$height_in)*25.4)
    accepted <- FALSE
    if (!is.null(route)) {
      replacement <- path[rep(1, nrow(route)-1), ]
      replacement$x0 <- head(route[, 1], -1); replacement$y0 <- head(route[, 2], -1)
      replacement$x1 <- tail(route[, 1], -1); replacement$y1 <- tail(route[, 2], -1)
      trial <- rbind(others, replacement)
      old_length <- sum(sqrt((path$x1-path$x0)^2+(path$y1-path$y0)^2))
      new_length <- sum(sqrt(rowSums((route[-1, , drop=FALSE]-route[-nrow(route), , drop=FALSE])^2)))
      accepted <- count(trial) < baseline_count && new_length <= old_length*3+5
      if (accepted) leaders <- trial
    }
    attempts[[as.character(source$source_row)]] <- list(route_found = !is.null(route), accepted = accepted)
    cat(name, source$source_row, 'accepted:', accepted, '\n'); flush.console()
  }
  output[[name]] <- list(before_count = count(before), after_count = count(leaders),
    before = before, after = leaders, attempts = attempts)
}
jsonlite::write_json(list(scope = 'Route-only prototype using recorded mm geometry; not a public layout result, no full render certification',
  cases = output), 'tools/long-term-audit/v0.27.0-leader-barrier-prototype.json',
  pretty = TRUE, auto_unbox = TRUE, na = 'null')
