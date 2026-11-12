# Scalar internal references from 665651c, used only by explicit geometry equivalence tools.
reference_route_shortest <- function (start, end, boxes, canvas, clearance = 0.6, barriers = NULL) 
{
    if (any(c(start, end) < 0) || any(c(start[1], end[1]) > canvas[1]) || 
        any(c(start[2], end[2]) > canvas[2])) 
        return(NULL)
    padded <- ggchord_route_inflate(boxes, clearance)
    if (any(ggchord_route_segment_hits(start, start, padded)) || 
        any(ggchord_route_segment_hits(end, end, padded))) 
        return(NULL)
    crosses <- function(a, b) {
        any(ggchord_route_segment_crossings(a, b, barriers))
    }
    if (!any(ggchord_route_segment_hits(start, end, padded)) && 
        !crosses(start, end)) 
        return(rbind(start, end))
    corners <- do.call(rbind, lapply(seq_len(nrow(padded)), function(i) as.matrix(expand.grid(x = c(padded$xmin[i], 
        padded$xmax[i]), y = c(padded$ymin[i], padded$ymax[i])))))
    if (!is.null(barriers) && nrow(barriers)) {
        endpoints <- rbind(as.matrix(barriers[c("x1", "y1")]), 
            as.matrix(barriers[c("x2", "y2")]))
        near <- do.call(rbind, lapply(seq_len(nrow(endpoints)), 
            function(i) sweep(as.matrix(expand.grid(x = c(-1, 
                1) * clearance, y = c(-1, 1) * clearance)), 2, 
                endpoints[i, ], "+")))
        corners <- rbind(corners, near)
    }
    corners <- corners[corners[, 1] >= clearance & corners[, 
        1] <= canvas[1] - clearance & corners[, 2] >= clearance & 
        corners[, 2] <= canvas[2] - clearance, , drop = FALSE]
    nodes <- unique(rbind(start, end, corners))
    distance <- rep(Inf, nrow(nodes))
    distance[1L] <- 0
    previous <- rep(NA_integer_, nrow(nodes))
    visited <- rep(FALSE, nrow(nodes))
    for (step in seq_len(nrow(nodes))) {
        available <- which(!visited)
        at <- available[which.min(distance[available])]
        if (!is.finite(distance[at])) 
            return(NULL)
        if (at == 2L) 
            break
        visited[at] <- TRUE
        for (to in which(!visited)) {
            cost <- sqrt(sum((nodes[to, ] - nodes[at, ])^2))
            if (distance[at] + cost >= distance[to] || any(ggchord_route_segment_hits(nodes[at, 
                ], nodes[to, ], padded)) || crosses(nodes[at, 
                ], nodes[to, ])) 
                next
            distance[to] <- distance[at] + cost
            previous[to] <- at
        }
    }
    if (!is.finite(distance[2L])) 
        return(NULL)
    route <- 2L
    while (route[1L] != 1L) route <- c(previous[route[1L]], route)
    nodes[route, , drop = FALSE]
}
reference_rail_crossings <- function (tasks, rows, edge_x, edge_y) 
{
    if (length(rows) < 2L) 
        return(0L)
    crossings <- 0L
    for (a in seq_len(length(rows) - 1L)) {
        i <- rows[a]
        first <- data.frame(x1 = tasks$route_x[i], y1 = tasks$route_y[i], 
            x2 = edge_x, y2 = edge_y[a])
        for (b in seq.int(a + 1L, length(rows))) {
            j <- rows[b]
            if (abs(tasks$route_x[i] - tasks$route_x[j]) < 1e-08 && 
                abs(tasks$route_y[i] - tasks$route_y[j]) < 1e-08) 
                next
            second <- data.frame(x1 = tasks$route_x[j], y1 = tasks$route_y[j], 
                x2 = edge_x, y2 = edge_y[b])
            crossings <- crossings + as.integer(ggchord_external_leader_crosses(first, 
                second))
        }
    }
    crossings
}
