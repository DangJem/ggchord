# Joint layout operates on trained guides, without retraining scales or data.
ggchord_guide_role <- function(aesthetic) {
  if (grepl("^track_", aesthetic)) return("track")
  switch(aesthetic, seq_colour="seq", ribbon_fill="ribbon", ribbon_colour="ribbon",
    ribbon_alpha="ribbon", ribbon_linetype="ribbon", gene_fill="gene",
    feature_fill="feature", feature_shape="feature", region_fill="region", "link")
}

ggchord_guide_constraints <- function(source, params) {
  settings <- ggchord_plot_settings(source)
  explicit <- source$ggchord$theme_explicit %||% list()
  constraints <- list(position=NULL, direction=NULL, nrow=NULL, ncol=NULL,
    width=FALSE, height=FALSE)
  for (aes in params$aesthetic) {
    role <- settings$legends[[ggchord_guide_role(aes)]]
    scale <- source$scales$get_scales(aes)
    guide <- source$guides$guides[[aes]] %||% scale$guide
    user <- if (inherits(guide, "Guide") && !isTRUE(attr(guide, "ggchord_fusion_default"))) guide$params else list()
    constraints$position <- constraints$position %||% user$position %||% role$position %||%
      settings$legend$position %||% explicit$legend.position
    constraints$direction <- constraints$direction %||% user$direction %||% role$direction %||%
      settings$legend$direction %||% explicit$legend.direction
    for (name in c("nrow", "ncol")) constraints[[name]] <- constraints[[name]] %||% user[[name]]
    for (axis in c("width", "height")) {
      field <- paste0("key.",axis)
      value <- user$theme[[paste0("legend.",field)]] %||% user$theme$legend.key.size %||%
        role[[field]] %||% role$key.size %||% settings$legend[[field]] %||% settings$legend$key.size %||%
        explicit[[paste0("legend.",field)]] %||% explicit$legend.key.size
      constraints[[axis]] <- constraints[[axis]] || !is.null(value)
    }
  }
  constraints
}

ggchord_guide_size <- function(grob) c(
  width=convertWidth(grobWidth(grob),"mm",TRUE),
  height=convertHeight(grobHeight(grob),"mm",TRUE))

ggchord_guide_candidates <- function(guide, params, theme, source) {
  generated <- all(vapply(params$aesthetic, function(aes) {
    g <- source$guides$guides[[aes]] %||% source$scales$get_scales(aes)$guide
    isTRUE(attr(g, "ggchord_fusion_default"))
  }, logical(1)))
  locks <- ggchord_guide_constraints(source, params)
  position <- params$position %||% theme$legend.position %||% "right"
  if (is.numeric(position)) position <- "inside"
  direction <- params$direction %||% theme$legend.direction %||%
    if (position %in% c("top","bottom")) "horizontal" else "vertical"
  known <- class(guide)[1] %in% c("GuideLegend", "GuideColourbar", "GuideColoursteps")
  discrete <- inherits(guide,"GuideLegend") && !inherits(guide,"GuideColourbar")
  positions <- if (!known || !is.null(locks$position) || position == "inside") position else
    unique(c(position,"right","bottom","left","top"))
  columns <- if (!known || !discrete || !is.null(locks$nrow) || !is.null(locks$ncol)) list(params$ncol) else
    c(list(params$ncol), as.list(setdiff(unique(c(seq_len(min(8L,nrow(params$key))),
      ceiling(nrow(params$key)/seq_len(min(8L,nrow(params$key)))))), params$ncol %||% integer())))
  candidates <- list()
  for (where in positions) {
    directions <- if (!known || !is.null(locks$direction)) direction else
      unique(c(if(where %in% c("top","bottom")) "horizontal" else "vertical", direction))
    for (orientation in directions) for (column in columns) {
      p <- params
      p$position <- where; p$direction <- orientation
      if (discrete && known && is.null(locks$nrow) && is.null(locks$ncol)) {
        p$ncol <- column; p$nrow <- NULL
      }
      if (!discrete && known && orientation != direction && !locks$width && !locks$height) {
        width <- p$theme$legend.key.width; height <- p$theme$legend.key.height
        if (!is.null(width) && !is.null(height)) {
          p$theme$legend.key.width <- height; p$theme$legend.key.height <- width
        }
      }
      if (generated) {
        default_guide <- ggchord_fusion_identity_guide(source, list(source$ggchord$theme_explicit), direction = orientation)
        p$theme <- default_guide$params$theme
      }
      grob <- guide$draw(theme=theme, position=where, direction=orientation, params=p)
      size <- ggchord_guide_size(grob)
      if (any(!is.finite(size))) next
      changed <- as.integer(where != position) * 2 + as.integer(orientation != direction) +
        as.integer(!isTRUE(all.equal(column,params$ncol)))
      candidates[[length(candidates)+1L]] <- list(params=p,grob=grob,size=size,
        position=where,direction=orientation,change=changed,locks=locks,
        requested_position=position,known=known)
    }
  }
  candidates
}

# Measure the actual guide-box viewport in the same gtable grid as the panel.
ggchord_guide_bounds <- function(table, canvas) {
  origin <- grid::deviceLoc(unit(0,"npc"),unit(0,"npc"),valueOnly=TRUE)
  origin <- c(origin$x,origin$y)*25.4
  pushViewport(viewport(layout=grid::grid.layout(length(table$heights),length(table$widths),
    widths=table$widths,heights=table$heights,respect=table$respect)))
  on.exit(popViewport(),add=TRUE)
  rows <- lapply(which(grepl("^guide-box",table$layout$name)),function(i) {
    grob <- table$grobs[[i]]
    if (inherits(grob,"zeroGrob")) return(NULL)
    entry <- table$layout[i,]
    pushViewport(viewport(layout.pos.row=seq.int(entry$t,entry$b),
      layout.pos.col=seq.int(entry$l,entry$r)))
    on.exit(popViewport(),add=TRUE)
    if (!is.null(grob$vp)) {
      pushViewport(grob$vp)
      on.exit(popViewport(),add=TRUE)
    }
    lower <- grid::deviceLoc(unit(0,"npc"),unit(0,"npc"),valueOnly=TRUE)
    upper <- grid::deviceLoc(unit(1,"npc"),unit(1,"npc"),valueOnly=TRUE)
    bounds <- c(lower$x,upper$x,lower$y,upper$y)*25.4-rep(origin,each=2L)
    data.frame(guide=table$layout$name[i],xmin_mm=bounds[1],xmax_mm=bounds[2],
      ymin_mm=bounds[3],ymax_mm=bounds[4],
      overflow=any(bounds[c(1,3)] < -.01 | bounds[c(2,4)] > canvas+.01))
  })
  ggchord_rbind_fill(rows) %||% data.frame()
}

ggchord_guide_joint_layout <- function(built, source, evaluate = NULL, positions = NULL) {
  canvas <- ggchord_render_dimensions()
  if (any(!is.finite(canvas))) return(list(built=built,report=list(status="unavailable")))
  if (identical(built$plot$theme$legend.position, "none")) return(list(built = built,
    table = ggchord_base_gtable(built), report = list(status = "no_guides", guides = data.frame())))
  collection <- ggproto(NULL,built$plot$guides)
  guides <- lapply(collection$guides,function(g) ggproto(NULL,g))
  params <- collection$params
  if (!length(guides)) return(list(built=built,report=list(status="no_guides",guides=data.frame())))
  theme <- built$plot$theme
  # Native complete themes can leave directional spacing to inheritance.
  # package_box consumes resolved units, whereas guide drawing resolves text.
  for (name in c("legend.spacing.x", "legend.spacing.y"))
    theme[[name]] <- ggplot2::calc_element(name, theme)
  candidates <- lapply(seq_along(guides),function(i)
    ggchord_guide_candidates(guides[[i]],params[[i]],theme,source))
  if (!is.null(positions)) candidates <- lapply(candidates, function(items)
    Filter(function(item) item$position %in% positions, items))
  if (any(lengths(candidates)==0L)) return(list(built=built,report=list(status="unmeasurable_guide")))
  score <- function(choices) {
    sizes <- matrix(0,4,2,dimnames=list(c("left","right","top","bottom"),c("width","height")))
    for (side in rownames(sizes)) {
      group <- Filter(function(x) x$position==side,choices)
      if (length(group)) sizes[side,] <- ggchord_guide_size(collection$package_box(
        lapply(group,`[[`,"grob"),side,theme))
    }
    remaining <- canvas-c(sum(sizes[c("left","right"),"width"]),sum(sizes[c("top","bottom"),"height"]))
    overflow <- sum(pmax(sizes[,1]/canvas[1]-1,0))+sum(pmax(sizes[,2]/canvas[2]-1,0))+
      sum(pmax(.4-remaining/canvas,0))
    extra <- if (is.null(evaluate)) list(penalty = 0, feasible = TRUE) else evaluate(sizes)
    list(score=overflow*1e4+extra$penalty+sum(vapply(choices,`[[`,numeric(1),"change"))*.01+
        sum(sizes[,1]*sizes[,2])/prod(canvas)*.001,
      budget_ok=overflow<=1e-10 && extra$feasible,remaining=remaining)
  }
  # Bounded deterministic beam: joint candidates compete across all visible
  # guides. Preserve trained ordering, labels, keys and scale semantics.
  beam <- list(list(choices=list(),index=integer(),score=0))
  for (i in seq_along(candidates)) {
    next_beam <- list()
    for (state in beam) for (j in seq_along(candidates[[i]])) {
      chosen <- c(state$choices,list(candidates[[i]][[j]]))
      evaluation <- score(chosen)
      next_beam[[length(next_beam)+1L]] <- list(choices=chosen,index=c(state$index,j),score=evaluation$score)
    }
    ordering <- order(vapply(next_beam,`[[`,numeric(1),"score"))
    beam <- next_beam[head(ordering,16L)]
  }
  finalists <- lapply(head(beam,4L),function(state) {
    trial <- built
    trial$plot$guides <- ggproto(NULL,collection)
    trial$plot$guides$guides <- guides
    trial$plot$guides$params <- setNames(lapply(state$choices,`[[`,"params"),names(params))
    table <- ggchord_base_gtable(trial)
    bounds <- ggchord_guide_bounds(table,canvas)
    panel <- ggchord_panel_context(table,trial)
    evaluation <- score(state$choices)
    feasible <- evaluation$budget_ok && !any(bounds$overflow) && identical(panel$status,"measured")
    list(built=trial,table=table,choices=state$choices,index=state$index,bounds=bounds,
      feasible=feasible,score=state$score+if(feasible) 0 else 1e6)
  })
  winner <- finalists[[which.min(vapply(finalists,`[[`,numeric(1),"score"))]]
  rows <- lapply(seq_along(winner$choices),function(i) {
    choice <- winner$choices[[i]]; p <- choice$params
    discrete <- choice$known && !inherits(guides[[i]], "GuideColourbar")
    grid <- if (discrete) guides[[i]]$setup_params(p) else list()
    data.frame(aesthetic=paste(p$aesthetic,collapse=","),title=paste(p$title,collapse=" "),
      type=if(inherits(guides[[i]],"GuideColourbar")) "continuous" else
        if(choice$known) "discrete" else "custom",
      requested_position=choice$requested_position,final_position=choice$position,direction=choice$direction,
      rows=grid$nrow %||% NA_integer_,columns=grid$ncol %||% NA_integer_,
      n_items=nrow(p$key),width_mm=unname(choice$size[1]),height_mm=unname(choice$size[2]),
      position_locked=!is.null(choice$locks$position),direction_locked=!is.null(choice$locks$direction),
      grid_locked=!is.null(choice$locks$nrow)||!is.null(choice$locks$ncol),
      candidate_count=length(candidates[[i]]),selected_candidate=winner$index[i],
      reason=if(!choice$known) "custom_guide_preserved" else if(!winner$feasible)
        "constraints_exceed_space" else if(choice$change==0) "requested_layout_fits" else "joint_space_budget")
  })
  list(built=winner$built,table=winner$table,report=list(
    status=if(winner$feasible) "feasible" else "insufficient_space",
    guides=ggchord_rbind_fill(rows),bounds=winner$bounds,
    canvas_mm=canvas,minimum_panel_fraction=.4,
    candidate_sizes=lapply(candidates,function(items) do.call(rbind,lapply(items,function(x)
      data.frame(position=x$position,direction=x$direction,
        width_mm=unname(x$size[1]),height_mm=unname(x$size[2])))))))
}
ggchord_geometry_guides <- function(parent) {
  force(parent)
  ggproto(NULL, parent,
    train = function(self, scales, labels) {
      ggproto_parent(parent, self)$train(scales, labels)
      params <- self$params
      groups <- split(seq_along(params), vapply(params, function(p)
        paste(p$order, p$hash), character(1)))
      for (rows in groups) {
        if (length(rows) < 2L || !any(vapply(params[rows], function(p)
          any(startsWith(p$aesthetic, "ggchord_scope_")), logical(1)))) next
        palettes <- lapply(params[rows], ggchord_legend_palette)
        colour <- which(lengths(palettes) > 0L)
        if (length(unique(palettes[colour])) < 2L) next
        for (i in colour) {
          at <- rows[i]
          params[[at]]$hash <- paste0(params[[at]]$hash, "_", rlang::hash(palettes[[i]]))
        }
      }
      self$params <- params
    },
    process_layers = function(self, layers, data = NULL, theme = NULL) {
      ggproto_parent(parent, self)$process_layers(layers, data, theme)
      for (i in seq_along(self$params)) {
        if (inherits(self$guides[[i]], "GuideLegend"))
          self$params[[i]]$decor <- ggchord_legend_separate_keys(self$params[[i]]$decor)
      }
    })
}

# The native legend hash includes title and labels but omits mapped colours.
ggchord_legend_palette <- function(params) {
  aesthetic <- grep("(^|_)(colou?r|fill(_from|_to)?)$", params$aesthetic, value = TRUE)
  if (!length(aesthetic)) return(NULL)
  colours <- params$key[[aesthetic[1L]]]
  if (!length(colours)) return(NULL)
  rgba <- function(x) if (length(x)) unname(grDevices::col2rgb(x, alpha = TRUE))
  list(breaks = rgba(colours), ramp = if (is.list(params$decor)) rgba(params$decor$colour))
}

ggchord_legend_separate_keys <- function(decor) {
  if (length(decor) < 2L || !all(vapply(decor, function(x)
    is.list(x) && is.function(x$draw_key) && !is.null(x$data$.id), logical(1)))) return(decor)
  # ggplot2 retains a row for each break in each layer, including hidden keys.
  # Give the visible glyphs distinct slots within that one labelled row.
  ids <- decor[[1L]]$data$.id
  active <- vapply(decor, function(x) {
    data <- x$data
    ids %in% data$.id[data$.draw %||% TRUE]
  }, logical(length(ids)))
  active <- matrix(active, nrow = length(ids))
  count <- rowSums(active)
  if (all(count < 2L)) return(decor)
  rank <- integer(length(ids))
  for (i in seq_along(decor)) {
    rank <- rank + active[, i]
    decor[[i]]$draw_key <- ggchord_legend_place_key(decor[[i]]$draw_key, ids, rank, count)
  }
  decor
}

ggchord_legend_place_key <- function(draw, ids, rank, count) {
  force(draw); force(ids); force(rank); force(count)
  function(data, params, size) {
    glyph <- draw(data, params, size)
    row <- match(data$.id, ids)
    if (length(row) != 1L || is.na(row) || count[row] < 2L) return(glyph)
    n <- count[row]
    key <- grobTree(glyph, vp = viewport(y = 1 - (rank[row] - .5) / n,
      height = unit(1 / n, "npc"), clip = "on"))
    attr(key, "width") <- attr(glyph, "width")
    attr(key, "height") <- max(.38, attr(glyph, "height"), rep_len(size, 2L)[2L] / 10) * n
    key
  }
}
