# Branch after ggplot2 maps scales (including after_scale), so equal raw
# categories are not confused with equal visible styles. The same geometry
# is cached for export, retaining membership without duplicating the trunk.
ggchord_branch_built <- function(built) {
  plot <- built$plot
  layout <- plot$ggchord$ref$layout
  if (is.null(layout)) return(built)
  changed <- FALSE
  for (i in seq_along(plot$layers)) {
    layer <- plot$layers[[i]]
    params <- layer$ggchord_params
    if (is.null(params) || !params$type %in% c("link", "ribbon") ||
        is.null(params$link_branch) || params$link_branch == "none") next
    data <- built$data[[i]]
    original <- layer$data
    if (!nrow(data)) next
    if (nrow(data) != nrow(original)) {
      warning("Shared links require a stat/position that preserves vertex rows; drawing separate links", call. = FALSE)
      next
    }
    # ggplot2's group numbering may change, but row order is preserved by the
    # supported identity and ggchord ribbon statistics.
    data$source_row <- original$source_row
    input <- layout$layer_inputs[[layer$ggchord_layer_id]][[layer$ggchord_type]]
    ribbon <- params$type == "ribbon"
    side <- params$link_branch
    cols <- if (ribbon) {
      if (side == "query") c("qaccver","qstart","qend") else c("saccver","sstart","send")
    } else if (side == "query") c("qaccver","qpos") else c("saccver","spos")
    style_cols <- intersect(c("PANEL", "ribbon_fill","ribbon_alpha","ribbon_colour",
      "ribbon_linetype","link_colour","link_alpha","link_linewidth","link_linetype",
      "fill","colour","alpha","linewidth","linetype"),names(data))
    row_groups <- split(seq_len(nrow(data)), data$source_row)
    items <- lapply(row_groups, function(idx) {
      row <- data$source_row[idx[1]]
      xy <- as.matrix(data[idx,c("x","y")])
      if (ribbon) {
        nq <- original$.q_n[idx[1]]; ns <- original$.s_n[idx[1]]
        if (is.null(nq) || is.null(ns)) return(NULL)
        q <- xy[seq_len(nq),,drop=FALSE]
        s <- xy[nq+50+seq_len(ns),,drop=FALSE][ns:1,,drop=FALSE]
        near <- if (side == "query") q else s
        far <- if (side == "query") s else q
      } else {
        near <- xy[if (side == "query") 1 else nrow(xy),,drop=FALSE]
        far <- xy[if (side == "query") nrow(xy) else 1,,drop=FALSE]
      }
      # Varying aesthetics along a path cannot be represented by one shared style.
      styles <- unique(data[idx,style_cols,drop=FALSE])
      if (nrow(styles) != 1L) return(list(idx=idx, key=paste0("single-",row)))
      control <- params$link_ctrl_point %||% params$ribbon_ctrl_point
      if (is.list(control)) control <- control[[min(row,length(control))]]
      signature <- list(unname(as.list(input[row,cols,drop=FALSE])), unname(as.list(styles)),
        unname(near), control)
      key <- paste(as.integer(serialize(signature,NULL,version=2)),collapse=".")
      path <- if (!ribbon) xy else {
        right <- xy[nq+seq_len(50),,drop=FALSE]
        left <- xy[nrow(xy)-seq_len(50)+1L,,drop=FALSE]
        (right+left)/2
      }
      list(idx=idx,row=row,near=near,far=far,xy=xy,key=key,
        control=control,path_length=sum(sqrt(rowSums(diff(path)^2))))
    })
    if (any(vapply(items,is.null,logical(1)))) next
    groups <- split(seq_along(items), vapply(items, `[[`, character(1), "key"))
    output <- list(); geometries <- list(); next_group <- 0L; failures <- 0L
    emit <- function(item, xy, component, members, first=TRUE,last=TRUE,preserve=FALSE) {
      next_group <<- next_group + 1L
      rows <- if (preserve) item$idx else rep(item$idx[1], nrow(xy))
      d <- data[rows,,drop=FALSE]
      d$x <- xy[,1]; d$y <- xy[,2]; d$group <- next_group
      d$.component <- component; d$.arrow_first <- first; d$.arrow_last <- last
      d$source_row <- if (component == "trunk") NA_integer_ else data$source_row[item$idx[1]]
      d$source_rows <- rep(list(as.integer(members)),nrow(d))
      output[[length(output)+1L]] <<- d
      raw <- original[rows,,drop=FALSE]
      raw$x <- xy[,1]; raw$y <- xy[,2]; raw$group <- next_group
      raw$.component <- component; raw$source_row <- d$source_row
      raw$source_rows <- d$source_rows
      # A shared trunk does not claim any one member's destination or statistics.
      if (component == "trunk") {
        protected <- c("x","y","group",".component","source_row","source_rows",cols)
        for (nm in intersect(setdiff(names(input),protected),names(raw))) raw[[nm]] <- raw[[nm]][NA_integer_]
      }
      geometries[[length(geometries)+1L]] <<- raw
    }
    for (g in groups) {
      parts <- items[g]
      if (length(parts) < 2L || is.null(parts[[1]]$near)) {
        for (item in parts) {
          idx <- item$idx
          emit(item,as.matrix(data[idx,c("x","y")]),"branch",data$source_row[idx[1]],preserve=TRUE)
        }
        next
      }
      anchor <- colMeans(parts[[1]]$near)
      targets <- do.call(rbind,lapply(parts,function(x) colMeans(x$far)))
      vectors <- targets-matrix(anchor,nrow(targets),2,byrow=TRUE)
      distances <- sqrt(rowSums(vectors^2)); direction <- colMeans(vectors)
      norm <- sqrt(sum(direction^2))
      if (!is.finite(norm) || norm < 1e-8 || min(distances) < 1e-8) {
        failures <- failures+1L
        for (item in parts) emit(item,item$xy,"branch",item$row,preserve=TRUE)
        next
      }
      delta <- direction/norm*min(vapply(parts, `[[`, numeric(1), "path_length"))*params$link_branch_fraction
      near <- parts[[1]]$near
      node <- sweep(near,2,delta,"+")
      members <- vapply(parts,`[[`,integer(1),"row")
      if (!ribbon) {
        trunk <- bezier_pts(as.numeric(near),as.numeric(node),
          as.numeric(near)+delta/3,as.numeric(near)+delta*2/3,n=24)
        if (side == "subject") trunk <- trunk[nrow(trunk):1,,drop=FALSE]
        emit(parts[[1]],trunk,"trunk",members,first=side=="query",last=side=="subject")
        for (item in parts) {
          from <- as.numeric(node); to <- as.numeric(item$far)
          far_control <- (from+to)/2
          if (!is.null(item$control)) {
            controls <- ggchord_link_controls(item$control,1,far_control)
            far_control <- as.numeric(ggchord_rotate_points(matrix(controls[[if(side=="query") 2 else 1]],1,2),layout$rotation))
          }
          xy <- bezier_pts(from,to,from+delta,far_control,n=60)
          if (side == "subject") xy <- xy[nrow(xy):1,,drop=FALSE]
          emit(item,xy,"branch",item$row,first=side=="subject",last=side=="query")
        }
      } else {
        polygon <- function(a,b,c1,c2) {
          right <- bezier_pts(a[nrow(a),],b[nrow(b),],c1[2,],c2[2,],n=50)
          left <- bezier_pts(a[1,],b[1,],c1[1,],c2[1,],n=50)
          rbind(a,right,b[nrow(b):1,,drop=FALSE],left[nrow(left):1,,drop=FALSE])
        }
        edges <- near[c(1,nrow(near)),,drop=FALSE]
        trunk <- polygon(near,node,sweep(edges,2,delta/3,"+"),sweep(edges,2,delta*2/3,"+"))
        # Reject invalid strip sides rather than hiding a fold with overdraw.
        if (ggchord_front_invalid(node,near)) {
          failures <- failures+1L
          for (item in parts) emit(item,item$xy,"branch",item$row,preserve=TRUE)
          next
        }
        emit(parts[[1]],trunk,"trunk",members)
        tangent <- node[nrow(node),]-node[1,]
        ranks <- order(as.numeric(targets %*% tangent))
        for (k in seq_along(ranks)) {
          item <- parts[[ranks[k]]]
          # Partition the artificial fork cross-section, not the genomic
          # endpoint, so branches do not overpaint the shared junction.
          f <- seq((k-1)/length(parts), k/length(parts), length.out=30)
          node_slice <- cbind(stats::approx(seq(0,1,length.out=nrow(node)),node[,1],f)$y,
                              stats::approx(seq(0,1,length.out=nrow(node)),node[,2],f)$y)
          a <- node_slice[c(1,nrow(node_slice)),,drop=FALSE]
          b <- item$far[c(1,nrow(item$far)),,drop=FALSE]
          xy <- polygon(node_slice,item$far,sweep(a,2,delta,"+"),(a+b)/2)
          emit(item,xy,"branch",item$row)
        }
      }
    }
    if (failures) warning(sprintf("Shared links used separate paths for %d group(s)",failures),call.=FALSE)
    built$data[[i]] <- ggchord_rbind_fill(output)
    geometry <- ggchord_rbind_fill(geometries)
    layout$layer_geometry[[layer$ggchord_layer_id]][[layer$ggchord_type]] <- geometry
    changed <- TRUE
  }
  if (changed) {
    for (pair in list(c("link","link_lines"),c("ribbon","ribbon_polys"))) {
      parts <- lapply(layout$layer_geometry,`[[`,pair[1]); parts <- Filter(function(x) !is.null(x)&&nrow(x)>0,parts)
      if (length(parts)) layout[[pair[2]]] <- ggchord_rbind_fill(parts)
    }
    plot$ggchord$ref$layout <- layout
  }
  built
}
