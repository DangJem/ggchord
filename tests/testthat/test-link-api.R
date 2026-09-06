local_edition(3)

test_that("accver migration and grid helpers preserve explicit inputs", {
  old_dir <- setwd(tempdir())
  on.exit(setwd(old_dir), add = TRUE)
  seq <- data.frame(accver = c("A", "B"), length = c(1000, 1000))
  old <- seq
  names(old)[1] <- "seq_id"
  expect_warning(p <- ggchord(old), "seq_id is deprecated")
  expect_named(p$ggchord$data$seq_data, c("accver", "length"))
  expect_error(ggchord(transform(seq, seq_id = accver)), "only one")
  expect_warning(geom_seq(aes(seq_id = id)), "mapping is deprecated")
  expect_identical(unit, grid::unit)
  expect_identical(arrow, grid::arrow)
})

test_that("six-column ribbons and unannotated features are usable", {
  old_dir <- setwd(tempdir())
  on.exit(setwd(old_dir), add = TRUE)
  seq <- data.frame(accver = c("A", "B"), length = c(1000, 1000))
  rib <- data.frame(qaccver = "A", saccver = "B", qstart = 100,
                    qend = 200, sstart = 300, send = 400)
  gene <- data.frame(accver = "A", start = 100, end = 200, strand = "+")
  expect_true(validate_ggchord_data(seq, rib, gene)$valid)
  expect_false(any(c("length", "pident") %in% names(clean_ggchord_data(seq, rib, gene)$ribbon_data)))
  p <- ggchord(seq, rib, gene) + geom_seq() + geom_link_ribbon() + geom_gene()
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(ggchord(seq) + geom_seq() + geom_feature(data = gene)), "ggplot_built")
  expect_equal(nrow(deduplicate_ggchord_ribbons(rib, keep = "first")$data), 1L)
  expect_equal(nrow(bundle_ggchord_ribbons(rib, seq, weight = "count")$data), 1L)
  expect_error(bundle_ggchord_ribbons(rib, seq, weight = "length"), "length")
  expect_error(filter_ggchord_ribbons(rib, min_pident = 90), "pident")
})

test_that("point links retain layer data, mapping and arrows", {
  old_dir <- setwd(tempdir())
  on.exit(setwd(old_dir), add = TRUE)
  seq <- data.frame(accver = c("A", "B"), length = c(1000, 1000))
  points <- data.frame(qaccver = "A", saccver = "B", qpos = c(100, 200), spos = c(300, 500), kind = c("a", "b"))
  p <- ggchord(seq) + geom_seq() + geom_link_line(aes(link_colour = kind), data = points,
    arrow = arrow(length = unit(2, "mm"))) +
    geom_link_line(data = points[1, ], link_type = "straight", colour = "red")
  built <- ggplot2::ggplot_build(p)
  expect_length(unique(built$data[[2]]$group), 2L)
  expect_equal(nrow(built$data[[3]]), 2L)
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
  exported <- export_ggchord_layout(p, include = "link", original_data = TRUE)
  expect_setequal(unique(exported$link$source_row), c(1L, 2L))
  expect_error(geom_link_line(link_type = "straight", link_branch = "query"), "Straight")
  expect_error(ggplot2::ggplot_build(ggchord(seq) + geom_link_line(data = transform(points, qpos = 2000))), "positions")
})

test_that("shared links preserve membership and mapped style boundaries", {
  old_dir <- setwd(tempdir())
  on.exit(setwd(old_dir), add = TRUE)
  seq <- data.frame(accver = c("A", "B", "C"), length = 1000)
  d <- data.frame(qaccver = "A", saccver = c("B", "C"), qpos = 200, spos = 500,
                  kind = c("one", "two"))
  base <- ggchord(seq) + geom_seq()
  p <- base + geom_link_line(data = d, link_branch = "query")
  out <- export_ggchord_layout(p, include = "link", original_data = TRUE)$link
  trunk <- out[out$.component == "trunk", ]
  expect_gt(nrow(trunk), 0)
  expect_equal(trunk$source_rows[[1]], c(1L, 2L))
  expect_true(all(is.na(trunk$source_row)))
  expect_setequal(out$source_row[!is.na(out$source_row)], c(1L, 2L))
  different <- ggplot2::ggplot_build(base + geom_link_line(aes(link_colour = kind), data = d, link_branch = "query"))
  expect_false(any(different$data[[2]]$.component == "trunk"))
  same <- ggplot2::ggplot_build(base + geom_link_line(aes(link_colour = kind), data = d, link_branch = "query") +
    scale_link_colour_manual(values = c(one = "red", two = "red")))
  expect_true(any(same$data[[2]]$.component == "trunk"))
  r <- transform(d, qstart = 150, qend = 250, sstart = 400, send = 600)
  p <- ggchord(seq, r) + geom_seq() + geom_link_ribbon(fill = "orange", link_branch = "query")
  out <- export_ggchord_layout(p, include = "ribbon")$ribbon
  expect_true(any(out$.component == "trunk"))
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
})

test_that("local entity clearance keeps unobstructed fronts close", {
  old_dir <- setwd(tempdir()); on.exit(setwd(old_dir), add=TRUE)
  seq <- data.frame(accver=c("A","B"),length=1000)
  rib <- data.frame(qaccver="A",saccver="B",qstart=100,qend=800,sstart=100,send=800)
  feature <- data.frame(accver="A",start=350,end=450,strand="-",type="block")
  base <- ggchord(seq,rib)+geom_seq()+geom_feature(data=feature,feature_shape="block")
  default <- get_chord_layout(base+geom_link_ribbon())
  smooth <- get_chord_layout(base+geom_link_ribbon(link_avoid="smooth"))
  uniform <- get_chord_layout(base+geom_link_ribbon(link_avoid="uniform"))
  expect_equal(nrow(default$obstacles), 0L)
  expect_equal(unique(default$ribbon_polys$q_gap), .035)
  expect_true(all(c("accver","side","normal_min","normal_max","source_layer") %in% names(smooth$obstacles)))
  expect_true(all(smooth$obstacles$normal_max > 0))
  front <- function(layout) {
    d <- layout$ribbon_polys; n <- d$.q_n[1]
    xy <- ggchord_rotate_points(as.matrix(d[seq_len(n),c("x","y")]),-layout$rotation)
    sqrt(rowSums(xy^2))
  }
  expect_gt(diff(range(front(smooth))), .03)
  expect_lt(diff(range(front(uniform))), .005)
  explicit <- get_chord_layout(base+geom_link_ribbon(ribbon_gap=.08))
  expect_equal(unique(explicit$ribbon_polys$q_gap), .08)
  expect_error(ggplot2::ggplot_build(ggchord(seq,rib)+geom_link_ribbon(aes(ribbon_fill=pident))), "pident")
})

test_that("feature legends draw the mapped shapes and a taller gene arrow", {
  old_dir <- setwd(tempdir()); on.exit(setwd(old_dir), add=TRUE)
  s <- data.frame(accver=c("A","B"),length=1000)
  d <- data.frame(accver="A",start=c(100,300,500,700),end=c(200,400,600,800),
    strand="+",type=c("CDS","tRNA","repeat","promoter"))
  p <- ggchord(s)+geom_seq(show.legend=FALSE)+geom_feature(aes(feature_shape=type),data=d)+
    scale_feature_shape_manual(values=c(CDS="arrow",tRNA="block","repeat"="chevron",promoter="lollipop"))
  keys <- list()
  visit <- function(g) {
    if (!inherits(g,"gtable")) return(invisible(NULL))
    idx <- which(grepl("^key-", g$layout$name))
    if (length(idx)) keys <<- c(keys,lapply(g$grobs[idx],function(x) tail(as.list(x$children),1)[[1]]))
    invisible(lapply(g$grobs,visit))
  }
  visit(ggplot2::ggplotGrob(p))
  expect_equal(vapply(keys,function(x) class(x)[1],character(1)),c("polygon","rect","polygon","gTree"))
  expect_gt(diff(range(as.numeric(keys[[1]]$y))),.5)
  expect_true(any(vapply(keys[[4]]$children,inherits,logical(1),"circle")))
})


test_that("unshared links preserve vertex-wise staged styles", {
  old_dir <- setwd(tempdir()); on.exit(setwd(old_dir), add = TRUE)
  s <- data.frame(accver = c("A", "B", "C"), length = 1000)
  d <- data.frame(qaccver = "A", saccver = c("B", "C"), qpos = 200, spos = 500)
  build <- function(branch) ggplot2::ggplot_build(
    ggchord(s) + geom_seq() + geom_link_line(
      aes(link_alpha = after_scale(seq_along(x) / length(x))),
      data = d, link_branch = branch))$data[[2]]
  separate <- build("none")
  shared <- build("query")
  expect_false(any(shared$.component == "trunk"))
  expect_equal(shared$link_alpha, separate$link_alpha)
  expect_equal(shared[c("x", "y")], separate[c("x", "y")])
})
