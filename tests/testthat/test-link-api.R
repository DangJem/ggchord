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
})


test_that("shared ribbons retain explicit edge control points", {
  old_dir <- setwd(tempdir()); on.exit(setwd(old_dir), add = TRUE)
  s <- data.frame(accver = c("A", "B", "C"), length = 1000)
  r <- data.frame(qaccver = "A", saccver = c("B", "C"),
    qstart = 100, qend = 300, sstart = 400, send = 600)
  p <- ggchord(s, r) + geom_seq() + geom_link_ribbon(
    fill = "orange", link_branch = "query", ribbon_ctrl_point = c(4, 4))
  d <- export_ggchord_layout(p, include = "ribbon")$ribbon
  expect_true(any(d$.component == "trunk"))
  # A deliberately external control must pull the branches beyond the tracks.
  branches <- d[d$.component == "branch", ]
  expect_gt(max(sqrt(branches$x^2 + branches$y^2)), 2)
})
