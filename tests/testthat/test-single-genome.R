test_that("v0.13 removes legacy sequence and ribbon entries", {
  old <- data.frame(seq_id="A",length=100)
  expect_error(ggchord(old), "removed in v0.13.0")
  expect_error(geom_seq(aes(seq_id=id)), "removed in v0.13.0")
  expect_false("geom_ribbon" %in% getNamespaceExports("ggchord"))
})

test_that("coord_genome owns a one-sequence circular contract", {
  seq <- data.frame(accver="circle",length=1000)
  closed <- ggchord(seq,validate="none") + geom_seq() + coord_genome()
  arc <- get_chord_layout(closed)$seq_arcs[[1L]]
  expect_s3_class(coord_genome(),"CoordGenome")
  expect_equal(unname(unlist(arc[1,c("x","y")])),unname(unlist(arc[nrow(arc),c("x","y")])),tolerance=1e-8)

  open <- ggchord(seq,validate="none") + geom_seq() + coord_genome(gap=20)
  open_arc <- get_chord_layout(open)$seq_arcs[[1L]]
  expect_gt(sum((open_arc[1,c("x","y")]-open_arc[nrow(open_arc),c("x","y")])^2),0.01)
  expect_equal(unname(unlist(open_arc[1,c("x","y")])),c(0,1),tolerance=.02)
  ccw <- ggchord(seq,validate="none") + geom_seq() +
    coord_genome(gap=20,direction="counterclockwise")
  ccw_arc <- get_chord_layout(ccw)$seq_arcs[[1L]]
  expect_equal(unname(unlist(ccw_arc[1,c("x","y")])),c(0,1),tolerance=.02)
  expect_error(get_chord_layout(ggchord(rbind(seq,transform(seq,accver="other")),validate="none") + geom_seq() + coord_genome()),"exactly one")
  expect_error(get_chord_layout(ggchord(seq,validate="none") + geom_seq(seq_gap=.1) + coord_genome()),"owns the opening")
})

test_that("restriction search reports reproducible biological fields", {
  sites <- find_restriction_sites(c(g="AAAAGAATTCTTTGGATCC"),c("EcoRI","BamHI"),circular=FALSE)
  expect_equal(sites$enzyme,c("EcoRI","BamHI"))
  expect_equal(sites$position,c(6,15))
  expect_true(all(c("recognition_sequence","cut_top","cut_bottom","end_type","source") %in% names(sites)))
  expect_identical(attr(sites,"enzyme_database_version"),"ggchord-builtin-2026.09")
  expect_equal(nrow(find_restriction_sites("AATTCG", "EcoRI", circular=TRUE)),1L)
  expect_equal(nrow(find_restriction_sites("AATTCG", "EcoRI", circular=FALSE)),0L)
  expect_equal(nrow(find_restriction_sites("GAATTCGAATTC", "EcoRI", max_cuts=1)),0L)
  empty <- find_restriction_sites("AAAAAAAA", "EcoRI")
  p <- ggchord(data.frame(accver="sequence",length=8),validate="none") +
    geom_seq() + coord_genome() + geom_restriction_site(data=empty)
  expect_s3_class(ggplot2::ggplot_build(p),"ggplot_built")
})

test_that("restriction geometry retains every site and builds", {
  seq <- data.frame(accver="g",length=1000)
  sites <- data.frame(accver="g",position=c(100,105,700),enzyme=c("A","B","C"))
  p <- ggchord(seq,validate="none") + geom_seq() + coord_genome(gap=8) +
    geom_restriction_site(data=sites,branch_threshold=.02)
  expect_s3_class(ggplot2::ggplot_build(p),"ggplot_built")
  expect_s3_class(ggplot2::ggplotGrob(p),"gtable")
  layout <- get_chord_layout(p)
  geometry <- layout$restriction_sites
  labels <- geometry[geometry$.component=="label",,drop=FALSE]
  expect_equal(sort(unique(labels$source_row)),seq_len(nrow(sites)))
  expect_equal(nrow(labels),nrow(sites))
  expect_equal(nrow(geometry),17L)
  exported <- export_ggchord_layout(p,include="restriction",original_data=TRUE)
  expect_equal(sort(unique(stats::na.omit(exported$restriction$source_row))),seq_len(nrow(sites)))
  expect_identical(exported$metadata$coordinate,"genome")
})
