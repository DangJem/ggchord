test_that("named tracks share domains and compose with selected axes and a bp reference", {
  data <- data.frame(accver = c("A", "B"), start = 100, end = 200, value = c(-1, 2))
  plot <- ggchord(smoke_sequences()) + geom_seq() +
    geom_track_col(aes(track_value = value), data = data,
      position = position_track("signal")) +
    scale_track_value_continuous("signal", limits = c(-2, 3),
      guide = ggplot2::guide_axis(theme = theme(axis.text = element_text(colour = "purple")))) +
    scale_seq_position_continuous(sequences = "B") + geom_seq_scale(length = 100) +
    coord_collinear()
  out <- export_ggchord_layout(plot, include = c("track", "axis"))
  expect_equal(out$track_registry$track_id, c("signal", "signal"))
  expect_true(all(vapply(out$track_registry$limits, function(x) all(x == c(-2, 3)), logical(1))))
  expect_true(any(out$axis$component == "seq_scale"))
  expect_s3_class(build_ggchord_smoke(plot), "ggplot_built")
})

test_that("track bars respect explicit source groups before projection", {
  events <- data.frame(accver=c("A","A","B"),pos=c(100,100,100),weight=c(2,3,7),
    category=c("a","b","a"))
  p <- ggchord(smoke_sequences())+geom_seq()+
    geom_track(aes(x=pos,weight=weight,fill=category,group=1),data=events,geom="bar",track=position_track("counts"))+
    coord_collinear()
  out <- suppressWarnings(export_ggchord_layout(p,include="track"))
  expect_equal(out$track_statistics[[1]]$count,c(5,7))
  expect_equal(out$track_statistics[[1]]$source_rows[[1]],1:2)
})

test_that("track bins retain the observations counted at outer boundaries", {
  events <- data.frame(accver="A",pos=c(100,150,200))
  p <- ggchord(smoke_sequences())+geom_seq()+
    geom_track(aes(x=pos),data=events,geom="col",stat="bin",breaks=c(100,200),
      track=position_track("bins"))+coord_collinear()
  out <- export_ggchord_layout(p,include="track")$track_statistics[[1]]
  expect_equal(out$count,3)
  expect_equal(out$source_rows[[1]],1:3)
})
