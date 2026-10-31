# Direct checks of the new public workflows, not a testthat/release run.
pkgload::load_all(quiet = TRUE)
local({
  old <- getwd(); setwd(tempdir()); on.exit(setwd(old), add = TRUE)
  sequences <- data.frame(accver = c("A", "B"), length = c(1000, 1000))
  values <- data.frame(accver = c("A", "B"), start = 100, end = 300, value = c(-2, 8))
  points <- data.frame(accver = "A", pos = seq(100, 600, 100), value = c(1, 2, NA, 4, 3, 2))
  plot <- ggchord(sequences) + geom_seq() +
    geom_track_bar(aes(track_value = value), data = values,
      position = position_track("coverage")) +
    geom_track_line(aes(track_value = value), data = points,
      position = position_track("signal", order = 2)) +
    scale_track_value_continuous("coverage", limits = c(-5, 10), unit = "reads") +
    coord_collinear()
  exported <- export_ggchord_layout(plot, include = c("track", "axis"))
  stopifnot(nrow(exported$track_registry) == 4L,
    all(vapply(exported$track_registry$limits[exported$track_registry$track_id == "coverage"],
      function(x) all(x == c(-5, 10)), logical(1))),
    any(exported$track_diagnostics$reason == "missing_value"),
    any(exported$axis$component == "track_axis"))
  stopifnot(inherits(ggplot2::ggplot_build(plot), "ggplot_built"))
  # Line chunks retain the NA gap; no single group may span both sides.
  line <- exported$track[exported$track$track_id == "signal" &
    exported$track$.component == "track_line", ]
  stopifnot(all(vapply(split(line$source_row, line$group), function(x)
    !(any(x < 3) && any(x > 3)), logical(1))))
  for (coordinate in list(coord_circular(), coord_chord())) {
    ring <- ggchord(sequences[1L, ]) + geom_seq() +
      geom_track_tile(aes(track_fill = category),
        data = transform(values[1L, ], category = "site"),
        position = position_track("kind")) +
      scale_track_fill_manual("kind", values = c(site = "#334455")) +
      geom_track_bar(aes(track_value = value), data = values[1L, ],
        position = position_track("coverage", order = 2)) +
      geom_track_line(aes(track_value = value), data = points,
        position = position_track("signal", order = 3)) + coordinate
    built <- ggplot2::ggplot_build(ring)
    key <- ggchord:::ggchord_track_aes("fill", "kind")
    stopifnot(all(built$data[[2L]][[key]] == "#334455"))
  }
  windows <- data.frame(accver = c("A", "A", "B", "B"),
    start = c(100, 400, 100, 400), end = c(200, 600, 200, 600),
    locus_id = c("A1", "A2", "B1", "B2"))
  view <- segment_ggchord_data(sequences, windows = windows)
  focused <- focus_ggchord_track(points, sequences, windows)
  segmented <- ggchord(view$seq_data) + geom_seq() +
    geom_track_line(aes(track_value = value), data = focused,
      position = position_track("signal")) +
    geom_seq_scale(length = 50) +
    scale_seq_position_continuous(sequences = "A2", coordinates = "source") +
    coord_collinear(segments = view$segments)
  out <- export_ggchord_layout(segmented)
  stopifnot(nrow(out$track) > 0L, all(out$axis_selection$selected == (out$axis_selection$accver == "A2")))
  stopifnot(inherits(ggplot2::ggplot_build(segmented), "ggplot_built"))
  report <- diagnose_ggchord(plot, width = 8, height = 5)
  if (isTRUE(report$measurement$converged)) {
    physical_width <- with(report$track_registry, (normal_max - normal_min) /
      report$measurement$units_per_mm[["y"]])
    stopifnot(all(abs(physical_width - 4) < .15))
  }
  invalid <- tryCatch(ggplot2::ggplot_build(plot +
    scale_track_value_continuous("coverage", limits = c(1, 10))), error = identity)
  stopifnot(inherits(invalid, "error"), grepl("baseline", conditionMessage(invalid)))
  empty <- points[0, ]
  invisible(ggplot2::ggplot_build(ggchord(sequences) + geom_seq() +
    geom_track_line(aes(track_value = value), data = empty,
      position = position_track("empty")) + coord_collinear()))
})
cat("Scoped quantitative tracks, coordinate mappings, gaps and physical bandwidth passed.\n")

# v0.16.3: shared declarations/domain and selected boundary semantics.
local({
  results <- list()
  record <- function(case, detail) results[[length(results)+1L]] <<-
    data.frame(case=case,result='PASS',detail=detail)
  fails <- function(expr, pattern) {
    e <- tryCatch(force(expr), error=identity)
    stopifnot(inherits(e,'error'),grepl(pattern,conditionMessage(e),ignore.case=TRUE))
  }
  s <- data.frame(accver=c('A','B'),length=1000)
  a <- data.frame(accver='A',start=100,end=200,value=10)
  b <- data.frame(accver='B',pos=c(100,200),value=c(50,100))
  position <- position_track('shared',width=unit(4,'mm'))
  scale <- scale_track_value_continuous('shared',unit='reads')
  p <- ggchord(s)+geom_seq()+geom_track_bar(aes(track_value=value),data=a,
    position=position)+geom_track_line(aes(track_value=value),data=b,
    position=position)+scale+coord_collinear()
  original_domain <- scale$range$range
  out <- export_ggchord_layout(p,include='track',original_data=TRUE)
  stopifnot(nrow(out$track_registry)==2L,
    all(vapply(out$track_registry$limits,function(x) isTRUE(all.equal(x,c(0,100))),logical(1))),
    all(lengths(out$track_registry$source_layers)==2L),
    identical(scale$range$range,original_domain),all(out$track_diagnostics$retained))
  record('same_id_domain','Bar baseline and two independent layers/sequences train one cloned numerical domain/band')
  conflict <- p+geom_track_line(aes(track_value=value),data=b,
    position=position_track('shared',width=unit(5,'mm')))
  fails(export_ggchord_layout(conflict),'Conflicting placement')
  record('same_id_conflict','Inconsistent physical declarations are rejected before drawing')
  values <- data.frame(accver='A',start=c(100,200,300),end=c(150,250,350),value=c(.1,10,1000))
  scaled <- ggchord(s[1,])+geom_seq()+geom_track_bar(aes(track_value=value),
    data=values,baseline=1,position=position_track('log'))+
    scale_track_value_continuous('log',limits=c(1,100),transform='log10')+coord_collinear()
  d <- export_ggchord_layout(scaled,include='track')
  stopifnot(identical(d$track_diagnostics$retained,c(FALSE,TRUE,FALSE)),
    all(d$track_diagnostics$reason[c(1,3)]=='outside_scale_limits'))
  fails(export_ggchord_layout(scaled+scale_track_value_continuous('log',
    limits=c(1,100),transform='log10',oob=scales::oob_keep)),'physical band')
  record('transform_oob_baseline','Transformed limits censor outliers; an OOB policy cannot escape its band; valid baseline remains explicit')
  bins <- data.frame(accver='A',start=c(-5,150,2000),end=c(30,200,2100),value=c(1,NA,3))
  clipped <- ggchord(s[1,])+geom_seq()+geom_track_bar(aes(track_value=value),
    data=bins,bounds='clip',position=position_track('bounds'))+coord_collinear()
  d <- export_ggchord_layout(clipped,include='track',original_data=TRUE)
  stopifnot(identical(d$track_diagnostics$retained,c(TRUE,FALSE,FALSE)),
    d$track_diagnostics$coordinate_clipped[1],
    d$track_diagnostics$reason[2]=='missing_value',
    identical(d$original_data[[unique(d$track$layer_id)[1]]]$track$start,bins$start))
  fails(export_ggchord_layout(ggchord(s[1,])+geom_seq()+geom_track_bar(aes(track_value=value),
    data=bins,position=position_track('bounds'))+coord_collinear()),'outside sequence')
  record('bounds_missing_source','Clipping/omission/NA are reported separately and source coordinates remain unchanged')
  points <- data.frame(accver='A',pos=c(100,300,700,900),value=c(1,2,4,3))
  ring <- function(data,closed=TRUE,gap=0) ggchord(s[1,])+geom_seq()+
    geom_track_line(aes(track_value=value),data=data,closed=closed,
      position=position_track('ring'))+coord_circular(gap=gap)
  complete <- export_ggchord_layout(ring(points),include='track')$track
  complete <- complete[complete$.component=='track_line',]
  stopifnot(isTRUE(all.equal(unname(unlist(complete[1,c('x','y')])),
    unname(unlist(complete[nrow(complete),c('x','y')])))))
  missing <- points;missing$value[2]<-NA
  open <- export_ggchord_layout(ring(missing),include='track')
  line <- open$track[open$track$.component=='track_line',]
  stopifnot(!2L %in% line$source_row,
    !isTRUE(all.equal(unname(unlist(line[1,c('x','y')])),
      unname(unlist(line[nrow(line),c('x','y')])))))
  fails(export_ggchord_layout(ring(points,gap=10)),'closed=TRUE')
  record('closed_na_seam','Complete circular series explicitly close; NA prevents seam closure; gapped circles reject closed')
  result <- do.call(rbind,results)
  write.csv(result,'design-roadmap/audit-generated/engine-track-validation-v0.16.3.csv',row.names=FALSE)
  print(result,row.names=FALSE)
})
