# Focused engine equivalence and public ownership checks, never a release run.
# Optional: capture <rds> before extraction, then compare <rds> after each batch.
# The RDS is temporary evidence, not a checked-in numerical snapshot contract.
pkgload::load_all(quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
local({
  path <- tempfile(fileext = '.pdf')
  grDevices::pdf(path, width = 8, height = 6, family = 'Helvetica')
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  set.seed(163)
  s <- data.frame(accver = c('A', 'B'), length = c(1000, 1200))
  g <- data.frame(accver = c('A', 'B'), start = c(100, 500), end = c(300, 800),
    strand = c('+', '-'), anno = c('one', 'two'))
  f <- data.frame(id = 'A', lo = c(400, 710), hi = c(560, 830),
    chain = c('+', '.'), kind = c('CDS', 'misc'), name = c('first', 'unknown'))
  r <- data.frame(qaccver = 'A', saccver = 'B', qstart = c(100, 650),
    qend = c(220, 820), sstart = c(500, 1000), send = c(700, 850), pident = c(95, 80))
  entities <- function(seqs = s) ggchord(seqs, ribbon_data = if(nrow(seqs)>1) r else NULL) +
    geom_seq() + geom_gene(data = g[g$accver %in% seqs$accver, ]) +
    geom_feature(aes(accver = id, start = lo, end = hi, strand = chain,
      type = kind, label = name), data = f, shape = 'block', width = .04,
      position = position_strand(.13))
  sites <- data.frame(accver = 'A', position = c(80, 350, 650), enzyme = c('e1','e2','e3'))
  windows <- data.frame(accver = c('A','A','B'), start = c(1,601,201),
    end = c(400,1000,1000), locus_id = c('A1','A2','B1'))
  view <- segment_ggchord_data(s, windows = windows)
  points <- data.frame(accver = 'A', pos = c(100,200,300,700,800,900),
    value = c(1,2,NA,3,4,2))
  focused <- focus_ggchord_track(points, s, windows)
  plots <- list(
    chord = entities() + geom_link_ribbon() + geom_gene_label(data = g),
    curved = ggchord(s, r) + geom_seq(radius = c(A=1.1,B=.9),
      curvature = c(A=-.6,B=.8), orientation = c(A=-1,B=1)) +
      geom_gene(data = g) + geom_gene_label_repel(data = g, label_layout='radial') +
      geom_link_ribbon(avoid='uniform') + coord_chord(rotation=135),
    auto = entities() + geom_gene_label_repel(data=g, label_layout='auto'),
    circular = entities(s[1, ]) + geom_feature_label_repel(aes(accver=id,
      start=lo,end=hi,strand=chain,type=kind,label=name),data=f) +
      geom_restriction_site(data=sites) + coord_circular(direction='counterclockwise', gap=20),
    collinear = entities() + geom_link_ribbon() + geom_gene_label_repel(data=g) +
      coord_collinear(),
    segmented = ggchord(view$seq_data) + geom_seq() +
      geom_track_line(aes(track_value=value),data=focused,position=position_track('signal')) +
      coord_collinear(segments=view$segments))
  include <- c('seq','ribbon','gene','feature','labels','link','restriction','axis','track')
  snapshots <- lapply(plots, function(p) {
    exported <- export_ggchord_layout(p, include=include, original_data=TRUE)
    b <- ggplot2::ggplot_build(p)
    stopifnot(inherits(b, 'ggplot_built'))
    exported
  })
  for (name in names(snapshots)) {
    x <- snapshots[[name]]
    if (!nrow(x$seq)) stop('Empty sequence export: ', name)
    stopifnot(nrow(x$seq)>0, identical(x$metadata$units, 'data'),
      all(c('layer_id','source_row','component') %in% names(x$seq)))
    if (name %in% c('chord','auto','circular','collinear')) {
      stopifnot(nrow(x$gene)>0, nrow(x$feature)>0,
        all(c('one','two')[seq_len(if(name=='circular') 1L else 2L)] %in% unique(x$gene$anno)),
        '.' %in% x$original_data[[unique(x$feature$layer_id)[1]]]$gene_poly$chain)
    }
  }
  if (length(args) && args[1]=='capture') {
    saveRDS(snapshots, args[2]); cat('Captured six public engine baselines.\n')
  } else if (length(args) && args[1]=='compare') {
    original <- readRDS(args[2])
    for (name in names(snapshots)) {
      delta <- all.equal(original[[name]], snapshots[[name]], check.attributes=TRUE)
      if (!isTRUE(delta)) stop(name, ': ', paste(delta, collapse='; '))
    }
    cat('Six live before/after public exports agree under the same dimensions/seed/font.\n')
  }
  out <- data.frame(case=names(snapshots), result='PASS',
    detail='Representative public build/export, data retention and source ownership')
  write.csv(out, 'design-roadmap/audit-generated/engine-baseline-v0.16.3.csv',row.names=FALSE)
  print(out,row.names=FALSE)
})
