# Directed geometry controls, not an exhaustive device or release matrix.
source('tools/audit-source-snapshot.R')
audit <- audit_source_snapshot(additional = 'tools/long-term-audit/check-circular-feature-topology.R')
pkgload::load_all(quiet = TRUE)
seqs <- data.frame(accver='circle',length=1000)
results <- list(); panels <- list()
# One wrapped input, four biological directions, closed/open and rotated/reversed.
for (shape in c('arrow','promoter_arrow','block','capped_line','marker','lollipop')) {
 for (strand in c('+','-','+/-','.')) for (gap in c(0,24)) {
  d <- data.frame(accver='circle',start=900,end=120,strand=strand,anno='wrapped')
  p <- ggchord(seqs) + geom_seq() + geom_feature(data=d,shape=shape,
    position=position_plasmid(),width=.07,fill='white',show.legend=FALSE) +
    coord_circular(gap=gap) + theme_ggchord_plasmid()
  l <- get_chord_layout(p)
  poly <- l$gene_polys[l$gene_polys$.component=='polygon',]
  groups <- unique(poly$group)
  # Markers/lollipops are intentionally multipart glyphs, not interval bodies.
  body <- shape %in% c('arrow','promoter_arrow','block')
  if (body) stopifnot(length(groups)==if(gap==0) 1L else 2L)
  ref <- l$sequence_reference$refs[[1]]
  stopifnot(identical(ggchord_closed_reference(ref),gap==0))
  # Inspect head tips in the shape factory, before curve sampling can blur them.
  tips <- numeric()
  if (shape %in% c('arrow','promoter_arrow')) {
    pieces <- ggchord_feature_intervals(900,120,1000,if(strand=='-') '-' else '+',
      circular=TRUE,continuous=gap==0)
    for (i in seq_along(pieces)) {
      interval <- ggchord_feature_angle_interval(pieces[[i]],1000,-1,-2*pi,0)
      parts <- ggchord_feature_geometry(shape,interval[1],interval[2],.9,.07,1,ref,
        draw_head=pieces[[i]]$draw_head && strand!='.',
        draw_start_head=i==1L && strand=='+/-',bidirectional=strand=='+/-')
      for (part in parts) if (!is.null(part$angle))
        tips <- c(tips,part$angle[abs(part$radius-.9)<1e-10])
    }
    expected <- if(strand=='+/-') 2L else if(strand=='.') 0L else 1L
    stopifnot(length(tips)==expected)
  }
  key <- paste(shape,strand,gap,sep='/')
  results[[key]] <- list(shape=shape,strand=strand,gap=gap,
    polygons=length(groups),head_tips=length(tips),finite=all(is.finite(c(poly$x,poly$y))))
  stopifnot(results[[key]]$finite)
  if(shape=='promoter_arrow' && strand %in% c('+/-','-')) panels[[key]] <- p
 }
}
# Short features at the seam must wrap their minimum-size glyph, not collapse.
short <- data.frame(accver='circle',start=999,end=2,strand='+',anno='short')
p <- ggchord(seqs) + geom_seq() + geom_feature(data=short,shape='promoter_arrow',
  position=position_plasmid(),width=.07,show.legend=FALSE) + coord_circular()
l <- get_chord_layout(p)
stopifnot(length(unique(l$gene_polys$group))==1L,
  all(is.finite(c(l$gene_polys$x,l$gene_polys$y))))
results$short <- list(polygons=length(unique(l$gene_polys$group)))
# Preserve truly disconnected source runs; they must never be seam-joined.
d <- data.frame(accver='circle',start=100,end=500,strand='+',anno='disconnected')
d$segments <- I(list(data.frame(segment_index=1:3,segment_type=c('standard','gap','standard'),
 start=c(100,201,400),end=c(200,399,500),color='#ffffff')))
p <- ggchord(seqs) + geom_seq() + geom_feature(data=d,shape='arrow',position=position_plasmid()) + coord_circular()
l <- get_chord_layout(p)
stopifnot(length(unique(l$gene_polys$group))==2L,
  length(unique(l$gene_polys$source_row))==1L)
results$disconnected <- list(polygons=2L,biological_sources=1L)
# A true wrapped gap must remain empty, while adjacent equal-colour segments
# may join through the origin using the known sequence length.
wrapped <- data.frame(accver='circle',start=900,end=120,strand='+/-',anno='wrapped segments')
for (kind in c('gap','adjacent','colour')) {
  segments <- if(kind=='gap') data.frame(segment_index=1:3,
    segment_type=c('standard','gap','standard'),start=c(900,951,20),end=c(950,19,120),color='#ffffff') else
    data.frame(segment_index=1:2,segment_type='standard',start=c(900,1),end=c(1000,120),
      color=if(kind=='colour') c('#ffffff','#31979c') else '#ffffff')
  wrapped$segments <- I(list(segments))
  p <- ggchord(seqs) + geom_seq() + geom_feature(data=wrapped,shape='promoter_arrow',
    position=position_plasmid()) + coord_circular()
  l <- get_chord_layout(p)
  poly <- l$gene_polys[l$gene_polys$.component=='polygon',]
  covered <- any(vapply(split(poly,poly$group),function(p)
    ggchord_point_in_polygon(0,.9,p),logical(1)))
  expected <- if(kind=='adjacent') 1L else 2L
  stopifnot(length(unique(poly$group))==expected,
    length(unique(poly$source_row))==1L)
  if(kind!='colour') stopifnot(covered==(kind!='gap'))
  # At a colour join the origin lies on polygon boundaries, where strict
  # interior tests need not return TRUE. Probe both sides away from that edge.
  if(kind=='colour') for(angle in c(-.02,.02)) stopifnot(any(vapply(
    split(poly,poly$group),function(p)ggchord_point_in_polygon(.9*sin(angle),.9*cos(angle),p),logical(1))))
  if(kind=='colour') stopifnot(length(unique(poly$feature_fill_explicit))==2L)
  results[[paste0('wrapped_',kind)]] <- list(polygons=expected,origin_covered=covered,
    colours=unique(poly$feature_fill_explicit),biological_sources=1L)
}
# Real catalogue inventory: report every wrap and its original segment count.
ids <- c('pUC19','pTRE_Tight_BI','pCAMBIA1300','pSB1C3','pBR322',
 'pBluescript_II_SK_plus','pET_28a_plus','pETDuet_1','pDONR221',
 'pcDNA3_1_plus','pTRIPZ','pEarleyGate_201','pSpCas9_BB_2A_GFP_PX458')
inventory <- list()
for(id in ids) {
 env <- new.env(); data(list=paste0('plasmid_example_',id),package='ggchord',envir=env)
 d <- find_common_features(env[[paste0('plasmid_example_',id)]])
 at <- which(d$start>d$end)
 if(length(at)) inventory[[id]] <- data.frame(label=d$anno[at],start=d$start[at],end=d$end[at],
  strand=d$strand[at],segments=vapply(d$segments[at],nrow,integer(1)))
}
path <- 'tools/long-term-audit/v0.27.0-circular-topology-controls.png'
ggchord_diagnostic_device(path,'png',10,8,110)
grid::grid.newpage(); grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,2)))
for(i in seq_along(panels)) {
 grid::pushViewport(grid::viewport(layout.pos.row=(i-1L)%/%2L+1L,layout.pos.col=(i-1L)%%2L+1L))
 grid::grid.draw(ggplot2::ggplotGrob(panels[[i]])); grid::grid.text(names(panels)[i],y=.03)
 grid::popViewport()
}
grid::popViewport(); grDevices::dev.off()
audit_assert_source_snapshot(audit)
jsonlite::write_json(list(cases=results,real_wrap_inventory=inventory,
 runtime_sha256=as.list(audit$engine),input_sha256=as.list(audit$inputs),
 not_checked=c('all devices/fonts','arbitrary custom grobs','arbitrary overlapping/invalid segment tables')),
 'tools/long-term-audit/v0.27.0-circular-topology-results.json',pretty=TRUE,auto_unbox=TRUE)
cat('Circular topology controls completed\n')
