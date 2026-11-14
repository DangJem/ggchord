pkgload::load_all(quiet = TRUE)
sequence <- data.frame(accver = LETTERS[1:4], length = 1000)
hsp <- data.frame(qaccver = 'A', saccver = c('B','C','D'),
  qstart = 100, qend = 200, sstart = c(400,300,600), send = c(300,400,500),
  kind = c('a','b','c'))
base <- ggchord(sequence, hsp, validate = 'none') + geom_seq()
coords <- list(Chord = coord_chord(), Collinear = coord_collinear(),
  Radial = coord_radial_synteny(), Sankey = coord_sankey())
styles <- function(d) {
  value <- d$ribbon_linetype %||% d$linetype
  unname(vapply(split(as.character(value), d$group), unique, character(1)))
}
result <- list()
for (name in names(coords)) {
  plot <- base + geom_link_ribbon(adjacent_only = FALSE, linewidth = .12) +
    geom_link_ribbon(adjacent_only = FALSE, linetype = 'dotted', linewidth = .2,
      nonadjacent_linewidth = 0) +
    geom_link_ribbon(aes(ribbon_linetype = kind), adjacent_only = FALSE) +
    scale_ribbon_linetype_manual(values = c(a='dotdash',b='longdash',c='dotted')) + coords[[name]]
  built <- ggplot2::ggplot_build(plot)
  result[[name]] <- list(default=styles(built$data[[2]]),fixed=styles(built$data[[3]]),
    mapped=styles(built$data[[4]]),default_widths=sort(unique(built$data[[2]]$linewidth)),
    geometry_sha256=digest::digest(lapply(built$data[2:4],function(d)d[c('x','y','group')]),algo='sha256'))
}
args <- commandArgs(TRUE)
if (!'--before' %in% args) {
  for (name in names(result)) {
    expected <- if(name=='Collinear') c('solid','dashed','dashed') else rep('solid',3)
    stopifnot(identical(result[[name]]$default,expected),
      identical(result[[name]]$fixed,rep('dotted',3)),
      identical(result[[name]]$mapped,c('dotdash','longdash','dotted')))
    if(name!='Collinear')stopifnot(identical(result[[name]]$default_widths,.12))
  }
}
path <- paste0('tools/site-audit/ribbon-styles-',if('--before' %in% args)'before' else 'after','.json')
runtime_files <- list.files('R', '[.]R$', full.names=TRUE)
runtime <- as.list(setNames(vapply(runtime_files,function(f) digest::digest(file=f,algo='sha256',serialize=FALSE),character(1)),runtime_files))
jsonlite::write_json(list(views=result,runtime_sha256=runtime,scope='Four public builds; automatic style scope and explicit fixed/mapped line types; no package test suite'),path,auto_unbox=TRUE,pretty=TRUE)
cat('Recorded',path,'\n')
