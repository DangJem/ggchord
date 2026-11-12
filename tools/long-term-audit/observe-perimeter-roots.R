# Observe how a real input chooses directed roots before shared placement.
# This collector returns the original task table unchanged.
source('tools/audit-source-snapshot.R')
snapshot <- audit_source_snapshot(c('tools/long-term-audit/observe-perimeter-roots.R',
  'tools/documentation-examples.R'))
pkgload::load_all(quiet=TRUE)
namespace <- asNamespace('ggchord')
original_tasks <- get('ggchord_perimeter_tasks',namespace)
observations <- list()
collector <- function(registry, units_per_inch) {
  tasks <- original_tasks(registry,units_per_inch)
  if (!nrow(tasks)) return(tasks)
  for (i in which(tasks$kind != 'restriction')) {
    d <- registry[[tasks$id[i]]]$gene_label_repel
    d <- d[d$.component %in% 'segment' & d$source_row == tasks$source[i],,drop=FALSE]
    starts <- vapply(seq_len(nrow(d)),function(j) !any(
      (d$xend-d$x[j])^2+(d$yend-d$y[j])^2 < 1e-18),logical(1))
    unique_root <- sum(starts)==1L
    same_vertex <- if(unique_root) isTRUE(all.equal(
      unname(c(tasks$root_x[i],tasks$root_y[i])),unname(c(d$x[starts],d$y[starts])))) else NA
    observations[[length(observations)+1L]] <<- list(layer=tasks$id[i],
      source_row=tasks$source[i],kind=tasks$kind[i],directed_root_unique=unique_root,
      task_uses_directed_root=same_vertex)
  }
  tasks
}
unlockBinding('ggchord_perimeter_tasks',namespace)
assign('ggchord_perimeter_tasks',collector,namespace)
lockBinding('ggchord_perimeter_tasks',namespace)
lines <- readLines('tools/documentation-examples.R',warn=FALSE)
markers <- grep('^# ---- [a-z]+ ----$',lines)
chunks <- setNames(lapply(seq_len(length(markers)-1L),function(i)
  lines[seq.int(markers[i]+1L,markers[i+1L]-1L)]),
  sub('^# ---- ([a-z]+) ----$','\\1',lines[markers[-length(markers)]]))
for(name in c('docstyle','showcaseplasmiddata','showcaseplasmid')) eval(parse(text=chunks[[name]]))
path <- tempfile(fileext='.png')
ggchord_diagnostic_device(path,'png',8,8,110)
built <- ggplot2::ggplot_build(p_show_plasmid)
grDevices::dev.off();unlink(path)
audit_assert_source_snapshot(snapshot)
jsonlite::write_json(list(scope='One real canonical flagship build before gtable refinement; collector returns original tasks, not a full quality or arbitrary-input certification',
  geometry_override=FALSE, observations=observations,
  runtime_sha256=as.list(snapshot$engine)),
  'tools/long-term-audit/v0.27.0-perimeter-root-observations.json',pretty=TRUE,auto_unbox=TRUE,na='null')
cat('Directed perimeter root observations recorded\n')
