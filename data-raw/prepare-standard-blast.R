# Preserve historical observations; change only the file column layout.
# This is NOT a new BLAST run. The original outfmt7 records remain untouched.
pkgload::load_all(quiet=TRUE)
source_dir <- 'examples/blastn'
destination <- 'examples/blastn-standard'
dir.create(destination,showWarnings=FALSE)
records <- list()
for(path in sort(list.files(source_dir,'\\.o7$',full.names=TRUE))) {
 lines <- readLines(path,warn=FALSE)
 rows <- lines[nzchar(lines) & !startsWith(lines,'#')]
 fields <- strsplit(rows,'\t',fixed=TRUE)
 stopifnot(all(lengths(fields)==17L))
 output <- file.path(destination,paste0(sub('\\.o7$','',basename(path)),'.tsv'))
 writeLines(vapply(fields,function(x)paste(x[1:12],collapse='\t'),''),output,useBytes=TRUE)
 records[[length(records)+1L]] <- data.frame(source=path,file=output,rows=length(rows),
  source_md5=unname(tools::md5sum(path)),file_md5=unname(tools::md5sum(output)),
  transformation='First 12 standard BLAST fields; source HSP rows unchanged; not a rerun')
}
write.csv(do.call(rbind,records),file.path(destination,'provenance.csv'),row.names=FALSE)
