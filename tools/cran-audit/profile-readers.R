# Explicit local reader measurements. Run in separate fresh processes:
# Rscript --vanilla tools/cran-audit/profile-readers.R vcf [small]
# Rscript --vanilla tools/cran-audit/profile-readers.R bigwig
args <- commandArgs(TRUE)
reader <- match.arg(args[1], c("vcf", "bigwig"))
dependency <- if(reader=='vcf') 'VariantAnnotation' else 'rtracklayer'
root <- getwd()
stopifnot(!dependency %in% loadedNamespaces())
measure <- function(expr) {
  elapsed <- system.time(value <- force(expr))
  list(value=value, timing=as.list(elapsed[c('user.self','sys.self','elapsed')]))
}
load <- measure(stopifnot(requireNamespace(dependency,quietly=TRUE)))
cat('dependency loaded:',dependency,load$timing$elapsed,'\n');flush.console()
package <- measure(pkgload::load_all(root,quiet=TRUE))
path <- if(reader=='vcf') system.file('extdata','chr22.vcf.gz',package='VariantAnnotation') else file.path(root,'inst/extdata/pUC19-gc-100.bw')
preparation <- NULL
if(length(args)>1 && args[2] %in% c("small", "stream")) {
  preparation <- measure({
    con <- gzfile(path, 'rt')
    if (args[2] == "stream") {
      lines <- character(); records <- 0L
      while (records < 3L) {
        line <- readLines(con, n = 1L)
        if (!length(line)) break
        lines <- c(lines, line)
        if (!startsWith(line, "#")) records <- records + 1L
      }
    } else lines <- readLines(con)
    close(con)
    header <- startsWith(lines, '#')
    path <- tempfile(fileext='.vcf')
    writeLines(c(lines[header],head(lines[!header],3L)),path)
  })
}
read <- if(reader=='vcf') function() ggchord::read_vcf(path) else function() ggchord::read_bigwig(path)
first <- measure(read())
cat('first read:',first$timing$elapsed,'\n');flush.console()
warm <- measure(read())
selected <- if(reader=='bigwig') measure(ggchord::read_bigwig(path,which=data.frame(accver='pUC19',start=50,end=150))) else NULL
preview <- if(reader=='vcf') first$value$features$variant[c('accver','start','end','feature_id')] else first$value
output <- measure(capture.output(print(utils::head(preview))))
result <- list(reader=reader,dependency=dependency,dependency_version=as.character(utils::packageVersion(dependency)),R=R.version.string,platform=R.version$platform,namespace_initially_unloaded=TRUE,source_package_load=package$timing,namespace_first_load=load$timing,example_preparation=preparation$timing,first_read=first$timing,warm_read=warm$timing,selected_read=selected$timing,preview_output=output$timing,output_lines=length(output$value),output_bytes=sum(nchar(output$value,type='bytes')),records=nrow(preview),input_bytes=file.info(path)$size,notes='Independent R process per decoder; cold namespace, filesystem caches not flushed. Reader timings exclude decoder namespace and source-package loading. Output captures only the displayed preview of at most six rows.')
output_dir <- Sys.getenv("GGCHORD_PROFILE_DIR", unset = tempdir())
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(result,file.path(output_dir,paste0('ggchord-',reader,if(length(args)>1) paste0('-', args[2]) else '', '-profile.json')),auto_unbox=TRUE,pretty=TRUE)
