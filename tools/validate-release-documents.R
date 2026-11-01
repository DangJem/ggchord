# Explicit release acceptance for source-backed documentation plots.
# Rscript tools/validate-release-documents.R [output-directory]
main <- function() {
  pkgload::load_all(quiet = TRUE)
  args <- commandArgs(trailingOnly = TRUE)
  out <- if(length(args)) args[1] else tempfile('ggchord-release-documents-')
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(NULL, width = 12, height = 10)
  on.exit(grDevices::dev.off(), add = TRUE)
  lines <- readLines('tools/documentation-examples.R', warn = FALSE)
  starts <- grep('^# ---- [a-z]+ ----$', lines)
  snippets <- setNames(lapply(seq_len(length(starts)-1L),function(i)
    lines[seq.int(starts[i]+1L,starts[i+1L]-1L)]),
    sub('^# ---- ([a-z]+) ----$','\\1',lines[starts[-length(starts)]]))
  env <- new.env(parent = globalenv())
  for(key in setdiff(names(snippets),c('save','preview','imports','firstsave')))
    eval(parse(text = snippets[[key]]), env)
  # Real bundled rows, independently assigned to two layers of each class.
  split_genes <- env$genes$accver %in% env$seqs$accver[1:2]
  split_links <- env$links$qaccver == env$seqs$accver[1]
  mixed <- ggchord(env$seqs) +
    geom_seq(seq_radius = c(1.1,1,.9,1.2), seq_curvature = c(.8,1.2,1,.7),
      seq_orientation = c(1,-1,1,-1), seq_gap = c(.03,.04,.06,.05)) +
    geom_gene(data = env$genes[split_genes, ], gene_fill = '#268DAB',
      position = position_strand(.1)) +
    geom_gene(data = env$genes[!split_genes, ], gene_fill = '#CA644D',
      position = position_strand(.1)) +
    geom_link_ribbon(data = env$links[split_links, ], fill = '#5D72BA', alpha = .6) +
    geom_link_ribbon(data = env$links[!split_links, ], fill = '#D68B39', alpha = .6) +
    geom_seq_label(size = 4) + coord_chord(rotation = 25) + env$doc_theme
  layers <- export_ggchord_layout(mixed, include = c('gene','ribbon'))
  stopifnot(length(unique(layers$gene$layer_id)) == 2L,
    length(unique(layers$ribbon$layer_id)) == 2L)
  plots <- list(chord = env$p_real_chord, circular = env$p_plasmid,
    collinear = env$p_real_linear, local = env$p_homologs,
    tracks = env$p_multitrack, labels_auto = env$label_plots$auto,
    labels_radial = env$label_plots$radial, labels_arc = env$label_plots$arc,
    independent_layers = mixed)
  records <- list()
  for(name in names(plots)) {
    # Only the plasmid labels need Unicode. Keep numerical colour bars on
    # standard PDF; Quartz can invert raster colour bars on some macOS builds.
    pdf <- if(name != 'circular') grDevices::pdf else if(isTRUE(capabilities('aqua'))) {
      function(filename, width, height, ...) grDevices::quartz(type = 'pdf',
        file = filename, width = width, height = height, ...)
    } else grDevices::cairo_pdf
    width <- if(startsWith(name,'labels_')) 12 else 10
    height <- if(startsWith(name,'labels_')) 10 else if(name=='local') 5 else 8
    for(format in c('png','pdf','svg')) {
      path <- file.path(out,paste0(name,'.',format))
      warnings <- character()
      withCallingHandlers(ggplot2::ggsave(path, plots[[name]], width = width,
        height = height, units = 'in', dpi = 144, bg = 'white',
        device = switch(format, pdf = pdf, svg = svglite::svglite, png = 'png')),
        warning = function(w) warnings <<- c(warnings, conditionMessage(w)))
      if(length(warnings)) stop(name,'/',format,': ',paste(warnings,collapse='; '))
      stopifnot(file.exists(path), file.info(path)$size > 0)
      if(format=='pdf') stopifnot(qpdf::pdf_length(path)==1L)
      if(format=='svg') stopifnot(xml2::xml_name(xml2::read_xml(path))=='svg')
      records[[length(records)+1L]] <- data.frame(case = name, format,
        width, height, units = 'in', file = basename(path), bytes = file.info(path)$size,
        sha256 = digest::digest(file=path,algo='sha256'), result = 'PASS')
      cat('PASS:',name,format,'\n')
    }
  }
  write.csv(do.call(rbind,records),file.path(out,'devices.csv'),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(out,'session.txt'))
  cat('PASS: 9 source-backed cases × PNG/PDF/SVG; independent layers and mixed sequence geometry.\n')
}
main()
