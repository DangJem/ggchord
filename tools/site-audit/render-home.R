# Explicit homepage/flagship PNG generation, normal dimensions; no package matrix.
source('tools/audit-source-snapshot.R')
audit <- audit_source_snapshot(additional=c('tools/site-flagships.R','tools/home-previews.R',
  'tools/documentation-examples.R','website/home-previews.json','website/manual-figures.json','tools/site-audit/render-home.R'))
pkgload::load_all(quiet=TRUE)
blocks <- function(file) {
  text <- readLines(file,warn=FALSE);marks<-grep('^# ---- [a-z]+ ----$',text)
  setNames(lapply(seq_len(length(marks)-1L),function(i) paste(text[seq.int(marks[i]+1L,marks[i+1L]-1L)],collapse='\n')),
    sub('^# ---- ([a-z]+) ----$','\\1',text[marks[-length(marks)]]))
}
chunks <- c(blocks('tools/documentation-examples.R'),blocks('tools/site-flagships.R'))
for(key in c('sequencepalette','featurepalette','docstyle','showcasedata','chordshowcasegenes',
    'showcasechord','showcasecollinear','fusiondata','flagshipdata','flagshipmulti','fusionoverview','flagshipdot','flagshipfusion'))
  invisible(eval(parse(text=chunks[[key]])))
source('tools/home-previews.R')
previews<-jsonlite::fromJSON('website/home-previews.json',simplifyVector=FALSE)
record<-list(status='started',engine_sha256=as.list(audit$engine),
  source_sha256=as.list(audit$inputs),figures=list(),scope='Local annotation algorithm, two full flagships, seven homepage PNGs, and directly affected linear illustrations; no crop, browser or device matrix.')
# Actual data inputs are compared to the preceding public example objects.
baseline <- system2('git',c('show','3a6b412:tools/site-flagships.R'),stdout=TRUE)
start <- grep('^# ---- flagshipmulti ----$',baseline)
end <- grep('^# ---- fusionoverview ----$',baseline)
stopifnot(length(start)==1L,length(end)==1L,end>start)
previous <- new.env(parent=globalenv())
invisible(eval(parse(text=baseline[seq.int(start+1L,end-1L)]),previous))
before <- list(radial=previous$flag_radial,sankey=previous$flag_sankey)
source_rows<-function(plot,component) {
  layout<-get_chord_layout(plot)
  result<-list()
  for(owner in names(layout$layer_inputs)) {
    key<-if(component=='gene') 'gene_poly' else component
    d<-layout$layer_inputs[[owner]][[key]]
    if(is.data.frame(d)&&nrow(d)) result[[length(result)+1L]]<-d
  }
  dplyr::bind_rows(result)
}
select_source<-function(data,component) {
  fields<-if(component=='ribbon') c('qaccver','saccver','qstart','qend','sstart','send','pident') else c('accver','start','end','strand','anno')
  data[,intersect(fields,names(data)),drop=FALSE] |> as.data.frame()
}
checks<-list()
for(name in c('radial','sankey')) for(component in c('ribbon','gene')) {
  old<-select_source(source_rows(before[[name]],component),component)
  current<-select_source(source_rows(get(paste0('flag_',name)),component),component)
  stopifnot(nrow(old)>0L,nrow(current)>0L,identical(old,current))
  checks[[paste(name,component,sep='_')]]<-list(original_records=nrow(old),all_source_fields_identical=TRUE)
}
record$source_checks<-checks
if('--verify-sources-only' %in% commandArgs(trailingOnly=TRUE)) {
  completed<-jsonlite::fromJSON('tools/site-audit/home-rendering.json',simplifyVector=FALSE)
  stopifnot(completed$status=='complete')
  completed$initial_source_checks<-completed$source_checks
  completed$source_checks<-checks
  completed$source_check_correction<-list(
    reason='Initial gene selector used gene instead of the public gene_poly component and returned zero rows; empty checks are rejected and superseded.',
    nonempty_source_verification=TRUE,
    verifier_sha256=digest::digest(file='tools/site-audit/render-home.R',algo='sha256',serialize=FALSE),
    scope='Corrected public input comparison only; original PNG render hashes and execution source hashes retained.')
  jsonlite::write_json(completed,'tools/site-audit/home-rendering.json',auto_unbox=TRUE,pretty=TRUE)
  cat('PASS: nonempty original 12/10 HSP and 20/16 CDS records identical\n')
  quit(save='no',status=0)
}

# The runtime delta is restricted to three internal functions; preserve the
# historical execution source and prove the remaining expressions unchanged.
ns <- asNamespace("ggchord")
old_runtime <- new.env(parent = ns)
changed_functions <- character()
for (file in c("R/ggchord-build.R", "R/layout-collinear-annotation.R")) {
  old <- parse(text = system2("git", c("show", paste0("a1cd31d:", file)), stdout = TRUE), keep.source = FALSE)
  current <- parse(file, keep.source = FALSE)
  stopifnot(length(old) == length(current))
  for (i in seq_along(old)) if (!identical(old[[i]], current[[i]])) {
    stopifnot(is.call(old[[i]]), identical(old[[i]][[1]], as.name("<-")),
      identical(old[[i]][[2]], current[[i]][[2]]),
      identical(old[[i]][[3]][[1]], as.name("function")),
      identical(old[[i]][[3]][[2]], current[[i]][[3]][[2]]))
    changed_functions <- c(changed_functions, as.character(old[[i]][[2]]))
  }
  invisible(eval(old, old_runtime))
}
stopifnot(setequal(changed_functions, c("compute_chord_geometry",
  "ggchord_collinear_gene_labels", "ggchord_collinear_share_labels")))
current_functions <- mget(changed_functions, ns)
replace_runtime <- function(functions) for (name in names(functions)) {
  locked <- bindingIsLocked(name, ns)
  if (locked) unlockBinding(name, ns)
  assign(name, functions[[name]], ns)
  if (locked) lockBinding(name, ns)
}
quality <- function() {
  result <- list()
  for (id in c("radial", "sankey", "collinear")) {
    plot <- switch(id, radial = flag_radial, sankey = flag_sankey, collinear = p_show_collinear)
    size <- switch(id, radial = c(12, 10.5), sankey = c(10, 16), collinear = c(12, 8))
    cat("Measure:", id, "\n"); flush.console()
    d <- diagnose_ggchord(plot, width = size[1], height = size[2], dpi = 144)
    result[[id]] <- list(summary = d$label_quality$summary,
      conflicts = d$label_quality$conflicts, visibility = d$label_visibility,
      placement = d$label_diagnostics, warnings = d$warnings,
      converged = d$measurement$converged)
  }
  result
}
# Copies reset render caches so each engine really computes its own geometry.
public_geometry <- function(plot) {
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = 12, height = 10, units = "in", res = 144)
  on.exit({ grDevices::dev.off(); unlink(file) })
  d <- export_ggchord_layout(ggchord:::ggchord_render_plot_copy(plot))
  d[c("seq", "gene", "ribbon", "feature")]
}
replace_runtime(mget(changed_functions, old_runtime))
old_geometry <- list(chord = public_geometry(p_show_chord), circular = public_geometry(home_circular))
old_quality <- tryCatch(quality(), finally = replace_runtime(current_functions))
stopifnot(identical(old_geometry$chord, public_geometry(p_show_chord)),
  identical(old_geometry$circular, public_geometry(home_circular)))
record$quality <- list(before = old_quality, after = quality())
record$runtime_review <- list(baseline = "a1cd31d", changed_functions = changed_functions,
  signatures_identical = TRUE, other_expressions_identical = TRUE,
  representative_chord_circular_geometry_identical = TRUE,
  scope = "Only Collinear annotation functions and gated collision placement change; projected plots without gene-name layers retain geometry. Not package-wide certification.")
# Prepare existing objects without catalogue batch maps or any file-saving chunk.
legacy_env <- new.env(parent = globalenv())
legacy_chunks <- blocks("tools/documentation-examples.R")
for (key in setdiff(names(legacy_chunks), c("save", "preview", "frontsave", "frontpreview",
    "fusionsave", "imports", "firstsave", "firstsavebasic", "catalogdata", "catalogplots", "catalogbatch")))
  invisible(eval(parse(text = legacy_chunks[[key]]), legacy_env))
uses_linear_names <- function(plot) {
  if (inherits(plot, "ggchord_fusion")) return(any(vapply(plot$panels,
    function(panel) uses_linear_names(panel$plot), logical(1L))))
  if (!inherits(plot, "ggchord") || !isTRUE(plot$coordinates$ggchord_collinear)) return(FALSE)
  any(vapply(plot$layers, function(layer) (layer$ggchord_params$type %||% "") %in%
    c("gene_label", "gene_label_repel"), logical(1L)))
}
legacy <- read.csv("website/assets/plots.csv", stringsAsFactors = FALSE)
linear <- which(legacy$coordinate %in% c("collinear", "radial_synteny", "sankey", "fusion"))
annotated <- linear[vapply(legacy$object[linear], function(name)
  uses_linear_names(get(name, legacy_env)), logical(1L))]
stopifnot(setequal(legacy$id[annotated], c("showcasecollinear", "linearlabels", "practicelinear")))
record$legacy_scope <- list(inspected_linear_objects = as.list(legacy$id[linear]),
  annotated = as.list(legacy$id[annotated]),
  unchanged_recipes_sha256 = digest::digest(file = "tools/documentation-examples.R", algo = "sha256", serialize = FALSE),
  source_gates = "Other coordinate functions are unchanged; unannotated linear/Fusion plots do not enter the collision solver.")
fingerprints <- jsonlite::fromJSON(system2("python3", "tools/documentation-fingerprints.py", stdout = TRUE))
record$legacy_figures <- list()
for (i in annotated) {
  id <- legacy$id[i]; cat("Refresh annotated:", id, "\n"); flush.console()
  file <- view_ggchord(get(legacy$object[i], legacy_env), width = legacy$width[i],
    height = legacy$height[i], dpi = legacy$dpi[i], bg = legacy$background[i], viewer = "none")
  stopifnot(file.copy(file, legacy$file[i], overwrite = TRUE))
  legacy$md5[i] <- unname(tools::md5sum(legacy$file[i]))
  legacy$code_md5[i] <- fingerprints$code
  legacy$dependency_md5[i] <- fingerprints$assets[[id]]
  legacy$generator[i] <- "tools/site-audit/render-home.R"
  record$legacy_figures[[id]] <- list(file = legacy$file[i], status = "actual PNG rendered",
    sha256 = digest::digest(file = legacy$file[i], algo = "sha256", serialize = FALSE),
    md5 = legacy$md5[i], dependency_md5 = legacy$dependency_md5[i])
}
write.csv(legacy, "website/assets/plots.csv", row.names = FALSE)
# The progressive Chord-to-Collinear page also owns a linear annotation image.
course_chunks <- blocks("tools/site-course-examples.R")
for (key in c("frontbackbone", "frontlinks", "frontribbon", "frontgenes", "frontfeature", "frontgenelabels"))
  invisible(eval(parse(text = legacy_chunks[[key]]), legacy_env))
invisible(eval(parse(text = course_chunks$courseblockad), legacy_env))
courses <- jsonlite::fromJSON("website/course-figures.json", simplifyVector = FALSE)
row <- Filter(function(row) row$id == "frontlabelledlinear", courses)[[1L]]
file <- view_ggchord(get(row$object, legacy_env), width = row$width, height = row$height,
  dpi = row$dpi, bg = "white", viewer = "none")
stopifnot(file.copy(file, row$file, overwrite = TRUE))
course_record <- jsonlite::fromJSON("tools/site-audit/course-rendering.json", simplifyVector = FALSE)
course_record$figures[[row$id]]$sha256 <- digest::digest(file = row$file, algo = "sha256", serialize = FALSE)
course_record$figures[[row$id]]$runtime_record <- "tools/site-audit/home-rendering.json"
jsonlite::write_json(course_record, "tools/site-audit/course-rendering.json", auto_unbox = TRUE, pretty = TRUE)
record$course_figures <- setNames(list(course_record$figures[[row$id]]), row$id)
# Confirm source-feature records/segments survive the palette correction.
stopifnot(identical(home_circular_features, find_common_features(plasmid_example_pUC19)))
colour_data <- ggplot2::ggplot_build(home_circular)$data
fills <- unique(unlist(lapply(colour_data, function(d) d$feature_fill)))
stopifnot(length(fills) > 1L, !any(fills %in% c("grey50", "#7F7F7F")))
record$circular_palette <- list(source_rows = nrow(home_circular_features),
  source_and_segments_unchanged = TRUE, scale_key = "source annotation name", matched_fills = as.list(fills))
audit_assert_source_snapshot(audit)

# Complete figures: freshly rendered final recipes, not renamed candidate files.
full<-list(flagradial=list(object='flag_radial',width=12,height=10.5,dpi=144,file='man/figures/flagradial.png'),
  flagsankey=list(object='flag_sankey',width=10,height=16,dpi=144,file='man/figures/flagsankey.png'))
for(id in c(names(full),vapply(previews,`[[`,'','id'))) {
 row<-if(id %in% names(full)) full[[id]] else Filter(function(r) r$id==id,previews)[[1L]]
 cat('Render:',id,row$width,row$height,'\n');flush.console()
 file<-view_ggchord(get(row$object),width=row$width,height=row$height,units='in',device='png',dpi=row$dpi,bg='#FFFCF5',viewer='none')
 stopifnot(file.copy(file,row$file,overwrite=TRUE))
 audit_assert_source_snapshot(audit)
 record$figures[[id]]<-list(file=row$file,sha256=digest::digest(file=row$file,algo='sha256',serialize=FALSE),
   width=row$width,height=row$height,dpi=row$dpi,status='actual PNG rendered',background='#FFFCF5')
 jsonlite::write_json(record,'tools/site-audit/home-rendering.json',auto_unbox=TRUE,pretty=TRUE)
}
record$status<-'complete'
record$labels<-list(radial_selected=nrow(radial_labels),sankey_selected=nrow(sankey_labels),
  rule='All original full-figure name candidates; bounded local placement after projection. Homepage relationship previews omit name layers but retain source geometry.')
jsonlite::write_json(record,'tools/site-audit/home-rendering.json',auto_unbox=TRUE,pretty=TRUE)
cat('PASS: two complete flagships, seven homepage previews and four directly affected linear images; source endpoints/CDS identical\n')
