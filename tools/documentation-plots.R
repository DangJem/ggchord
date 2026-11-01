# Shared documentation assets only: one PNG device, explicit sizes, no release matrix.
pkgload::load_all(quiet = TRUE)
text <- readLines('tools/documentation-examples.R', warn = FALSE)
markers <- grep('^# ---- [a-z]+ ----$', text)
chunks <- setNames(lapply(seq_len(length(markers)-1L), function(i)
  text[seq.int(markers[i]+1L,markers[i+1L]-1L)]),
  sub('^# ---- ([a-z]+) ----$','\\1',text[markers[-length(markers)]]))
env <- new.env(parent=globalenv())
for(name in setdiff(names(chunks),c('save','preview','imports'))) {
  cat('Example:',name,'\n')
  eval(parse(text=chunks[[name]]),env)
}
env$p_labelauto <- env$label_plots$auto
env$p_labelradial <- env$label_plots$radial
env$p_labelarc <- env$label_plots$arc
assets <- data.frame(
 id=c('chord','circular','collinear','tracks','focus','labels','realchord','reallinear','plasmid',
   'density','pointlinks','featurestack','multitrack','homologs','decorations','labelauto','labelradial','labelarc','styling'),
 object=c('p_chord','p_circular','p_collinear','p_tracks','p_focus','p_labels','p_real_chord','p_real_linear','p_plasmid',
   'p_density','p_points','p_stack','p_multitrack','p_homologs','p_decorated','p_labelauto','p_labelradial','p_labelarc','p_styled'),
 width=c(7,7,8,7,8,8,10,10,9,10,7,9,7,8,8,8,8,8,8),
 height=c(6,6,5,6,5,5,8,6,8,8,6,4,6,5,5,5,5,5,5),
 units='in',dpi=144,device='png',background='white',generator='tools/documentation-plots.R')
assets$input <- ifelse(assets$id %in% c('realchord','reallinear','density'),
 'Original examples/fasta, examples/gff3 and examples/blastn; explicit filters in shared code',
 ifelse(assets$id=='plasmid','plasmid_example_pUC19 copied from examples/plasmid/pUC19.fna',
 'Explicit synthetic teaching inputs; not biological observations'))
requested <- commandArgs(trailingOnly=TRUE)
if(length(requested) && any(!requested %in% assets$id)) stop('Unknown documentation figure ID')
if(length(requested)) {
 previous <- read.csv('website/assets/plots.csv')
 if(any(previous$code_md5 != unname(tools::md5sum('tools/documentation-examples.R'))))
  stop('Shared example code changed; regenerate all figures before selecting individual exports')
}
selected <- if(length(requested)) which(assets$id %in% requested) else seq_len(nrow(assets))
for(i in selected) {
 cat('Figure:',assets$id[i],'\n')
 ggplot2::ggsave(paste0('man/figures/',assets$id[i],'.png'),get(assets$object[i],env),
   width=assets$width[i],height=assets$height[i],units='in',dpi=144,device='png',bg='white')
}
for(name in c('p_mapped','p_sites','p_segments','p_migrated','p_offsets','p_bundle','p_manual_labels','p_marked'))
 invisible(ggplot2::ggplot_build(get(name,env)))
assets$file<-paste0('man/figures/',assets$id,'.png')
assets$md5<-unname(tools::md5sum(assets$file))
assets$code_md5<-unname(tools::md5sum('tools/documentation-examples.R'))
write.csv(assets,'website/assets/plots.csv',row.names=FALSE)
results<-list(
 real=list(sequences=nrow(env$seq_real),features=nrow(env$feature_real),displayed_chord_features=nrow(env$real_display_features),raw_links=nrow(env$link_real),
   displayed_links=nrow(env$real_filtered$data),labelled_features=nrow(env$real_labels)),
 cleaning=list(input_rows=nrow(env$problem_links),retained_rows=nrow(env$repaired$ribbon_data)),
 ribbons=list(filtered=nrow(env$real_filtered$data),deduplicated=nrow(env$unique_links$data),
   merged=nrow(env$merged_links$data),bundled=nrow(env$bundled_links$data)),
 plasmid=list(sequence_length=get("plasmid_example_pUC19", envir = env)$length,features=nrow(env$plasmid_features),
   primers=nrow(env$plasmid_primers),all_sites=nrow(env$plasmid_sites),display_sites=nrow(env$display_sites)),
 projection=list(genomic_features=nrow(env$genomic),projected_links=nrow(env$projected_links)),
 optional=list(genbank=if(is.null(env$genbank_result))'SKIP: Biopython unavailable' else 'PASS',
   ranges=if(is.null(env$range_inputs))'SKIP: GenomicRanges unavailable' else 'PASS'))
writeLines(jsonlite::toJSON(results,auto_unbox=TRUE,pretty=TRUE),'website/assets/example-results.json')
stopifnot(nrow(env$genomic)>0,nrow(env$primers)>0,nrow(env$sites)>0,nrow(env$rows)>0,
 nrow(env$homologs$data)>0,nrow(env$projected_links)>0,env$hierarchy$valid)
cat('Shared examples and',length(selected),'documentation PNG exports completed;',nrow(assets),'assets recorded.\n')
