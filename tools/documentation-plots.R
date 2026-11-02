# Shared documentation assets only: one PNG device, explicit sizes, no release matrix.
pkgload::load_all(quiet = TRUE)
text <- readLines('tools/documentation-examples.R', warn = FALSE)
fingerprints <- jsonlite::fromJSON(system2('python3', 'tools/documentation-fingerprints.py', stdout = TRUE))
markers <- grep('^# ---- [a-z]+ ----$', text)
chunks <- setNames(lapply(seq_len(length(markers)-1L), function(i)
  text[seq.int(markers[i]+1L,markers[i+1L]-1L)]),
  sub('^# ---- ([a-z]+) ----$','\\1',text[markers[-length(markers)]]))
env <- new.env(parent=globalenv())
for(name in setdiff(names(chunks),c('save','preview','imports','firstsave'))) {
  cat('Example:',name,'\n')
  eval(parse(text=chunks[[name]]),env)
}
env$p_labelauto <- env$label_plots$auto
env$p_labelradial <- env$label_plots$radial
env$p_labelarc <- env$label_plots$arc
assets <- data.frame(
 id=c('chord','circular','collinear','tracks','focus','labels','realchord','reallinear','plasmid',
   'density','pointlinks','featurestack','multitrack','homologs','decorations','labelauto','labelradial','labelarc','styling','paf'),
 object=c('p_chord','p_circular','p_collinear','p_tracks','p_focus','p_labels','p_real_chord','p_real_linear','p_plasmid',
   'p_density','p_points','p_stack','p_multitrack','p_homologs','p_decorated','p_labelauto','p_labelradial','p_labelarc','p_styled','p_paf'),
 width=c(10,9,10,9,10,10,10,10,10,10,10,10,9,10,10,10,10,10,10,10),
 height=c(8,8,6,8,6,6,8,6,9,8,8,5,8,5,6,8,8,8,6,5),
 units='in',dpi=144,device='png',background='white',generator='tools/documentation-plots.R')
assets <- rbind(assets, transform(assets[1, ], id="first", object="p_first", width=10, height=5))
assets$input <- ifelse(assets$id=='paf','Public source-backed EMALE minimap2 PAF pair; original row values and tags retained',ifelse(assets$id %in% c('realchord','reallinear','density'),
 'Original repository FASTA/GFF3; unchanged HSPs converted to the standard BLAST 12-column layout',
 ifelse(assets$id %in% c('plasmid','circular','tracks','multitrack','featurestack'),
 'Bundled pUC19; curated annotations or exact catalogue matches; GC computed from its DNA',
 ifelse(assets$id=='homologs','Bundled pUC19 and pETDuet-1 AmpR; cluster membership from exact oriented CDS DNA equality',
 'Bundled seq_data_example/ribbon_data_example/gene_data_example; explicit minimal columns and restored biological source_strand'))))
assets$width[assets$id %in% c("labelauto", "labelradial", "labelarc")] <- 12
assets$height[assets$id %in% c("labelauto", "labelradial", "labelarc")] <- 10
assets <- rbind(assets, data.frame(
 id=c("showcasechord","showcaseplasmid","showcasecollinear"),
 object=c("p_show_chord","p_show_plasmid","p_show_collinear"),
 width=c(12,16.5,12),height=c(10,8,8),units="in",dpi=144,device="png",
 background="#FFFCF5",generator="tools/documentation-plots.R",
 input=c(
  "Four original FASTA/GFF3 genomes; at most three observed BLAST HSPs per pair, length >=900 bp and identity >=75%; explicit backbone curvature",
  "Bundled pETDuet-1; 20 curated features, five exact reference primers, commercial_all unique_6plus restriction subset",
  "Same four original genomes and observed HSP subset; all CDS, selected long-CDS highlights and GC computed from original FASTA in 2 kb windows"
 )
))
learning_ids <- c("learnseq", "learnlinks", "learngenes", "learnlabels", "linearseq", "linearlinks", "linearlabels", "plasmidbackbone", "plasmidfeatures", "plasmidtext")
learning_objects <- c("p_learn_seq", "p_learn_links", "p_learn_genes", "p_learn_labels", "p_linear_seq", "p_linear_links", "p_linear_labels", "p_plasmid_backbone", "p_plasmid_features", "p_plasmid_text")
catalog_ids <- c("catalogpsb", "catalogpuc", "catalogpbr", "catalogblue", "catalogpet", "catalogduet", "catalogdonr", "catalogpcdna", "catalogtre", "catalogtripz", "catalogcambia", "catalogearley", "catalogcas")
for (i in seq_along(catalog_ids)) assign(paste0("p_", catalog_ids[i]), env$catalog_plots[[i]], env)
extra_ids <- c(learning_ids, catalog_ids)
extra_objects <- c(learning_objects, paste0("p_", catalog_ids))
extra <- assets[rep(1L, length(extra_ids)), ]
extra$id <- extra_ids
extra$object <- extra_objects
extra$width <- 10
extra$height <- c(rep(8, 4), rep(6, 3), rep(8, 3), rep(8, 13))
extra$input <- c(rep("Progressive layers using the same four repository FASTA/GFF3 records and explicitly selected observed HSPs/CDS", 7), rep("Bundled pUC19; stepwise backbone, source annotations and labels", 3), paste("Bundled source-backed reference", env$catalog_summary$accver, "with current feature recognition; see plasmid-catalog.csv"))
assets <- rbind(assets, extra)
assets$background <- "#FFFCF5"
requested <- commandArgs(trailingOnly=TRUE)
if(length(requested) && any(!requested %in% assets$id)) stop('Unknown documentation figure ID')
if(length(requested)) {
 previous <- read.csv('website/assets/plots.csv')
 current <- unname(unlist(fingerprints$assets)[previous$id])
 if(!'dependency_md5' %in% names(previous) ||
    any(!previous$id[previous$dependency_md5 != current] %in% requested))
  stop('Changed plot dependencies extend beyond requested figures; regenerate all affected assets')
}
selected <- if(length(requested)) which(assets$id %in% requested) else seq_len(nrow(assets))
for(i in selected) {
 cat('Figure:',assets$id[i],'\n')
 plot <- get(assets$object[i],env)
 stopifnot(inherits(plot$theme$panel.background, 'element_blank'),
   identical(toupper(plot$theme$plot.background$fill), '#FFFCF5'))
 ggplot2::ggsave(paste0('man/figures/',assets$id[i],'.png'),plot,
   width=assets$width[i],height=assets$height[i],units='in',dpi=144,device='png',bg=assets$background[i])
}
for(name in c('p_mapped','p_sites','p_segments','p_migrated','p_offsets','p_bundle','p_manual_labels','p_marked','p_parameters'))
 invisible(ggplot2::ggplot_build(get(name,env)))
assets$file<-paste0('man/figures/',assets$id,'.png')
assets$md5<-unname(tools::md5sum(assets$file))
stopifnot(identical(fingerprints$code, unname(tools::md5sum('tools/documentation-examples.R'))))
assets$code_md5 <- if(length(requested)) {
 ifelse(assets$id %in% requested, fingerprints$code, previous$code_md5[match(assets$id, previous$id)])
} else rep(fingerprints$code, nrow(assets))
assets$dependency_md5 <- unname(unlist(fingerprints$assets)[assets$id])
write.csv(assets,'website/assets/plots.csv',row.names=FALSE)
catalog <- env$catalog_summary
catalog$figure <- catalog_ids
write.csv(catalog, 'website/assets/plasmid-catalog.csv', row.names=FALSE)
results<-list(
 showcase=list(sequences=nrow(env$show_sequences),CDS=nrow(env$show_genes),
   observed_HSPs=nrow(env$show_alignments),displayed_HSPs=nrow(env$show_links),
   selected_labels=nrow(env$show_labels),displayed_CDS=nrow(env$show_features),
   collinear_HSPs=nrow(env$show_adjacent)+nrow(env$show_cross),
   collinear_CDS=nrow(env$show_linear_features),
   plasmid_features=nrow(env$show_plasmid_features),plasmid_primers=nrow(env$show_plasmid_primers),
   plasmid_cut_positions=length(unique(env$show_plasmid_sites$position))),
 real=list(sequences=nrow(env$seq_real),features=nrow(env$feature_real),displayed_chord_features=nrow(env$real_display_features),raw_links=nrow(env$link_real),
   displayed_links=nrow(env$real_filtered$data),labelled_features=nrow(env$real_labels)),
 cleaning=list(input_rows=nrow(env$problem_links),retained_rows=nrow(env$repaired$ribbon_data)),
 ribbons=list(filtered=nrow(env$real_filtered$data),deduplicated=nrow(env$unique_links$data),
   merged=nrow(env$merged_links$data),bundled=nrow(env$bundled_links$data)),
 plasmid=list(sequence_length=get("plasmid_example_pUC19", envir = env)$length,features=nrow(env$plasmid_features),
   primers=nrow(env$plasmid_primers),all_sites=nrow(env$plasmid_sites),display_sites=nrow(env$display_sites)),
 projection=list(genomic_features=nrow(env$genomic),projected_links=nrow(env$projected_links)),
 minimal=list(sequence_columns=ncol(env$seqs),link_columns=ncol(env$links),gene_columns=ncol(env$genes)),
 protein=list(input_HSPs=nrow(env$protein_source),complete_query_CDS=nrow(env$emale_coding),retained_query_HSPs=nrow(env$protein_kept),query_aa_start=env$selected_hsp$qstart,query_aa_end=env$selected_hsp$qend,projected_regions=nrow(env$projected_query)),
 optional=list(genbank=if(is.null(env$genbank_result))'SKIP: Biopython unavailable' else 'PASS',
   ranges=if(is.null(env$range_inputs))'SKIP: GenomicRanges unavailable' else 'PASS',
   interactions=if(is.null(env$interaction_inputs))'SKIP: InteractionSet or GenomicRanges unavailable' else 'PASS'))
writeLines(jsonlite::toJSON(results,auto_unbox=TRUE,pretty=TRUE),'website/assets/example-results.json')
stopifnot(nrow(env$genomic)>0,nrow(env$primers)>0,nrow(env$sites)>0,nrow(env$rows)>0,
 nrow(env$homologs$data)>0,nrow(env$projected_links)>0,env$hierarchy$valid)
cat('Shared examples and',length(selected),'documentation PNG exports completed;',nrow(assets),'assets recorded.\n')
