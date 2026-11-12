# Direct v0.21.0 public behavior checks, not a package/release test suite.
pkgload::load_all(quiet = TRUE)
library(dplyr)
out <- 'tools/long-term-audit'
started <- proc.time()[['elapsed']]
expect_error <- function(expr, pattern) {
  e <- tryCatch({force(expr);NULL},error=identity)
  stopifnot(inherits(e,'error'),grepl(pattern,conditionMessage(e)))
}
dna <- read_aligned_fasta('examples/alignment/opuntia-aligned.fasta','DNA',source_file=TRUE)
aa <- read_aligned_fasta('examples/alignment/protein-aligned.fasta','AA',source_file=TRUE)
for (d in list(dna,aa)) {
  stopifnot(all(is.na(d$source_position)==d$is_gap),any(d$is_gap),all(file.exists(d$.source_file)),
    all(d$.source_line>0),all(!d$source_verified))
  for (id in unique(d$accver)) {
    row <- filter(d,accver==id,!is_gap)
    stopifnot(identical(as.numeric(row$source_position),as.numeric(seq_len(nrow(row)))))
    forth <- map_alignment_position(d,id,row$source_position,inverse=TRUE)
    back <- map_alignment_position(d,id,forth$alignment_column)
    stopifnot(identical(back$source_position,row$source_position),all(forth$matched),all(back$matched))
  }
}
# Original aligned characters and row identities must survive adaptation exactly.
original <- read_fasta_lengths('examples/alignment/opuntia-aligned.fasta',keep_sequence=TRUE)
for (id in original$accver) stopifnot(paste0(dna$symbol[dna$accver==id],collapse='')==original$sequence[original$accver==id])
# A declared reverse orientation of observed aligned fragments: a coordinate
# transformation of upstream data, not new alignment/strand evidence.
reverse <- vapply(original$sequence,function(s) {
  text <- paste0(rev(strsplit(s,'',fixed=TRUE)[[1]]),collapse='')
  chartr('ACGTRYSWKMBDHVNacgtryswkmbdhvn','TGCAYRSWMKVHDBNtgcayrswmkvhdbn',text)
},'')
names(reverse) <- original$accver
source <- data.frame(accver=original$accver,source_accver=original$accver,
  start=nchar(gsub('[-.]','',original$sequence)),strand='-',sequence=gsub('[-.]','',original$sequence))
r <- as_ggchord_alignment(reverse,'DNA',source)
stopifnot(all(r$source_verified))
for(id in original$accver) {
  actual <- filter(r,accver==id,!is_gap)
  stopifnot(all(diff(actual$source_position)== -1),
    identical(map_alignment_position(r,id,actual$source_position,inverse=TRUE)$alignment_column,actual$alignment_column))
}
# Explicit synthetic format edge: unknown/out-of-window, RNA and missing cells.
a <- as_ggchord_alignment(c(A='aC-G',B='ACTG'),'DNA')
m <- map_alignment_position(a,'A',c(3,5,NA))
stopifnot(m$matched[1],m$is_gap[1],is.na(m$source_position[1]),all(!m$matched[2:3]))
stopifnot(map_alignment_position(filter(a,alignment_column<=2),'A',3,inverse=TRUE)$matched==FALSE,
  as_ggchord_alignment(c(A='AU-G',B='AUGG'),'RNA')$alphabet[1]=='RNA')
expect_error(as_ggchord_alignment(c(A='AC',B='A-G'),'DNA'),'same positive width')
expect_error(as_ggchord_alignment(c(A='AC-G',B='ACTG'),'RNA'),'alphabet')
expect_error(as_ggchord_alignment(c(A='AC'), 'DNA',data.frame(accver='A',start=1,strand='+',sequence='TT')),'disagree')
# Real local windows; retain original alignment columns and fragment positions.
dna_ids <- unique(dna$accver)[c(1,3,7)]
aa_ids <- unique(aa$accver)[c(2,3,4)]
dna_window <- filter(dna,accver %in% dna_ids,between(alignment_column,321,380))
aa_window <- filter(aa,accver %in% aa_ids,between(alignment_column,221,260))
stopifnot(any(dna_window$is_gap),any(aa_window$is_gap))
# Dataset selection is explicit; adjacent rows do not create new alignments.
plot_window <- function(d,title) {
  d$accver <- factor(d$accver,levels=rev(unique(d$accver)))
  ggplot2::ggplot(d,aes(alignment_column,accver,label=symbol))+
    geom_residue(aes(fill=symbol),overflow='error',size=3)+
    ggplot2::scale_fill_manual(values=c(A='#B8D8BA',C='#A6CBE3',G='#F2D28C',T='#EAA8A8',
      N='#D8D8D8','-'='#F5F5F5','.'='#F5F5F5'),na.value='#E5EDF2',guide='none')+
    ggplot2::labs(x='Original alignment column (not source position)',y=NULL,title=title)+
    ggplot2::theme_minimal(base_size=11)+
    ggplot2::theme(panel.grid=ggplot2::element_blank(),axis.text.y=ggplot2::element_text(size=8))
}
p_dna <- plot_window(dna_window,'Opuntia DNA: original gapped alignment columns 321–380')
p_aa <- plot_window(aa_window,'Protein alignment: original columns 221–260 (AA)')+
  ggplot2::scale_fill_hue(l=85,c=35,guide='none')
before <- dna_window
# Actual device draw triggers deferred metrics, shared with ordinary ggsave.
draw_at <- function(plot,width,height) {
  close <- ggchord:::ggchord_measurement_device(text="",width=width,height=height)
  on.exit(close())
  result <- tryCatch({grid::grid.newpage();grid::grid.draw(ggplot2::ggplotGrob(plot));NULL},error=identity)
  close()
  if (inherits(result,'error')) stop(result)
}
render_devices <- grDevices::dev.list()
draw_at(p_dna,12,3.4)
small <- p_dna+ggplot2::theme(axis.text.y=ggplot2::element_blank())
expect_error(draw_at(small,2.5,2),'cannot fit')
# A replacement layer uses the same table and finite font, explicit tile policy.
p_tile <- ggplot2::ggplot(dna_window,aes(alignment_column,accver,label=symbol))+
  geom_residue(aes(fill=symbol),size=3,overflow='tile')+
  ggplot2::scale_fill_manual(values=c(A='#B8D8BA',C='#A6CBE3',G='#F2D28C',T='#EAA8A8',N='#D8D8D8','-'='#F5F5F5','.'='#F5F5F5'))+
  ggplot2::theme_void()+ggplot2::theme(legend.position='none')
draw_at(p_tile,2.5,2)
draw_at(p_dna+ggplot2::scale_x_reverse(),12,3.4)
# Independent native layers, including a standard data callback and user scale.
first_row <- filter(dna_window,accver==dna_ids[1])
other_rows <- filter(dna_window,accver!=dna_ids[1])
multi <- ggplot2::ggplot(dna_window,aes(alignment_column,accver,label=symbol))+
  geom_residue(aes(fill=symbol),data=function(d) d[d$accver==dna_ids[1],])+
  geom_residue(data=other_rows,fill="#123456")+
  ggplot2::scale_fill_manual(values=stats::setNames(rep("#B8D8BA",length(unique(dna_window$symbol))),unique(dna_window$symbol)))
built <- ggplot2::ggplot_build(multi)
stopifnot(nrow(built$data[[1]])==nrow(first_row),nrow(built$data[[2]])==nrow(other_rows),
  all(built$data[[1]]$fill=="#B8D8BA"),all(built$data[[2]]$fill=="#123456"))
expect_error(draw_at(p_dna+ggplot2::coord_polar(),12,3.4),'Cartesian')
expect_error(draw_at(p_dna+ggplot2::scale_x_log10(),12,3.4),'identity or reverse')
stopifnot(identical(before,dna_window))
ggplot2::ggsave(file.path(out,'v0.21.0-alignment-dna.png'),p_dna,width=12,height=3.4,dpi=120,bg='white')
ggplot2::ggsave(file.path(out,'v0.21.0-alignment-protein.png'),p_aa,width=12,height=3.4,dpi=120,bg='white')
ggplot2::ggsave(file.path(out,'v0.21.0-alignment-tiles.png'),p_tile,width=2.5,height=2,dpi=120,bg='white')
stopifnot(identical(grDevices::dev.list(),render_devices))
utils::write.csv(dna_window,file.path(out,'v0.21.0-dna-columns.csv'),row.names=FALSE)
utils::write.csv(aa_window,file.path(out,'v0.21.0-protein-columns.csv'),row.names=FALSE)
jsonlite::write_json(list(version='0.21.0.9000',source_md5=as.list(tools::md5sum(c(
  'examples/alignment/opuntia-aligned.fasta','examples/alignment/protein-aligned.fasta'))),
  dna_records=length(unique(dna$accver)),aa_records=length(unique(aa$accver)),
  dna_cells=nrow(dna),aa_cells=nrow(aa),gap_cells=sum(dna$is_gap)+sum(aa$is_gap),
  exact_source_column_roundtrip=TRUE,negative_orientation_and_source_verification=TRUE,
  note='Source positions refer to ungapped input fragments, not independently verified whole genomes',
  window_columns=list(DNA=c(321,380),AA=c(221,260)),
  physical_font_budget=list(wide='letters at 12 x 3.4 in',narrow='explicit error or tiles at 2.5 x 2 in',size_mm=3),
  source_table_unchanged=TRUE,graphics_devices_restored=TRUE,independent_layers_and_user_scale=TRUE,elapsed_seconds=unname(proc.time()[['elapsed']]-started),
  not_in_scope=c('MSA computation','CDS translation/phase','CIGAR decoding','macro view projection','release device matrix')),
  file.path(out,'v0.21.0-alignment-results.json'),pretty=TRUE,auto_unbox=TRUE)
message('Alignment: real DNA/AA gaps, exact maps, declared reverse, source verification and physical text/tile drawing passed')
