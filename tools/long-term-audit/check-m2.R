# Direct M2 checks; no package test/build/check or publication device matrix.
pkgload::load_all(quiet = TRUE)
library(dplyr)
out <- 'tools/long-term-audit'
fail <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(e, 'error'), grepl(pattern, conditionMessage(e)))
}
e <- jsonlite::read_json('examples/cds/expected.json', simplifyVector = TRUE)
stopifnot(digest::digest(file='examples/cds/NC_000932.gb',algo='sha256',serialize=FALSE)==e$source_sha256,
  digest::digest(file='examples/cds/NC_000932.fna',algo='sha256',serialize=FALSE)==e$fasta_sha256)
source <- read_fasta_lengths('examples/cds/NC_000932.fna', keep_sequence = TRUE, source_file = TRUE)
cds_rows <- lapply(seq_len(nrow(e$models)), function(i) {
  m <- e$models[i, ]; s <- m$segments[[1]]
  d <- data.frame(feature_id=m$feature_id,accver=m$accver,strand=m$strand,
    start=min(s$start),end=max(s$end),codon_start=m$codon_start,genetic_code=m$genetic_code)
  d$segments <- I(list(s)); d
})
cds_input <- bind_rows(cds_rows)
cds <- build_ggchord_cds(cds_input)
translated <- translate_ggchord_cds(cds, source, initiation='start')
expected <- e$models$expected_protein[match(translated$models$feature_id,e$models$feature_id)]
stopifnot(identical(sub('[*]$','',translated$models$protein_sequence),expected))
# Independently audit all genomic letters and codon positions, not endpoint interpolation.
for (id in cds$models$feature_id) {
  d <- filter(translated$nucleotides,accver==id)
  a <- filter(translated$residues,accver==id)
  observed <- substring(source$sequence,d$source_position,d$source_position)
  if (d$strand[1]=='-') observed <- as.character(Biostrings::complement(Biostrings::DNAStringSet(observed)))
  stopifnot(identical(toupper(observed),toupper(d$symbol)),
    identical(unname(unlist(a$source_positions)),d$source_position),
    all(lengths(a$source_positions)==3L))
  for (j in seq_len(nrow(a))) {
    projected <- project_ggchord_features(data.frame(feature_id='codon',parent_id=id,start=j,end=j),cds)
    pieces <- projected$segments[[1]]
    positions <- unlist(lapply(seq_len(nrow(pieces)),function(k) seq.int(pieces$start[k],pieces$end[k])))
    stopifnot(setequal(positions,a$source_positions[[j]]))
  }
}
spliced <- filter(translated$residues,accver %in% e$models$feature_id[e$models$gene %in% c('rps16','petD')])
stopifnot(any(vapply(spliced$source_cds_segments,function(x) length(unique(x))>1L,logical(1))))
# A second independently sourced negative CDS (EMALE) matches the original protein FASTA.
g <- read_gff3('examples/format-inputs/emales-BVI_023A.gff3') |>
  filter(feature_id=='BVI_023A_0009')
emale_source <- read_fasta_lengths('examples/format-inputs/emales-BVI_023A.fna',keep_sequence=TRUE)
emale <- translate_ggchord_cds(build_ggchord_cds(g,genetic_code=11),emale_source,initiation='start')
expected_em <- read_fasta_lengths('examples/format-inputs/protein-pair.faa',keep_sequence=TRUE)$sequence[1]
stopifnot(sub('[*]$','',emale$models$protein_sequence)==expected_em)
# Explicit synthetic boundary controls, not claimed as biological annotation evidence.
syn <- data.frame(feature_id='control',accver='source',start=1,end=9,strand='+',phase=0)
model <- build_ggchord_cds(syn,genetic_code=1)
stopifnot(translate_ggchord_cds(model,c(source='GTGTGATGA'))$models$protein_sequence=='V**',
  translate_ggchord_cds(build_ggchord_cds(syn,genetic_code=2),c(source='GTGTGATGA'))$models$protein_sequence=='VWW',
  translate_ggchord_cds(build_ggchord_cds(syn,genetic_code=11),c(source='GTGTGATGA'),initiation='start')$models$protein_sequence=='M**',
  translate_ggchord_cds(model,c(source='GUGUGAUGA'),alphabet='RNA')$models$protein_sequence=='V**')
fail(translate_ggchord_cds(model,c(source='ATGNNNTAA')),'translation failed')
stopifnot(translate_ggchord_cds(model,c(source='ATGNNNTAA'),ambiguous='solve')$models$protein_sequence=='MX*')
partial <- syn; partial$phase <- 1; partial$end <- 8
fail(build_ggchord_cds(partial,genetic_code=11),'partial CDS')
partial_model <- build_ggchord_cds(partial,genetic_code=11,partial='trim')
p <- translate_ggchord_cds(partial_model,c(source='AATGTAAC'))
stopifnot(p$models$protein_sequence=='M*',p$models$terminal_bases==1,
  identical(as.numeric(p$nucleotides$source_position),as.numeric(2:7)),p$models$translation_status=='observed_partial')
fail(translate_ggchord_cds(partial_model,c(source='AATGTAAC'),initiation='start'),'complete 5-prime')
fail(translate_ggchord_cds(build_ggchord_cds(syn),c(source='ATGTATTAA')),'genetic_code')
fail(translate_ggchord_cds(model,c(other='ATGTATTAA')),'missing source')
fail(translate_ggchord_cds(model,c(source='ATG')),'exceeds source')
fail(translate_ggchord_cds(model,c(source='ATG---TAA')),'ungapped IUPAC')
fail(translate_ggchord_cds(build_ggchord_cds(syn,genetic_code=999),c(source='ATGTATTAA')),'unsupported genetic_code')
exception <- syn; exception$translation_exception <- 'Sec'
fail(build_ggchord_cds(exception,genetic_code=11),'translation_exception')
bad_phase <- cds_input[1,,drop=FALSE]; bad_phase$segments[[1]]$phase <- c(0,0)
fail(build_ggchord_cds(bad_phase),'phase continuity')
# Real gapped 60-column DNA and 40-column AA excerpts; fragment origins are explicit.
dna <- read_aligned_fasta('examples/alignment/opuntia-aligned.fasta','DNA',source_file=TRUE)
aa <- read_aligned_fasta('examples/alignment/protein-aligned.fasta','AA',source_file=TRUE)
windows <- list(dna=filter(dna,accver %in% unique(accver)[c(1,3)],between(alignment_column,321,380)),
  aa=filter(aa,accver %in% unique(accver)[c(2,3)],between(alignment_column,221,260)))
make_plot <- function(d, title) {
  ggplot2::ggplot(d,aes(alignment_column,accver,label=symbol)) +
    geom_residue(aes(fill=symbol),size=2.7,overflow='error',show.legend=FALSE) +
    ggplot2::scale_fill_hue(l=85,c=40,guide='none') +
    ggplot2::labs(x='Original alignment / coding column',y=NULL,title=title) +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(panel.grid=ggplot2::element_blank(),axis.text.y=ggplot2::element_text(size=7))
}
figures <- list(); exports <- list()
for (name in names(windows)) {
  d <- windows[[name]]; original <- d
  stopifnot(any(d$is_gap),any(vapply(split(d$symbol,d$alignment_column),function(x) length(unique(x))>1,logical(1))))
  for (id in unique(d$accver)) {
    row <- filter(d,accver==id,!is_gap)
    stopifnot(all(map_alignment_position(d,id,row$source_position,inverse=TRUE)$alignment_column==row$alignment_column))
  }
  whole_data <- if(name=='dna') dna else aa
  seqs <- whole_data |> filter(!is_gap) |> group_by(accver) |>
    summarise(length=max(source_position),sequence_unit=first(sequence_unit),.groups='drop') |>
    filter(accver %in% unique(d$accver))
  region <- d |> filter(!is_gap) |> group_by(source_accver) |>
    summarise(start=min(source_position),end=max(source_position),.groups='drop') |>
    rename(accver=source_accver)
  macro <- ggchord(seqs)+geom_seq()+coord_collinear()
  plot <- register_ggchord_alignment(make_plot(d,paste(name,'real aligned excerpt')),d)
  # Reverse the native x scale; registered exact endpoints must follow it.
  if(name=='aa') plot <- plot+ggplot2::scale_x_reverse()
  figures[[name]] <- magnify_ggchord(macro,region,plot,center_size=c(.22,.7),detail_size=c(.5,.8))
  exports[[name]] <- export_ggchord_fusion(figures[[name]],width=16,height=5)
  stopifnot(nrow(exports[[name]]$leaders)==4,nrow(exports[[name]]$connections)==0,
    identical(original,d),identical(exports[[name]]$source_data$detail,d))
  # Independent trained-scale/Coord projection using exported final panel bounds.
  built <- ggplot2::ggplot_build(plot)
  for(i in seq_len(nrow(region))) {
    row <- filter(d,source_accver==region$accver[i],!is_gap)
    cols <- row$alignment_column[match(c(region$start[i],region$end[i]),row$source_position)]+c(-.5,.5)
    xy <- built$layout$coord$transform(data.frame(
      x=built$layout$panel_scales_x[[1]]$transform(cols),
      y=built$layout$panel_scales_y[[1]]$map(unique(row$accver))),built$layout$panel_params[[1]])
    panel <- filter(exports[[name]]$panels,panel_id=='detail')
    target <- filter(exports[[name]]$focus_fronts,panel_id=='detail',region_row==i)
    stopifnot(all(abs(target$x-(panel$panel_xmin+xy$x*(panel$panel_xmax-panel$panel_xmin)))<.01),
      all(abs(target$y-(panel$panel_ymin+xy$y*(panel$panel_ymax-panel$panel_ymin)))<.01))
  }
  ggplot2::ggsave(file.path(out,paste0('v0.26.0-m2-',name,'.png')),figures[[name]],width=16,height=5,dpi=130)
}
# Source-verified reverse complement of the original alignment is a direction control.
fasta <- read_fasta_lengths('examples/alignment/opuntia-aligned.fasta',keep_sequence=TRUE)
forward <- gsub('[-.]','',fasta$sequence)
reversed <- vapply(fasta$sequence,function(s) ggchord:::ggchord_alignment_complement(
  paste0(rev(strsplit(s,'',fixed=TRUE)[[1]]),collapse=''),'DNA'),'')
names(reversed) <- fasta$accver
revdata <- as_ggchord_alignment(reversed,'DNA',data.frame(accver=fasta$accver,
  start=nchar(forward),strand='-',sequence=forward)) |>
  filter(accver %in% fasta$accver[1:2],between(alignment_column,321,380))
revplot <- register_ggchord_alignment(make_plot(revdata,'Declared reverse-complement direction control'),revdata)
revseq <- data.frame(accver=fasta$accver[1:2],length=nchar(forward)[1:2])
revregion <- revdata |> filter(!is_gap) |> group_by(source_accver) |>
  summarise(start=min(source_position),end=max(source_position),.groups='drop') |> rename(accver=source_accver)
revfusion <- magnify_ggchord(ggchord(revseq)+geom_seq()+coord_collinear(),revregion,revplot,center_size=c(.22,.7),detail_size=c(.5,.8))
revexport <- export_ggchord_fusion(revfusion,16,5)
stopifnot(all(revdata$source_verified),all(diff(filter(revexport$focus_fronts,panel_id=='detail',region_row==1)$local_x)<0))
# Real spliced negative CDS window: 60 coding bases, with translations on codon centres.
id <- e$models$feature_id[e$models$gene=='rps16']
nt <- filter(translated$nucleotides,accver==id,alignment_column<=60)
res <- filter(translated$residues,accver==id,alignment_column<=20)
cds_plot <- make_plot(nt,'rps16: spliced negative CDS, 60 bp / 20 aa') +
  ggplot2::geom_text(data=res,aes(x=alignment_column*3-1,y='translation',label=symbol),
    inherit.aes=FALSE,family='mono',size=2.7) +
  ggplot2::scale_y_discrete(limits=c(id,'translation'))
cds_plot <- register_ggchord_alignment(cds_plot,nt)
regions <- nt |> group_by(source_cds_segment) |>
  summarise(accver=first(source_accver),start=min(source_position),end=max(source_position),.groups='drop')
whole <- ggchord(source |> select(accver,length))+geom_seq()+coord_circular()
cds_fusion <- magnify_ggchord(whole,regions,cds_plot,center_size=c(.18,.7),detail_size=c(.52,.8))
cds_export <- export_ggchord_fusion(cds_fusion,16,5)
stopifnot(nrow(cds_export$leaders)==2*nrow(regions),nrow(cds_export$focus_regions)==2)
# Intronic envelope is not fully represented in this spliced character view.
fail(export_ggchord_fusion(magnify_ggchord(whole,data.frame(accver=source$accver,
  start=min(nt$source_position),end=max(nt$source_position)),cds_plot),13,5),'no fully visible')
ggplot2::ggsave(file.path(out,'v0.26.0-m2-cds.png'),cds_fusion,width=16,height=5,dpi=130)
# Missing cells and modified registered positions fail explicitly.
d <- windows$dna
bad <- make_plot(d,'bad')+ggplot2::aes(x=alignment_column+1)
fail(register_ggchord_alignment(bad,d),'no longer match')
fail(register_ggchord_alignment(make_plot(d,'bad')+ggplot2::scale_x_log10(),d),'identity/reverse')
# Multiple alignment rows sharing one source require explicit detail_instance.
small <- as_ggchord_alignment(c(one='ACG',two='ACG'),'DNA',data.frame(
  accver=c('one','two'),source_accver='shared',start=1,strand='+',sequence='ACG'))
sp <- register_ggchord_alignment(make_plot(small,'shared'),small)
main <- ggchord(data.frame(accver='shared',length=3))+geom_seq()+coord_collinear()
fail(export_ggchord_fusion(magnify_ggchord(main,data.frame(accver='shared',start=1,end=3),sp)),'ambiguous')
selected <- export_ggchord_fusion(magnify_ggchord(main,data.frame(accver='shared',start=1,end=3,detail_instance='two'),sp))
stopifnot(selected$focus_regions$detail_instance=='two')
# Draw-time physical overflow preserves the existing explicit error or tile policy.
close <- ggchord:::ggchord_measurement_device(width=2,height=2)
err <- tryCatch({grid::grid.newpage();grid::grid.draw(ggplot2::ggplotGrob(make_plot(d,'narrow')));NULL},error=identity)
close()
stopifnot(inherits(err,'error'),grepl('cannot fit',conditionMessage(err)))
tile <- ggplot2::ggplot(d,aes(alignment_column,accver,label=symbol))+geom_residue(overflow='tile')
ggplot2::ggsave(tempfile(fileext='.png'),tile,width=2,height=2)
local({
a <- as_ggchord_alignment(c(A='AC-G',B='ACTG'),'DNA')
p <- register_ggchord_alignment(ggplot2::ggplot(a,aes(alignment_column,accver,label=symbol))+geom_residue(),a)
w <- ggchord(data.frame(accver=c('A','B'),length=4))+geom_seq()+coord_collinear()
r <- data.frame(accver='A',start=1,end=3,detail_panel=c('letters','skeleton'))
f <- magnify_ggchord(w,r,list(letters=p,skeleton=w),side=c(letters='top',skeleton='bottom'))
e <- export_ggchord_fusion(f)
stopifnot(nrow(e$leaders)==4,length(e$source_data)==3)
f2 <- ggchord_fusion(w)+fusion_panel(p,'letters')+fusion_ribbon(
  data.frame(qaccver='A',saccver='A',qstart=1,qend=3,sstart=1,send=3),
  qpanel='center',spanel='letters',fill='steelblue')
e2 <- export_ggchord_fusion(f2)
stopifnot(nrow(e2$connections)==1,all(e2$connections$retained))
# Existing Collinear-only path remains available.
e3 <- export_ggchord_fusion(magnify_ggchord(w,r[1,1:3],w))
stopifnot(nrow(e3$leaders)==2)
cat('PASS mixed character/Collinear details, explicit-fill Fusion link, old single-detail path\n')
})
jsonlite::write_json(list(status='passed direct M2 checks; not release acceptance',
  real_cds=e$models[,c('gene','feature_id','strand')],codons=nrow(translated$residues),
  source_translation_matches=4,spliced_codon_count=sum(vapply(spliced$source_cds_segments,
    function(x) length(unique(x))>1L,logical(1))),
  alignment_windows=lapply(windows,function(d) list(cells=nrow(d),gaps=sum(d$is_gap),columns=range(d$alignment_column))),
  focus_regions=lapply(exports,`[[`,'focus_regions'),cds_regions=cds_export$focus_regions,
  reverse_regions=revexport$focus_regions,biostrings_version=translated$provenance$engine_version),
  file.path(out,'v0.26.0-m2-results.json'),auto_unbox=TRUE,pretty=TRUE,na='null')
cat('PASS: real DNA/AA, exact source anchors, source-backed CDS/phase/translation and failures\n')
