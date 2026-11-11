# Local development evidence only; no package-wide test/check or device matrix.
pkgload::load_all(quiet = TRUE)
library(dplyr)
source("tools/long-term-audit/semantic-prototypes.R")
out <- "tools/long-term-audit"
started <- proc.time()[["elapsed"]]
expect_error <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}
paths <- c("examples/fasta/MT108731.1.fna", "examples/fasta/MT118296.1.fna")
dna <- read_fasta_lengths(files = paths, sequence_unit = "bp", keep_sequence = TRUE, source_file = TRUE)
aa <- read_fasta_lengths("examples/protein/four-proteins.faa", sequence_unit = "aa", keep_sequence = TRUE)
plain <- read_fasta_lengths(paths[1])
stopifnot(!"sequence" %in% names(plain), all(nchar(dna$sequence) == dna$length),
  all(nchar(aa$sequence) == aa$length), all(file.exists(dna$.source_file)))
letters <- bind_rows(lapply(dna$accver, function(id) source_letters(dna, id, 101, 160, "DNA")),
  lapply(aa$accver[1:2], function(id) source_letters(aa, id, 21, 60, "AA")))
for (id in dna$accver) {
  forward <- source_letters(dna, id, 101, 160, "DNA")
  backward <- source_letters(dna, id, 101, 160, "DNA", reverse = TRUE)
  stopifnot(identical(backward$source_position, rev(forward$source_position)),
    identical(backward$symbol, strsplit(chartr("ACGTacgt", "TGCAtgca", paste(rev(forward$symbol), collapse="")), "")[[1]]))
}
# Explicitly synthetic alphabet edge, not new biological evidence.
synthetic <- data.frame(accver="alphabet", length=15, sequence="aCgTrYsWkMbDhVn")
ambiguous <- source_letters(synthetic,"alphabet",1,15,"DNA",reverse=TRUE)
stopifnot(paste(ambiguous$symbol,collapse="")=="nBdHvKmWsRyAcGt")
expect_error(source_letters(plain, plain$accver[1], 1, 10), "actual sequence")
expect_error(source_letters(dna, dna$accver[1], 1, 10, "RNA"), "alphabet")
local <- focus_ggchord_data(dna, loci=data.frame(accver=dna$accver,start=101,end=160))
again <- focus_ggchord_data(local$seq_data, loci=data.frame(accver=dna$accver,start=11,end=30))
stopifnot(identical(again$seq_data$sequence, substr(dna$sequence,111,130)),
  all(again$seq_data$.source_origin == 111),
  isTRUE(all.equal(source_letters(again$seq_data,dna$accver[1],1,20)$source_position,111:130)))
expect_error(focus_ggchord_data(transform(dna, sequence="A"),
  loci=data.frame(accver=dna$accver,start=1,end=10)), "match input lengths")
expect_error(focus_ggchord_data(dna,
  loci=data.frame(accver=dna$accver,start=1.5,end=10)), "integer locus")
# Actual measured per-position width; do not shrink characters indefinitely.
close <- ggchord:::ggchord_measurement_device(width=10,height=6)
wide <- letter_budget(letters[letters$accver==dna$accver[1],], 180)
narrow <- letter_budget(letters[letters$accver==dna$accver[1],], 12)
close()
stopifnot(wide$representation=="letters",narrow$representation=="tiles")
utils::write.csv(letters,file.path(out,"v0.20.0-source-letters.csv"),row.names=FALSE)
# M0: true pETDuet-1 annotations, 1800 bp, independent linear re-layout.
# Origin-crossing feature/primer spans [5404,5420]+[1,2/3] do not overlap
# this 101:1900 window; exclude explicitly before the linear focus projection.
plasmid <- plasmid_example_pETDuet_1
plasmid_file <- "examples/plasmid/pETDuet-1.fna"
stopifnot(identical(read_fasta_lengths(plasmid_file,keep_sequence=TRUE)$sequence,plasmid$sequence))
features <- find_common_features(plasmid) |> filter(type != "primer_bind")
primers <- find_primer_bindings(plasmid,set="universal")
sites <- find_restriction_sites(plasmid) |>
  filter_restriction_sites(enzymes=c("EcoRI","BamHI","HindIII"))
region <- data.frame(accver=plasmid$accver,start=101,end=1900)
selected <- focus_ggchord_data(plasmid,feature_data=filter(features,start<=end),loci=region)
local_primers <- focus_ggchord_data(plasmid,feature_data=filter(primers,start<=end),loci=region)$feature_data
local_sites <- focus_ggchord_data(plasmid,feature_data=sites,loci=region)$feature_data
stopifnot(nchar(selected$seq_data$sequence)==1800,
  identical(selected$seq_data$sequence,substr(plasmid$sequence,101,1900)),
  nrow(selected$feature_data)>1,nrow(local_primers)>0,nrow(local_sites)>0,
  all(selected$feature_data$start>=1),all(selected$feature_data$end<=1800),
  all(selected$feature_data$start+100==pmax(selected$feature_data$.source_start,101)),
  all(selected$feature_data$end+100==pmin(selected$feature_data$.source_end,1900)))
whole <- ggchord(plasmid,feature_data=features) + geom_seq(show.legend=FALSE) +
  geom_seq_region(data=region) +
  geom_feature(position=position_feature_stack(base_position=position_plasmid()),show.legend=FALSE) +
  geom_seq_center_label() + coord_circular(rotation=-30) +
  ggplot2::theme(plot.background=ggplot2::element_rect(fill=NA,colour=NA),
    panel.background=ggplot2::element_rect(fill=NA,colour=NA))
# Three nearby recognition intervals retain three geometry records and share
# one explicit display label; no site or cut evidence is removed or invented.
site_label <- local_sites |>
  summarise(accver=first(accver),start=min(start),end=max(end),strand=".",
    anno=paste(unique(enzyme),collapse=" / "))
detail_labels <- selected$feature_data |> filter(type %in% c("CDS","rep_origin"),anno!="ATG")
detail <- ggchord(selected$seq_data,feature_data=selected$feature_data) +
  geom_seq(seq_orientation=-1,show.legend=FALSE) +
  geom_feature(position=position_feature_stack(base_position=position_strand(.1)),show.legend=FALSE) +
  geom_feature_label(data=detail_labels,label_orientation="horizontal",label_overlap="allow",colour="#202020",size=3,show.legend=FALSE) +
  geom_feature(data=local_primers,shape="arrow",colour="#7C3A94",width=.035,
    position=position_strand(.22),show.legend=FALSE) +
  geom_feature_label(ggplot2::aes(feature_label=name),data=local_primers,label_orientation="horizontal",label_overlap="allow",colour="#202020",
    position=position_strand(.22),size=3,show.legend=FALSE) +
  geom_feature(data=local_sites,shape="block",colour="#B45309",width=.035,
    position=position_strand(.34),show.legend=FALSE) +
  geom_feature_label(ggplot2::aes(feature_label=anno),data=site_label,label_orientation="horizontal",label_overlap="allow",colour="#202020",
    label_radial_offset=.3,
    position=position_strand(.34),size=3,show.legend=FALSE) + coord_collinear(seq_axis="vertical") +
  scale_seq_position_continuous(coordinates="source",breaks=c(101,500,1000,1500,1900)) +
  ggplot2::theme(plot.background=ggplot2::element_rect(fill=NA,colour=NA),
    panel.background=ggplot2::element_rect(fill=NA,colour=NA))
composition <- ggchord_fusion(whole,center_size=c(.4,.85)) +
  fusion_panel(detail,"detail","right",size=c(.42,.88))
e <- export_ggchord_fusion(composition,12,8)
arc <- prototype_front(e,"center",101,1900)
line <- prototype_front(e,"detail",101,1900)
stopifnot(nrow(e$connections)==0,identical(arc$source_boundary[c(1,80)],c(100,1900)),
  identical(line$source_boundary[c(1,80)],c(100,1900)),all(is.finite(c(arc$x,arc$y,line$x,line$y))))
# Same source window rebuilt at a second physical size, no additional file.
e2 <- export_ggchord_fusion(composition,10,7)
stopifnot(!isTRUE(all.equal(prototype_front(e2,"center",101,1900)$x,arc$x)))
prototype <- grid::gTree(composition=composition,start=101,end=1900,cl="ggchord_m0_prototype")
ggplot2::ggsave(file.path(out,"v0.20.0-m0-detail.png"),plot=prototype,width=12,height=8,dpi=120,bg="white")
utils::write.csv(select(selected$feature_data,accver,start,end,strand,type,anno,starts_with(".source")),
  file.path(out,"v0.20.0-m0-features.csv"),row.names=FALSE)
annotation_rows <- bind_rows(
  transmute(selected$feature_data,kind="feature",accver,start,end,strand,label=anno,.source_start,.source_end),
  transmute(local_primers,kind="primer",accver,start,end,strand,label=name,.source_start,.source_end),
  transmute(local_sites,kind="recognition",accver,start,end,strand,label=enzyme,.source_start,.source_end))
utils::write.csv(annotation_rows,file.path(out,"v0.20.0-m0-annotations.csv"),row.names=FALSE)
utils::write.csv(bind_rows(mutate(arc,panel="center"),mutate(line,panel="detail")),
  file.path(out,"v0.20.0-m0-fronts.csv"),row.names=FALSE)
jsonlite::write_json(list(version="0.20.0.9000",elapsed_seconds=unname(proc.time()[["elapsed"]]-started),source_md5=as.list(tools::md5sum(c(paths,"examples/protein/four-proteins.faa",plasmid_file))),
  character_rows=nrow(letters),dna_examples=2,aa_examples=2,source_mapping=TRUE,reverse_complement=TRUE,
  repeated_focus=TRUE,wide=wide,narrow=narrow,
  M0=list(source=plasmid$accver,start=101,end=1900,features=nrow(selected$feature_data),primers=nrow(local_primers),
    restriction_sites=nrow(local_sites),site_label_groups=1,statistics="none: no recomputed or implied signal",
    feature_labels=nrow(detail_labels),label_selection="CDS excluding literal ATG and replication origins; primer names and selected enzyme names",
    primer_display="generic arrow at real binding interval",restriction_display="three original recognition blocks; one shared enzyme-name label; not cut-site ticks",
    feature_source_overlap_verified=TRUE,
    source_boundaries=c(100,1900),decorative_leaders=2,biological_connections=nrow(e$connections),
    coordinates="identity position scales only; shared continuous source frame and measured panel context",
    sizes_inches=list(c(12,8),c(10,7)),PNG="v0.20.0-m0-detail.png"),
  not_verified=c("PDF/SVG release matrix","true alignment columns/mismatch/indel","CDS translation/phase detail", "M1-M3","arbitrary inset positions/scales")),
  file.path(out,"v0.20.0-semantic-results.json"),pretty=TRUE,auto_unbox=TRUE)
message("Semantic detail: true DNA/AA source mapping, repeated focus, width budget and M0 passed")
