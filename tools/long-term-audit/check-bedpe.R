# Direct v0.20.0 importer checks; no package-wide test/check pipeline.
pkgload::load_all(quiet = TRUE)
hits <- read_paf("examples/format-inputs/emales-pair.paf")
path <- tempfile(fileext = ".bedpe")
six <- tempfile(fileext = ".bedpe")
unknown <- tempfile(fileext = ".bedpe")
bad <- tempfile(fileext = ".bedpe")
# This is an explicit format re-encoding of observed HSPs. Query '+' is the
# chosen alignment reference axis, not a gene strand or paired-read claim.
rows <- data.frame(chrom1=hits$qaccver,start1=hits$qstart-1,end1=hits$qend,
  chrom2=hits$saccver,start2=pmin(hits$sstart,hits$send)-1,end2=pmax(hits$sstart,hits$send),
  name=paste0("paf-line:",hits$.source_line),score=".",strand1="+",strand2=hits$relative_strand,
  identity_evidence=as.character(hits$pident),record=as.character(hits$.source_line))
writeLines(c("# Derived from original EMALE PAF; not paired-end sequencing evidence",
  "track name=coordinate_roundtrip", ""), path)
utils::write.table(rows,path,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
imported <- read_bedpe(path,source_file=TRUE)
cols <- c("qaccver","saccver","qstart","qend","sstart","send")
stopifnot(nrow(imported)==nrow(hits), all(imported$q_mapped & imported$s_mapped),
  isTRUE(all.equal(imported[cols],hits[cols],check.attributes=FALSE)),
  identical(imported$.source_line,seq_len(nrow(hits))+3L),
  all(imported$.source_file==path),all(imported$score=="."),
  !any(c("pident","length","qlen","slen","source_hsp_id") %in% names(imported)),
  imported$bedpe_extra[[1]][["column_11"]]==as.character(hits$pident[1]),
  attr(imported,"ggchord_import")$format=="BEDPE")
# One representative full-interval build. Explicit fill does not invent Identity.
seqs <- unique(rbind(data.frame(accver=hits$qaccver,length=hits$qlen),
  data.frame(accver=hits$saccver,length=hits$slen)))
p <- ggchord(seqs,imported,validate="none")+geom_seq()+
  geom_link_ribbon(fill="steelblue")+coord_collinear()
invisible(ggplot2::ggplotGrob(p))
# Synthetic format cases: one-position six-column records, unknown anchors,
# arbitrary score strings, and an empty trailing custom field.
writeLines("track\t0\t1\tbrowser\t4\t5",six)
a <- read_bedpe(six)
stopifnot(a$qaccver=="track",a$saccver=="browser",a$qstart==1,a$qend==1,
  a$qstrand==".",!a$orientation_known,is.na(a$direction))
combined <- read_bedpe(files=c(path,six))
stopifnot(nrow(combined)==12L,length(combined$bedpe_extra[[12]])==0L)
writeLines("A\t0\t1\t.\t-1\t-1\tmissing\tp=1e-5\t+\t.\tmeta\t",unknown)
expect_error <- function(expr,pattern) {
  error <- tryCatch({force(expr);NULL},error=identity)
  stopifnot(inherits(error,"error"),grepl(pattern,conditionMessage(error)))
}
expect_error(read_bedpe(unknown),"incomplete anchor")
u <- read_bedpe(unknown,unmapped="keep")
stopifnot(u$q_mapped,!u$s_mapped,is.na(u$saccver),is.na(u$sstart),
  !u$orientation_known,is.na(u$direction),u$score=="p=1e-5",
  identical(unname(u$bedpe_extra[[1]]),c("meta","")))
writeLines("A\t0\t0\tB\t1\t2",bad)
expect_error(read_bedpe(bad),"end1|zero-width")
writeLines(c("A\t0\t1\tB\t1\t2","A\t0\t1\tB\t1"),bad)
expect_error(read_bedpe(bad),"consistent field count")
jsonlite::write_json(list(version=as.character(utils::packageVersion("ggchord")),
  source="examples/format-inputs/emales-pair.paf",source_rows=nrow(hits),
  canonical_endpoint_roundtrip=TRUE,physical_lines_preserved=TRUE,
  extra_fields_preserved=TRUE,unknown_anchor_policy="error by default; explicit audit keep",
  representative_build="Collinear full-interval ribbon with explicit constant fill",
  evidence_scope="Observed HSP re-encoding plus explicitly synthetic format cases; no real paired-read claim",
  specification="https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format"),
  "tools/long-term-audit/v0.20.0-bedpe-results.json",auto_unbox=TRUE,pretty=TRUE)
unlink(c(path,six,unknown,bad))
message("BEDPE: 11 source HSPs round-trip; raw fields, unknowns, provenance and representative build passed")
