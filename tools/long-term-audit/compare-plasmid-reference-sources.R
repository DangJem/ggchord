# Read-only source comparison; do not source the database generator itself.
source('tools/audit-source-snapshot.R')
audit <- audit_source_snapshot(additional=c('data-raw/generate_common_feature_database.R',
 'tools/long-term-audit/compare-plasmid-reference-sources.R'))
pkgload::load_all(quiet=TRUE)
reader <- new.env(parent=globalenv())
for(expr in parse('data-raw/generate_common_feature_database.R')) {
 if(is.call(expr) && identical(expr[[1]],as.name('<-')) && is.symbol(expr[[2]]) &&
   as.character(expr[[2]]) %in% c('read_uint32_be','value_or','read_binary_plasmid_reference')) eval(expr,reader)
}
files <- list.files('examples/plasmid',pattern='[.]dna$',full.names=TRUE)
rows <- list(); details <- list()
for(file in files) {
 source <- reader$read_binary_plasmid_reference(file)
 label <- source$sequence$label
 seq <- data.frame(accver=label,length=nchar(source$sequence$sequence),sequence=source$sequence$sequence)
 f <- find_common_features(seq)
 s <- find_restriction_sites(seq)
 selected <- filter_restriction_sites(s,set='reference',parent_set='commercial_nonredundant')
 target <- unique(vapply(seq$length*c(1,3,5)/6,function(x)which.min(abs(s$position-x)),integer(1)))
 # Recover source segment attributes independently: the generator currently
 # discards some of these fields, which is itself part of this comparison.
 bytes <- readBin(file,'raw',n=file.info(file)$size); pos <- 1L; xml <- NULL
 while(pos+4L<=length(bytes)) {
  type<-as.integer(bytes[pos]);n<-reader$read_uint32_be(bytes[(pos+1L):(pos+4L)])
  if(type==10L)xml<-rawToChar(bytes[(pos+5L):(pos+4L+n)])
  pos<-pos+5L+n
 }
 doc <- xml2::read_xml(xml)
 named <- lapply(xml2::xml_find_all(doc,'.//Feature'),function(node) {
  segments<-xml2::xml_find_all(node,'./Segment')
  attrs<-lapply(segments,xml2::xml_attrs)
  at<-which(vapply(attrs,function(x)!is.na(x['name']) && nzchar(x['name']),logical(1)))
  if(!length(at))return(NULL)
  data.frame(parent=xml2::xml_attr(node,'name'),name=vapply(attrs[at],function(x)unname(x['name']),character(1)),
   range=vapply(attrs[at],function(x)unname(x['range']),character(1)))
 })
 named <- ggchord_rbind_fill(named) %||% data.frame()
 ids<-source$features$common_feature_id
 missing<-source$features[!ids %in% f$common_feature_id,]
 extra<-f[!f$common_feature_id %in% ids,]
 coordinate_changes<-lapply(seq_len(nrow(source$features)),function(i) {
  original<-source$features[i,];at<-match(original$common_feature_id,f$common_feature_id)
  if(is.na(at))return(NULL)
  detected<-f[at,]
  if(any(c(original$start!=detected$start,original$end!=detected$end,original$strand!=detected$strand,
    original$name!=detected$anno,original$type!=detected$type)))return(original$common_feature_id)
  NULL
 })
 rows[[label]]<-data.frame(plasmid=label,length=seq$length,source_features=nrow(source$features),
  recognized_features=nrow(f),missing_features=nrow(missing),extra_features=nrow(extra),
  changed_identity_coordinates=length(Filter(Negate(is.null),coordinate_changes)),
  source_primer_bindings=nrow(source$primers),named_segments=nrow(named),
  all_enzyme_site_rows=nrow(s),all_cut_positions=length(unique(s$position)),
  old_sampled_site_rows=length(target),reference_profile=source$enzyme_profile$set_name,
  reference_selected_site_rows=nrow(selected),reference_selected_positions=length(unique(selected$position)))
 details[[label]]<-list(source_path=file,source_sha256=source$sequence$source_sha256,
  image_path=sub('[.]dna$','.png',file),image_sha256=digest::digest(sub('[.]dna$','.png',file),algo='sha256',file=TRUE),
  source_features=source$features,recognized_features=f,missing=missing,extra=extra,
  source_primers=source$primers[c('name','start','end','strand')],named_segments=named,
  source_profile=source$enzyme_profile,selected_sites=selected,
  observed_methylation=unique(s$methylation_status))
 cat(label,nrow(source$features),nrow(f),'reference sites',nrow(selected),'\n')
}
summary<-ggchord_rbind_fill(rows)
write.csv(summary,'tools/long-term-audit/v0.27.0-reference-source-summary.csv',row.names=FALSE)
audit_assert_source_snapshot(audit)
jsonlite::write_json(list(summary=summary,cases=details,runtime_sha256=as.list(audit$engine),
 input_sha256=as.list(audit$inputs),scope='Original 13 .dna features/primers/profiles versus complete current recognition; no PNG label OCR or exact export equivalence claimed'),
 'tools/long-term-audit/v0.27.0-reference-source-results.json',auto_unbox=TRUE,pretty=TRUE,na='null')
