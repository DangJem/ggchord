# Restriction pattern parsing, searching, filtering and rendering.

ggchord_iupac_regex <- function(motif) {
  code <- c(A="A", C="C", G="G", T="T", R="[AG]", Y="[CT]",
    S="[GC]", W="[AT]", K="[GT]", M="[AC]", B="[CGT]",
    D="[AGT]", H="[ACT]", V="[ACG]", N="[ACGT]")
  chars <- strsplit(toupper(motif), "", fixed = TRUE)[[1L]]
  if (any(!chars %in% names(code))) ggchord_stop("Restriction motifs may contain only IUPAC DNA symbols")
  paste0(unname(code[chars]), collapse = "")
}

ggchord_reverse_complement <- function(sequence) {
  code <- c(A="T", C="G", G="C", T="A", R="Y", Y="R", S="S", W="W",
            K="M", M="K", B="V", D="H", H="D", V="B", N="N")
  chars <- strsplit(toupper(sequence), "", fixed = TRUE)[[1L]]
  if (any(!chars %in% names(code))) ggchord_stop("DNA text may contain only IUPAC DNA symbols")
  paste0(rev(unname(code[chars])), collapse = "")
}

ggchord_restriction_match_starts <- function(search, query, max_start) {
  m <- nchar(query)
  runs <- gregexpr("[ACGT]+", query, perl=TRUE)[[1L]]
  run_lengths <- attr(runs,"match.length")
  if (length(runs)==1L && runs[1L]<0L) {
    hits <- gregexpr(paste0("(?=",ggchord_iupac_regex(query),")"),search,perl=TRUE)[[1L]]
    if (length(hits)==1L && hits[1L]<0L) return(integer())
    return(hits[hits<=max_start])
  }
  chosen <- which.max(run_lengths)
  anchor_start <- runs[chosen]; anchor <- substr(
    query,anchor_start,anchor_start+run_lengths[chosen]-1L
  )
  hits <- gregexpr(anchor,search,fixed=TRUE)[[1L]]
  if (length(hits)==1L && hits[1L]<0L) return(integer())
  starts <- unique(hits-anchor_start+1L)
  starts <- starts[starts>=1L & starts<=max_start]
  if (!length(starts)) return(integer())
  candidates <- substring(search,starts,starts+m-1L)
  starts[grepl(paste0("^",ggchord_iupac_regex(query),"$"),candidates,perl=TRUE)]
}

ggchord_parse_rebase <- function(path) {
  required <- c("VERSION", "embossa_e.txt", "embossa_r.txt", "embossa_s.txt")
  missing <- required[!file.exists(file.path(path, required))]
  if (length(missing)) ggchord_stop("REBASE directory is missing: ", paste(missing, collapse = ", "))
  version_text <- trimws(readLines(file.path(path, "VERSION"), warn = FALSE)[1L])
  version <- sub("^.*?([0-9]+).*$", "\\1", version_text)
  if (!nzchar(version)) ggchord_stop("Cannot parse REBASE VERSION")

  enzyme_lines <- readLines(file.path(path, "embossa_e.txt"), warn = FALSE)
  source_rows <- which(nzchar(trimws(enzyme_lines)) & !grepl("^#", enzyme_lines))
  fields <- strsplit(trimws(enzyme_lines[source_rows]), "[[:space:]]+")
  if (any(lengths(fields) != 9L)) ggchord_stop("Invalid embossa_e.txt row width")
  mat <- do.call(rbind, fields)
  defs <- data.frame(
    pattern_id = sprintf("rebase%s:e:%06d", version, source_rows),
    pattern_source_row = source_rows,
    enzyme = mat[,1L], motif = toupper(mat[,2L]),
    motif_length = as.integer(mat[,3L]), ncuts = as.integer(mat[,4L]),
    blunt = as.logical(as.integer(mat[,5L])),
    cut_offset_1 = as.integer(mat[,6L]), cut_offset_2 = as.integer(mat[,7L]),
    cut_offset_3 = as.integer(mat[,8L]), cut_offset_4 = as.integer(mat[,9L]),
    stringsAsFactors = FALSE
  )
  if (any(nchar(defs$motif) != defs$motif_length) ||
      any(!defs$ncuts %in% c(0L,1L,2L,4L))) {
    ggchord_stop("Invalid REBASE pattern lengths or cut counts")
  }

  supplier_lines <- readLines(file.path(path, "embossa_s.txt"), warn = FALSE)
  supplier_lines <- supplier_lines[nzchar(trimws(supplier_lines)) & !grepl("^#", supplier_lines)]
  supplier_code <- sub("[[:space:]].*$", "", supplier_lines)
  supplier_name <- sub("^[^[:space:]]+[[:space:]]+", "", supplier_lines)
  supplier_map <- stats::setNames(supplier_name, supplier_code)

  ref_lines <- readLines(file.path(path, "embossa_r.txt"), warn = FALSE)
  ref_lines <- ref_lines[!grepl("^#", ref_lines)]
  while (length(ref_lines) && !nzchar(trimws(ref_lines[1L]))) ref_lines <- ref_lines[-1L]
  blocks <- list(); current <- character()
  for (line in ref_lines) {
    if (identical(line, "//")) {
      blocks[[length(blocks)+1L]] <- current; current <- character()
    } else current <- c(current, line)
  }
  info <- lapply(blocks, function(x) {
    codes <- if (length(x) >= 6L) strsplit(trimws(x[6L]), "", fixed=TRUE)[[1L]] else character()
    codes <- codes[nzchar(codes)]
    data.frame(
      enzyme=x[1L], organism=if (length(x)>=2L) x[2L] else "",
      supplier_codes=paste(codes, collapse=""),
      suppliers=paste(unname(supplier_map[codes]), collapse="; "),
      commercial=length(codes)>0L, stringsAsFactors=FALSE
    )
  })
  info <- if (length(info)) do.call(rbind, info) else data.frame()
  matched <- match(defs$enzyme, info$enzyme)
  defs$organism <- info$organism[matched]
  defs$supplier_codes <- info$supplier_codes[matched]
  defs$suppliers <- info$suppliers[matched]
  defs$commercial <- info$commercial[matched]
  defs$commercial[is.na(defs$commercial)] <- FALSE
  defs$database_version <- version
  defs$source <- paste0("REBASE ", version)
  rownames(defs) <- NULL
  defs
}

ggchord_builtin_rebase <- function() {
  if (exists("ggchord_rebase_database", inherits=TRUE)) {
    return(get("ggchord_rebase_database", inherits=TRUE))
  }
  configured <- getOption("ggchord.rebase.path", NULL)
  candidates <- unique(c(configured, file.path(getwd(), "examples", "rebase"),
                         file.path(dirname(getwd()), "examples", "rebase")))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[vapply(candidates, function(x) all(file.exists(file.path(
    x, c("VERSION", "embossa_e.txt", "embossa_r.txt", "embossa_s.txt")))), logical(1))]
  if (length(hit)) return(ggchord_parse_rebase(hit[1L]))
  data.frame(
    pattern_id=paste0("ggchord:fallback:",seq_len(5L)),
    pattern_source_row=seq_len(5L),
    enzyme=c("EcoRI","BamHI","HindIII","PstI","SmaI"),
    motif=c("GAATTC","GGATCC","AAGCTT","CTGCAG","CCCGGG"),
    motif_length=6L,ncuts=2L,blunt=c(FALSE,FALSE,FALSE,FALSE,TRUE),
    cut_offset_1=c(1L,1L,1L,5L,3L),cut_offset_2=c(5L,5L,5L,1L,3L),
    cut_offset_3=0L,cut_offset_4=0L,organism=NA_character_,
    supplier_codes=NA_character_,suppliers=NA_character_,commercial=FALSE,
    database_version="fallback",source="ggchord fallback",stringsAsFactors=FALSE
  )
}

ggchord_normalize_restriction_patterns <- function(patterns) {
  if (is.character(patterns) && length(patterns)) {
    if (is.null(names(patterns)) || anyNA(names(patterns)) ||
        any(!nzchar(names(patterns)))) {
      ggchord_stop("Named custom motifs must use enzyme names")
    }
    patterns <- data.frame(
      enzyme = names(patterns), motif = unname(patterns),
      stringsAsFactors = FALSE
    )
  }
  ggchord_require_columns(patterns,c("enzyme","motif"),"restriction patterns")
  out <- as.data.frame(patterns,stringsAsFactors=FALSE)
  if (!"pattern_id"%in%names(out)) out$pattern_id <- sprintf("custom:%06d",seq_len(nrow(out)))
  if (!"pattern_source_row"%in%names(out)) out$pattern_source_row <- seq_len(nrow(out))
  if (!"motif_length"%in%names(out)) out$motif_length <- nchar(out$motif)
  if (!"ncuts"%in%names(out)) out$ncuts <- 0L
  if (!"blunt"%in%names(out)) out$blunt <- FALSE
  if (all(c("cut_top","cut_bottom")%in%names(out))) {
    out$cut_offset_1 <- out$cut_top; out$cut_offset_2 <- out$cut_bottom
    out$ncuts[out$ncuts==0L] <- 2L
  }
  for(i in 1:4){nm<-paste0("cut_offset_",i);if(!nm%in%names(out))out[[nm]]<-0L}
  defaults <- list(organism=NA_character_,supplier_codes=NA_character_,
    suppliers=NA_character_,commercial=FALSE,database_version="custom",source="custom")
  for(nm in names(defaults))if(!nm%in%names(out))out[[nm]]<-defaults[[nm]]
  out$enzyme<-as.character(out$enzyme);out$motif<-toupper(as.character(out$motif))
  if(anyNA(out$enzyme)||any(!nzchar(out$enzyme))||anyNA(out$motif)||
     any(!nzchar(out$motif))||anyDuplicated(out$pattern_id))
    ggchord_stop("Restriction patterns require non-empty enzymes, motifs, and unique pattern_id values")
  out$pattern_id <- as.character(out$pattern_id)
  out$pattern_source_row <- as.integer(out$pattern_source_row)
  out$motif_length <- as.integer(out$motif_length)
  out$ncuts <- as.integer(out$ncuts)
  out$blunt <- as.logical(out$blunt)
  if (anyNA(out$pattern_source_row) || any(out$pattern_source_row < 1L) ||
      anyNA(out$motif_length) || any(out$motif_length != nchar(out$motif)) ||
      anyNA(out$ncuts) || any(!out$ncuts %in% c(0L, 1L, 2L, 4L)) ||
      anyNA(out$blunt)) {
    ggchord_stop(
      "Restriction patterns require valid source rows, motif lengths, ",
      "ncuts (0, 1, 2, or 4), and blunt flags"
    )
  }
  for (i in 1:4) {
    nm <- paste0("cut_offset_", i)
    out[[nm]] <- as.numeric(out[[nm]])
    if (any(!is.finite(out[[nm]]))) {
      ggchord_stop("Restriction cut offsets must contain finite numbers")
    }
  }
  invisible(lapply(out$motif,ggchord_iupac_regex));out
}

#' Find restriction pattern matches and cleavage coordinates
#'
#' @param sequence DNA text, a named character vector, or a data frame with
#' `accver` and `sequence`.
#' @param enzymes Optional enzyme names. A data frame is a compatibility
#' spelling for `patterns`.
#' @param patterns Optional custom definitions.
#' @param circular Allow recognition and cleavage coordinates to wrap origin.
#' @param database Optional parsed pattern database.
#' @return One row per pattern match without display filtering or deduplication.
#' @export
find_restriction_sites <- function(sequence,enzymes=NULL,patterns=NULL,
                                   circular=TRUE,database=NULL){
  if(is.data.frame(sequence)){
    ggchord_require_columns(sequence,c("accver","sequence"),"find_restriction_sites()")
    ids<-as.character(sequence$accver);seqs<-as.character(sequence$sequence)
  }else if(is.character(sequence)&&length(sequence)){
    seqs<-as.character(sequence);ids<-names(sequence)
    if(is.null(ids))ids<-if(length(seqs)==1L)"sequence" else paste0("sequence_",seq_along(seqs))
    bad<-is.na(ids)|!nzchar(ids);ids[bad]<-paste0("sequence_",which(bad))
  }else ggchord_stop("find_restriction_sites(): invalid sequence input")
  if(anyDuplicated(ids))ggchord_stop("find_restriction_sites(): sequence IDs must be unique")
  seqs<-toupper(gsub("[[:space:]]","",seqs))
  if(any(!grepl("^[ACGTRYSWKMBDHVN]+$",seqs)))ggchord_stop("find_restriction_sites(): invalid DNA symbols")
  if(!is.logical(circular)||length(circular)!=1L||is.na(circular))ggchord_stop("find_restriction_sites(): circular must be TRUE or FALSE")
  if(is.data.frame(enzymes)){if(!is.null(patterns))ggchord_stop("Supply custom definitions once");patterns<-enzymes;enzymes<-NULL}
  defs<-if(!is.null(patterns))patterns else database%||%ggchord_builtin_rebase()
  defs<-ggchord_normalize_restriction_patterns(defs)
  if(!is.null(enzymes)){
    if(!is.character(enzymes)||anyNA(enzymes))ggchord_stop("find_restriction_sites(): enzymes must be character")
    unknown<-setdiff(enzymes,unique(defs$enzyme));if(length(unknown))ggchord_stop("Unknown restriction enzyme(s): ",paste(unknown,collapse=", "))
    defs<-defs[defs$enzyme%in%enzymes,,drop=FALSE]
    defs<-defs[order(match(defs$enzyme,enzymes),seq_len(nrow(defs))),,drop=FALSE]
  }
  empty_template <- function(){
    out<-data.frame(accver=character(),match_id=character(),pattern_id=character(),pattern_source_row=integer(),enzyme=character(),motif=character(),motif_length=integer(),recognition_sequence=character(),start=integer(),end=integer(),crosses_origin=logical(),strand=character(),ncuts=integer(),blunt=logical(),end_type=character(),position=numeric(),display_position=numeric(),anchor_kind=character(),organism=character(),supplier_codes=character(),suppliers=character(),commercial=logical(),database_version=character(),source=character(),stringsAsFactors=FALSE)
    for(k in 1:4){out[[paste0("cut_offset_",k)]]<-numeric();out[[paste0("cut_",k,"_unwrapped")]]<-numeric();out[[paste0("cut_",k)]]<-numeric()};out
  }
  rows<-list();match_number<-0L
  for(s in seq_along(seqs)){
    target<-seqs[s];len<-nchar(target);match_cache<-new.env(hash=TRUE,parent=emptyenv())
    for(d in seq_len(nrow(defs))){
    motif<-defs$motif[d];m<-nchar(motif)
    search<-if(circular&&m>1L)paste0(target,substr(target,1L,m-1L))else target
    rc<-ggchord_reverse_complement(motif);strands<-if(identical(rc,motif))"+" else c("+","-")
    for(site_strand in strands){
      query<-if(site_strand=="+")motif else rc
      cache_key<-paste0(query,"\r",m)
      if(exists(cache_key,envir=match_cache,inherits=FALSE))starts<-get(cache_key,envir=match_cache,inherits=FALSE) else {starts<-ggchord_restriction_match_starts(search,query,len);assign(cache_key,starts,envir=match_cache)}
      if(length(starts)){
        count<-length(starts);numbers<-seq.int(match_number+1L,match_number+count);match_number<-match_number+count
        offsets<-as.numeric(unlist(defs[d,paste0("cut_offset_",1:4)],use.names=FALSE))
        known<-defs$ncuts[d]>0L & offsets!=0
        # REBASE numbers residues as ... -2, -1, 1, 2, ... and defines each
        # cut immediately to the right of that residue. Convert to a
        # zero-based boundary displacement from just before motif base 1.
        boundary_offsets<-ifelse(offsets>0,offsets,offsets+1)
        if(site_strand=="-"){
          boundary_offsets<-m-c(boundary_offsets[2L],boundary_offsets[1L],boundary_offsets[4L],boundary_offsets[3L])
          known<-known[c(2L,1L,4L,3L)]
        }
        unwrapped<-matrix(NA_real_,nrow=count,ncol=4L)
        for(k in which(known))unwrapped[,k]<-starts+boundary_offsets[k]-1
        cuts<-unwrapped
        if(circular)cuts[,known]<-((unwrapped[,known,drop=FALSE]-1)%%len)+1 else cuts[cuts<1|cuts>len]<-NA_real_
        finish_unwrapped<-starts+m-1L;finish<-if(circular)((finish_unwrapped-1L)%%len)+1L else finish_unwrapped
        end_type<-if(defs$ncuts[d]==0L)"unknown" else if(defs$ncuts[d]==1L)"nick" else if(defs$ncuts[d]==4L)"complex" else if(isTRUE(defs$blunt[d]))"blunt" else if(defs$cut_offset_1[d]<defs$cut_offset_2[d])"5_prime_overhang" else "3_prime_overhang"
        known_indices<-which(known);has_cut<-length(known_indices)>0L
        first_known<-if(has_cut)known_indices[1L] else NA_integer_
        position<-if(has_cut)cuts[,first_known] else starts
        row<-data.frame(accver=rep(ids[s],count),match_id=sprintf("%s:%s:%s:%d:%06d",ids[s],defs$pattern_id[d],site_strand,starts,numbers),pattern_id=rep(defs$pattern_id[d],count),pattern_source_row=rep(defs$pattern_source_row[d],count),enzyme=rep(defs$enzyme[d],count),motif=rep(motif,count),motif_length=rep(m,count),recognition_sequence=substring(search,starts,starts+m-1L),start=starts,end=finish,crosses_origin=finish_unwrapped>len,strand=rep(site_strand,count),ncuts=rep(defs$ncuts[d],count),blunt=rep(defs$blunt[d],count),end_type=rep(end_type,count),position=position,display_position=position,anchor_kind=rep(if(has_cut)"cut" else "recognition",count),organism=rep(defs$organism[d],count),supplier_codes=rep(defs$supplier_codes[d],count),suppliers=rep(defs$suppliers[d],count),commercial=rep(defs$commercial[d],count),database_version=rep(defs$database_version[d],count),source=rep(defs$source[d],count),stringsAsFactors=FALSE)
        for(k in 1:4){row[[paste0("cut_offset_",k)]]<-rep(as.numeric(defs[[paste0("cut_offset_",k)]][d]),count);row[[paste0("cut_",k,"_unwrapped")]]<-unwrapped[,k];row[[paste0("cut_",k)]]<-cuts[,k]}
        rows[[length(rows)+1L]]<-row
      }
    }
  }}
  out<-if(length(rows))ggchord_rbind_fill(rows)else empty_template();rownames(out)<-NULL
  attr(out,"enzyme_database_version")<-unique(defs$database_version);out
}

#' Filter restriction sites for display
#' @param sites Result from [find_restriction_sites()].
#' @param set Display preset.
#' @param enzymes,min_site_length,cuts,window,commercial_only Additional filters.
#' @return A row subset whose biological coordinates are unchanged.
#' @export
filter_restriction_sites<-function(sites,set=c("all","unique","unique_dual","six_plus","unique_6plus","commercial"),enzymes=NULL,min_site_length=NULL,cuts=NULL,window=NULL,commercial_only=FALSE){
  if(!is.data.frame(sites))ggchord_stop("filter_restriction_sites(): sites must be a data frame")
  set<-match.arg(set);if(!nrow(sites))return(sites)
  ggchord_require_columns(sites,c("accver","enzyme","motif_length","position"),"filter_restriction_sites()")
  key<-paste(sites$accver,sites$enzyme,sep="\r");site_count<-as.integer(table(key)[key]);keep<-rep(TRUE,nrow(sites))
  if(set=="unique")keep<-keep&site_count==1L
  if(set=="unique_dual")keep<-keep&site_count%in%c(1L,2L)
  if(set=="six_plus")keep<-keep&sites$motif_length>=6L
  if(set=="unique_6plus")keep<-keep&site_count==1L&sites$motif_length>=6L
  if(set=="commercial"){
    if(!"commercial"%in%names(sites))ggchord_stop("filter_restriction_sites(): sites is missing commercial")
    keep<-keep&!is.na(sites$commercial)&sites$commercial
  }
  if(!is.null(enzymes)){
    if(!is.character(enzymes)||anyNA(enzymes)||any(!nzchar(enzymes)))ggchord_stop("filter_restriction_sites(): enzymes must be non-empty character values")
    keep<-keep&sites$enzyme%in%enzymes
  }
  if(!is.null(min_site_length)){
    if(!is.numeric(min_site_length)||length(min_site_length)!=1L||!is.finite(min_site_length)||min_site_length<0)ggchord_stop("filter_restriction_sites(): min_site_length must be one non-negative number")
    keep<-keep&sites$motif_length>=min_site_length
  }
  if(!is.null(cuts)){
    if(!is.numeric(cuts)||anyNA(cuts)||any(!is.finite(cuts))||any(cuts<0)||any(cuts!=as.integer(cuts)))ggchord_stop("filter_restriction_sites(): cuts must contain non-negative integers")
    keep<-keep&site_count%in%as.integer(cuts)
  }
  if(!is.null(window)){
    if(!is.numeric(window)||length(window)!=2L||any(!is.finite(window)))ggchord_stop("filter_restriction_sites(): invalid window")
    keep<-keep&if(window[1L]<=window[2L])sites$position>=window[1L]&sites$position<=window[2L] else sites$position>=window[1L]|sites$position<=window[2L]
  }
  if(!is.logical(commercial_only)||length(commercial_only)!=1L||is.na(commercial_only))ggchord_stop("filter_restriction_sites(): commercial_only must be TRUE or FALSE")
  if(commercial_only){
    if(!"commercial"%in%names(sites))ggchord_stop("filter_restriction_sites(): sites is missing commercial")
    keep<-keep&!is.na(sites$commercial)&sites$commercial
  }
  attrs<-attributes(sites);out<-sites[keep,,drop=FALSE];rownames(out)<-NULL
  for(nm in setdiff(names(attrs),c("names","row.names","class")))attr(out,nm)<-attrs[[nm]];out
}

GeomRestrictionSite<-ggplot2::ggproto("GeomRestrictionSite",ggplot2::Geom,required_aes=c("x","y"),default_aes=ggplot2::aes(xend=NA_real_,yend=NA_real_,label=NA_character_,.component=NA_character_,group=NA_integer_,colour="#B42318",alpha=1,linewidth=.35,linetype=1,size=2.6,angle=0,hjust=.5,vjust=.5,family="",fontface=1,lineheight=1.2),draw_key=ggplot2::draw_key_path,draw_panel=function(data,panel_params,coord,na.rm=FALSE,segment_params=list(),text_params=list()){paths<-data[data$.component=="path",,drop=FALSE];labels<-data[data$.component=="label",,drop=FALSE];for(nm in names(segment_params))if(nm%in%names(paths))paths[[nm]]<-segment_params[[nm]];for(nm in names(text_params))if(nm%in%names(labels))labels[[nm]]<-text_params[[nm]];grobs<-list();if(nrow(paths))grobs[[length(grobs)+1L]]<-ggplot2::GeomPath$draw_panel(paths,panel_params,coord,lineend="round",linejoin="round",na.rm=na.rm);if(nrow(labels))grobs[[length(grobs)+1L]]<-ggplot2::GeomText$draw_panel(labels,panel_params,coord,parse=FALSE,check_overlap=FALSE,na.rm=na.rm);do.call(grid::grobTree,grobs)})

#' Draw restriction sites, labels and deterministic leaders
#' @param mapping,data Standard layer inputs.
#' @param label Draw labels.
#' @param label_style Label enzyme and position, enzyme only, or position only.
#' @param label_side Draw labels outside, inside, or automatically.
#' @param leader Leader routing style.
#' @param min_label_gap Minimum genomic fraction between label slots.
#' @param tick_length,label_offset Local-normal distances.
#' @param colour,linewidth,label_size Fixed appearance.
#' @param position,show.legend,inherit.aes Standard layer arguments.
#' @param ... Additional geom parameters.
#' @return A ggplot2 layer.
#' @export
geom_restriction_site<-function(mapping=NULL,data=NULL,label=TRUE,label_style=c("enzyme_position","enzyme","position"),label_side=c("outside","inside","auto"),leader=c("trunk","elbow","straight"),min_label_gap=.02,tick_length=.045,label_offset=.28,colour="#B42318",linewidth=.35,label_size=2.6,position="identity",show.legend=FALSE,inherit.aes=FALSE,...){
  colour_supplied<-!missing(colour);linewidth_supplied<-!missing(linewidth);label_size_supplied<-!missing(label_size)
  if(is.null(data))ggchord_stop("geom_restriction_site(): supply restriction-site data")
  label_style<-match.arg(label_style);label_side<-match.arg(label_side);leader<-match.arg(leader)
  vals<-c(min_label_gap,tick_length,label_offset);if(!is.numeric(vals)||any(!is.finite(vals))||any(vals<0))ggchord_stop("geom_restriction_site(): invalid gaps or offsets")
  segment_params<-list();text_params<-list()
  if(colour_supplied){segment_params$colour<-colour;text_params$colour<-colour}
  if(linewidth_supplied)segment_params$linewidth<-linewidth
  if(label_size_supplied)text_params$size<-label_size
  lyr<-ggplot2::layer(data=data.frame(x=numeric(),y=numeric()),mapping=ggplot2::aes(x=x,y=y,group=group,label=label,.component=I(.component)),stat="identity",geom=GeomRestrictionSite,position=position,show.legend=show.legend,inherit.aes=inherit.aes,check.aes=FALSE,check.param=FALSE,params=c(list(na.rm=FALSE,colour=colour,linewidth=linewidth,size=label_size,segment_params=segment_params,text_params=text_params),list(...)))
  lyr$ggchord_type<-"restriction_site";lyr$ggchord_theme_components<-c(segment_params="ggchord.restriction.label.segment",text_params="ggchord.restriction.label")
  lyr$ggchord_params<-list(type="restriction_site",label=label,label_style=label_style,label_side=label_side,leader=leader,min_label_gap=min_label_gap,tick_length=tick_length,label_offset=label_offset)
  ggchord_capture_layer_input(lyr,data,mapping,c("accver","position","enzyme"))
}

ggchord_restriction_geometry <- function(data, params, layout, seq_data) {
  ggchord_require_columns(
    data, c("accver", "position", "enzyme"), "geom_restriction_site()"
  )
  lens <- stats::setNames(seq_data$length, seq_data$accver)
  if (any(!data$accver %in% names(lens))) {
    ggchord_stop("geom_restriction_site(): unknown accver")
  }
  if (!is.numeric(data$position) || any(!is.finite(data$position))) {
    ggchord_stop("geom_restriction_site(): position must contain finite numbers")
  }

  output <- list()
  gid <- 0L
  append_path <- function(points, component, source_row = NA_integer_,
                          anchor_position = NA_real_, slot_position = NA_real_,
                          member_rows = "") {
    gid <<- gid + 1L
    points$.component <- "path"
    points$restriction_component <- component
    points$group <- gid
    points$label <- NA_character_
    points$source_row <- source_row
    points$anchor_position <- anchor_position
    points$slot_position <- slot_position
    points$member_rows <- member_rows
    output[[length(output) + 1L]] <<- points
  }

  for (id in unique(as.character(data$accver))) {
    idx <- which(as.character(data$accver) == id)
    pattern_order <- if ("pattern_id" %in% names(data)) {
      as.character(data$pattern_id[idx])
    } else rep("", length(idx))
    source_order <- if ("pattern_source_row" %in% names(data)) {
      data$pattern_source_row[idx]
    } else idx
    idx <- idx[order(data$position[idx], pattern_order, source_order, idx)]
    arc <- layout$seq_arcs[[id]]
    n <- nrow(arc)
    frac <- data$position[idx] / lens[id]

    # Open the circular ordering at its largest empty gap so labels next to
    # genomic origin remain neighbours during deterministic packing.
    if (isTRUE(layout$circular) && length(frac) > 1L) {
      gaps <- diff(c(frac, frac[1L] + 1))
      cut <- which.max(gaps)
      after <- if (cut < length(frac)) seq.int(cut + 1L, length(frac)) else integer()
      ord <- c(after, seq_len(cut))
      idx <- idx[ord]
      frac <- data$position[idx] / lens[id]
      frac[frac < frac[1L]] <- frac[frac < frac[1L]] + 1
    }

    point_at <- function(fraction, offset = 0) {
      fraction <- fraction %% 1
      k <- pmax(1L, pmin(n, round(1 + fraction * (n - 1))))
      kp <- pmax(1L, k - 1L)
      kn <- pmin(n, k + 1L)
      tx <- arc$x[kn] - arc$x[kp]
      ty <- arc$y[kn] - arc$y[kp]
      tangent_length <- sqrt(tx^2 + ty^2)
      tangent_length[tangent_length <= 1e-12] <- 1
      nx <- -ty / tangent_length
      ny <- tx / tangent_length
      flip <- nx * arc$x[k] + ny * arc$y[k] < 0
      nx[flip] <- -nx[flip]
      ny[flip] <- -ny[flip]
      data.frame(x = arc$x[k] + offset * nx, y = arc$y[k] + offset * ny)
    }

    side_sign <- if (params$label_side == "inside") -1 else 1
    gap <- params$min_label_gap
    slots <- frac
    if (length(slots) > 1L) {
      for (j in 2:length(slots)) slots[j] <- max(slots[j], slots[j - 1L] + gap)
      if (max(slots) > 1) slots <- slots - (max(slots) - 1) / 2
    }
    cluster <- cumsum(c(TRUE, diff(frac) > gap))

    for (members in split(seq_along(idx), cluster)) {
      shared <- isTRUE(params$label) && params$leader == "trunk" &&
        length(members) > 1L
      trunk_offset <- side_sign * params$label_offset * .55
      member_text <- paste(idx[members], collapse = ",")
      if (shared) {
        trunk_range <- range(c(frac[members], slots[members]))
        trunk_fractions <- seq(
          trunk_range[1L], trunk_range[2L],
          length.out = max(12L, ceiling(diff(trunk_range) * 360))
        )
        append_path(
          point_at(trunk_fractions, trunk_offset), "trunk",
          member_rows = member_text
        )
      }

      for (member in members) {
        row <- idx[member]
        biological_position <- data$position[row]
        base <- point_at(frac[member], 0)
        tip <- point_at(frac[member], side_sign * params$tick_length)
        append_path(
          rbind(base, tip), "tick", row, biological_position,
          biological_position, as.character(row)
        )

        if (isTRUE(params$label)) {
          lp <- point_at(slots[member], side_sign * params$label_offset)
          if (shared) {
            site_trunk <- point_at(frac[member], trunk_offset)
            label_trunk <- point_at(slots[member], trunk_offset)
            append_path(
              rbind(tip, site_trunk), "site_branch", row,
              biological_position, slots[member] * lens[id], member_text
            )
            append_path(
              rbind(label_trunk, lp), "label_branch", row,
              biological_position, slots[member] * lens[id], member_text
            )
          } else {
            branch <- if (params$leader == "elbow") {
              elbow <- data.frame(x = lp$x, y = tip$y)
              rbind(tip, elbow, lp)
            } else rbind(tip, lp)
            append_path(
              branch, "leader", row, biological_position,
              slots[member] * lens[id], as.character(row)
            )
          }

          position_text <- format(
            biological_position, big.mark = ",", scientific = FALSE, trim = TRUE
          )
          label_text <- switch(
            params$label_style,
            enzyme_position = paste0(data$enzyme[row], " (", position_text, ")"),
            enzyme = as.character(data$enzyme[row]),
            position = position_text
          )
          gid <- gid + 1L
          hjust <- if (lp$x > .12) 0 else if (lp$x < -.12) 1 else .5
          vjust <- if (abs(lp$x) <= .12 && lp$y > 0) 0 else if (
            abs(lp$x) <= .12 && lp$y < 0
          ) 1 else .5
          output[[length(output) + 1L]] <- data.frame(
            x = lp$x, y = lp$y, label = label_text,
            .component = "label", restriction_component = "label",
            group = gid, source_row = row,
            anchor_position = biological_position,
            slot_position = slots[member] * lens[id],
            member_rows = as.character(row), size = 2.6,
            angle = 0, hjust = hjust, vjust = vjust,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (!length(output)) {
    return(data.frame(
      x = numeric(), y = numeric(), label = character(),
      .component = character(), restriction_component = character(),
      group = integer(), source_row = integer(),
      anchor_position = numeric(), slot_position = numeric(),
      member_rows = character(), stringsAsFactors = FALSE
    ))
  }
  ggchord_rbind_fill(output)
}
