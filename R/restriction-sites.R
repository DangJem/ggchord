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
  if (!exists("ggchord_rebase_database", inherits = TRUE)) {
    ggchord_stop("The internal REBASE database is missing; reinstall ggchord")
  }
  get("ggchord_rebase_database", inherits = TRUE)
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
  normalized <- ggchord_normalize_sequence_input(
    sequence, "find_restriction_sites()"
  )
  ids <- normalized$ids
  seqs <- normalized$sequences
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

GeomRestrictionSite <- ggplot2::ggproto(
  "GeomRestrictionSite", ggplot2::Geom,
  required_aes = c("x", "y"),
  default_aes = ggplot2::aes(
    xend = NA_real_, yend = NA_real_, label = NA_character_,
    .component = NA_character_, group = NA_integer_, colour = "#202020",
    alpha = 1, linewidth = .28, linetype = 1, size = 2.9, angle = 0,
    hjust = .5, vjust = .5, family = "", fontface = 1, lineheight = 1.2
  ),
  draw_key = ggplot2::draw_key_path,
  draw_panel = function(data, panel_params, coord, na.rm = FALSE,
                        segment_params = list(), text_params = list(),
                        composite_labels = TRUE) {
    paths <- data[data$.component == "path", , drop = FALSE]
    labels <- data[data$.component == "label", , drop = FALSE]
    for (nm in names(segment_params)) {
      if (nm %in% names(paths)) paths[[nm]] <- segment_params[[nm]]
    }
    for (nm in names(text_params)) {
      if (nm %in% names(labels)) labels[[nm]] <- text_params[[nm]]
    }
    grobs <- list()
    if (nrow(paths)) {
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPath$draw_panel(
        paths, panel_params, coord, lineend = "round", linejoin = "round",
        na.rm = na.rm
      )
    }
    if (nrow(labels)) {
      composite <- isTRUE(composite_labels) &&
        "plotmath_label" %in% names(labels)
      if (composite) labels$label <- labels$plotmath_label
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomText$draw_panel(
        labels, panel_params, coord, parse = composite,
        check_overlap = FALSE, na.rm = na.rm
      )
    }
    do.call(grid::grobTree, grobs)
  }
)

#' Draw restriction sites, labels and deterministic leaders
#' @param mapping,data Standard layer inputs.
#' @param label Draw labels.
#' @param label_style Label enzyme and position, enzyme only, or position only.
#' @param label_side Draw labels outside, inside, or automatically.
#' @param leader Leader routing style. The default `"radial"` measures actual
#'   text boxes, separates site anchors from final radial label positions,
#'   and routes either a direct connector or an independent radial stub plus
#'   fan segment.
#'   Most labels follow one common radial contour; collision-bound clusters may
#'   move to a nearby outer contour. `"trunk"` is retained only as a
#'   compatibility alias.
#' @param min_label_gap Minimum genomic fraction between label slots. `NULL`
#'   derives a compact device-aware default from the rendered labels.
#' @param tick_length,label_offset Local-normal distances.
#' @param colour,linewidth,label_size,fontface,family Fixed appearance.
#' @param label_order Automatic mirrored label order, enzyme first, or
#'   position first.
#' @param position,show.legend,inherit.aes Standard layer arguments.
#' @param ... Additional geom parameters.
#' @return A ggplot2 layer.
#' @export
geom_restriction_site<-function(mapping=NULL,data=NULL,label=TRUE,label_style=c("enzyme_position","enzyme","position"),label_order=c("auto","enzyme_position","position_enzyme"),label_side=c("outside","inside","auto"),leader=c("radial","elbow","straight","trunk"),min_label_gap=NULL,tick_length=.04,label_offset=.18,colour="#202020",linewidth=.28,label_size=3.2,fontface=NULL,family=NULL,position="identity",show.legend=FALSE,inherit.aes=FALSE,...){
  colour_supplied<-!missing(colour);linewidth_supplied<-!missing(linewidth);label_size_supplied<-!missing(label_size);fontface_supplied<-!missing(fontface)&&!is.null(fontface);family_supplied<-!missing(family)&&!is.null(family)
  if(is.null(data))ggchord_stop("geom_restriction_site(): supply restriction-site data")
  label_style<-match.arg(label_style);label_order<-match.arg(label_order);label_side<-match.arg(label_side);leader<-match.arg(leader)
  if (identical(leader, "trunk")) leader <- "radial"
  vals<-c(tick_length,label_offset);if(!is.numeric(vals)||any(!is.finite(vals))||any(vals<0)||(!is.null(min_label_gap)&&(!is.numeric(min_label_gap)||length(min_label_gap)!=1L||!is.finite(min_label_gap)||min_label_gap<0)))ggchord_stop("geom_restriction_site(): invalid gaps or offsets")
  segment_params<-list();text_params<-list()
  if(colour_supplied){segment_params$colour<-colour;text_params$colour<-colour}
  if(linewidth_supplied)segment_params$linewidth<-linewidth
  if(label_size_supplied)text_params$size<-label_size
  if(fontface_supplied)text_params$fontface<-fontface
  if(family_supplied)text_params$family<-family
  lyr<-ggplot2::layer(data=data.frame(x=numeric(),y=numeric()),mapping=ggplot2::aes(x=x,y=y,group=group,label=label,.component=I(.component)),stat="identity",geom=GeomRestrictionSite,position=position,show.legend=show.legend,inherit.aes=inherit.aes,check.aes=FALSE,check.param=FALSE,params=c(list(na.rm=FALSE,colour=colour,linewidth=linewidth,size=label_size,segment_params=segment_params,text_params=text_params,composite_labels=!fontface_supplied),list(...)))
  lyr$ggchord_type<-"restriction_site";lyr$ggchord_theme_components<-c(segment_params="ggchord.restriction.label.segment",text_params="ggchord.restriction.label")
  lyr$ggchord_params<-list(type="restriction_site",label=label,label_style=label_style,label_order=label_order,label_side=label_side,leader=leader,min_label_gap=min_label_gap,tick_length=tick_length,label_offset=label_offset,label_size=label_size)
  ggchord_capture_layer_input(lyr,data,mapping,c("accver","position","enzyme"))
}

# Restriction labels use the same radial principles as gene-label repel, with
# local site clusters so a dense MCS may rise one level without moving every
# sparse callout. Most inner label edges therefore follow one common offset
# contour; only collision-bound clusters use a nearby outer contour.
ggchord_restriction_label_lanes <- function(
    gl, seq_arcs, side = "outside", units_per_inch = .3,
    label_offset = .13, min_label_gap = NULL) {
  frame <- ggchord_label_curve_frame(gl, seq_arcs)
  anchor_clearance <- max(abs(frame$signed_distance), na.rm = TRUE)
  if (!is.finite(anchor_clearance)) anchor_clearance <- 0
  point_padding <- max(.01,
    label_offset - anchor_clearance - .18 * units_per_inch)

  # min_label_gap is expressed as a sequence fraction. Convert it to the
  # solver's physical box padding, shared equally by adjacent labels.
  box_padding <- .05
  if (!is.null(min_label_gap)) {
    arc_lengths <- vapply(seq_arcs, function(arc) sum(sqrt(
      diff(arc$x)^2 + diff(arc$y)^2)), numeric(1))
    box_padding <- max(box_padding,
      min_label_gap * max(arc_lengths, na.rm = TRUE) /
        max(2 * units_per_inch, 1e-8))
  }
  labels <- gl
  labels$.radial_bend_x <- NA_real_
  labels$.radial_bend_y <- NA_real_
  labels$.radial_parameter <- NA_real_
  directions <- rep(NA_character_, nrow(gl))
  tracks <- rep(NA_integer_, nrow(gl))
  occupied <- ggchord_text_boxes(data.frame())
  fixed_paths <- data.frame(
    x0 = numeric(), y0 = numeric(), x1 = numeric(), y1 = numeric(),
    group = integer()
  )
  arc_ids <- vapply(seq_arcs, function(arc) as.character(arc$accver[1L]),
    character(1))

  # A single gene-label group intentionally shares one clearance level. A
  # plasmid restriction map is different: several independent site fans can
  # occur around the same circle. Solve each locally dense fan with the same
  # radial engine so an MCS does not push every other label onto a huge ring.
  groups <- list()
  totals <- numeric()
  for (id in unique(as.character(gl$accver))) {
    rows <- which(as.character(gl$accver) == id)
    arc <- seq_arcs[[match(id, arc_ids)]]
    distance <- c(0, cumsum(sqrt(diff(arc$x)^2 + diff(arc$y)^2)))
    keep <- !duplicated(distance)
    arc <- arc[keep, , drop = FALSE]
    distance <- distance[keep]
    total <- tail(distance, 1L)
    preferred <- vapply(rows, function(i) distance[which.min(
      (arc$x - gl$anchor_x[i])^2 + (arc$y - gl$anchor_y[i])^2
    )], numeric(1))
    ord <- order(preferred, rows)
    if (length(ord) > 1L) {
      sorted <- preferred[ord]
      gaps <- c(diff(sorted), total - tail(sorted, 1L) + sorted[1L])
      cut <- which.max(gaps)
      ord <- c(if (cut < length(ord)) ord[seq.int(cut + 1L, length(ord))]
        else integer(), ord[seq_len(cut)])
    }
    opened <- preferred[ord]
    if (length(opened) > 1L) {
      for (j in 2:length(opened)) if (opened[j] < opened[j - 1L]) {
        opened[j:length(opened)] <- opened[j:length(opened)] + total
      }
    }
    fan_gap <- max(.025, min_label_gap %||% .014) * total
    fan <- cumsum(c(TRUE, diff(opened) > fan_gap))
    local <- split(rows[ord], fan)
    groups <- c(groups, local)
    totals <- c(totals, rep(total, length(local)))
  }

  for (g in seq_along(groups)) {
    rows <- groups[[g]]
    measured <- gl[rows, , drop = FALSE]
    # Plotmath renders the enzyme fragment in bold, which is wider than the
    # plain composite string available to the base-device measurer. Invisible
    # side bearings reserve that horizontal difference without inflating the
    # vertical spacing of dense right/left fans.
    measured$text <- paste0("  ", gl$text[rows], "  ")
    solved <- ggchord_radial_label_lanes(
      measured, seq_arcs, side = side,
      units_per_inch = units_per_inch, box_padding = box_padding,
      point_padding = point_padding, repel_boxes = occupied,
      .fixed_paths = fixed_paths, .balance = FALSE,
      .tangent_limit = .055 * totals[g]
    )
    labels[rows, names(solved$labels)] <- solved$labels
    directions[rows] <- solved$directions
    tracks[rows] <- solved$tracks
    labels$.radial_parameter[rows] <-
      (solved$labels$.radial_parameter %% totals[g]) / totals[g]
    occupied <- rbind(occupied, ggchord_text_boxes(
      solved$labels, units_per_inch = units_per_inch,
      box_padding = box_padding
    ))
    for (i in seq_along(rows)) {
      bend <- c(solved$labels$.radial_bend_x[i],
        solved$labels$.radial_bend_y[i])
      direct <- sqrt(sum((bend - c(gl$anchor_x[rows[i]],
        gl$anchor_y[rows[i]]))^2)) < 1e-8
      if (direct) {
        fixed_paths <- rbind(fixed_paths, data.frame(
          x0 = gl$anchor_x[rows[i]], y0 = gl$anchor_y[rows[i]],
          x1 = solved$labels$text_x[i], y1 = solved$labels$text_y[i],
          group = rows[i]
        ))
      } else {
        fixed_paths <- rbind(fixed_paths, data.frame(
          x0 = c(gl$anchor_x[rows[i]], bend[1]),
          y0 = c(gl$anchor_y[rows[i]], bend[2]),
          x1 = c(bend[1], solved$labels$text_x[i]),
          y1 = c(bend[2], solved$labels$text_y[i]),
          group = rows[i]
        ))
      }
    }
  }
  labels$text <- gl$text
  boxes <- ggchord_text_boxes(labels, units_per_inch = units_per_inch)
  # Near twelve and six o'clock, solve all local clusters together. Ordering
  # by the real anchor x coordinate preserves genomic order across the origin;
  # projecting the packed centres back to one radius keeps the SnapGene-like
  # circular label contour instead of producing a flat Cartesian row.
  anchor_radius <- sqrt(gl$anchor_x^2 + gl$anchor_y^2)
  anchor_radius[anchor_radius <= 1e-10] <- 1
  for (sector in c("top", "bottom")) {
    rows <- if (sector == "top") {
      which(gl$anchor_y > 0 & abs(gl$anchor_x) / anchor_radius < .45)
    } else {
      which(gl$anchor_y < 0 & abs(gl$anchor_x) / anchor_radius < .68)
    }
    if (length(rows) < 2L) next
    half_width <- boxes$w[rows] * .62
    packed <- ggchord_pack_label_axis(
      boxes$cx[rows], gl$anchor_x[rows], half_width, half_width,
      gap = .030
    )
    # Do not let packing across the 0/180-degree seam swap a label to the
    # other side of the vertical axis.  SnapGene keeps pre-seam sites in the
    # left queue and post-seam sites in the right queue; this is also what
    # determines which enzyme-bearing text edge receives the connector.
    seam_gap <- .015
    left <- gl$anchor_x[rows] < 0
    if (any(left)) {
      overflow <- max(packed[left] + half_width[left] + seam_gap)
      if (overflow > 0) packed[left] <- packed[left] - overflow
    }
    if (any(!left)) {
      underflow <- min(packed[!left] - half_width[!left] - seam_gap)
      if (underflow < 0) packed[!left] <- packed[!left] - underflow
    }
    labels$text_x[rows] <- labels$text_x[rows] +
      packed - boxes$cx[rows]
    target_radius <- max(sqrt(
      packed^2 + boxes$cy[rows]^2
    ), max(abs(packed)) + .02)
    projected_y <- sqrt(pmax(0, target_radius^2 - packed^2))
    labels$text_y[rows] <- labels$text_y[rows] +
      if (sector == "top") projected_y - boxes$cy[rows] else
        -projected_y - boxes$cy[rows]
    directions[rows] <- sector
  }
  boxes <- ggchord_text_boxes(labels, units_per_inch = units_per_inch)
  labels$.text_width <- boxes$w
  labels$.text_height <- boxes$h
  labels$.text_center_x <- boxes$cx
  labels$.text_center_y <- boxes$cy
  list(labels = labels, directions = directions, tracks = tracks,
    lanes = paste(gl$accver, side, sep = "\r"),
    draw_segment = !is.na(gl$text) & nzchar(gl$text))
}

ggchord_restriction_text_metrics <- function(plotmath, size,
                                             units_per_inch) {
  n <- length(plotmath)
  size <- rep_len(size, n)
  width <- height <- numeric(n)
  valid <- !is.na(plotmath) & nzchar(plotmath)
  if (!any(valid)) return(data.frame(width = width, height = height))
  close_device <- ggchord_measurement_device()
  on.exit(close_device())
  for (i in which(valid)) {
    label <- tryCatch(parse(text = plotmath[i])[[1L]],
      error = function(e) plotmath[i])
    grob <- grid::textGrob(
      label, gp = grid::gpar(fontsize = size[i] * (72.27 / 25.4))
    )
    width[i] <- grid::convertWidth(
      grid::grobWidth(grob), "inches", valueOnly = TRUE
    ) * units_per_inch
    height[i] <- grid::convertHeight(
      grid::grobHeight(grob), "inches", valueOnly = TRUE
    ) * units_per_inch
  }
  data.frame(width = width, height = height)
}

# Retained experimental Cartesian solver. The radial layout above is preferred
# because plasmid labels should follow a common offset contour.
ggchord_restriction_label_lanes_cartesian <- function(
    gl, seq_arcs, side = "outside", units_per_inch = .3,
    label_offset = .18, min_label_gap = NULL) {
  frame <- ggchord_label_curve_frame(gl, seq_arcs)
  labels <- gl
  boxes <- ggchord_text_boxes(
    gl, units_per_inch = units_per_inch, box_padding = .035
  )
  # The rendered plotmath label bolds only the enzyme token. A small allowance
  # covers that weight change while keeping the measured and visible edges in
  # the same place; the former 1.5 multiplier displaced long-label anchors.
  width <- boxes$w * 1.08
  height <- boxes$h
  side_sign <- if (identical(side, "inside")) -1 else 1
  extra <- max(.04, label_offset - .04)
  labels$text_x <- gl$anchor_x + frame$outward_x * side_sign * extra
  labels$text_y <- gl$anchor_y + frame$outward_y * side_sign * extra
  labels$.radial_bend_x <- gl$anchor_x
  labels$.radial_bend_y <- gl$anchor_y
  labels$.radial_parameter <- NA_real_
  directions <- rep(NA_character_, nrow(gl))
  tracks <- integer(nrow(gl))
  arc_ids <- vapply(seq_arcs, function(arc) as.character(arc$accver[1L]),
    character(1))

  pack_ordered <- function(preferred, genomic, before, after, gap) {
    if (length(preferred) < 2L) return(preferred)
    correlation <- suppressWarnings(stats::cor(genomic, preferred))
    direction <- if (is.finite(correlation) && correlation < 0) -1 else 1
    ggchord_pack_label_axis(
      preferred, direction * genomic, before, after, gap = gap
    )
  }

  for (id in unique(as.character(gl$accver))) {
    rows <- which(as.character(gl$accver) == id)
    arc <- seq_arcs[[match(id, arc_ids)]]
    distance <- c(0, cumsum(sqrt(diff(arc$x)^2 + diff(arc$y)^2)))
    keep <- !duplicated(distance)
    arc <- arc[keep, , drop = FALSE]
    distance <- distance[keep]
    total <- tail(distance, 1L)
    parameter <- vapply(rows, function(i) distance[which.min(
      (arc$x - gl$anchor_x[i])^2 + (arc$y - gl$anchor_y[i])^2
    )], numeric(1))
    ord <- order(parameter, rows)
    if (length(ord) > 1L) {
      sorted <- parameter[ord]
      gaps <- c(diff(sorted), total - tail(sorted, 1L) + sorted[1L])
      cut <- which.max(gaps)
      ord <- c(if (cut < length(ord)) ord[seq.int(cut + 1L, length(ord))]
        else integer(), ord[seq_len(cut)])
    }
    opened <- parameter[ord]
    if (length(opened) > 1L) {
      for (j in 2:length(opened)) if (opened[j] < opened[j - 1L]) {
        opened[j:length(opened)] <- opened[j:length(opened)] + total
      }
    }
    fan_gap <- max(.04, min_label_gap %||% .014) * total
    fan <- cumsum(c(TRUE, diff(opened) > fan_gap))
    cluster <- integer(length(rows))
    cluster[ord] <- fan

    radius <- sqrt(gl$anchor_x[rows]^2 + gl$anchor_y[rows]^2)
    radius[radius <= 1e-8] <- 1
    top_bottom <- abs(gl$anchor_x[rows]) / radius < .38
    directions[rows] <- ifelse(
      top_bottom,
      ifelse(gl$anchor_y[rows] >= 0, "top", "bottom"),
      ifelse(gl$anchor_x[rows] >= 0, "right", "left")
    )

    # Dense lateral clusters share their inner text edge and use an ordered
    # vertical axis. This gives long MCS callouts enough room without allowing
    # labels to exchange genomic order.
    for (cluster_id in unique(cluster)) {
      members <- rows[cluster == cluster_id]
      lateral <- unique(directions[members])
      if (length(members) < 3L || length(lateral) != 1L ||
          !lateral %in% c("left", "right")) next
      local_parameter <- if ("genomic_position" %in% names(gl)) {
        gl$genomic_position[members]
      } else parameter[match(members, rows)]
      labels$text_y[members] <- pack_ordered(
        labels$text_y[members], local_parameter,
        height[members] * .58, height[members] * .58,
        gap = .018
      )
      if (lateral == "right") {
        inner_edge <- max(
          labels$text_x[members] - width[members] / 2,
          gl$anchor_x[members] + .10
        )
        labels$text_x[members] <- inner_edge + width[members] / 2
      } else {
        inner_edge <- min(
          labels$text_x[members] + width[members] / 2,
          gl$anchor_x[members] - .10
        )
        labels$text_x[members] <- inner_edge - width[members] / 2
      }
      tracks[members] <- 1L
    }

    # Split the top and bottom around x = 0, then pack each queue along x.
    # Central labels receive a little more normal clearance, which forms a
    # restrained staircase instead of a flat row.
    for (sector in c("top", "bottom")) {
      sector_rows <- rows[directions[rows] == sector]
      if (!length(sector_rows)) next
      for (half in c("left", "right")) {
        members <- sector_rows[
          if (half == "left") gl$anchor_x[sector_rows] < 0 else
            gl$anchor_x[sector_rows] >= 0
        ]
        if (!length(members)) next
        local_parameter <- if ("genomic_position" %in% names(gl)) {
          gl$genomic_position[members]
        } else parameter[match(members, rows)]
        labels$text_x[members] <- pack_ordered(
          labels$text_x[members], local_parameter,
          width[members] * .30, width[members] * .30,
          gap = .018
        )
        if (half == "left") {
          intrusion <- max(labels$text_x[members] + width[members] / 2) + .035
          if (intrusion > 0) {
            labels$text_x[members] <- labels$text_x[members] - intrusion
          }
        } else {
          intrusion <- .035 - min(
            labels$text_x[members] - width[members] / 2
          )
          if (intrusion > 0) {
            labels$text_x[members] <- labels$text_x[members] + intrusion
          }
        }
        centre_order <- order(abs(gl$anchor_x[members]), members)
        stair <- integer(length(members))
        stair[centre_order] <- rev(seq_along(members)) - 1L
        stair_step <- max(height[members] * .90 + .012)
        labels$text_y[members] <- labels$text_y[members] +
          if (sector == "top") stair * stair_step else -stair * stair_step
        tracks[members] <- as.integer(length(members) > 1L)
      }
    }
    labels$.radial_parameter[rows] <- parameter / total
  }

  labels$text <- gl$text
  boxes <- ggchord_text_boxes(labels, units_per_inch = units_per_inch)
  labels$.text_width <- boxes$w
  labels$.text_height <- boxes$h
  labels$.text_center_x <- boxes$cx
  labels$.text_center_y <- boxes$cy
  list(labels = labels, directions = directions, tracks = tracks,
    lanes = paste(gl$accver, side, sep = "\r"),
    draw_segment = !is.na(gl$text) & nzchar(gl$text))
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

  output <- list(); gid <- 0L
  append_path <- function(points, component, source_row = NA_integer_,
                          anchor_position = NA_real_, slot_position = NA_real_,
                          source_rows = integer(), cluster_id = NA_character_,
                          junction_id = NA_character_,
                          label_direction = NA_character_,
                          label_order = NA_character_,
                          label_connection_side = NA_character_) {
    gid <<- gid + 1L
    points$.component <- "path"
    points$restriction_component <- component
    points$group <- gid
    points$label <- NA_character_
    points$source_row <- source_row
    points$anchor_position <- anchor_position
    points$slot_position <- slot_position
    points$source_rows <- I(rep(list(as.integer(source_rows)), nrow(points)))
    points$cluster_id <- cluster_id
    points$junction_id <- junction_id
    points$label_direction <- label_direction
    points$label_order <- label_order
    points$label_connection_side <- label_connection_side
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
    # Rendering collapses enzymes at exactly the same cleavage coordinate into
    # one deterministic callout. Biological search rows remain untouched and
    # are preserved in source_rows for export/provenance.
    site_members <- unname(split(idx,
      factor(data$position[idx], levels = unique(data$position[idx]))))
    idx <- vapply(site_members, `[`, integer(1), 1L)
    site_enzyme <- vapply(site_members, function(rows) paste(
      unique(as.character(data$enzyme[rows])), collapse = " - "
    ), character(1))
    arc <- layout$seq_arcs[[id]]
    n <- nrow(arc)
    frac <- data$position[idx] / lens[id]

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
    position_text <- format(data$position[idx], big.mark = "",
      scientific = FALSE, trim = TRUE)
    raw_labels <- switch(params$label_style,
      enzyme_position = paste0(site_enzyme, " (", position_text, ")"),
      enzyme = site_enzyme, position = position_text)
    # Cluster IDs remain useful export metadata, but they do not imply shared
    # geometry. Open at the largest circular gap so origin-adjacent sites keep
    # stable neighbouring IDs.
    cluster <- seq_along(idx)
    if (length(frac) > 1L) {
      gap_limit <- max(.025, params$min_label_gap %||% .014)
      opened <- frac
      if (isTRUE(layout$circular)) {
        gaps <- diff(c(frac, frac[1L] + 1))
        cut <- which.max(gaps)
        ord <- c(if (cut < length(frac)) seq.int(cut + 1L, length(frac)) else
          integer(), seq_len(cut))
        opened <- frac[ord]
        opened[opened < opened[1L]] <- opened[opened < opened[1L]] + 1
        opened_cluster <- cumsum(c(TRUE, diff(opened) > gap_limit))
        cluster[ord] <- opened_cluster
      } else {
        cluster <- cumsum(c(TRUE, diff(frac) > gap_limit))
      }
    }
    cluster_names <- paste0(id, ":cluster:", sprintf("%03d", cluster))
    junction_names <- paste0(id, ":site:", sprintf("%06d", idx))
    bases <- point_at(frac, 0)
    tips <- point_at(frac, side_sign * params$tick_length)

    for (member in seq_along(idx)) {
      row <- idx[member]
      append_path(
        rbind(bases[member, , drop = FALSE], tips[member, , drop = FALSE]),
        "tick", row, data$position[row], data$position[row],
        site_members[[member]],
        cluster_names[member], junction_names[member]
      )
    }
    if (!isTRUE(params$label)) next

    gl <- data.frame(
      text = raw_labels, text_x = tips$x, text_y = tips$y,
      text_angle = 0, hjust = .5, vjust = .5,
      size = params$label_size %||% 2.9, accver = id,
      group = seq_along(idx), source_row = idx,
      genomic_position = data$position[idx],
      anchor_x = tips$x, anchor_y = tips$y,
      stringsAsFactors = FALSE
    )
    device <- if (grDevices::dev.cur() == 1L) c(8, 8) else
      grDevices::dev.size("in")
    span <- c(diff(range(arc$x)), diff(range(arc$y)))
    units_per_inch <- max(span / pmax(device - c(2.2, 1.5), device * .45))
    units_per_inch <- max(.20, min(.42, units_per_inch))
    solved <- ggchord_restriction_label_lanes(
      gl, layout$seq_arcs, side = params$label_side,
      units_per_inch = units_per_inch,
      label_offset = params$label_offset %||% .13,
      min_label_gap = params$min_label_gap
    )
    labels <- solved$labels
    directions <- solved$directions
    tracks <- solved$tracks
    bend_x <- labels$.radial_bend_x
    bend_y <- labels$.radial_bend_y
    invalid_bend <- !is.finite(bend_x) | !is.finite(bend_y)
    bend_x[invalid_bend] <- tips$x[invalid_bend]
    bend_y[invalid_bend] <- tips$y[invalid_bend]
    centre_x <- labels$.text_center_x
    centre_y <- labels$.text_center_y
    # Text order and attachment edge follow the label's visual half of the
    # plasmid, not merely the centre-to-bend vector.  The latter is unstable
    # near twelve and six o'clock because a wide label can straddle its bend.
    visual_x <- centre_x
    visual_x[abs(visual_x) < .010] <- gl$anchor_x[abs(visual_x) < .010]
    connection_side <- ifelse(visual_x < 0, "right", "left")
    order_values <- rep(params$label_order, length(idx))
    if (params$label_order == "auto") {
      order_values <- ifelse(connection_side == "right",
        "position_enzyme", "enzyme_position")
    }
    quote_text <- function(x) paste0("'", gsub("'", "\\\\'", x,
      fixed = TRUE), "'")
    coordinate <- format(data$position[idx], big.mark = "",
      scientific = FALSE, trim = TRUE)
    label_text <- vapply(seq_along(idx), function(member) {
      switch(params$label_style,
        enzyme_position = if (order_values[member] == "position_enzyme")
          paste0("(", coordinate[member], ") ", site_enzyme[member]) else
          paste0(site_enzyme[member], " (", coordinate[member], ")"),
        enzyme = site_enzyme[member], position = coordinate[member])
    }, character(1))
    plotmath <- vapply(seq_along(idx), function(member) {
      switch(params$label_style,
        enzyme_position = if (order_values[member] == "position_enzyme")
          paste0(quote_text(paste0("(", coordinate[member], ")")),
            "~bold(", quote_text(site_enzyme[member]), ")") else
          paste0("bold(", quote_text(site_enzyme[member]), ")~",
            quote_text(paste0("(", coordinate[member], ")"))),
        enzyme = paste0("bold(", quote_text(site_enzyme[member]), ")"),
        position = quote_text(coordinate[member]))
    }, character(1))
    # Measure the expression that GeomText will actually draw. Measuring the
    # unparsed plain string makes bold enzyme names several pixels wider than
    # the routing box, which visually puts the connector inside the first or
    # last glyph even when the nominal side is correct.
    metrics <- ggchord_restriction_text_metrics(
      plotmath, params$label_size %||% 2.9, units_per_inch
    )
    routing_half_width <- metrics$width / 2 + .010
    routing_half_height <- metrics$height / 2 + .010
    endpoint_x <- centre_x + ifelse(
      connection_side == "left", -routing_half_width, routing_half_width
    )
    # A selected left edge must lie to the right of its bend, and a selected
    # right edge to the left.  Bottom/top packing can otherwise leave a wide
    # bbox straddling the bend, so the nominally correct edge is approached
    # backwards. Move only as far as needed to restore the actual approach.
    route_clearance <- .012
    shift_x <- ifelse(
      connection_side == "left",
      pmax(0, bend_x + route_clearance - endpoint_x),
      pmin(0, bend_x - route_clearance - endpoint_x)
    )
    centre_x <- centre_x + shift_x
    endpoint_x <- endpoint_x + shift_x
    labels$.text_center_x <- centre_x
    labels$text_x <- labels$text_x + shift_x
    # SnapGene-style site labels always connect at the enzyme-bearing left or
    # right edge, never at the middle of the top/bottom edge. Along that chosen
    # vertical edge, use the point nearest the bend; a top label therefore
    # naturally connects at its lower-left/lower-right corner.
    endpoint_y <- pmax(
      centre_y - routing_half_height,
      pmin(bend_y, centre_y + routing_half_height)
    )
    # Keep the rendered label centred on the box used by the routing solver.
    # The leader endpoint is a separate bbox intersection, so justification no
    # longer shifts the visible edge after routing has been computed.
    labels$hjust <- .5
    labels$vjust <- .5
    labels$text_x <- centre_x
    labels$text_y <- centre_y
    slot_positions <- if (".radial_parameter" %in% names(labels)) {
      labels$.radial_parameter * lens[id]
    } else vapply(seq_along(idx), function(member) {
      k <- which.min((arc$x - labels$text_x[member])^2 +
        (arc$y - labels$text_y[member])^2)
      (k - 1L) / max(1L, nrow(arc) - 1L) * lens[id]
    }, numeric(1))

    leader_frame <- ggchord_label_curve_frame(gl, layout$seq_arcs)
    segments <- lapply(seq_along(idx), function(member) {
      vector <- c(endpoint_x[member] - tips$x[member],
        endpoint_y[member] - tips$y[member])
      radial <- c(leader_frame$outward_x[member] * side_sign,
        leader_frame$outward_y[member] * side_sign)
      denominator <- sqrt(sum(vector^2) * sum(radial^2))
      centre_vector <- c(labels$.text_center_x[member] - tips$x[member],
        labels$.text_center_y[member] - tips$y[member])
      centre_length <- sqrt(sum(centre_vector^2))
      aligned <- is.finite(denominator) && denominator > 1e-12 &&
        is.finite(centre_length) && centre_length > 1e-12 &&
        sum(centre_vector * radial) > 0 &&
        abs(centre_vector[1] * radial[2] - centre_vector[2] * radial[1]) /
          centre_length < .02
      # A label still on its natural radial ray uses one segment. Tangentially
      # displaced labels use exactly the bend selected by the radial solver.
      if (identical(params$leader, "straight") || aligned) {
        return(data.frame(x0 = tips$x[member], y0 = tips$y[member],
          x1 = endpoint_x[member], y1 = endpoint_y[member], group = member))
      }
      bend <- c(labels$.radial_bend_x[member],
        labels$.radial_bend_y[member])
      if (any(!is.finite(bend)) || sqrt(sum((bend - unlist(tips[member, ]))^2)) <
          1e-8) {
        distance <- sqrt(sum(vector^2))
        shoulder <- min(.085, max(.045, distance * .30))
        bend <- c(
          tips$x[member] + radial[1] * shoulder,
          tips$y[member] + radial[2] * shoulder
        )
      }
      data.frame(
        x0 = c(tips$x[member], bend[1]),
        y0 = c(tips$y[member], bend[2]),
        x1 = c(bend[1], endpoint_x[member]),
        y1 = c(bend[2], endpoint_y[member]), group = member
      )
    })
    segments <- do.call(rbind, segments)
    for (segment in seq_len(nrow(segments))) {
      member <- segments$group[segment]
      row <- idx[member]
      append_path(
        data.frame(x = c(segments$x0[segment], segments$x1[segment]),
          y = c(segments$y0[segment], segments$y1[segment])),
        "leader", row, data$position[row], slot_positions[member],
        site_members[[member]],
        cluster_names[member], junction_names[member], directions[member],
        order_values[member], connection_side[member]
      )
    }

    for (member in seq_along(idx)) {
      row <- idx[member]
      direction <- directions[member]
      order_value <- order_values[member]
      gid <- gid + 1L
      output[[length(output) + 1L]] <- data.frame(
        x = labels$text_x[member], y = labels$text_y[member],
        label = label_text[member], plotmath_label = plotmath[member],
        .component = "label", restriction_component = "label",
        group = gid, source_row = row,
        anchor_position = data$position[row],
        slot_position = slot_positions[member],
        source_rows = I(list(as.integer(site_members[[member]]))),
        cluster_id = cluster_names[member],
        junction_id = junction_names[member],
        label_direction = direction, label_order = order_value,
        label_connection_side = connection_side[member],
        size = params$label_size %||% 2.9, angle = 0,
        hjust = labels$hjust[member], vjust = labels$vjust[member],
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(output)) {
    return(data.frame(
      x = numeric(), y = numeric(), label = character(),
      .component = character(), restriction_component = character(),
      group = integer(), source_row = integer(),
      anchor_position = numeric(), slot_position = numeric(),
      source_rows = I(list()), cluster_id = character(),
      junction_id = character(), label_direction = character(),
      label_order = character(), label_connection_side = character(),
      stringsAsFactors = FALSE
    ))
  }
  ggchord_rbind_fill(output)
}
