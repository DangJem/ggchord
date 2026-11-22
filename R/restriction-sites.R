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
  # gregexpr(..., fixed = TRUE) consumes a match, so overlapping anchors such
  # as GGG in GGGG were previously skipped.  A zero-width look-ahead retains
  # every possible anchor start while keeping the inexpensive fixed anchor.
  hits <- gregexpr(paste0("(?=", anchor, ")"), search, perl=TRUE)[[1L]]
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
  source_motif <- mat[,2L]
  defs <- data.frame(
    pattern_id = sprintf("rebase%s:e:%06d", version, source_rows),
    pattern_source_row = source_rows,
    enzyme = mat[,1L], motif = toupper(source_motif),
    source_motif = source_motif,
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
  # EMBOSS uses a lower-case recognition sequence when an equivalent enzyme's
  # cleavage orientation is reversed relative to the normalized motif.  Keep
  # the source spelling, but orient the top/bottom cut pairs consistently for
  # downstream coordinate calculation.
  reverse_cleavage <- source_motif == tolower(source_motif) &
    source_motif != toupper(source_motif)
  if (any(reverse_cleavage)) {
    first <- defs$cut_offset_1[reverse_cleavage]
    defs$cut_offset_1[reverse_cleavage] <-
      defs$cut_offset_2[reverse_cleavage]
    defs$cut_offset_2[reverse_cleavage] <- first
    third <- defs$cut_offset_3[reverse_cleavage]
    defs$cut_offset_3[reverse_cleavage] <-
      defs$cut_offset_4[reverse_cleavage]
    defs$cut_offset_4[reverse_cleavage] <- third
  }

  supplier_lines <- readLines(file.path(path, "embossa_s.txt"), warn = FALSE)
  supplier_lines <- supplier_lines[nzchar(trimws(supplier_lines)) & !grepl("^#", supplier_lines)]
  supplier_code <- sub("[[:space:]].*$", "", supplier_lines)
  supplier_name <- sub("^[^[:space:]]+[[:space:]]+", "", supplier_lines)
  supplier_map <- setNames(supplier_name, supplier_code)

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

  equivalence_path <- file.path(path, "misc", "embossre.equ")
  preferred <- setNames(character(), character())
  if (file.exists(equivalence_path)) {
    equivalence_lines <- readLines(equivalence_path, warn = FALSE)
    equivalence_lines <- equivalence_lines[
      nzchar(trimws(equivalence_lines)) & !grepl("^#", equivalence_lines)
    ]
    equivalence_fields <- strsplit(trimws(equivalence_lines), "[[:space:]]+")
    if (any(lengths(equivalence_fields) != 2L)) {
      ggchord_stop("Invalid misc/embossre.equ row width")
    }
    equivalence <- do.call(rbind, equivalence_fields)
    preferred <- setNames(equivalence[, 2L], equivalence[, 1L])
  }
  mapped <- unname(preferred[defs$enzyme])
  defs$preferred_enzyme <- ifelse(is.na(mapped), defs$enzyme, mapped)
  defs$is_preferred_enzyme <- defs$enzyme == defs$preferred_enzyme
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

# Source-backed methylation sensitivities used by the generic site-status
# evaluator.  This is enzyme metadata, never a plasmid/name rendering special
# case. Unknown enzymes remain unknown rather than being guessed.
ggchord_restriction_methylation_sensitivity <- function() {
  data.frame(
    enzyme = "PflMI", blocked_by = "Dcm", methylation_motif = "CCWGG",
    sensitivity_source = paste0(
      "https://www.snapgene.com/plasmids/fluorescent_protein_genes_and_plasmids/",
      "DsRed-Express"
    ),
    stringsAsFactors = FALSE
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
  defaults <- list(source_motif=out$motif,organism=NA_character_,
    supplier_codes=NA_character_,suppliers=NA_character_,commercial=FALSE,
    preferred_enzyme=out$enzyme,is_preferred_enzyme=TRUE,
    database_version="custom",source="custom", blocked_by=NA_character_,
    methylation_motif=NA_character_, sensitivity_source=NA_character_,
    methylation_status=NA_character_, display_warning=NA_character_)
  for(nm in names(defaults))if(!nm%in%names(out))out[[nm]]<-defaults[[nm]]
  sensitivity <- ggchord_restriction_methylation_sensitivity()
  sensitivity_row <- match(out$enzyme, sensitivity$enzyme)
  for (nm in c("blocked_by", "methylation_motif", "sensitivity_source")) {
    missing_value <- is.na(out[[nm]]) | !nzchar(as.character(out[[nm]]))
    fill <- sensitivity[[nm]][sensitivity_row]
    use <- missing_value & !is.na(fill)
    out[[nm]][use] <- fill[use]
  }
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
#' @param methylation Methylation context. `"auto"` uses the audited context
#'   for recognized reference plasmids and otherwise makes no blocked-site
#'   claim; `"none"` disables blocking; `"dam_dcm"` evaluates source-backed
#'   Dam/Dcm sensitivities against the supplied sequence.
#' @return One row per pattern match without display filtering or deduplication.
#'   Methylation-aware columns include `methylation_status`, `blocked_by`,
#'   `display_warning`, `methylation_motif`, and `sensitivity_source`.
#'   `enzyme_site_count` records the number of distinct cleavage/display
#'   positions for that enzyme on the complete supplied sequence, so it remains
#'   stable after display filtering. For an exactly recognized bundled
#'   reference sequence, the `reference_enzyme_profiles` attribute records its
#'   saved display set for use by [filter_restriction_sites()].
#' @seealso \code{\link{ggchord_data_contract}}
#' @export
#' @examples
#' target <- c(A="ACGGAATTCGTT")
#' found <- find_restriction_sites(target, enzymes="EcoRI", circular=FALSE)
#' ref <- as_ggchord_seq(data.frame(accver="A", length=nchar(target)), sequence_unit="bp")
#' ggchord(ref) + geom_restriction_site(data=found) + coord_collinear()
find_restriction_sites <- function(sequence, enzymes = NULL, patterns = NULL,
    circular = NULL, database = NULL, methylation = "unknown") {
  ggchord_local_console_errors()
  caller <- "find_restriction_sites()"
  methylation <- ggchord_match_arg(methylation, c("unknown", "none", "dam_dcm"))
  normalized <- ggchord_normalize_sequence_input(
    sequence, "find_restriction_sites()"
  )
  ids <- normalized$ids
  seqs <- normalized$sequences
  topology <- ggchord_search_topology(sequence, circular, ids, caller)
  if (!is.null(patterns) && !is.null(database)) ggchord_stop(caller, ": supply patterns or database")
  defs<-if(!is.null(patterns))patterns else database%||%ggchord_builtin_rebase()
  defs<-ggchord_normalize_restriction_patterns(defs)
  if(!is.null(enzymes)){
    if(!is.character(enzymes)||anyNA(enzymes))ggchord_stop("find_restriction_sites(): enzymes must be character")
    unknown<-setdiff(enzymes,unique(defs$enzyme));if(length(unknown))ggchord_stop("Unknown restriction enzyme(s): ",paste(unknown,collapse=", "))
    defs<-defs[defs$enzyme%in%enzymes,,drop=FALSE]
    defs<-defs[order(match(defs$enzyme,enzymes),seq_len(nrow(defs))),,drop=FALSE]
  }
  empty_template <- function(){
    out<-data.frame(accver=character(),match_id=character(),pattern_id=character(),pattern_source_row=integer(),enzyme=character(),preferred_enzyme=character(),is_preferred_enzyme=logical(),motif=character(),source_motif=character(),motif_length=integer(),recognition_sequence=character(),start=integer(),end=integer(),crosses_origin=logical(),strand=character(),ncuts=integer(),blunt=logical(),end_type=character(),position=numeric(),display_position=numeric(),enzyme_site_count=integer(),anchor_kind=character(),organism=character(),supplier_codes=character(),suppliers=character(),commercial=logical(),database_version=character(),source=character(),blocked_by=character(),methylation_motif=character(),sensitivity_source=character(),methylation_status=character(),display_warning=character(),stringsAsFactors=FALSE)
    for(k in 1:4){out[[paste0("cut_offset_",k)]]<-numeric();out[[paste0("cut_",k,"_unwrapped")]]<-numeric();out[[paste0("cut_",k)]]<-numeric()};out
  }
  rows<-list();match_number<-0L
  for(s in seq_along(seqs)){
    target <- seqs[s]
    circular <- topology[s]
    len <- nchar(target)
    match_cache <- new.env(hash = TRUE, parent = emptyenv())
    for(d in seq_len(nrow(defs))){
    motif<-defs$motif[d];m<-nchar(motif)
    if (m > len) next
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
        if (circular) cuts[, known] <- ((unwrapped[, known, drop = FALSE] - 1) %% len) + 1
        else cuts[cuts < 0 | cuts > len] <- NA_real_
        finish_unwrapped<-starts+m-1L;finish<-if(circular)((finish_unwrapped-1L)%%len)+1L else finish_unwrapped
        end_type<-if(defs$ncuts[d]==0L)"unknown" else if(defs$ncuts[d]==1L)"nick" else if(defs$ncuts[d]==4L)"complex" else if(isTRUE(defs$blunt[d]))"blunt" else if(defs$cut_offset_1[d]<defs$cut_offset_2[d])"5_prime_overhang" else "3_prime_overhang"
        known_indices<-which(known);has_cut<-length(known_indices)>0L
        first_known<-if(has_cut)known_indices[1L] else NA_integer_
        position<-if(has_cut)cuts[,first_known] else starts
        row<-data.frame(accver=rep(ids[s],count),match_id=sprintf("%s:%s:%s:%d:%06d",ids[s],defs$pattern_id[d],site_strand,starts,numbers),pattern_id=rep(defs$pattern_id[d],count),pattern_source_row=rep(defs$pattern_source_row[d],count),enzyme=rep(defs$enzyme[d],count),preferred_enzyme=rep(defs$preferred_enzyme[d],count),is_preferred_enzyme=rep(defs$is_preferred_enzyme[d],count),motif=rep(motif,count),source_motif=rep(defs$source_motif[d],count),motif_length=rep(m,count),recognition_sequence=substring(search,starts,starts+m-1L),start=starts,end=finish,crosses_origin=finish_unwrapped>len,strand=rep(site_strand,count),ncuts=rep(defs$ncuts[d],count),blunt=rep(defs$blunt[d],count),end_type=rep(end_type,count),position=position,display_position=position,anchor_kind=rep(if(has_cut)"cut" else "recognition",count),organism=rep(defs$organism[d],count),supplier_codes=rep(defs$supplier_codes[d],count),suppliers=rep(defs$suppliers[d],count),commercial=rep(defs$commercial[d],count),database_version=rep(defs$database_version[d],count),source=rep(defs$source[d],count),blocked_by=rep(defs$blocked_by[d],count),methylation_motif=rep(defs$methylation_motif[d],count),sensitivity_source=rep(defs$sensitivity_source[d],count),methylation_status=rep(defs$methylation_status[d],count),display_warning=rep(defs$display_warning[d],count),stringsAsFactors=FALSE)
        for(k in 1:4){row[[paste0("cut_offset_",k)]]<-rep(as.numeric(defs[[paste0("cut_offset_",k)]][d]),count);row[[paste0("cut_",k,"_unwrapped")]]<-unwrapped[,k];row[[paste0("cut_",k)]]<-cuts[,k]}
        rows[[length(rows)+1L]]<-row
      }
    }
  }}
  out<-if(length(rows))ggchord_rbind_fill(rows)else empty_template();rownames(out)<-NULL
  if(nrow(out)){
    enzyme_key<-paste(out$accver,out$enzyme,sep="\r")
    count_by_enzyme<-vapply(split(out$position,enzyme_key),function(x)
      length(unique(x[is.finite(x)])),integer(1L))
    out$enzyme_site_count<-unname(count_by_enzyme[enzyme_key])
  }
  reference_profiles <- data.frame()
  common_db <- tryCatch(ggchord_builtin_common_features(),
    error = function(e) NULL)
  if (!is.null(common_db) && all(c("reference_sequences",
      "reference_enzyme_profiles") %in% names(common_db))) {
    profile_rows <- lapply(seq_along(seqs), function(i) {
      matched <- common_db$reference_sequences$reference_id[
        common_db$reference_sequences$sequence == seqs[i]]
      if (!length(matched)) return(NULL)
      profile_index <- match(matched[1L],
        common_db$reference_enzyme_profiles$reference_id)
      if (is.na(profile_index)) return(NULL)
      profile <- common_db$reference_enzyme_profiles[
        profile_index, , drop = FALSE]
      profile$accver <- ids[i]
      profile
    })
    profile_rows <- Filter(Negate(is.null), profile_rows)
    if (length(profile_rows)) reference_profiles <-
      ggchord_rbind_fill(profile_rows)
  }
  if (nrow(out) && methylation == "unknown") {
    out$methylation_status <- rep("unknown", nrow(out))
    formerly_blocked <- out$display_warning %in% "methylation_blocked"
    out$display_warning[formerly_blocked] <- NA_character_
  }
  if (nrow(out) && methylation != "unknown") {
    explicit_status <- !is.na(out$methylation_status) &
      nzchar(out$methylation_status)
    if (methylation == "none") {
      inferred <- ifelse(is.na(out$blocked_by),
        "unknown", "unblocked")
      out$methylation_status[!explicit_status] <- inferred[!explicit_status]
    } else {
      target_by_id <- setNames(seqs, ids)
      context_matches <- vapply(seq_len(nrow(out)), function(i) {
        if (is.na(out$blocked_by[i]) || is.na(out$methylation_motif[i])) {
          return(FALSE)
        }
        blocker <- tolower(out$blocked_by[i])
        if (!blocker %in% c("dam", "dcm")) return(FALSE)
        target <- target_by_id[[as.character(out$accver[i])]]
        len <- nchar(target)
        positions <- seq.int(out$start[i] - 6L,
          out$start[i] + out$motif_length[i] + 5L)
        positions <- if (topology[[as.character(out$accver[i])]]) ((positions - 1L) %% len) + 1L else
          positions[positions >= 1L & positions <= len]
        window <- paste0(strsplit(target, "", fixed = TRUE)[[1L]][positions],
          collapse = "")
        grepl(ggchord_iupac_regex(out$methylation_motif[i]), window,
          perl = TRUE)
      }, logical(1L))
      known <- !is.na(out$blocked_by)
      out$methylation_status[known & !explicit_status] <- "unblocked"
      out$methylation_status[context_matches & !explicit_status] <- "blocked"
      inferred_warning <- context_matches & !explicit_status
      out$display_warning[inferred_warning] <- "methylation_blocked"
    }
  }
  out$site_id <- out$match_id
  out$label <- out$enzyme
  out$pos <- out$recognition_start <- out$start
  out$recognition_end <- out$end
  out$cut_top <- out$cut_1
  out$cut_bottom <- out$cut_2
  out$segments <- I(ggchord_segment_search_matches(out$start, out$end, nchar(seqs[match(out$accver, ids)])))
  attr(out, "reference_enzyme_profiles") <- reference_profiles
  ggchord_attach_data_report(out, "find_restriction_sites", inputs = ids,
    metadata = list(database_version = unique(defs$database_version),
      topology = topology, methylation = methylation, reference_profiles = reference_profiles))
}

ggchord_filter_restriction_profile<-function(sites,set=c("all","unique","unique_dual","six_plus","unique_6plus","commercial","reference"),enzymes=NULL,min_site_length=NULL,cuts=NULL,window=NULL,commercial_only=FALSE,parent_set=c("all","commercial_all","commercial_nonredundant")){
  if(!is.data.frame(sites))ggchord_stop("filter_restriction_sites(): sites must be a data frame")
  set<-ggchord_match_arg(set);parent_set<-ggchord_match_arg(parent_set);if(!nrow(sites))return(sites)
  ggchord_require_columns(sites,c("accver","enzyme","motif_length","position"),"filter_restriction_sites()")
  if (set == "reference") {
    profiles <- attr(sites, "reference_enzyme_profiles")
    if (!is.data.frame(profiles) || !nrow(profiles)) {
      ggchord_stop(
        "filter_restriction_sites(set = \"reference\"): no exact reference ",
        "enzyme profile is available"
      )
    }
    missing_profiles <- setdiff(unique(as.character(sites$accver)),
      as.character(profiles$accver))
    if (length(missing_profiles)) {
      ggchord_stop(
        "filter_restriction_sites(set = \"reference\"): no profile for ",
        paste(missing_profiles, collapse = ", ")
      )
    }
    filtered <- lapply(unique(as.character(sites$accver)), function(id) {
      part <- sites[as.character(sites$accver) == id, , drop = FALSE]
      profile <- profiles[match(id, profiles$accver), , drop = FALSE]
      if (profile$profile_type == "none") return(part[FALSE, , drop = FALSE])
      if (profile$profile_type == "custom") {
        return(ggchord_filter_restriction_profile(part, set = "all",
          enzymes = profile$enzymes[[1L]], min_site_length = min_site_length,
          cuts = cuts, window = window, commercial_only = commercial_only,
          parent_set = "all"))
      }
      preset <- c(
        "Unique Cutters" = "unique",
        "Unique & Dual Cutters" = "unique_dual",
        "6+ Cutters" = "six_plus",
        "Unique 6+ Cutters" = "unique_6plus"
      )[[profile$set_name]]
      if (is.null(preset)) {
        ggchord_stop("Unsupported reference enzyme set: ", profile$set_name)
      }
      ggchord_filter_restriction_profile(part, set = preset,
        min_site_length = min_site_length, cuts = cuts, window = window,
        commercial_only = commercial_only, parent_set = parent_set)
    })
    attrs<-attributes(sites);out<-ggchord_rbind_fill(filtered);rownames(out)<-NULL
    for(nm in setdiff(names(attrs),c("names","row.names","class")))
      attr(out,nm)<-attrs[[nm]]
    return(out)
  }
  key<-paste(sites$accver,sites$enzyme,sep="\r")
  site_count<-if("enzyme_site_count"%in%names(sites))
    as.integer(sites$enzyme_site_count) else as.integer(table(key)[key])
  keep<-rep(TRUE,nrow(sites))
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
  if(parent_set!="all"){
    if(!"commercial"%in%names(sites))ggchord_stop("filter_restriction_sites(): sites is missing commercial")
    keep<-keep&!is.na(sites$commercial)&sites$commercial
  }
  if(parent_set=="commercial_nonredundant"){
    ggchord_require_columns(sites,c("preferred_enzyme","is_preferred_enzyme"),"filter_restriction_sites()")
    candidate<-which(keep)
    if(length(candidate)){
      group<-paste(sites$accver[candidate],sites$preferred_enzyme[candidate],
        sites$start[candidate],sep="\r")
      source_order<-if("pattern_source_row"%in%names(sites))
        sites$pattern_source_row[candidate] else candidate
      rank<-order(group,!sites$is_preferred_enzyme[candidate],source_order,
        sites$enzyme[candidate])
      chosen<-candidate[rank][!duplicated(group[rank])]
      keep[candidate]<-FALSE;keep[chosen]<-TRUE
    }
  }
  attrs<-attributes(sites);out<-sites[keep,,drop=FALSE];rownames(out)<-NULL
  for(nm in setdiff(names(attrs),c("names","row.names","class")))attr(out,nm)<-attrs[[nm]];out
}
