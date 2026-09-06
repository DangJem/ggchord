# Restriction-site calculation and geometry are intentionally separate.

ggchord_restriction_enzymes <- data.frame(
  enzyme = c("EcoRI", "BamHI", "HindIII", "PstI", "SmaI"),
  motif = c("GAATTC", "GGATCC", "AAGCTT", "CTGCAG", "CCCGGG"),
  cut_top = c(1, 1, 1, 5, 3), cut_bottom = c(5, 5, 5, 1, 3),
  stringsAsFactors = FALSE
)

ggchord_iupac_regex <- function(motif) {
  code <- c(A="A",C="C",G="G",T="T",R="[AG]",Y="[CT]",S="[GC]",
    W="[AT]",K="[GT]",M="[AC]",B="[CGT]",D="[AGT]",H="[ACT]",V="[ACG]",N="[ACGT]")
  chars <- strsplit(toupper(motif), "", fixed=TRUE)[[1L]]
  if (any(!chars %in% names(code))) ggchord_stop("Restriction motifs may contain only IUPAC DNA symbols")
  paste0(unname(code[chars]), collapse="")
}

ggchord_reverse_complement <- function(sequence) {
  code <- c(A="T",C="G",G="C",T="A",R="Y",Y="R",S="S",W="W",
            K="M",M="K",B="V",D="H",H="D",V="B",N="N")
  chars <- strsplit(toupper(sequence), "", fixed=TRUE)[[1L]]
  paste0(rev(unname(code[chars])), collapse="")
}

#' Find restriction-enzyme recognition and cut sites
#'
#' Searches DNA without requiring Biostrings. Built-in definitions are
#' versioned with ggchord; custom definitions can be supplied as a data frame.
#'
#' @param sequence A DNA string, named character vector, or data frame with
#'   `accver` and `sequence`.
#' @param enzymes Built-in enzyme names, or a data frame with `enzyme`, `motif`,
#'   `cut_top`, and `cut_bottom`.
#' @param circular Whether matches may cross the sequence origin.
#' @param min_cuts,max_cuts Optional inclusive cut-count filters, evaluated per
#'   enzyme and sequence before applying `window`.
#' @param window Optional inclusive genomic interval of length two.
#' @return A data frame containing enzyme, recognition, cut and source fields.
#' @details Built-in definitions currently include EcoRI, BamHI, HindIII,
#'   PstI, and SmaI. `cut_top` and `cut_bottom` identify the nucleotide
#'   immediately following each strand's cleavage boundary; `position` is the
#'   top-strand cut used by the plotting layer. The database version is stored
#'   both in `source` and in the `enzyme_database_version` attribute.
#' @export
#' @examples
#' find_restriction_sites("AAAAGAATTCTTTT", enzymes = "EcoRI")
find_restriction_sites <- function(sequence, enzymes=c("EcoRI","BamHI"), circular=TRUE,
                                   min_cuts=NULL, max_cuts=NULL, window=NULL) {
  if (is.data.frame(sequence)) {
    ggchord_require_columns(sequence, c("accver","sequence"), "find_restriction_sites()")
    ids <- as.character(sequence$accver); seqs <- as.character(sequence$sequence)
  } else if (is.character(sequence) && length(sequence)) {
    seqs <- as.character(sequence); ids <- names(sequence)
    if (is.null(ids)) ids <- if (length(seqs)==1L) "sequence" else paste0("sequence_",seq_along(seqs))
    ids[!nzchar(ids)] <- paste0("sequence_",which(!nzchar(ids)))
  } else ggchord_stop("find_restriction_sites(): sequence must be DNA text, a named character vector, or a data frame")
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) ggchord_stop("find_restriction_sites(): sequence IDs must be non-missing and unique")
  seqs <- toupper(gsub("[[:space:]]","",seqs))
  if (any(!grepl("^[ACGTRYSWKMBDHVN]+$",seqs))) ggchord_stop("find_restriction_sites(): invalid IUPAC DNA symbols")
  if (!is.logical(circular) || length(circular)!=1L || is.na(circular)) ggchord_stop("find_restriction_sites(): circular must be TRUE or FALSE")
  custom <- is.data.frame(enzymes)
  if (custom) {
    ggchord_require_columns(enzymes,c("enzyme","motif","cut_top","cut_bottom"),"find_restriction_sites()")
    defs <- enzymes[,c("enzyme","motif","cut_top","cut_bottom"),drop=FALSE]; source <- "custom"
  } else {
    if (!is.character(enzymes) || !length(enzymes)) ggchord_stop("find_restriction_sites(): enzymes must name built-ins or provide definitions")
    unknown <- setdiff(enzymes,ggchord_restriction_enzymes$enzyme)
    if (length(unknown)) ggchord_stop("Unknown built-in restriction enzyme(s): ",paste(unknown,collapse=", "))
    defs <- ggchord_restriction_enzymes[match(enzymes,ggchord_restriction_enzymes$enzyme),,drop=FALSE]
    source <- "ggchord-builtin-2026.09"
  }
  defs$enzyme <- as.character(defs$enzyme); defs$motif <- toupper(as.character(defs$motif))
  if (anyNA(defs$enzyme)||any(!nzchar(defs$enzyme))||anyDuplicated(defs$enzyme)||
      anyNA(defs$motif)||any(!nzchar(defs$motif))||!is.numeric(defs$cut_top)||
      !is.numeric(defs$cut_bottom)||any(!is.finite(defs$cut_top))||any(!is.finite(defs$cut_bottom)))
    ggchord_stop("find_restriction_sites(): invalid enzyme definitions")
  for (item in list(min_cuts=min_cuts,max_cuts=max_cuts)) {
    if (!is.null(item) && (!is.numeric(item)||length(item)!=1L||!is.finite(item)||item<0))
      ggchord_stop("find_restriction_sites(): cut-count filters must be non-negative finite numbers")
  }
  if (!is.null(window) && (!is.numeric(window)||length(window)!=2L||any(!is.finite(window))||window[1]>window[2]))
    ggchord_stop("find_restriction_sites(): window must be two increasing finite coordinates")

  rows <- list()
  for (s in seq_along(seqs)) for (d in seq_len(nrow(defs))) {
    motif <- defs$motif[d]; m <- nchar(motif); target <- seqs[s]
    search <- if (circular && m>1L) paste0(target,substr(target,1L,m-1L)) else target
    strands <- c("+", if (ggchord_reverse_complement(motif)!=motif) "-")
    for (site_strand in strands) {
      search_motif <- if (site_strand=="+") motif else ggchord_reverse_complement(motif)
      starts <- gregexpr(paste0("(?=",ggchord_iupac_regex(search_motif),")"),search,perl=TRUE)[[1L]]
      if (length(starts)==1L && starts[1L] < 0L) next
      for (start in starts[starts<=nchar(target)]) {
      finish <- ((start+m-2L)%%nchar(target))+1L
      offsets <- if (site_strand=="+") c(defs$cut_top[d],defs$cut_bottom[d]) else c(m-defs$cut_bottom[d],m-defs$cut_top[d])
      top <- ((start-1L+offsets[1])%%nchar(target))+1L
      bottom <- ((start-1L+offsets[2])%%nchar(target))+1L
      end_type <- if (defs$cut_top[d]==defs$cut_bottom[d]) "blunt" else if (defs$cut_top[d]<defs$cut_bottom[d]) "5_prime_overhang" else "3_prime_overhang"
      rows[[length(rows)+1L]] <- data.frame(accver=ids[s],enzyme=defs$enzyme[d],motif=motif,
        recognition_sequence=substr(search,start,start+m-1L),start=start,end=finish,strand=site_strand,
        cut_top=top,cut_bottom=bottom,position=top,end_type=end_type,source=source,stringsAsFactors=FALSE)
      }
    }
  }
  template <- data.frame(accver=character(),enzyme=character(),motif=character(),recognition_sequence=character(),start=integer(),end=integer(),strand=character(),cut_top=integer(),cut_bottom=integer(),position=integer(),end_type=character(),source=character(),stringsAsFactors=FALSE)
  out <- if (length(rows)) do.call(rbind,rows) else template
  if (nrow(out) && (!is.null(min_cuts)||!is.null(max_cuts))) {
    counts <- table(paste(out$accver,out$enzyme,sep="\r"))
    n <- unname(counts[paste(out$accver,out$enzyme,sep="\r")])
    keep_min <- if (is.null(min_cuts)) rep(TRUE,length(n)) else n>=min_cuts
    keep_max <- if (is.null(max_cuts)) rep(TRUE,length(n)) else n<=max_cuts
    keep <- keep_min & keep_max
    out <- out[keep,,drop=FALSE]
  }
  if (nrow(out) && !is.null(window)) out <- out[out$position>=window[1]&out$position<=window[2],,drop=FALSE]
  rownames(out) <- NULL
  attr(out,"enzyme_database_version") <- if (custom) "custom" else "ggchord-builtin-2026.09"
  out
}

GeomRestrictionSite <- ggplot2::ggproto("GeomRestrictionSite",ggplot2::Geom,
  required_aes=c("x","y"),
  default_aes=ggplot2::aes(xend=NA_real_,yend=NA_real_,label=NA_character_,.component=NA_character_,group=NA_integer_,colour="#B42318",alpha=1,linewidth=0.35,linetype=1,size=2.6,angle=0,hjust=0.5,vjust=0.5,family="",fontface=1,lineheight=1.2),
  draw_key=ggplot2::draw_key_path,
  draw_panel=function(data,panel_params,coord,na.rm=FALSE) {
    paths <- data[data$.component=="path",,drop=FALSE]; labels <- data[data$.component=="label",,drop=FALSE]; grobs <- list()
    if (nrow(paths)) grobs[[length(grobs)+1L]] <- ggplot2::GeomPath$draw_panel(paths,panel_params,coord,lineend="round",linejoin="round",na.rm=na.rm)
    if (nrow(labels)) grobs[[length(grobs)+1L]] <- ggplot2::GeomText$draw_panel(labels,panel_params,coord,parse=FALSE,check_overlap=FALSE,na.rm=na.rm)
    do.call(grid::grobTree,grobs)
  })

#' Add restriction sites to a sequence track
#'
#' Draws cut ticks and optional labels. Nearby labels share a deterministic
#' trunk while each cut remains a separate observation.
#' @param mapping,data Standard layer inputs. Data require `accver`, `position`, and `enzyme`.
#' @param label Whether to draw enzyme names.
#' @param tick_length,label_offset Coordinate-unit offsets from the sequence.
#' @param branch_threshold Maximum genomic fraction for labels to share a trunk.
#' @param colour,linewidth,label_size Fixed appearance.
#' @inheritParams geom_seq
#' @return A ggplot2 layer.
#' @export
geom_restriction_site <- function(mapping=NULL,data=NULL,label=TRUE,tick_length=0.045,
    label_offset=0.16,branch_threshold=0.02,colour="#B42318",linewidth=0.35,
    label_size=2.6,position="identity",show.legend=FALSE,inherit.aes=FALSE,...) {
  if (is.null(data)) ggchord_stop("geom_restriction_site(): supply restriction-site data")
  vals <- c(tick_length,label_offset)
  if (!is.numeric(vals)||any(!is.finite(vals))||any(vals<0)) ggchord_stop("geom_restriction_site(): offsets must be finite non-negative numbers")
  if (!is.numeric(branch_threshold)||length(branch_threshold)!=1L||!is.finite(branch_threshold)||branch_threshold<0||branch_threshold>0.5) ggchord_stop("geom_restriction_site(): branch_threshold must be in [0, 0.5]")
  lyr <- ggplot2::layer(data=data.frame(x=numeric(),y=numeric()),
    mapping=ggplot2::aes(x=x,y=y,group=group,label=label,
                         .component=I(.component)),
    stat="identity",geom=GeomRestrictionSite,position=position,
    show.legend=show.legend,inherit.aes=inherit.aes,check.aes=FALSE,
    check.param=FALSE,
    params=c(list(na.rm=FALSE,colour=colour,linewidth=linewidth,size=label_size),list(...)))
  lyr$ggchord_type <- "restriction_site"
  lyr$ggchord_params <- list(type="restriction_site",label=label,tick_length=tick_length,label_offset=label_offset,branch_threshold=branch_threshold)
  ggchord_capture_layer_input(lyr,data,mapping,c("accver","position","enzyme"))
}

ggchord_restriction_geometry <- function(data,params,layout,seq_data) {
  ggchord_require_columns(data,c("accver","position","enzyme"),"geom_restriction_site()")
  lens <- stats::setNames(seq_data$length,seq_data$accver)
  unknown <- setdiff(unique(as.character(data$accver)),names(lens))
  if (length(unknown)) ggchord_stop("geom_restriction_site(): unknown accver: ",paste(unknown,collapse=", "))
  if (!is.numeric(data$position)||any(!is.finite(data$position))||any(data$position<1|data$position>lens[as.character(data$accver)])) ggchord_stop("geom_restriction_site(): positions must be finite and within sequence lengths")
  output <- list(); gid <- 0L
  for (id in unique(as.character(data$accver))) {
    idx <- which(as.character(data$accver)==id); idx <- idx[order(data$position[idx],data$enzyme[idx])]
    arc <- layout$seq_arcs[[id]]; n <- nrow(arc); frac <- data$position[idx]/lens[id]
    point_at <- function(f,offset=0) {
      k <- pmax(1L,pmin(n,round(1+f*(n-1)))); kp <- pmax(1L,k-1L); kn <- pmin(n,k+1L)
      tx <- arc$x[kn]-arc$x[kp]; ty <- arc$y[kn]-arc$y[kp]; norm <- sqrt(tx^2+ty^2); tx <- tx/norm; ty <- ty/norm
      nx <- -ty; ny <- tx; flip <- nx*arc$x[k]+ny*arc$y[k]<0; nx[flip] <- -nx[flip]; ny[flip] <- -ny[flip]
      data.frame(x=arc$x[k]+offset*nx,y=arc$y[k]+offset*ny)
    }
    split_id <- cumsum(c(TRUE,diff(frac)>params$branch_threshold))
    for (cluster in split(idx,split_id)) {
      cf <- data$position[cluster]/lens[id]; trunk <- point_at(mean(range(cf)),params$label_offset*0.55)
      base <- point_at(cf,0); tip <- point_at(cf,params$tick_length)
      shared <- isTRUE(params$label) && length(cluster)>1L
      if (shared) {
        gid <- gid+1L
        trunk_path <- rbind(point_at(mean(range(cf)),params$tick_length),trunk)
        trunk_path$.component <- "path"; trunk_path$group <- gid
        # Attribute the shared connector to the first deterministic member so
        # layer-input enrichment preserves it without inventing a data point.
        trunk_path$label <- NA_character_; trunk_path$source_row <- cluster[1L]
        output[[length(output)+1L]] <- trunk_path
      }
      for (j in seq_along(cluster)) {
        gid <- gid+1L; path <- rbind(base[j,],tip[j,])
        path$.component <- "path"; path$group <- gid
        path$label <- NA_character_; path$source_row <- cluster[j]
        output[[length(output)+1L]] <- path
        if (isTRUE(params$label)) {
          lf <- max(0,min(1,cf[j]+(j-(length(cluster)+1)/2)*params$branch_threshold*0.7)); lp <- point_at(lf,params$label_offset)
          gid <- gid+1L
          branch <- if (shared) rbind(trunk,lp) else rbind(tip[j,],lp)
          branch$.component <- "path"; branch$group <- gid
          branch$label <- NA_character_; branch$source_row <- cluster[j]
          output[[length(output)+1L]] <- branch
          output[[length(output)+1L]] <- data.frame(x=lp$x,y=lp$y,label=as.character(data$enzyme[cluster[j]]),.component="label",group=gid,source_row=cluster[j],stringsAsFactors=FALSE)
        }
      }
    }
  }
  if (!length(output)) {
    return(data.frame(x=numeric(),y=numeric(),label=character(),
      .component=character(),group=integer(),source_row=integer()))
  }
  ggchord_rbind_fill(output)
}
