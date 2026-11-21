# Source-only bridges; mature parsers supply coordinates, never an aligner.
ggchord_genome_alignment_script <- function() r"---(
import json,sys,gzip
try:
    import Bio
    from Bio import Align
except ImportError:
    raise SystemExit('The selected Python needs optional biopython with Bio.Align MAF/chain support. Select an existing environment; no package is installed automatically.')
def preserve(x):
    if x is None or isinstance(x,(str,int,float,bool)): return x
    if isinstance(x,dict): return {str(k):preserve(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [preserve(v) for v in x]
    if hasattr(x,'item'): return preserve(x.item())
    if hasattr(x,'id') and hasattr(x,'seq'):
        return {'source_class':type(x).__name__,'id':x.id,'length':len(x.seq),'annotations':preserve(x.annotations)}
    return {'source_class':type(x).__name__,'representation':str(x)}
with open(sys.argv[1],encoding='utf-8') as handle: config=json.load(handle)
files=config['files']
if isinstance(files,str): files=[files]
parsed=[]
for path in files:
    blocks=[]
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path,'rt',encoding='utf-8-sig') as source:
        for number,alignment in enumerate(Align.parse(source,config['format']),1):
            rows=[]
            for index,record in enumerate(alignment.sequences):
                row={'accver':record.id,'length':len(record.seq),'annotations':preserve(record.annotations)}
                if config['format']=='maf': row['symbols']=alignment[index]
                rows.append(row)
            blocks.append({'source_block':number,'sequences':rows,
                'coordinates':alignment.coordinates.tolist(),
                'annotations':preserve(getattr(alignment,'annotations',{})),
                'score':getattr(alignment,'score',None)})
    parsed.append({'file':path,'blocks':blocks})
with open(sys.argv[2],'w',encoding='utf-8') as handle:
    json.dump({'decoder':{'package':'biopython','version':Bio.__version__},'files':parsed},handle,allow_nan=False)
)---"

ggchord_alignment_source <- function(file, files, python, format, caller) {
  paths <- ggchord_resolve_files(file, files, sub('()', '', caller, fixed = TRUE))
  parsed <- ggchord_parse_external_source(ggchord_genome_alignment_script(),
    list(files = unname(normalizePath(paths)), format = format), python, caller)
  list(paths = paths, parsed = parsed$data, messages = parsed$messages,
    raw_files = setNames(lapply(paths, ggchord_read_text, caller = caller), paths))
}

ggchord_alignment_dictionary <- function(rows, sequence_unit, caller) {
  if (!length(rows)) return(data.frame(accver = character(), length = numeric(), sequence_unit = character()))
  out <- do.call(rbind, rows)
  if (anyNA(out$accver) || any(!nzchar(out$accver))) ggchord_stop(caller, ': source IDs must be nonempty')
  if (any(!is.finite(out$length) | out$length < 1 | out$length > 2^53-1 | out$length != floor(out$length)))
    ggchord_stop(caller, ': source lengths must be exact positive integers')
  for (id in unique(out$accver)) if (length(unique(out$length[out$accver == id])) != 1L)
    ggchord_stop(caller, ': conflicting full lengths for accession ', id, '; supply distinct source identifiers explicitly')
  out <- out[!duplicated(out$accver), , drop = FALSE]
  out$sequence_unit <- sequence_unit; rownames(out) <- NULL
  out
}

ggchord_alignment_report <- function(source, convention) list(
  decoder = source$parsed$decoder, messages = source$messages,
  parsed_files = source$parsed$files, raw_files = source$raw_files,
  input_coordinates = convention, output_coordinates = '1-based-inclusive',
  interpretation = 'existing alignment hypotheses, not newly validated homology; source block indices are provenance, not independent observation IDs')

#' Read true alignment characters and relationships from MAF
#'
#' Keep block context and real character columns. A block's first component
#' is its query anchor. An empty anchor yields no derived relation.
#' Select a single source file and block using `source_input` and
#' `source_block` when displaying a character alignment.
#'
#' @param file Local MAF paths.
#' @param alphabet Explicit DNA, RNA, or AA.
#' @param python Python executable with Biopython available.
#' @param source_file Whether to include row-level source paths.
#' @return A ggchord_import collection with sequences, links$maf,
#'   alignments$maf, and block provenance.
#' @export
#' @examples
#' \dontrun{
#' f <- system.file("extdata", "tutorial", "ucsc_mm9_chr10.maf", package = "ggchord")
#' x <- read_maf(f, python = "python3")
#' head(x$alignments$maf)
#' }
read_maf <- function(file, alphabet = "DNA", python = "python3", source_file = FALSE) {
  caller <- 'read_maf()'
  alphabet <- match.arg(alphabet, c('DNA', 'RNA', 'AA'))
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop(caller, ': source_file must be TRUE or FALSE')
  source <- ggchord_alignment_source(file, NULL, python, 'maf', caller)
  dictionary <- list(); links <- list(); cells <- list(); empty <- list()
  for (i in seq_along(source$parsed$files)) {
    input <- source$parsed$files[[i]]
    for (block in input$blocks) {
      rows <- block$sequences
      coordinate <- do.call(rbind, lapply(block$coordinates, unlist))
      if (!length(rows)) next
      ids <- vapply(rows, `[[`, character(1), 'accver')
      lengths <- vapply(rows, `[[`, numeric(1), 'length')
      if (anyDuplicated(ids)) ggchord_stop(caller, ': repeated component IDs within one block cannot identify a unique source row')
      dictionary[[length(dictionary)+1L]] <- data.frame(accver = ids, length = lengths)
      first <- coordinate[, 1L]; last <- coordinate[, ncol(coordinate)]
      plus <- last > first
      start <- ifelse(plus, first+1, first); end <- ifelse(plus, last, last+1)
      nonempty <- first != last
      if (any(!nonempty)) empty[[length(empty)+1L]] <- data.frame(.source_input = i,
        source_block = block$source_block, source_component = which(!nonempty), accver = ids[!nonempty],
        reason = 'zero-length source component; no placed characters or anchor relation')
      if (!any(nonempty)) next
      if (any((pmin(start, end) < 1 | pmax(start, end) > lengths)[nonempty]))
        ggchord_stop(caller, ': decoded source interval exceeds its declared full length')
      strings <- data.frame(accver = ids[nonempty], sequence = vapply(rows, `[[`, character(1), 'symbols')[nonempty])
      positions <- data.frame(accver = ids[nonempty], source_accver = ids[nonempty],
        start = start[nonempty], strand = ifelse(plus[nonempty], '+', '-'))
      aligned <- as_ggchord_alignment(strings, alphabet, source_data = positions)
      aligned$source_length <- lengths[match(aligned$accver, ids)]
      aligned$source_component <- match(aligned$accver, ids)
      aligned$source_input <- rep(i, nrow(aligned)); aligned$source_block <- rep(block$source_block, nrow(aligned))
      cells[[length(cells)+1L]] <- ggchord_add_source_file(aligned, source$paths[i], source_file)
      other <- which(nonempty & seq_along(rows) != 1L)
      if (!nonempty[1L] || !length(other)) next
      out <- data.frame(qaccver = ids[1L], saccver = ids[other],
        qstart = start[1L], qend = end[1L], sstart = start[other], send = end[other],
        qstrand = if (plus[1L]) "+" else "-", sstrand = ifelse(plus[other], "+", "-"),
        qlen = lengths[1L], slen = lengths[other], .source_input = i,
        source_block = block$source_block, source_component = other, stringsAsFactors = FALSE)
      links[[length(links)+1L]] <- ggchord_add_source_file(out, source$paths[i], source_file)
    }
  }
  report <- ggchord_alignment_report(source, 'MAF forward 0-based boundary coordinates from Bio.Align')
  report$unplaced_components <- ggchord_rbind_fill(empty) %||% data.frame()
  alignment <- ggchord_rbind_fill(cells) %||% data.frame()
  # rbind must not retain only the first block's single-alignment provenance.
  attr(alignment, 'ggchord_conversion') <- list(format = 'MAF', alphabet = alphabet,
    files = source$paths, context_columns = c('source_input', 'source_block'),
    decoder = source$parsed$decoder,
    coordinates = 'per-block alignment columns and nullable one-based source positions')
  ggchord_import_collection(
    ggchord_alignment_dictionary(dictionary, if (alphabet == "AA") "aa" else "bp", caller),
    links = list(maf = ggchord_canonical_read_links(links, caller)),
    alignments = list(maf = alignment), report = report)
}

ggchord_alignment_rename <- function(ids, map, caller) {
  if (is.null(map)) return(ids)
  if (!is.character(map) || is.null(names(map)) || anyNA(map) || any(!nzchar(map)) ||
      anyNA(names(map)) || any(!nzchar(names(map))) || anyDuplicated(names(map)) || any(!names(map) %in% ids))
    ggchord_stop(caller, ': ID maps must be named character vectors of existing original accessions')
  at <- match(ids, names(map)); ids[!is.na(at)] <- unname(map[at[!is.na(at)]])
  ids
}

ggchord_alignment_identity_maps <- function(query, subject, query_ids, subject_ids, caller) {
  q <- setNames(ggchord_alignment_rename(query, query_ids, caller), query)
  s <- setNames(ggchord_alignment_rename(subject, subject_ids, caller), subject)
  overlap <- intersect(unname(q), unname(s))
  if (length(overlap) && (any(!names(q)[q %in% overlap] %in% names(query_ids)) ||
      any(!names(s)[s %in% overlap] %in% names(subject_ids)))) ggchord_stop(caller,
        ': query and target IDs collide; explicitly declare query_ids/subject_ids for distinct assemblies or shared sources')
  list(query = q, subject = s)
}

#' Read block relationships from a chain file
#'
#' Use an existing Biopython environment to decode chain block coordinates.
#' Store lengths as reference information and retain the original chain
#' fields and directional alignment evidence. Chain scores preserve their
#' source meaning.
#'
#' @param file Local chain paths.
#' @param query_ids,subject_ids Optional named source ID translations.
#'   Renamed IDs must refer to distinct, consistent reference sequences.
#' @param sequence_unit Native coordinate unit.
#' @param python Python executable with Biopython available.
#' @param source_file Whether to include row-level source paths.
#' @return A ggchord_import collection with sequences, links$chain, and the
#'   source and decoder report.
#' @export
#' @examples
#' \dontrun{
#' f <- system.file("extdata", "tutorial", "dna_rna.chain", package = "ggchord")
#' x <- read_chain(f, python = "python3")
#' head(x$links$chain)
#' }
read_chain <- function(file, query_ids = NULL, subject_ids = NULL,
                       sequence_unit = "bp", python = "python3", source_file = FALSE) {
  caller <- "read_chain()"
  sequence_unit <- match.arg(sequence_unit, c("bp", "aa"))
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop(caller, ': source_file must be TRUE or FALSE')
  source <- ggchord_alignment_source(file, NULL, python, 'chain', caller)
  blocks <- unlist(lapply(source$parsed$files, `[[`, 'blocks'), recursive = FALSE)
  q_original <- unique(vapply(blocks, function(b) b$sequences[[2L]]$accver, character(1)))
  s_original <- unique(vapply(blocks, function(b) b$sequences[[1L]]$accver, character(1)))
  maps <- ggchord_alignment_identity_maps(q_original, s_original, query_ids, subject_ids, caller)
  q_names <- maps$query; s_names <- maps$subject
  dictionary <- list(); links <- list()
  for (i in seq_along(source$parsed$files)) for (block in source$parsed$files[[i]]$blocks) {
    if (length(block$sequences) != 2L) ggchord_stop(caller, ': chain decoder must return exactly two source sequences')
    coordinate <- do.call(rbind, lapply(block$coordinates, unlist))
    first <- coordinate[, -ncol(coordinate), drop = FALSE]
    last <- coordinate[, -1L, drop = FALSE]
    keep <- which(first[1L, ] != last[1L, ] & first[2L, ] != last[2L, ])
    original <- vapply(block$sequences, `[[`, character(1), 'accver')
    ids <- c(s_names[[original[1L]]], q_names[[original[2L]]])
    lengths <- vapply(block$sequences, `[[`, numeric(1), 'length')
    dictionary[[length(dictionary)+1L]] <- data.frame(accver = ids, length = lengths)
    if (!length(keep)) next
    start <- ifelse(last > first, first+1, first)
    end <- ifelse(last > first, last, last+1)
    if (any(pmin(start[, keep], end[, keep]) < 1) ||
        any(pmax(start[1L, keep], end[1L, keep]) > lengths[1L]) ||
        any(pmax(start[2L, keep], end[2L, keep]) > lengths[2L]))
      ggchord_stop(caller, ': chain block lies outside its declared source length')
    out <- data.frame(qaccver = ids[2L], saccver = ids[1L],
      qstart = start[2L, keep], qend = end[2L, keep],
      sstart = start[1L, keep], send = end[1L, keep], qlen = lengths[2L], slen = lengths[1L],
      qstrand = ifelse(last[2L, keep] > first[2L, keep], "+", "-"),
      sstrand = ifelse(last[1L, keep] > first[1L, keep], "+", "-"),
      q_source_accver = original[2L], s_source_accver = original[1L],
      .source_input = i, source_block = block$source_block, source_segment = keep, stringsAsFactors = FALSE)
    links[[length(links)+1L]] <- ggchord_add_source_file(out, source$paths[i], source_file)
  }
  report <- ggchord_alignment_report(source, 'chain 0-based boundary coordinates from Bio.Align')
  report$explicit_id_maps <- list(query = query_ids, subject = subject_ids)
  ggchord_import_collection(ggchord_alignment_dictionary(dictionary, sequence_unit, caller),
    links = list(chain = ggchord_canonical_read_links(links, caller)), report = report)
}
