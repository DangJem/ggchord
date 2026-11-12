# Source-only bridges; mature parsers supply coordinates, never an aligner.
ggchord_parse_external_source <- function(script_text, configuration, python, caller) {
  ggchord_conversion_require('jsonlite', caller)
  if (!is.character(python) || length(python) != 1L || is.na(python) || !nzchar(python))
    ggchord_stop(caller, ': python must name one existing interpreter')
  executable <- if (file.exists(python)) file.path(normalizePath(dirname(python)), basename(python)) else unname(Sys.which(python))
  if (!nzchar(executable)) ggchord_stop(caller, ': selected Python interpreter was not found')
  temporary <- tempfile('ggchord-source-'); dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  script <- file.path(temporary, 'parse.py'); config <- file.path(temporary, 'input.json')
  output <- file.path(temporary, 'output.json'); log <- file.path(temporary, 'parser.log')
  writeLines(script_text, script, useBytes = TRUE)
  jsonlite::write_json(configuration, config, auto_unbox = TRUE)
  status <- suppressWarnings(system2(executable, c(shQuote(script), shQuote(config), shQuote(output)), stdout = log, stderr = log))
  messages <- readLines(log, warn = FALSE)
  if (!identical(status, 0L) || !file.exists(output)) ggchord_stop(caller,
    ': optional source parser failed: ', paste(messages, collapse = '\n'))
  list(data = jsonlite::read_json(output, simplifyVector = FALSE), messages = messages)
}

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
    raw_files = stats::setNames(lapply(paths, ggchord_read_text, caller = caller), paths))
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

#' Read existing MAF blocks with exact source characters
#'
#' Optional Biopython Bio.Align parses local MAF through an explicitly selected
#' Python environment. No alignment is computed. Full source IDs, lengths,
#' block scores, coordinate matrices and original text remain traceable. Empty
#' components and quality/context lines remain in source text/parser metadata;
#' no characters or relations are invented for absent components.
#' Zero-length s components are listed in report$unplaced_components. They do
#' not acquire source positions or links; an empty first anchor is not replaced
#' by another component implicitly.
#'
#' Each block's first component is the query anchor, paired with every other
#' nonempty component. No all-pairs clique or Identity is inferred. Source
#' positions come from the decoder's forward zero-based boundary matrix and
#' are converted once to directional one-based inclusive endpoints. Negative
#' rows decrease from their true forward source coordinate, not the offset on
#' the reverse-complement string. Blocks and input files are separate contexts.
#'
#' Alignment cells reuse [as_ggchord_alignment()]. Original symbols and gaps
#' are retained; gaps have missing source positions. Full lengths do not verify
#' source characters independently, so source_verified stays FALSE. Filter
#' .source_input and source_block before column lookup or drawing a character
#' panel: alignment columns restart at one in each block. No blocks are joined,
#' renumbered or promoted to complete source sequences.
#'
#' @param file Optional local path vector or glob patterns, before files.
#' @param files Optional additional MAF paths or patterns.
#' @param alphabet Declared DNA, RNA or AA alphabet, default "DNA".
#' @param python Existing interpreter with optional Biopython, default "python3".
#' @param source_file Add .source_file to alignment/ribbon rows, default FALSE.
#' @return A list with seq_data (declared full lengths), ribbon_data, alignment
#'   cells and report. Source indices and original decoder metadata remain in
#'   the report. No pident or read support is manufactured.
#' @seealso [read_chain()], [map_alignment_position()]
#' @export
#' @examples
#' \dontrun{
#' # Requires an existing Python environment with optional Biopython.
#' blocks <- read_maf("examples/genome-alignments/ucsc_mm9_chr10.maf")
#' local <- blocks$alignment[blocks$alignment$source_block == 1, ]
#' }
read_maf <- function(file = NULL, files = NULL, alphabet = 'DNA',
                     python = 'python3', source_file = FALSE) {
  caller <- 'read_maf()'
  alphabet <- match.arg(alphabet, c('DNA', 'RNA', 'AA'))
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop(caller, ': source_file must be TRUE or FALSE')
  source <- ggchord_alignment_source(file, files, python, 'maf', caller)
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
      aligned$.source_input <- rep(i, nrow(aligned)); aligned$source_block <- rep(block$source_block, nrow(aligned))
      cells[[length(cells)+1L]] <- ggchord_add_source_file(aligned, source$paths[i], source_file)
      other <- which(nonempty & seq_along(rows) != 1L)
      if (!nonempty[1L] || !length(other)) next
      out <- data.frame(qaccver = ids[1L], saccver = ids[other],
        qstart = start[1L], qend = end[1L], sstart = start[other], send = end[other],
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
    files = source$paths, context_columns = c('.source_input', 'source_block'),
    decoder = source$parsed$decoder,
    coordinates = 'per-block alignment columns and nullable one-based source positions')
  list(seq_data = ggchord_alignment_dictionary(dictionary, if (alphabet == 'AA') 'aa' else 'bp', caller),
    ribbon_data = ggchord_rbind_fill(links) %||% data.frame(),
    alignment = alignment,
    report = report)
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
  q <- stats::setNames(ggchord_alignment_rename(query, query_ids, caller), query)
  s <- stats::setNames(ggchord_alignment_rename(subject, subject_ids, caller), subject)
  overlap <- intersect(unname(q), unname(s))
  if (length(overlap) && (any(!names(q)[q %in% overlap] %in% names(query_ids)) ||
      any(!names(s)[s %in% overlap] %in% names(subject_ids)))) ggchord_stop(caller,
        ': query and target IDs collide; explicitly declare query_ids/subject_ids for distinct assemblies or shared sources')
  list(query = q, subject = s)
}

#' Read chain alignment blocks without inventing sequence characters
#'
#' Optional Biopython Bio.Align parses chain coordinates. Each block aligned
#' on both sides becomes one directional interval pair; intervening gaps stay
#' in the source coordinate matrix rather than becoming a false uninterrupted
#' match. Source chain scores and IDs remain in report metadata; neither is
#' an Identity or read-support count. Chain contains no aligned characters.
#'
#' The chain query becomes qaccver and its target becomes saccver. Original
#' IDs stay in q_source_accver/s_source_accver. If both assemblies use the same
#' IDs, supply explicit query_ids/subject_ids rather than implicitly merging
#' distinct molecules or auto-prefixing them. Shared final IDs require explicit
#' declaration on both sides and consistent full lengths. This is data identity
#' mapping, not automatic homology discovery or display-instance creation.
#'
#' @param file Optional local path vector or glob patterns, before files.
#' @param files Optional additional chain paths or patterns.
#' @param query_ids,subject_ids Optional named character vectors mapping
#'   original query/target IDs to explicit source IDs. Unspecified IDs stay
#'   unchanged. Unknown map keys and conflicting full lengths fail.
#' @param sequence_unit Declared coordinate unit, "bp" (default) or "aa".
#'   The format itself cannot establish the biological unit.
#' @param python Existing interpreter with optional Biopython, default "python3".
#' @param source_file Add .source_file to ribbon rows, default FALSE.
#' @return A list with seq_data, ribbon_data and report containing original
#'   text, chain scores/IDs, full lengths and exact zero-based boundary matrices.
#'   Each link retains source_block and source_segment provenance; no source
#'   HSP ID, pident, base mismatch script or character alignment is fabricated.
#' @seealso [read_maf()], [geom_link_ribbon()]
#' @export
#' @examples
#' \dontrun{
#' # Requires an existing Python environment with optional Biopython.
#' chains <- read_chain("examples/genome-alignments/dna_rna.chain")
#' head(chains$ribbon_data)
#' }
read_chain <- function(file = NULL, files = NULL, query_ids = NULL, subject_ids = NULL,
                       sequence_unit = c('bp', 'aa'), python = 'python3', source_file = FALSE) {
  caller <- 'read_chain()'; sequence_unit <- match.arg(sequence_unit)
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop(caller, ': source_file must be TRUE or FALSE')
  source <- ggchord_alignment_source(file, files, python, 'chain', caller)
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
      q_source_accver = original[2L], s_source_accver = original[1L],
      .source_input = i, source_block = block$source_block, source_segment = keep, stringsAsFactors = FALSE)
    links[[length(links)+1L]] <- ggchord_add_source_file(out, source$paths[i], source_file)
  }
  report <- ggchord_alignment_report(source, 'chain 0-based boundary coordinates from Bio.Align')
  report$explicit_id_maps <- list(query = query_ids, subject = subject_ids)
  list(seq_data = ggchord_alignment_dictionary(dictionary, sequence_unit, caller),
    ribbon_data = ggchord_rbind_fill(links) %||% data.frame(), report = report)
}
