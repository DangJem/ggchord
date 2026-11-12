# Mature show-coords reads existing delta; no aligner or analysis is invoked.
# https://mummer4.github.io/manual/manual.html

ggchord_delta_source <- function(text, caller) {
  if (length(text) < 2L || !text[2L] %in% c('NUCMER', 'PROMER'))
    ggchord_stop(caller, ': expected an existing NUCMER or PROMER delta file')
  headers <- list(); alignments <- list(); current <- NULL; active <- NULL
  for (line in if (length(text) > 2L) seq.int(3L, length(text)) else integer()) {
    value <- trimws(text[line])
    if (!nzchar(value)) next
    fields <- strsplit(value, '[[:space:]]+')[[1L]]
    if (startsWith(value, '>')) {
      if (length(fields) != 4L) ggchord_stop(caller, ': malformed delta sequence header at line ', line)
      current <- list(saccver = substring(fields[1L], 2L), qaccver = fields[2L],
        slen = ggchord_format_integer(fields[3L], 'reference length', caller, 1),
        qlen = ggchord_format_integer(fields[4L], 'query length', caller, 1), .source_line = line)
      headers[[length(headers)+1L]] <- current; active <- NULL
    } else if (length(fields) == 7L) {
      if (is.null(current)) ggchord_stop(caller, ': delta alignment precedes its sequence header')
      coordinates <- ggchord_format_integer(fields, 'delta alignment fields', caller)
      active <- length(alignments)+1L
      alignments[[active]] <- c(current, list(sstart = coordinates[1L], send = coordinates[2L],
        qstart = coordinates[3L], qend = coordinates[4L], errors = coordinates[5L],
        similarity_errors = coordinates[6L], stops = coordinates[7L],
        alignment_line = line, delta_script = numeric()))
    } else if (length(fields) == 1L && !is.null(active) && grepl('^-?[0-9]+$', value)) {
      number <- as.numeric(value)
      if (!is.finite(number) || abs(number) > 2^53-1)
        ggchord_stop(caller, ': delta script exceeds the exact numeric range at line ', line)
      alignments[[active]]$delta_script <- c(alignments[[active]]$delta_script, number)
      if (number == 0) active <- NULL
    } else ggchord_stop(caller, ': malformed delta data at line ', line)
  }
  list(mode = text[2L], original_preamble = text[1L], headers = headers, alignments = alignments)
}

#' Read existing delta alignments through an optional command-line decoder
#'
#' Invokes an existing show-coords executable with fixed -THld arguments to
#' read NUCMER/PROMER delta. No alignment, filtering, sorting, merging or
#' missing-character reconstruction is performed. Source FASTA paths in the
#' delta preamble are retained as metadata, never opened or downloaded.
#'
#' Query becomes qaccver and reference becomes saccver. Original one-based
#' inclusive directional endpoints and full lengths are checked against the
#' tool output without shifting or inferring lengths. Original numeric-looking
#' IDs stay strings. Source errors and signed delta scripts (including terminal
#' zero) remain attached to each row, with physical alignment/header lines.
#' No sequence characters are reconstructed from these scripts alone.
#'
#' pident is the decoder's reported percent identity, with its printed precision;
#' it is not inferred from endpoints. For NUCMER this describes nucleotide
#' alignment; for PROMER it describes amino-acid alignment, while endpoint and
#' full-length coordinates still use bp. Similarity, stop percentages and frames
#' remain separate native fields. qspan/sspan are source interval widths, not
#' alignment-column counts. No BLAST length or read-support field is invented.
#'
#' @param file Optional local path vector or glob patterns, before files.
#' @param files Optional additional delta paths or patterns.
#' @param executable Existing show-coords executable name or path, default
#'   "show-coords". It must already be installed; no environment is modified.
#' @param query_ids,subject_ids Optional explicit named source-ID maps,
#'   following [read_chain()]. Shared final IDs need explicit declaration on
#'   both sides and consistent full lengths. Defaults never auto-prefix IDs.
#' @param source_file Add .source_file to ribbon rows, default FALSE.
#' @return A list with seq_data, ribbon_data and report. The report retains
#'   original text/headers, native decoder output/messages, executable and
#'   exact arguments. No executable version is guessed when the tool cannot
#'   report one; record the prepared environment's version separately. All
#'   source rows remain distinct; no source HSP IDs or observation counts are
#'   manufactured. Parsing is explicit data preparation, never a draw-time step.
#' @seealso [read_chain()], [read_maf()], [geom_link_ribbon()]
#' @export
#' @examples
#' \dontrun{
#' # Requires an existing show-coords executable and the source input bundle.
#' observed <- read_delta("examples/genome-alignments/nucmer.delta")
#' head(observed$ribbon_data)
#' }
read_delta <- function(file = NULL, files = NULL, executable = 'show-coords',
                       query_ids = NULL, subject_ids = NULL, source_file = FALSE) {
  caller <- 'read_delta()'
  if (!is.character(executable) || length(executable) != 1L || is.na(executable) || !nzchar(executable))
    ggchord_stop(caller, ': executable must name one existing show-coords program')
  program <- if (file.exists(executable)) file.path(normalizePath(dirname(executable)), basename(executable)) else unname(Sys.which(executable))
  if (!nzchar(program)) ggchord_stop(caller, ': optional show-coords executable not found; select an existing installation or supply a standard pair table')
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop(caller, ': source_file must be TRUE or FALSE')
  paths <- ggchord_resolve_files(file, files, 'read_delta')
  text <- stats::setNames(lapply(paths, ggchord_read_text, caller = caller), paths)
  parsed <- lapply(text, ggchord_delta_source, caller = caller)
  headers <- unlist(lapply(parsed, `[[`, 'headers'), recursive = FALSE)
  maps <- ggchord_alignment_identity_maps(
    unique(vapply(headers, `[[`, character(1), 'qaccver')),
    unique(vapply(headers, `[[`, character(1), 'saccver')), query_ids, subject_ids, caller)
  temporary <- tempfile('ggchord-delta-'); dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  links <- list(); dictionary <- list(); reports <- list()
  for (i in seq_along(paths)) {
    input <- parsed[[i]]; output <- file.path(temporary, paste0(i, '.tsv'))
    log <- file.path(temporary, paste0(i, '.log'))
    arguments <- c('-THld', shQuote(normalizePath(paths[i])))
    status <- suppressWarnings(system2(program, arguments, stdout = output, stderr = log))
    messages <- readLines(log, warn = FALSE)
    if (!identical(status, 0L)) ggchord_stop(caller, ': optional decoder failed: ', paste(messages, collapse = '\n'))
    native_text <- readLines(output, warn = FALSE)
    columns <- c('sstart', 'send', 'qstart', 'qend', 'sspan', 'qspan', 'pident',
      if (input$mode == 'PROMER') c('similarity_percent', 'stop_percent'),
      'slen', 'qlen', 'subject_frame', 'query_frame', 'saccver', 'qaccver')
    if (length(native_text)) {
      native <- utils::read.delim(text = paste(native_text, collapse = '\n'), header = FALSE,
        quote = '', comment.char = '', colClasses = 'character', na.strings = character(), check.names = FALSE)
      if (ncol(native) != length(columns)) ggchord_stop(caller, ': unexpected show-coords output columns; require -THld compatible output')
      names(native) <- columns
      for (name in setdiff(columns, c('qaccver', 'saccver'))) {
        native[[name]] <- suppressWarnings(as.numeric(native[[name]]))
        if (any(!is.finite(native[[name]]))) ggchord_stop(caller, ': decoder returned nonfinite numeric ', name)
      }
    } else native <- stats::setNames(as.data.frame(replicate(length(columns), character(), simplify = FALSE)), columns)
    if (nrow(native) != length(input$alignments)) ggchord_stop(caller,
      ': decoder row count differs from source alignments; no filtered output can be used')
    for (header in input$headers) dictionary[[length(dictionary)+1L]] <- data.frame(
      accver = c(maps$query[[header$qaccver]], maps$subject[[header$saccver]]),
      length = c(header$qlen, header$slen))
    if (nrow(native)) {
      for (name in c('qaccver', 'saccver', 'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen')) {
        original <- vapply(input$alignments, `[[`, if (name %in% c('qaccver', 'saccver')) character(1) else numeric(1), name)
        if (!identical(native[[name]], original)) ggchord_stop(caller,
          ': decoder reordered or altered source endpoints/identifiers; cannot associate source evidence safely')
      }
      if (any(pmin(native$qstart, native$qend) < 1 | pmax(native$qstart, native$qend) > native$qlen |
          pmin(native$sstart, native$send) < 1 | pmax(native$sstart, native$send) > native$slen) ||
          any(native$pident < 0 | native$pident > 100)) ggchord_stop(caller, ': invalid decoded bounds or percent identity')
      out <- native
      out$q_source_accver <- native$qaccver; out$s_source_accver <- native$saccver
      out$qaccver <- unname(maps$query[native$qaccver]); out$saccver <- unname(maps$subject[native$saccver])
      out$.source_input <- i; out$source_alignment <- seq_len(nrow(out))
      out$.source_line <- vapply(input$alignments, `[[`, integer(1), 'alignment_line')
      out$source_header_line <- vapply(input$alignments, `[[`, integer(1), '.source_line')
      for (name in c('errors', 'similarity_errors', 'stops')) out[[paste0('delta_', name)]] <- vapply(input$alignments, `[[`, numeric(1), name)
      out$delta_script <- I(lapply(input$alignments, `[[`, 'delta_script'))
      out$alignment_unit <- if (input$mode == 'PROMER') 'aa' else 'bp'
      links[[i]] <- ggchord_add_source_file(out, paths[i], source_file)
    }
    reports[[i]] <- list(file = paths[i], mode = input$mode, source = input,
      native_output = native_text, messages = messages, arguments = c('-THld', normalizePath(paths[i])))
  }
  list(seq_data = ggchord_alignment_dictionary(dictionary, 'bp', caller),
    ribbon_data = ggchord_rbind_fill(links) %||% data.frame(qaccver = character(),
      saccver = character(), qstart = numeric(), qend = numeric(), sstart = numeric(), send = numeric()),
    report = list(raw_files = text, records = reports,
      decoder = list(executable = program, version = NA_character_,
        version_note = 'No version is inferred from a command path; record prepared environment provenance separately'),
      explicit_id_maps = list(query = query_ids, subject = subject_ids),
      coordinates = 'original 1-based-inclusive bp; PROMER identity is an amino-acid alignment metric'))
}
