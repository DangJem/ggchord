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

#' Read an existing NUCmer or PROmer delta file
#'
#' Decode ranges with an already installed show-coords executable.
#' Retain source blocks, the delta edit script and original position
#' evidence. PROmer locations remain in bp and its reported identity
#' continues to describe a protein alignment.
#'
#' @param file Local delta paths.
#' @param executable Path or command for an existing show-coords.
#' @param query_ids,subject_ids Optional named mappings to unambiguous
#'   reference IDs.
#' @param source_file Whether to include row-level source paths.
#' @return A ggchord_import collection with sequences, links$delta and a
#'   source/decoder report.
#' @export
#' @examples
#' \dontrun{
#' f <- system.file("extdata", "tutorial", "arabidopsis-first-pair.delta",
#'   package = "ggchord")
#' x <- read_delta(f, executable = "show-coords")
#' head(x$links$delta)
#' }
read_delta <- function(file, executable = "show-coords",
                       query_ids = NULL, subject_ids = NULL, source_file = FALSE) {
  ggchord_local_console_errors()
  caller <- 'read_delta()'
  if (!is.character(executable) || length(executable) != 1L || is.na(executable) || !nzchar(executable))
    ggchord_stop(caller, ': executable must name one existing show-coords program')
  program <- if (file.exists(executable)) file.path(normalizePath(dirname(executable)), basename(executable)) else unname(Sys.which(executable))
  if (!nzchar(program)) ggchord_stop(caller, ': optional show-coords executable not found; select an existing installation or supply a standard pair table')
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop(caller, ': source_file must be TRUE or FALSE')
  paths <- ggchord_resolve_files(file, NULL, 'read_delta')
  text <- setNames(lapply(paths, ggchord_read_text, caller = caller), paths)
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
    } else native <- setNames(as.data.frame(replicate(length(columns), character(), simplify = FALSE)), columns)
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
  relations <- ggchord_canonical_read_links(links, caller, directional_order = TRUE)
  ggchord_import_collection(ggchord_alignment_dictionary(dictionary, 'bp', caller),
    links = list(delta = relations),
    report = list(raw_files = text, records = reports,
      decoder = list(executable = program, version = NA_character_,
        version_note = 'No version is inferred from a command path; record prepared environment provenance separately'),
      explicit_id_maps = list(query = query_ids, subject = subject_ids),
      coordinates = 'original 1-based-inclusive bp; PROMER identity is an amino-acid alignment metric'))
}
