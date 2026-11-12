# Source-backed coordinate/character validation, no aligner or package checks.
pkgload::load_all(quiet = TRUE)
python <- Sys.getenv('ALIGNMENT_PYTHON', 'python3')
fail <- function(expr, pattern) {
  error <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(error, 'error'), grepl(pattern, conditionMessage(error)))
}
path <- 'examples/genome-alignments/ucsc_mm9_chr10.maf'
maf <- read_maf(path, python = python, source_file = TRUE)
raw <- readLines(path)
block <- 0L; records <- 0L; negatives <- 0L; gaps <- 0L
# Inspect original standard s rows independently of the Python decoder.
for (line in raw) {
  if (startsWith(line, 'a ')) block <- block+1L
  if (!startsWith(line, 's ')) next
  fields <- strsplit(trimws(line), '[[:space:]]+')[[1L]]
  id <- fields[2L]; start0 <- as.numeric(fields[3L]); width <- as.numeric(fields[4L])
  strand <- fields[5L]; full <- as.numeric(fields[6L]); symbols <- fields[7L]
  local <- maf$alignment[maf$alignment$source_block == block & maf$alignment$accver == id, ]
  stopifnot(identical(paste0(local$symbol, collapse = ''), symbols),
    all(local$source_length == full), all(!local$source_verified),
    sum(!local$is_gap) == width, all(is.na(local$source_position[local$is_gap])),
    all(local$.source_file == path))
  first <- if (strand == '+') start0+1 else full-start0
  expected <- first + (if (strand == '+') 1 else -1)*seq.int(0, width-1)
  stopifnot(identical(local$source_position[!local$is_gap], expected))
  records <- records+1L; negatives <- negatives+as.integer(strand == '-')
  gaps <- gaps+sum(local$is_gap)
}
stopifnot(negatives > 0, gaps > 0,
  nrow(maf$ribbon_data) == records-block,
  identical(maf$report$raw_files[[1L]], raw),
  !any(c('pident', 'read_support') %in% names(maf$ribbon_data)))
# Select one source block before existing column/source inverse lookup.
first <- maf$alignment[maf$alignment$source_block == 1L, ]
id <- unique(first$accver)[1L]
selected <- first[first$accver == id & !first$is_gap, ][1L, ]
lookup <- map_alignment_position(first, id, selected$alignment_column)
stopifnot(lookup$source_position == selected$source_position)
reverse <- map_alignment_position(first, id, selected$source_position, inverse = TRUE)
stopifnot(reverse$alignment_column == selected$alignment_column)
chain_path <- 'examples/genome-alignments/dna_rna.chain'
chain <- read_chain(chain_path, python = python, source_file = TRUE)
# Original first chain: actual reverse transcript, three source blocks.
observed <- chain$ribbon_data[chain$ribbon_data$source_block == 1L, ]
stopifnot(nrow(observed) == 3, all(observed$qaccver == 'NR_046654.1'),
  all(observed$saccver == 'chr3'), all(observed$qstart > observed$qend),
  sum(abs(observed$qend-observed$qstart)+1) == 181,
  all(diff(observed$sstart) > 0),
  all(observed$sstart[-1L]-head(observed$send, -1L) > 1),
  chain$seq_data$length[chain$seq_data$accver == 'chr3'] == 198295559,
  identical(chain$report$raw_files[[1L]], readLines(chain_path)),
  !any(c('pident', 'alignment', 'sequence') %in% names(chain$ribbon_data)))
# Declared assembly identities, not implicit accession renaming.
fixture <- tempfile(fileext = '.chain')
writeLines(c('chain 10 chr 100 + 0 10 chr 200 - 0 10 1', '10'), fixture)
fail(read_chain(fixture, python = python), 'IDs collide')
renamed <- read_chain(fixture, query_ids = c(chr = 'query.chr'),
  subject_ids = c(chr = 'reference.chr'), python = python)
stopifnot(renamed$ribbon_data$qaccver == 'query.chr', renamed$ribbon_data$saccver == 'reference.chr',
  renamed$ribbon_data$q_source_accver == 'chr', renamed$ribbon_data$s_source_accver == 'chr',
  renamed$ribbon_data$qstart == 200, renamed$ribbon_data$qend == 191,
  renamed$ribbon_data$sstart == 1, renamed$ribbon_data$send == 10)
unlink(fixture)
# Empty first anchor stays empty; no implicit anchor replacement.
fixture <- tempfile(fileext = '.maf')
writeLines(c('##maf version=1', 'a score=1',
  's empty 0 0 + 100 ---', 's actual 0 3 + 100 ACG', ''), fixture)
empty_anchor <- read_maf(fixture, python = python)
stopifnot(nrow(empty_anchor$ribbon_data) == 0,
  nrow(empty_anchor$report$unplaced_components) == 1,
  all(empty_anchor$alignment$accver == 'actual'))
unlink(fixture)
# Actual MAF gap columns rendered as characters, not inferred from chain spans.
local <- first[first$alignment_column >= 61 & first$alignment_column <= 100, ]
plot <- ggplot2::ggplot(local, aes(alignment_column, accver, label = symbol)) +
  geom_residue(aes(fill = symbol), size = 3, overflow = 'error') +
  ggplot2::scale_fill_manual(values = c(A = '#B8D8BA', C = '#A6CBE3', G = '#F2D28C',
    T = '#EAA8A8', N = '#D8D8D8', '-' = '#F5F5F5'), guide = 'none') +
  labs(x = 'MAF alignment column', y = NULL, title = 'Existing MAF source characters') +
  ggplot2::theme_minimal()
ggplot2::ggsave('tools/long-term-audit/v0.27.0-maf-characters.png', plot,
  width = 12, height = 2.8, dpi = 110, bg = 'white')
jsonlite::write_json(list(status = 'passed', decoder = maf$report$decoder,
  maf_blocks = block, maf_components = records, negative_components = negatives,
  actual_gap_cells = gaps, anchor_pairs = nrow(maf$ribbon_data), character_cells = nrow(maf$alignment),
  raw_s_rows_checked = TRUE, inverse_lookup = TRUE, empty_anchor_not_replaced = TRUE,
  real_reverse_chain_blocks = nrow(observed), chain_gaps_not_bridged = TRUE,
  explicit_assembly_identity = TRUE, no_fake_identity_or_characters = TRUE,
  sources = list(maf = digest::digest(file = path, algo = 'sha256'),
    chain = digest::digest(file = chain_path, algo = 'sha256'))),
  'tools/long-term-audit/v0.27.0-genome-alignment-results.json', auto_unbox = TRUE, pretty = TRUE)
cat('MAF/chain adapters: passed\n')
