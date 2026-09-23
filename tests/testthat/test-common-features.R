test_that("internal common-feature database is complete and self-contained", {
  db <- ggchord:::ggchord_builtin_common_features()
  expect_equal(nrow(db$features), 1454L)
  expect_equal(nrow(db$segments), 1970L)
  expect_equal(nrow(db$qualifiers), 5216L)
  expect_equal(nrow(db$qualifier_links), 548L)
  expect_identical(db$metadata$source_workbook_sha256,
    "7d9294df1ba89e6d9437f54d2aed5681d86ff4663f098fcc3d866b9ec58db69b")
  expect_identical(db$metadata$database_version,
    "snapgene-8.2.3-export15")
  expect_named(db, c("metadata", "source_releases", "source_sequences",
    "feature_identities", "features", "segments", "qualifiers",
    "qualifier_links", "feature_type_summary", "legacy_feature_records",
    "feature_correspondence", "source_audit", "reference_sequences",
    "reference_features", "reference_primers", "search_indexes"))
  expect_equal(nrow(db$source_releases), 2L)
  expect_equal(nrow(db$source_sequences), 1L)
  expect_equal(db$source_sequences$length_bp, 903397L)
  expect_equal(nchar(db$source_sequences$sequence), 903397L)
  expect_false(any(c("reference_dna_top_strand",
    "reference_dna_feature_5to3", "reference_protein", "segment_ranges",
    "segment_colors", "translated_any") %in% names(db$features)))
  expect_false(any(startsWith(names(db$features), "q_")))
  expect_true(all(c("common_feature_id", "source_feature_id",
    "source_recent_id", "prioritize") %in% names(db$features)))
  expect_true(all(c("start_1based_inclusive", "end_1based_inclusive") %in%
    names(db$segments)))
  expect_false("dna_sequence_top_strand" %in% names(db$segments))
  expect_true(all(c("text", "int", "predef", "bool", "value_display") %in%
    names(db$qualifiers)))
  expect_true(all(c("link_index", "anchor_text", "url") %in%
    names(db$qualifier_links)))
  expect_equal(sum(db$features$detectionMode == "exactProteinMatch",
    na.rm = TRUE), 250L)
  exact_types <- table(db$features$type[
    db$features$detectionMode == "exactProteinMatch"])
  expect_identical(as.integer(exact_types[c("CDS", "sig_peptide")]),
    c(242L, 8L))
  expect_equal(sum(db$features$prioritize %in% TRUE, na.rm = TRUE), 1L)
  expect_equal(db$features$name[db$features$prioritize %in% TRUE],
    "mini-white")
  expect_equal(sum(db$segments$segment_type == "gap"), 40L)
  expect_true(any(!is.na(db$qualifiers$text) &
    !is.na(db$qualifiers$predef)))
  repeated_url <- paste(db$qualifier_links$common_feature_id,
    db$qualifier_links$qualifier_name, db$qualifier_links$value_index,
    db$qualifier_links$url, sep = "\r")
  expect_equal(max(tabulate(match(repeated_url, unique(repeated_url)))), 3L)
  expect_equal(nrow(db$reference_sequences), 13L)
  expect_equal(nrow(db$reference_features), 201L)
  expect_equal(sum(vapply(db$reference_features$segments, nrow, integer(1L))),
    218L)
  expect_equal(sum(lengths(db$reference_features$cleavage_arrows)), 11L)
  expect_equal(nrow(db$reference_primers), 7L)
})

test_that("common-feature migration provenance distinguishes source identities", {
  db <- ggchord:::ggchord_builtin_common_features()
  expect_false(anyDuplicated(db$features$common_feature_id))
  expect_false(any(db$features$common_feature_id ==
    db$features$source_feature_id))
  expect_equal(length(unique(na.omit(
    db$feature_correspondence$old_source_feature_id))), 1273L)
  expect_equal(length(unique(na.omit(
    db$feature_correspondence$new_source_feature_id))), 1454L)
  expect_equal(sum(db$feature_correspondence$relationship == "unresolved"),
    268L)
  expect_false(any(db$feature_correspondence$relationship == "removed"))
  expect_equal(table(db$feature_identities$created_in_release),
    structure(c(448L, 1006L), names = c(
      "snapgene-common-features-8.2.3-export15",
      "snapgene-common-features-legacy-export13")))
  tpa <- subset(db$feature_correspondence,
    old_name == "tPA signal/pro sequence")
  expect_equal(tpa$relationship, rep("split", 2L))
  expect_setequal(tpa$new_type, c("sig_peptide", "propeptide"))
  csy4 <- subset(db$feature_correspondence, old_name == "Csy4 site")
  expect_equal(csy4$relationship, "renamed")
  expect_equal(csy4$new_name, "Csy4 Site")
  beta <- subset(db$feature_correspondence,
    old_type == "polyA_site")
  expect_equal(beta$relationship, rep("reclassified", 2L))
  expect_equal(beta$new_type, rep("polyA_signal", 2L))
})

test_that("built-in matching derives DNA and protein from normalized relations", {
  db <- ggchord:::ggchord_builtin_common_features()
  source <- stats::setNames(db$source_sequences$sequence,
    db$source_sequences$source_sequence_id)
  feature_dna <- function(id) {
    rows <- db$segments[db$segments$common_feature_id == id, , drop = FALSE]
    rows <- rows[order(rows$segment_index), , drop = FALSE]
    rows <- rows[rows$segment_type == "standard", , drop = FALSE]
    paste0(mapply(function(source_id, start, end) {
      substr(source[[source_id]], start, end)
    }, rows$source_sequence_id, rows$start_1based_inclusive,
      rows$end_1based_inclusive, USE.NAMES = FALSE), collapse = "")
  }

  promoter <- db$features[db$features$name == "AmpR promoter" &
    is.na(db$features$detectionMode), , drop = FALSE][1L, ]
  dna_hit <- find_common_features(feature_dna(promoter$common_feature_id),
    database = db, features = promoter$common_feature_id, mode = "dna",
    circular = FALSE, resolve = "all")
  expect_true(any(dna_hit$match_method == "dna_exact"))

  signal <- db$features[db$features$type == "sig_peptide" &
    db$features$detectionMode == "exactProteinMatch", , drop = FALSE][1L, ]
  protein_hit <- find_common_features(feature_dna(signal$common_feature_id),
    database = db, features = signal$common_feature_id, mode = "auto",
    circular = FALSE, resolve = "all")
  expect_true(any(protein_hit$type == "sig_peptide" &
    protein_hit$match_method == "protein_exact"))
})

test_that("common-feature DNA matching handles orientation and origin", {
  db <- data.frame(
    common_feature_id = c("forward", "nondirectional", "bidirectional"),
    name = c("forward", "neutral", "both"),
    type = c("CDS", "misc_feature", "promoter"),
    directionality_label = c("forward", "nondirectional", "bidirectional"),
    sequence = c("ATGAAA", "CCCGGG", "AATTCC"),
    stringsAsFactors = FALSE
  )
  plus <- find_common_features("TTATGAAAGGCCCGGG", database = db,
    circular = FALSE, resolve = "all")
  expect_true(any(plus$common_feature_id == "forward" & plus$strand == "+"))
  expect_true(any(plus$common_feature_id == "nondirectional" & plus$strand == "."))
  reverse <- find_common_features("TTTCAT", database = db,
    circular = FALSE, resolve = "all")
  expect_true(any(reverse$common_feature_id == "forward" & reverse$strand == "-"))
  wrap <- find_common_features("CCAATT", database = db,
    features = "both", circular = TRUE, resolve = "all")
  expect_true(any(wrap$cross_origin & wrap$strand == "+/-"))
})

test_that("common-feature matching preserves segments and gap constraints", {
  features <- data.frame(
    common_feature_id = "multi", name = "multi", type = "misc_feature",
    directionality_label = "forward", reference_protein = NA_character_,
    detectionMode = NA_character_, geneticCode = NA_integer_,
    stringsAsFactors = FALSE
  )
  segments <- data.frame(
    feature_id = "multi", segment_index = 1:3,
    segment_type = c("standard", "gap", "standard"),
    length_bp = c(3L, 2L, 3L),
    dna_sequence_top_strand = c("AAA", NA, "CCC"),
    segment_name = c("left", NA, "right"), translated = FALSE,
    stringsAsFactors = FALSE
  )
  db <- list(metadata = list(database_id = "test", database_version = "1",
    source_workbook_sha256 = "sha"), features = features, segments = segments)
  hit <- find_common_features("GGAAATTCCCGG", database = db,
    circular = FALSE)
  expect_equal(nrow(hit), 1L)
  expect_equal(hit$segment_count, 3L)
  expect_equal(hit$segments[[1]]$segment_type,
    c("standard", "gap", "standard"))
  miss <- find_common_features("GGAAATCCCGG", database = db,
    circular = FALSE)
  expect_equal(nrow(miss), 0L)
})

test_that("protein matching supports exact and approximate six-frame hits", {
  protein <- paste(rep("ACDEFGHIKL", 5), collapse = "")
  # Encode an independent protein-rich target with the standard code.
  codon <- c(A="GCT", C="TGT", D="GAT", E="GAA", F="TTT", G="GGT",
    H="CAT", I="ATT", K="AAA", L="CTG")
  dna <- paste0(unname(codon[strsplit(protein, "", fixed = TRUE)[[1L]]]),
    collapse = "")
  db <- data.frame(common_feature_id = "protein", name = "protein",
    type = "CDS", directionality_label = "forward", sequence = "",
    reference_protein = protein, translated_any = TRUE,
    stringsAsFactors = FALSE)
  exact <- find_common_features(dna, database = db, mode = "protein")
  expect_true(any(exact$match_method == "protein_exact"))
  mutated <- paste0(substr(dna, 1, 29), "GCT", substr(dna, 33, nchar(dna)))
  approximate <- find_common_features(mutated, database = db, mode = "protein")
  expect_true(any(approximate$match_method %in%
    c("protein_exact", "protein_approx")))
})

test_that("reference FASTA-derived plasmids receive exact annotations", {
  data(plasmid_example_pBR322)
  data(plasmid_example_pUC19)
  data(plasmid_example_pBluescript_II_SK_plus)
  br <- find_common_features(plasmid_example_pBR322)
  uc <- find_common_features(plasmid_example_pUC19)
  blue <- find_common_features(plasmid_example_pBluescript_II_SK_plus)
  expect_true(all(c("AmpR", "rop", "bom", "ori") %in% br$anno))
  expect_true(any(br$anno %in% c("TetR", "TcR")))
  expect_true(all(c("AmpR", "ori", "MCS", "lacZα") %in% uc$anno))
  expect_true(all(c("AmpR", "f1 ori", "MCS", "T7 promoter", "T3 promoter",
    "KS primer", "SK primer") %in% blue$anno))
  expect_true(all(br$match_method == "reference_exact"))
  expect_true(all(uc$match_method == "reference_exact"))
  expect_true(all(blue$match_method == "reference_exact"))
  observed <- blue[match(c("lacZα", "AmpR"), blue$anno),
    c("start", "end", "strand")]
  rownames(observed) <- NULL
  expect_equal(observed, data.frame(start = c(241L, 1973L),
    end = c(816L, 2833L), strand = c("-", "-")))
})

test_that("auto matching follows stored exact-protein detection mode", {
  protein <- paste(rep("ACDEFGHIKL", 5), collapse = "")
  codon <- c(A="GCT", C="TGT", D="GAT", E="GAA", F="TTT", G="GGT",
    H="CAT", I="ATT", K="AAA", L="CTG")
  dna <- paste0(unname(codon[strsplit(protein, "", fixed = TRUE)[[1L]]]),
    collapse = "")
  db <- data.frame(common_feature_id = "protein", name = "protein",
    type = "CDS", directionality_label = "forward", sequence = dna,
    reference_protein = protein, translated_any = TRUE,
    detectionMode = "exactProteinMatch", stringsAsFactors = FALSE)
  hit <- find_common_features(dna, database = db, mode = "auto")
  expect_equal(hit$match_method, "protein_exact")
})

test_that("plasmid data objects exactly reproduce every reference FASTA", {
  paths <- sort(list.files(
    testthat::test_path("..", "..", "examples", "plasmid"),
    pattern = "\\.fna$", full.names = TRUE
  ))
  skip_if_not(all(file.exists(paths)))
  labels <- tools::file_path_sans_ext(basename(paths))
  suffixes <- gsub("\\(\\+\\)", "_plus", labels)
  suffixes <- gsub("[^[:alnum:]]+", "_", suffixes)
  suffixes <- gsub("^_+|_+$", "", suffixes)
  objects <- paste0("plasmid_example_", suffixes)
  expect_length(objects, 13L)
  expect_false(anyDuplicated(objects))
  for (i in seq_along(paths)) {
    env <- new.env(parent = emptyenv())
    utils::data(list = objects[i], envir = env)
    object <- env[[objects[i]]]
    lines <- readLines(paths[[i]], warn = FALSE)
    sequence <- paste0(lines[!grepl("^>", lines)], collapse = "")
    expect_equal(nrow(object), 1L)
    expect_identical(names(object), c("accver", "label", "length", "sequence"))
    expect_identical(object$accver, labels[i])
    expect_identical(object$label, labels[i])
    expect_identical(object$sequence, toupper(sequence))
    expect_identical(object$length, nchar(sequence))
  }
  expect_false(any(c("plasmid_sequence_example", "plasmid_feature_example",
    "plasmid_restriction_example") %in% utils::data(package = "ggchord")$results[, 3]))
})

test_that("built-in REBASE never depends on the working directory", {
  old <- setwd(tempdir())
  on.exit(setwd(old), add = TRUE)
  db <- ggchord:::ggchord_builtin_rebase()
  expect_equal(nrow(db), 4909L)
  expect_equal(length(unique(db$enzyme)), 4907L)
  sites <- find_restriction_sites("AAAAGAATTCTTT", enzymes = "EcoRI")
  expect_equal(sites$database_version, "609")
  expect_true(grepl("^rebase609:e:", sites$pattern_id))
})
