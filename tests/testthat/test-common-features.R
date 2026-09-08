test_that("internal common-feature database is complete and self-contained", {
  db <- ggchord:::ggchord_builtin_common_features()
  expect_equal(nrow(db$features), 1273L)
  expect_equal(nrow(db$segments), 1730L)
  expect_equal(nrow(db$qualifiers), 4551L)
  expect_equal(nrow(db$qualifier_links), 384L)
  expect_identical(db$metadata$source_workbook_sha256,
    "3d9e7a78c705de2beeb9d8af3c74dede714faa468cd6eee61f78e0f784b65312")
  expect_named(db, c("metadata", "features", "segments", "qualifiers",
    "qualifier_links", "feature_type_summary", "search_indexes"))
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

test_that("three FASTA-derived plasmids receive key dynamic annotations", {
  data(plasmid_example_pBR322)
  data(plasmid_example_pUC19c)
  data(plasmid_example_pBluescript_II_SK_plus)
  br <- find_common_features(plasmid_example_pBR322)
  uc <- find_common_features(plasmid_example_pUC19c)
  blue <- find_common_features(plasmid_example_pBluescript_II_SK_plus)
  expect_true(all(c("AmpR", "rop", "bom", "ori") %in% br$anno))
  expect_true(any(br$anno %in% c("TetR", "TcR")))
  expect_true(all(c("AmpR", "ori", "MCS", "lacZα") %in% uc$anno))
  expect_true(all(c("AmpR", "f1 ori", "MCS", "T7 promoter", "T3 promoter",
    "KS primer", "SK primer") %in% blue$anno))
})

test_that("new plasmid data objects exactly reproduce FASTA", {
  objects <- c("plasmid_example_pUC19c", "plasmid_example_pBR322",
    "plasmid_example_pBluescript_II_SK_plus")
  files <- c("pUC19c.fna", "pBR322.fna", "pBluescript II SK(+).fna")
  paths <- testthat::test_path("..", "..", "examples", "plasmid", files)
  skip_if_not(all(file.exists(paths)))
  for (i in seq_along(objects)) {
    data(list = objects[i])
    object <- get(objects[i])
    lines <- readLines(paths[i], warn = FALSE)
    sequence <- paste0(lines[!grepl("^>", lines)], collapse = "")
    expect_equal(nrow(object), 1L)
    expect_identical(names(object), c("accver", "label", "length", "sequence"))
    expect_identical(object$sequence, sequence)
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
