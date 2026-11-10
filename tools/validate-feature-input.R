# 定向公开行为验证；合成区间明确用于 API 教学，不是科研观测。
pkgload::load_all(".", quiet = TRUE)
seqs <- data.frame(accver = "A", length = 1000L)
features <- data.frame(accver = "A", start = c(10, 100, 200, 300, 400),
  end = c(50, 150, 250, 350, 450), strand = c("+", "+", "-", ".", "+/-"),
  type = c("CDS", "tRNA", "rRNA", "promoter", "repeat"), anno = letters[1:5])
features$evidence <- I(as.list(letters[1:5]))
original <- features
results <- list()
check <- function(name, expr) {
  value <- tryCatch({
    eval(substitute(expr), envir = new.env(parent = parent.frame()))
    "passed"
  }, error = function(e) {
    stop(name, ": ", conditionMessage(e), call. = FALSE)
  })
  results[[name]] <<- value
}
expect_error <- function(expr, pattern) {
  value <- tryCatch({ force(expr); NULL }, error = conditionMessage)
  stopifnot(is.character(value), grepl(pattern, value))
}
export <- function(p) export_ggchord_layout(p, include = c("gene", "feature", "labels"), original_data = TRUE)
base <- ggchord(seqs, feature_data = features) + geom_seq()
check("mixed_disjoint_layers", {
  out <- export(base + geom_gene() + geom_feature(
    data = function(x) x[x$type != "CDS", ], shape = "block") + coord_collinear())
  stopifnot(setequal(unique(out$gene$start), 10),
    setequal(unique(out$feature$start), features$start[-1]),
    setequal(unique(out$feature$type), features$type[-1]),
    identical(features, original), !nrow(out$annotation_duplicates))
})
check("full_feature_and_labels", {
  out <- export(base + geom_feature(shape = "block") + geom_feature_label() + coord_collinear())
  stopifnot(setequal(unique(out$feature$start), features$start),
    setequal(out$layer_sources$source_id, paste0("feature:", 1:5)),
    all(out$layer_sources$source_scope == "constructor"), !nrow(out$annotation_duplicates))
  input <- out$original_data[[out$layer_sources$layer_id[1]]]$gene_poly
  stopifnot(identical(input$evidence, features$evidence),
    identical(input$.feature_biological_strand, features$strand))
})
check("default_gene_and_duplicate_diagnostic", {
  out <- export(base + geom_gene() + geom_feature(shape = "block") + coord_collinear())
  stopifnot(setequal(unique(out$gene$type), "CDS"),
    identical(out$annotation_duplicates$source_id, "feature:1"))
})
check("gene_labels_default_selection", {
  out <- export(base + geom_gene_label() + coord_collinear())
  stopifnot(all(out$layer_sources$source_id == "feature:1"))
})
check("legacy_position_and_message", {
  legacy <- features[1, c("accver", "start", "end", "strand", "anno")]
  messages <- character()
  p <- withCallingHandlers(ggchord(seqs, NULL, legacy), message = function(m) {
    messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage")
  })
  stopifnot(length(messages) == 1L, grepl("feature_data", messages))
  build_messages <- character()
  out <- withCallingHandlers(export(p + geom_seq() + geom_gene()), message = function(m) {
    build_messages <<- c(build_messages, conditionMessage(m)); invokeRestart("muffleMessage")
  })
  stopifnot(!length(build_messages), setequal(unique(out$gene$start), legacy$start))
})
check("input_errors_and_validation", {
  expect_error(ggchord(seqs, gene_data = features, feature_data = features), "only one")
  expect_error(ggchord(seqs, feature_data = transform(features, strand = "unknown")), "strand")
  v <- validate_ggchord_data(seqs, feature_data = transform(features, accver = "missing"))
  stopifnot(!v$valid, all(v$errors$table == "feature"), identical(v$data_summary$n_feature, 5L))
})
check("unclassified_features_and_missing_strand", {
  untyped <- features[, c("accver", "start", "end", "anno")]
  p <- ggchord(seqs, feature_data = untyped) + geom_seq()
  out <- export(p + geom_feature() + coord_circular())
  stopifnot(all(out$original_data[["layer-0002"]]$gene_poly$.feature_biological_strand == "."))
  expect_error(export(p + geom_gene()), "explicit gene layer data")
})
check("aa_generic_features", {
  p <- ggchord(transform(seqs, sequence_unit = "aa"), feature_data = features) + geom_seq()
  out <- export(p + geom_feature(shape = "block") + coord_circular())
  stopifnot(identical(out$metadata$sequence_unit, "aa"))
  expect_error(export(p + geom_gene()), "require bp")
})
check("source_duplicates_are_not_deduplicated", {
  duplicates <- features[c(1, 1), ]
  p <- suppressWarnings(ggchord(seqs, feature_data = duplicates)) + geom_seq() + geom_feature()
  out <- export(p)
  stopifnot(setequal(out$layer_sources$source_id, c("feature:1", "feature:2")))
})
check("empty_gene_subset", {
  out <- export(ggchord(seqs, feature_data = features[-1, ]) + geom_seq() + geom_gene())
  stopifnot(!nrow(out$gene))
})
check("source_metadata_cannot_overwrite_geometry_identity", {
  data <- features[1, ]
  data$source_row <- 999L
  data$source_id <- "user-supplied-id"
  out <- export(ggchord(seqs, feature_data = data) + geom_seq() + geom_feature())
  original <- out$original_data[["layer-0002"]]$gene_poly
  stopifnot(setequal(out$feature$source_row, 1L),
    identical(original$source_row, 999L), identical(original$source_id, "user-supplied-id"),
    identical(out$layer_sources$source_id, "feature:1"))
})
check("source_backed_palette_inheritance", {
  data("plasmid_example_pUC19", package = "ggchord")
  annotations <- find_common_features(plasmid_example_pUC19)
  seqs <- plasmid_example_pUC19[, c("accver", "length")]
  inherited <- ggchord(seqs, feature_data = annotations) + geom_seq() +
    geom_feature_plasmid() + coord_circular()
  explicit <- ggchord(seqs) + geom_seq() + geom_feature_plasmid(data = annotations) + coord_circular()
  stopifnot(identical(ggplot2::ggplot_build(inherited)$data[[2]]$feature_fill,
    ggplot2::ggplot_build(explicit)$data[[2]]$feature_fill))
})
if (requireNamespace("GenomicRanges", quietly = TRUE)) check("granges_unknown_strand", {
  ranges <- GenomicRanges::GRanges("A", IRanges::IRanges(10, 30), strand = "*")
  out <- export(ggchord(seqs, feature_data = ranges) + geom_seq() + geom_feature())
  stopifnot(nrow(out$feature) > 0)
})
cat(paste(names(results), unlist(results), sep = ": "), sep = "\n")
jsonlite::write_json(list(version = as.character(utils::packageVersion("ggchord")),
  scope = "targeted feature-input checks; no package release checks", results = results),
  "tools/long-term-audit/feature-input-results-v0.17.0.json", pretty = TRUE, auto_unbox = TRUE)
