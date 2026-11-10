# 首轮源码清单与既有能力探针；从仓库根目录运行，不调用包级测试。
out <- "tools/long-term-audit"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
ns <- asNamespace("ggchord")
write_csv <- function(x, name) utils::write.csv(x, file.path(out, name), row.names = FALSE, na = "")
definitions <- list()
for (path in sort(list.files("R", "[.]R$", full.names = TRUE))) {
  expressions <- parse(path, keep.source = TRUE)
  refs <- attr(expressions, "srcref")
  for (i in seq_along(expressions)) {
    e <- expressions[[i]]
    if (!is.call(e) || !identical(e[[1]], as.name("<-")) || !is.symbol(e[[2]])) next
    name <- as.character(e[[2]])
    definitions[[name]] <- list(path = path, line = refs[[i]][1],
      alias = if (is.symbol(e[[3]])) as.character(e[[3]]) else "")
  }
}
exports <- sort(getNamespaceExports("ggchord"))
imports <- getNamespaceImports("ggchord")
inventory <- lapply(exports, function(name) {
  d <- definitions[[name]]
  object <- getExportedValue("ggchord", name)
  origin <- setdiff(names(Filter(function(x) name %in% x, imports)), "")
  data.frame(name = name, kind = if (!is.null(d)) {
    if (nzchar(d$alias)) "alias" else "native"
  } else if (length(origin)) "reexport" else "other",
  origin = paste(origin, collapse = ";"), path = if (is.null(d)) "" else d$path,
  line = if (is.null(d)) NA else d$line, alias = if (is.null(d)) "" else d$alias,
  signature = if (is.function(object)) paste(deparse(args(object)), collapse = " ") else class(object)[1])
})
write_csv(do.call(rbind, inventory), "api.csv")
parameters <- lapply(exports, function(name) {
  object <- getExportedValue("ggchord", name)
  if (!is.function(object) || is.null(formals(object))) return(NULL)
  f <- formals(object)
  data.frame(name = name, parameter = names(f),
    default = vapply(f, function(x) paste(deparse(x), collapse = " "), character(1)))
})
write_csv(do.call(rbind, parameters), "parameters.csv")
internal <- lapply(names(definitions), function(name) {
  object <- get0(name, ns, inherits = FALSE)
  if (!is.function(object)) return(NULL)
  d <- definitions[[name]]
  data.frame(name = name, exported = name %in% exports, path = d$path, line = d$line,
    signature = paste(deparse(args(object)), collapse = " "))
})
write_csv(do.call(rbind, internal), "functions.csv")
environment <- list(baseline_sha = system2("git", c("rev-parse", "main"), stdout = TRUE),
  package_version = as.character(utils::packageVersion("ggchord")), R = R.version.string,
  dependencies = lapply(c("ggplot2", "dplyr", "pkgload", "jsonlite"), function(p)
    list(name = p, version = as.character(utils::packageVersion(p)))))
jsonlite::write_json(environment, file.path(out, "environment.json"), auto_unbox = TRUE, pretty = TRUE)

# 随包真实记录只做三种已有坐标的代表 build，不跑 testthat。
data("seq_data_example", package = "ggchord")
data("ribbon_data_example", package = "ggchord")
seqs <- seq_data_example[, c("accver", "length")]
links <- ribbon_data_example[1:2, ]
checks <- list()
probe <- function(name, expr) {
  warnings <- character()
  start <- proc.time()[["elapsed"]]
  value <- tryCatch(withCallingHandlers(force(expr), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) list(error = conditionMessage(e)))
  checks[[name]] <<- list(result = value, warnings = warnings,
    elapsed_seconds = unname(proc.time()[["elapsed"]] - start))
}
for (mode in c("chord", "collinear", "circular")) {
  probe(paste0("baseline_", mode), {
    p <- if (mode == "circular") ggchord(seqs[1, ]) + geom_seq() + coord_circular() else
      ggchord(seqs, links) + geom_seq() + geom_link_ribbon() +
        if (mode == "chord") coord_chord() else coord_collinear()
    b <- ggplot2::ggplot_build(p)
    list(layers = length(b$data), nonempty = sum(vapply(b$data, nrow, integer(1)) > 0))
  })
}
probe("duplicate_seq_order_rejected", {
  p <- ggchord(seqs) + geom_seq(seq_order = rep(seqs$accver[1], 2))
  ggplot2::ggplot_build(p)
  "unexpected success"
})
probe("track_bar_is_identity", {
  l <- geom_track_bar(data = data.frame(accver = seqs$accver[1], start = 1L,
    end = 10L, depth = 17), aes(track_value = depth), position = position_track("depth"))
  p <- ggchord(seqs[1, ]) + geom_seq() + l + coord_collinear()
  e <- export_ggchord_layout(p, include = "track", original_data = TRUE)
  list(stat_class = class(l$stat)[1], registry = e$track_registry,
    original_data = e$original_data)
})
probe("mapping_and_export", {
  p <- ggchord(seqs[1:2, ]) + geom_seq() + coord_collinear()
  e <- export_ggchord_layout(p, include = "seq")
  list(metadata = e$metadata, sequence_layout = e$sequence_layout,
    positions = map_ggchord_position(p, seqs$accver[1], c(1, seqs$length[1])))
})
probe("guide_measurement", {
  p <- ggchord(seqs) + geom_seq() + geom_seq_label()
  d <- diagnose_ggchord(p, width = 8, height = 6, device = "png")
  list(measurement = d$measurement, guides = d$guides)
})
probe("missing_target_apis", {
  targets <- c("geom_track", "geom_track_col", "coord_dotplot", "coord_sankey", "coord_radial_synteny")
  list(exported = stats::setNames(targets %in% exports, targets),
    feature_data_formal = "feature_data" %in% names(formals(ggchord)),
    seq_axis_formal = "seq_axis" %in% names(formals(coord_collinear)))
})
probe("real_paf_endpoints", {
  paf <- read_paf("examples/format-inputs/emales-pair.paf")
  stopifnot(all(paf$qstart == paf$query_start0 + 1),
    all(paf$qend == paf$query_end0),
    all(paf$sstart[paf$relative_strand == "-"] == paf$target_end0[paf$relative_strand == "-"]),
    all(paf$send[paf$relative_strand == "-"] == paf$target_start0[paf$relative_strand == "-"] + 1))
  coordinates <- paf[, c("qaccver", "saccver", "qstart", "qend", "sstart", "send", ".source_line", "direction")]
  write_csv(coordinates, "real-paf-dotplot-endpoints.csv")
  list(rows = nrow(paf), reverse = sum(paf$relative_strand == "-"),
    source = "examples/format-inputs/emales-pair.paf", conversion_verified = TRUE)
})
jsonlite::write_json(checks, file.path(out, "baseline-checks.json"), auto_unbox = TRUE,
  pretty = TRUE, null = "null", na = "null", digits = NA)
cat("API inventory and baseline probes written\n")
