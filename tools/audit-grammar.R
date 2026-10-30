# Current plotting signatures and compatibility evidence, separate from the
# immutable v0.16.0 full-source audit. No plots or package tests are executed.
pkgload::load_all(quiet = TRUE)
out <- "design-roadmap/audit-generated"
files <- sort(c(list.files("R", pattern = "[.]R$", full.names = TRUE),
  list.files("man", pattern = "[.]Rd$", full.names = TRUE), "NAMESPACE",
  "DESCRIPTION", "tools/audit-grammar.R"))
before <- tools::md5sum(files)
exports <- sort(getNamespaceExports("ggchord"))
own <- read.csv(file.path(out, "api_exports.csv"), na.strings = character())
exports <- intersect(exports, own$symbol[!own$is_third_party_reexport])
functions <- exports[grepl("^(geom_|stat_|position_|coord_|scale_|guide_|theme_)", exports) |
  exports == "ggchord"]
owner <- function(fn, parameter) {
  if (parameter %in% c("mapping", "data", "show.legend", "inherit.aes")) return("layer input / mapping / guide participation")
  if (parameter == "position") return("local position adjustment")
  if (parameter == "...") return("documented forwarding / retired-argument guard")
  if (fn == "ggchord") return("constructor data / validation")
  if (grepl("^theme_", fn)) return("non-data appearance")
  if (grepl("^scale_", fn)) return("mapping rule / domain / breaks / guide")
  if (grepl("^guide_", fn)) return("guide construction")
  if (grepl("^coord_", fn)) return("global coordinate contract")
  if (grepl("^position_", fn)) return("local geometry placement")
  if (grepl("^stat_", fn)) return("stat computation or forwarded geom parameter; see help")
  if (grepl("^(label_|max_overlaps$)", parameter)) return("label layout")
  if (parameter %in% c("colour", "fill", "alpha", "linewidth", "linetype", "size", "fontface", "family")) return("fixed visual style; explicit value wins")
  "geom geometry / input selection; see entry help"
}
rows <- lapply(functions, function(fn) {
  args <- formals(getExportedValue("ggchord", fn))
  do.call(rbind, lapply(seq_along(args), function(i) {
    nm <- names(args)[i]
    default <- if (nm == "...") "<ellipsis>" else
      paste(deparse(args[[i]], width.cutoff = 400), collapse = " ")
    if (!nzchar(default)) default <- "<required>"
    data.frame(function_name = fn, parameter = nm, default_expression = default,
      responsibility = owner(fn, nm), stringsAsFactors = FALSE)
  }))
})
matrix <- read.csv("design-roadmap/audit-reviewed-details/coordinate-compatibility.csv",
  na.strings = character())
matrix$runtime_validation <- "NOT_RUN_IN_v0.16.1; SOURCE_CONTRACT_ONLY"
core <- matrix$layer %in% c("geom_seq", "geom_gene", "geom_feature",
  "geom_gene_label", "geom_feature_label", "geom_feature_label_repel")
matrix$runtime_validation[core] <- "LIMITED_COMBINATION: tools/validate-grammar.R"
density <- matrix$layer == "stat_ribbon_density" & matrix$coordinate != "coord_circular"
matrix$runtime_validation[density] <- "LIMITED_STAGED_COMBINATION: tools/validate-grammar.R"
routes <- matrix$layer == "geom_restriction_site" & matrix$coordinate == "coord_circular"
matrix$runtime_validation[routes] <- "LIMITED_ROUTE_COMPARISON: tools/validate-grammar.R"
links <- matrix$layer %in% c("geom_link_line", "geom_link_ribbon") &
  matrix$coordinate == "coord_collinear"
matrix$runtime_validation[links] <- "LIMITED_STYLE_COMBINATION: tools/validate-fixed-styles.R"
matrix$review_status <- "CURRENT_SUPPORT_CONTRACT_REVIEWED; NOT_PARAMETER_EXHAUSTIVE"
matrix$evidence <- paste("Historical locator (not current line):", matrix$evidence,
  "Current review: API-GRAMMAR-v0.16.1.md; same mode restrictions")
matrix$source_revision <- "v0.16.1 working snapshot; see grammar-source-hashes-v0.16.1.csv"
role_sets <- list(
  input_identifiers = c("accver", "qaccver", "saccver"),
  input_bp = c("start", "end", "qstart", "qend", "sstart", "send", "qpos", "spos"),
  input_biology = c("strand", "anno", "feature_type", "label"),
  geometry_mapping = c("seq_ring", "feature_shape"),
  scaled_visual_style = c("seq_colour", "gene_fill", "feature_fill", "region_fill",
    "ribbon_fill", "ribbon_colour", "ribbon_alpha", "ribbon_linetype",
    "link_colour", "link_alpha", "link_linewidth", "link_linetype"),
  computed_stat = c("bundle_n", "bundle_weight", "density"))
roles <- do.call(rbind, lapply(names(role_sets), function(kind) data.frame(
  aesthetic_or_role = role_sets[[kind]], responsibility = kind,
  evaluation = if (kind %in% c("input_identifiers", "input_bp", "input_biology", "geometry_mapping")) {
    "Input/geometry stage; late stat/scale stages forbidden for geometry roles"
  } else if (kind == "computed_stat") "Exposed by ribbon stat; accessible via after_stat" else
    "Visual mapping / scale / after_scale; explicit fixed style takes precedence",
  stringsAsFactors = FALSE)))
if (!identical(before, tools::md5sum(files))) stop("Source changed during inventory")
write.csv(do.call(rbind, rows), file.path(out, "grammar-parameters-v0.16.1.csv"), row.names = FALSE)
write.csv(matrix, file.path(out, "grammar-compatibility-v0.16.1.csv"), row.names = FALSE)
write.csv(roles, file.path(out, "grammar-aesthetics-v0.16.1.csv"), row.names = FALSE)
write.csv(data.frame(path = names(before), md5 = unname(before)),
  file.path(out, "grammar-source-hashes-v0.16.1.csv"), row.names = FALSE)
cat(length(functions), "plotting exports inventoried; compatibility runtime evidence remains limited\n")
