# Pair each documented argument with the exported functions that accept it.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[[1L]])
pkgload::load_all(root, quiet = TRUE)
exports <- getNamespaceExports("ggchord")
tagged <- function(nodes, tag) {
  Filter(function(node) identical(attr(node, "Rd_tag"), tag), nodes)
}
rows <- list()
for (path in list.files(file.path(root, "man"), "\\.Rd$", full.names = TRUE)) {
  rd <- tools::parse_Rd(path)
  section <- tagged(rd, "\\arguments")
  if (!length(section)) next
  items <- tagged(section[[1L]], "\\item")
  terms <- vapply(items, function(item) paste(unlist(item[[1L]]), collapse = ""), "")
  documented <- trimws(unlist(strsplit(terms, ",", fixed = TRUE)))
  aliases <- vapply(tagged(rd, "\\alias"), function(node) paste(unlist(node), collapse = ""), "")
  for (api in intersect(aliases, exports)) {
    object <- getExportedValue("ggchord", api)
    if (!is.function(object)) next
    parameters <- intersect(names(formals(object)), documented)
    if (!length(parameters)) next
    defaults <- vapply(parameters, function(name) deparse1(formals(object)[[name]]), "")
    rows[[api]] <- data.frame(api = api, parameter = parameters, default = unname(defaults),
                              topic = tools::file_path_sans_ext(basename(path)),
                              anchor = paste0("parameter-", parameters))
  }
}
result <- do.call(rbind, rows)
result <- result[order(result$parameter, result$api), ]
jsonlite::write_json(result, args[[2L]], dataframe = "rows", auto_unbox = TRUE)
