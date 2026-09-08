# Generate the internal restriction-pattern database from the four permitted
# REBASE/EMBOSS source files. examples/rebase/misc.zip is intentionally never
# inspected or referenced.
pkgload::load_all(".", quiet = TRUE)
ggchord_rebase_database <- ggchord_parse_rebase("examples/rebase")
stopifnot(
  nrow(ggchord_rebase_database) == 4909L,
  length(unique(ggchord_rebase_database$enzyme)) == 4907L
)
old <- new.env(parent = emptyenv())
if (file.exists("R/sysdata.rda")) load("R/sysdata.rda", envir = old)
objects <- list(ggchord_rebase_database = ggchord_rebase_database)
if (exists("ggchord_common_feature_database", envir = old, inherits = FALSE)) {
  objects$ggchord_common_feature_database <- get(
    "ggchord_common_feature_database", envir = old
  )
}
list2env(objects, envir = environment())
save(list = names(objects), file = "R/sysdata.rda", compress = "xz")
