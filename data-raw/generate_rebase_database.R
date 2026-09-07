# Generate the internal restriction-pattern database from the four permitted
# REBASE/EMBOSS source files. examples/rebase/misc.zip is intentionally never
# inspected or referenced.
ggchord_rebase_database <- ggchord_parse_rebase("examples/rebase")
stopifnot(
  nrow(ggchord_rebase_database) == 4909L,
  length(unique(ggchord_rebase_database$enzyme)) == 4907L
)
if (!identical(Sys.getenv("GGCHORD_REBASE_REDISTRIBUTION_CONFIRMED"), "true")) {
  stop(
    "REBASE files declare all rights reserved. Confirm redistribution terms ",
    "before setting GGCHORD_REBASE_REDISTRIBUTION_CONFIRMED=true and bundling ",
    "the generated database.", call. = FALSE
  )
}
save(
  ggchord_rebase_database,
  file = "R/sysdata.rda",
  compress = "xz"
)
