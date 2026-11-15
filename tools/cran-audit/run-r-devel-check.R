# Explicit two-platform CRAN candidate check, not part of default tests.
stopifnot(grepl("development", R.version$status, ignore.case = TRUE))
if (identical(Sys.getenv("GGCHORD_EXPECT_DEBIAN"), "true")) {
  stopifnot(any(grepl('^ID=debian$', readLines('/etc/os-release'))))
}
out <- Sys.getenv("GGCHORD_PROFILE_DIR", unset = file.path(getwd(), "cran-results"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(GGCHORD_PROFILE_DIR = normalizePath(out, winslash = "/"))
writeLines(capture.output(sessionInfo()), file.path(out, "session.txt"))
if (file.exists('/etc/os-release')) file.copy('/etc/os-release', file.path(out,'os-release'))
rscript <- file.path(R.home("bin"), "Rscript")
r <- file.path(R.home("bin"), "R")
run <- function(command, args, log) {
  status <- system2(command, args, stdout = file.path(out, log), stderr = file.path(out, log))
  if (status != 0L) stop(log, " exited with status ", status)
}
# Each process starts without a loaded decoder namespace.
run(rscript, c("--vanilla", "tools/cran-audit/profile-readers.R", "vcf", "small"), "vcf-before.log")
run(rscript, c("--vanilla", "tools/cran-audit/profile-readers.R", "vcf", "stream"), "vcf-stream.log")
if (.Platform$OS.type != "windows") {
  run(rscript, c("--vanilla", "tools/cran-audit/profile-readers.R", "bigwig"), "bigwig.log")
}
run(r, c("CMD", "build", "."), "build.log")
archive <- paste0("ggchord_", read.dcf("DESCRIPTION")[1,"Version"], ".tar.gz")
stopifnot(file.exists(archive))
run(r, c("CMD", "check", "--as-cran", archive), "check.log")
