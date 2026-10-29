# Public signature introspection

Rebuild with `Rscript tools/validate-api-signatures.R` after the static inventory.
This loads the local development namespace and compares local exported formals only.
The APIs themselves are not called. No data, geometry, parameter meaning, validation,
return value, side effect or coordinate compatibility is certified by this check.
Upstream ggplot2/grid reexports remain outside this repository's signature check.
Source revision, input fingerprint and runtime versions are stored on each row.
