## Why are labels different after resizing? {#resize}
The available physical space changed. Use [explicit-size preview and saving](output.html#size).
Check hidden-label diagnostics; do not infer missing input rows from absent text.

## Why did an old parameter stop working? {#removed}
The API cleanup removed execution of old aliases. Use the [migration table](migration.html).
Do not mix old and new spellings for one semantic setting. The colour/color
function spelling aliases remain, with `colour` preferred in documentation.

## Why is a feature rejected as a gene? {#strand}
The gene shortcut requires known `+`/`-` strand. Keep unknown strand in a general
feature layer; do not guess `+`. Amino-acid intervals require explicit projection.

## Why is a ribbon missing in linear mode? {#ribbons}
By default ribbons connect adjacent displayed sequences. Inspect pair selection,
clipping and validation. Use explicit `pairs` or `adjacent_only = FALSE` deliberately.
This default does not remove original input records.

## Can the package compute alignments or validate primer experiments? {#biology}
It visualizes supplied alignments and computational search results. It does not
run BLAST, infer clusters from names or validate wet-lab outcomes. Search confidence
is not a probability. See [recognition limitations](recognition.html).

## Is a measured diagnostic a release certification? {#quality}
No. It describes a requested render. Platform/device/font checks, full vignettes
and package checks belong to an explicitly authorized release. This development
version has only scoped documentation verification.

## Where should I report a problem? {#report}
Use the [issue tracker](https://github.com/DangJem/ggchord/issues) with a small input,
package/ggplot2 versions, code, intended width/height, device and diagnostic output.
Avoid attaching private sequence data; a synthetic reproduction is usually enough.
