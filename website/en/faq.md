## Why are labels different after resizing? {#resize}
The available physical space changed. Use [explicit-size preview and saving](output.html#size).
Check hidden-label diagnostics; do not infer missing input rows from absent text.

## How do I find supported parameters? {#removed}
Use the [function reference](reference.html) to check argument meanings and defaults. Put mappings in `aes()`, fixed styles on the layer, and colour rules in `scale_*()`; see [plotting concepts](modes.html). Functions accept both `colour` and `color` spellings; examples use `colour`.

## Why is a feature rejected as a gene? {#strand}
The gene shortcut requires known `+`/`-` strand. Keep unknown strand in a general
feature layer; do not guess `+`. Protein backbones accept native aa coordinates; CDS projection is needed to place those intervals on nucleotide backbones.

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
and package checks belong to an explicitly authorized release. See the release acceptance record for the actual scope of verification.

## Where should I report a problem? {#report}
Use the [issue tracker](https://github.com/DangJem/ggchord/issues) with a small input,
package/ggplot2 versions, code, intended width/height, device and diagnostic output.
Avoid attaching private sequence data; a synthetic reproduction is usually enough.
