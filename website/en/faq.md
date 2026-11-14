<a id="resize"></a>

## Why are labels different after resizing? {#resize}
The available physical space changed. Use [explicit-size preview and saving](output.md).
Check hidden-label diagnostics; do not infer missing input rows from absent text.

<a id="removed"></a>

## How do I find supported parameters? {#removed}
Use the [function reference](../../man/) to check argument meanings and defaults. Put mappings in `aes()`, fixed styles on the layer, and colour rules in `scale_*()`; see [plotting concepts](modes.md). Functions accept both `colour` and `color` spellings; examples use `colour`.

<a id="strand"></a>

## Why is a feature rejected as a gene? {#strand}
The gene shortcut requires known `+`/`-` strand. Keep unknown strand in a general
feature layer; do not guess `+`. Protein backbones accept native aa coordinates; CDS projection is needed to place those intervals on nucleotide backbones.

<a id="ribbons"></a>

## Why is a ribbon missing in linear mode? {#ribbons}
By default ribbons connect adjacent displayed sequences. Inspect pair selection,
clipping and validation. Use explicit `pairs` or `adjacent_only = FALSE` deliberately.
This default does not remove original input records.

<a id="biology"></a>

## Can the package compute alignments or validate primer experiments? {#biology}
It visualizes supplied alignments and computational search results. It does not
run BLAST, infer clusters from names or validate wet-lab outcomes. Search confidence
is not a probability. See [recognition limitations](recognition.md).

<a id="quality"></a>

## Is a measured diagnostic a release certification? {#quality}
No. It describes a requested render. Platform/device/font checks, full vignettes
and package checks belong to an explicitly authorized release. See the release acceptance record for the actual scope of verification.

<a id="report"></a>

## Where should I report a problem? {#report}
Use the [issue tracker](https://github.com/DangJem/ggchord/issues) with a small input,
package/ggplot2 versions, code, intended width/height, device and diagnostic output.
Avoid attaching private sequence data; a synthetic reproduction is usually enough.
