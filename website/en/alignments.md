<div class="task-intro"><p>Use standard BLAST tabular output without the historical custom 17-column layout. Distinguish file fields, minimal plotting roles and nucleotide versus protein units.</p></div>

## Prepare this page {#prepare}
In a fresh R session, run the code below, then follow this page in order. Colours and layout are optional display choices; inputs select only the fields needed. For file-based examples, first extract the [input bundle](../assets/tutorial-inputs.zip) so the working directory contains `examples/`.

```r
library(ggchord)
```

## Use the default tabular columns {#standard}
The historical custom list of 17 fields is unnecessary. Use standard tabular output:

```sh
blastn -query query.fna -subject subject.fna -outfmt 6 -out result.tsv
blastp -query proteins.faa -subject reference.faa -outfmt 7 -out result.tsv
```

Format `6` is a table; `7` adds comments. Without a custom field list both use
12 default fields. Omitting `-outfmt` entirely produces pairwise text, which this
tabular reader does not parse. Fields are query ID, subject ID, pident, length,
mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore. Their identifier
roles become `qaccver`/`saccver`; values are preserved and need not be accession versions.
See the [NCBI format documentation](https://www.ncbi.nlm.nih.gov/sites/books/NBK279684/table/appendices.T.options_common_to_all_blast/).

## File layout and minimal plotting fields differ {#minimal}
Twelve file columns do not require twelve plotting columns. After the
preparation code above, try the source-backed records:

{{code:standardblast}}

`blast_source` retains statistics and provenance; `blast_links` has six endpoint roles.
Identity requires `pident` only when mapped/filtered; alignment-length filters or
weights require `length`. `qlen`/`slen` support the corresponding coverage calculations.
They are not general plotting requirements. A custom headerless layout requires
explicit `col_names`; arbitrary six-column meanings are not guessed.
Use `as_ggchord_link(x, mapping = ...)` for an already-read non-BLAST table;
see [format conversion](formats.html#adapters).

Files in `examples/blastn-standard/` retain historical HSP rows and their first
12 values. They are column conversions, not new BLAST runs. Old `.o7` inputs are
unchanged; the conversion script and hashes record that distinction.

## Read BLASTp with explicit coordinate units {#protein}
This real BLASTP 2.9.0+ output and query CDS annotation come from the public
EMALE tutorial. Retain the original table, then select sufficiently long,
low-evalue query HSPs with complete known CDS mappings.

{{code:protein}}
{{results:protein}}

BLASTp endpoints are amino-acid coordinates. This example projects the query HSP
span into the known CDS's genomic coordinates, retaining strand and phase.
A gapped HSP still has a projectable query span; that envelope is not a base-by-base
nucleotide correspondence. Weak hits, gap counts and statistics remain in the source.

Projecting both ends requires the subject's own CDS map. It cannot be inferred
from an accession, maximum endpoint or annotation name. The `ungapped = TRUE`
assertion for `project_ggchord_links()` needs evidence. The [projection tutorial](projection.html)
uses two real plasmid CDS with identical oriented DNA to demonstrate a justified
ungapped two-ended case.

Native protein backbones with aa axes are outside this version. For `blastx` and
`tblastn`, check each endpoint's units independently. Parsing does not establish
that a table can be drawn directly in genomic bp coordinates.

## Sources and limits {#source}
Original EMALE records come from the [public gggenomes data/tutorial](https://thackl.github.io/gggenomes/articles/emales.html).
Only data records are reused; all figures are drawn with ggchord, without copying
upstream figures or plotting code. The input bundle's `sources.json` records URLs,
archive members, line ranges, hashes and selection. Parsing neither reruns
alignments nor establishes biological function.

## Check your result {#check-result}
A basic link table needs six fields. BLASTp reading retains aa coordinates; use a supported CDS mapping before drawing on this version’s bp axes.

<div class="next-step">Continue: <a href="projection.html">CDS coordinate projection →</a></div>
