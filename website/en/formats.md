<a id="choose-format"></a>

## How do I read the files I have? {#choose-format}

| Files | Start here |
| --- | --- |
| FASTA sequences, GFF3 annotations and BLAST alignments | [Draw a comparison from source files](real-data.md#read) |
| BED intervals or PAF (Pairwise mApping Format) alignments | [Check zero-based coordinates](import-intervals.md#zero-based) |
| GenBank sequence with annotations | [Preserve record IDs and segments](import-sequences.md#genbank) |
| bedGraph or BigWig continuous signals | [Keep values and missing ranges explicit](import-signals.md#signals) |
| VCF/BCF variants and structural breakpoints | [Retain original variant evidence](import-relations.md#variants) |
| Existing protein signature annotations | [Retain AA locations and source scores](import-relations.md#protein-signatures) |
| Existing feature pairs and groups | [Retain original relationships and memberships](import-relations.md#feature-relations) |
| BAM/CRAM aligned reads | [Use a mature coverage bridge](import-signals.md#aligned-reads) |
| Existing delta results | [Keep native metrics and coordinate units distinct](import-relations.md#delta) |
| Existing MAF or chain genomic alignments | [Preserve blocks and source coordinates](import-relations.md#genome-alignments) |
| Existing R tables | [Map required fields](import-sequences.md#adapters) |
| GRanges or GInteractions | [Convert interval objects](import-intervals.md#ranges) |

<a id="files"></a>

## From files to a plot {#files}

[Read FASTA, GFF3 and BLAST](import-sequences.md#files) or see the [complete sequence comparison](real-data.md#read). If your data are already in R, go straight to the [required fields and coordinates](data.md#tables).

## Connect reader results to layers {#reader-results}

| Format and result | Coordinates and next step | Additional decoder |
| --- | --- | --- |
| FASTA length table, GFF3 annotation table, BLAST link table | Declare bp or aa for FASTA. GFF3 uses one-based inclusive intervals. Check the source strand before drawing arrows. BLAST retains endpoint order and input units. Supply the length table as the reference to `ggchord()`. Pass the annotations and alignments as `data` to their respective layers. | Built in |
| BED, BEDPE, PAF and bedGraph each return a table | Zero-based half-open intervals become one-based inclusive intervals for features, paired links, local alignments and quantitative tracks, respectively. | Built in |
| BigWig signal table | Decoded intervals are already one-based and inclusive. Map `value` in a quantitative track. | `rtracklayer`, with existing bedGraph as an option on Windows |
| GenBank list | Pass the full-length table and `features$annotation` to the preparation function. One-based features retain their source parts. | Python and Biopython |
| Aligned FASTA long table, or MAF and chain lists | The alignment table records columns and source positions. MAF and chain provide full lengths and paired intervals. MAF also includes characters. Reports retain decoded information and original text, which coexist in memory for large imports. | Aligned FASTA is built in. MAF and chain use Python and Biopython |
| Delta list | Pass full lengths and links to preparation using decoder coordinates and the declared units. | Installed `show-coords` executable |
| VCF or BCF list | Use one-based reference allele spans in `features$variant`. Extract explicitly decoded breakends and links separately. | `VariantAnnotation`, plus pysam for BCF. See [variant inputs](import-relations.md#variants) for the structural variant extension |
| Signature or feature relation lists | Signatures provide `features$signature` in protein aa coordinates. Relations retain original pairs or group members. Only complete mapped original pairs produce drawable links. | Signature JSON requires `jsonlite`. Local tabular files are read directly |

## What do you need to do after import? {#data-operations}

- [Check identifiers, complete lengths and units](data.md#tables)
- [Filter, focus or organize comparison data](workflows.md)
- [Project between nucleotide and protein coordinates](projection.md)
- [Inspect source records and layout results](output.md#layout-results)
- [Find original files and build steps for the showcase figures](flagship-data.md)
