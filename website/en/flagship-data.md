Flagship plotting pages directly use the bundled datasets listed here. Loading requires only an installed ggchord; original files from the checkout are needed only to rebuild the datasets. Complete source tables and display subsets remain separate.

## Dataset entries {#datasets}

| Dataset | Flagships | Ready-to-use contents |
| --- | --- | --- |
| `flagship_comparison_example` | Chord, Collinear, Radial, Sankey | Per-view sequence tables, selected HSPs, CDS, label order and display orientation; Collinear also includes route controls for the fixed layout |
| `flagship_dotplot_example` | Four-sequence Dotplot | All 32 HSPs, complete sequence lengths, count labels for 16 cells and the fixed axis gap |
| `flagship_fusion_example` | Six-panel Fusion | Overview HSP/CDS, cross-panel relationships, source windows, local tables, 60 real coding bases with source positions, and magnification regions |
| `flagship_circular_example` | Dense circular flagship and 13 complete reference maps | Complete sequence lengths, features, independent reference primers, selected enzyme sets, candidate-site counts and selection policies |

{{code:courseblockg}}

```r
data(flagship_comparison_example)
chord_data <- flagship_comparison_example$chord
names(chord_data)
```

The four flagship entries in `flagship_comparison_example` are `chord`, `collinear`, `radial` and `sankey`. `flagship_circular_example$dense` corresponds to the existing dense flagship; `flagship_circular_example$references` is named by original sequence identifier. An additional `comparison` entry retains the original example subset for the strict four-coordinate comparison. Loading these entries requires no further discovery, filtering, focusing or character extraction.

## Original sources and units {#sources}

The four-genome inputs are the original FASTA, GFF3 and standard BLAST tables for MT108731.1, MT118296.1, OQ646790.1 and OR222515.1, retained under `examples/fasta/`, `examples/gff3/` and `examples/blastn-standard/`. Complete inputs remain traceable through `genome_comparison_example`: 4 complete sequences, 32 observed HSPs and 266 CDS. See the [real-input workflow](real-data.md) for original files and upstream records.

Circular maps use the 13 complete bundled plasmid sequences. Provenance and licenses for recognition databases, reference primers and enzyme sets are recorded in `inst/NOTICE` and the corresponding `data-raw/` generators. Original identifiers, sequence characters, coordinates, strands and record order are retained; biological data are not altered to improve appearance.

Nucleotide lengths and coordinates use bp and one-based inclusive intervals; Identity is a percentage. Reverse-HSP endpoints may descend and must not be sorted merely to remove a diagnostic. Display orientation and biological strand are stored separately.

## Display selections for comparison views {#comparison}

Chord and Radial use the same 12 HSPs and 20 CDS. HSPs have length at least 900 bp and Identity at least 75%, with up to three longest records per source pair. CDS are sampled by actual strand and source position; positive-strand sampling prioritizes overlap with displayed endpoints. Original short CDS and biological strands are retained.

Collinear and Sankey retain 16 CDS and 10 HSPs. Two crowded candidates are omitted only from this display, remaining in complete source tables. Adjacent rows retain selected HSPs; non-neighbour pairs retain only their longest selected record. Label input order is also saved.

Display directions and Collinear route controls are explicit presets for the existing order, gaps and alignment configuration. Recompute the relevant presets when changing those settings. They are not biological annotations, device measurements or label-placement solutions.

## Derived Dotplot and Fusion tables {#composition}

Dotplot retains all 32 HSPs without inventing reverse, self or missing records. Counts cover all 16 combinations of four sequences; zero means no record in this input. Label positions use native sequence lengths and a 4,000 bp axis gap; rebuild them when changing axis concatenation.

Fusion selects the observed HSP closest to 4,000 bp for each source pair and three CDS per sequence by source annotation order. Local windows follow the selected relationship endpoints, extend 800 bp on each side and clip at source boundaries. Complete sources, clipped local tables and context features are saved separately.

The character panel uses the first 60 real coding bases of the second displayed CDS on the lower local sequence, reverse-complemented for a negative strand. Each character retains its source position for character registration and magnification. Per-base alignment is not inferred from macro HSPs.

## Circular-map discovery and enzyme sets {#circular}

Feature discovery, reference-primer matching and restriction searches run only during dataset construction. The dense flagship retains unique sites with recognition length at least six from the commercial-all collection. Complete reference maps prefer saved enzyme profiles; an explicit None profile uses unique six-base-or-longer sites from the commercial nonredundant collection, with the choice documented in each record’s `policy`.

`sites` is the complete set actually drawn; `all_site_count` records the number of unfiltered search results. Large unselected candidate tables are not duplicated in the package and can be rebuilt from the original sequences and frozen recognition database. The older three-cut gallery selection is saved separately as `gallery_sites`; it must not be described as the complete enzyme set. Existing recognition and selection rules are retained.

For example, rebuild the dense circular map’s site selection as follows. Search results are used for construction and are not experimental cleavage measurements.

```r
all_sites <- find_restriction_sites(plasmid_example_pETDuet_1)
dense_sites <- filter_restriction_sites(
  all_sites, set = "unique_6plus", parent_set = "commercial_all"
)
dense_sites |>
  select(accver, enzyme, position) |>
  head()
```

## Rebuild and verify {#rebuild}

From the source root, after installing pkgload, dplyr (≥1.1.0), purrr and ggplot2, run:

```sh
Rscript data-raw/generate_website_data.R
Rscript data-raw/generate_flagship_data.R
```

Complete code is available in the [source and example-data generator](../../data-raw/generate_website_data.R) and the [flagship-data generator](../../data-raw/generate_flagship_data.R). Original input checksums are in the [source-file manifest](../../data-raw/website-data-sources.tsv); flagship generator dependencies are in the [flagship-input manifest](../../data-raw/flagship-data-sources.tsv).

Bundled objects store source/display tables and explicit fixed presets, not ggplot objects, grobs, device measurements or final label positions. Drawing still computes geometry and label placement at the physical output size. Verify source records, display order and rendered appearance after rebuilding, before updating the bundled data.

## Return to plotting {#plots}

[Four comparison views](multi-sequence-views.md#flagship) · [Dotplot](dotplot-start.md#matrix) · [Fusion](fusion-layout.md#flagship) · [Dense circular map](circular-dense-case.md) · [Complete reference maps](plasmid-catalog.md)
