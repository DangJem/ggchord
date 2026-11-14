Compare four complete proteins with two bundled tables, without downloading inputs or running BLASTp. Both backbones and alignment endpoints use amino-acid coordinates (aa).

## Plot BLASTp protein alignments directly {#protein}

`protein_seq_data_example` contains three EMALE proteins and one Mavirus protein. Inspect identifiers, full lengths and units:

{{code:nativeproteininput}}
{{output:nativeproteininput}}

The lengths are 302, 302, 306 and 310 aa, obtained from complete proteins rather than alignment endpoints. The `sequence` column retains complete residues; only the three plotting columns are shown here.

`protein_ribbon_data_example` contains six observed BLASTp comparisons between all four proteins. Inspect both interval endpoints and identity:

{{code:nativeproteinlinks}}
{{output:nativeproteinlinks}}

{{code:nativeproteinplot}}
{{figure:nativeprotein}}

The collinear layout connects adjacent rows by default, displaying three ribbons. Colours show identity and the scale bar spans 50 aa. The first two proteins are highly similar; the remaining adjacent comparisons are weaker. Intervals are inclusive and retain both endpoint spans of gapped alignments.

Switching to a chord layout displays all six pairs by default:

{{code:nativeproteinchord}}
{{figure:nativeproteinchord}}

Pairwise identity ranges from about 35.7% to 99.0%. The plot uses the default Identity scale without additional colour settings.

Use generic layers such as `geom_feature()` for annotations at amino-acid positions, declaring `directionality = "none"` for nondirectional intervals. One plot uses either aa or bp. The endpoints of `blastx` and `tblastn` have different units and require conversion to one coordinate domain or separate displays. Gene arrows, primers and restriction sites require nucleotide input. A circular protein layout describes a visual arrangement.

## Place protein hits on a nucleotide genome {#projection}

To locate hits on a nucleotide backbone, supply known CDS mappings and project aa intervals to bp. This separate example uses a source HSP for BVI_023A_0009 with a complete CDS, then projects its query end. This subsection reads original files from the [manual input bundle](formats.html).

```r
library(dplyr)
```

{{code:protein}}
{{output:protein}}
{{code:proteinproject}}
{{output:proteinproject}}
{{results:protein}}

Residues 50–289 project to 7,181–7,900 bp on BVI_023A, on the reverse strand. A single-end projection can retain the span of a gapped HSP. Connecting two genomes also requires the subject CDS model. See [CDS projection](projection.html) for a supported ungapped two-ended example.

## Data source {#source}

The bundled identifiers are BVI_029_0020, BVI_058_0011, E4-10_023_0015 and ADZ16414.1. Complete FASTA records come from public EMALE/Mavirus data. NCBI BLASTp 2.17.0+ computed all comparisons, retaining one HSP per pair. The raw output retains self hits and both directions, for 16 rows. The teaching table retains one direction per unordered pair in sequence-table order, without changing endpoints or identity.

The input bundle's `examples/protein/` directory retains complete FASTA records, raw BLASTp output and `sources.json` hashes. Rebuild rules are in `data-raw/prepare-bundled-proteins.py` and `data-raw/generate_protein_data.R`. Plotting uses bundled objects and does not require BLASTp.

## Result checks {#check-result}

Check complete lengths, endpoints and units. Native protein plots use aa. Locating hits on a genome requires CDS projection to bp. Use `view_ggchord()` to inspect the target output size.
