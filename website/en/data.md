Start with one row per reference sequence. Add a feature when you know a location on one reference, a link when you know locations on two references, or a track when you have a measurement or want to count observations. The first chord diagram uses just two reference rows and one alignment row.

<a id="fields"></a>

## Reference sequences {#fields}

The two rows selected from `seq_data_example` give the full genome lengths, including the parts outside the alignment. The ID in every other table must match a reference exactly, including any accession version. A local identifier such as a contig name works equally well.

<div class="field-definitions">

| Column | When needed | Meaning and unit |
| --- | --- | --- |
| `accver` | Required | Unique, nonempty sequence identifier, stored as text. |
| `length` | Required | Positive, full sequence length in bp for nucleotides or residues for proteins. |
| `sequence_unit` | Required for native protein coordinates, optional for nucleotide coordinates | Use `"aa"` for residues or `"bp"` for base pairs. Omission means bp. All references in one plot use the same unit. |

</div>

| accver | length (bp) |
| --- | ---: |
| MT108731.1 | 64,323 |
| MT118296.1 | 32,090 |

Coordinates start at 1 and include both ends. Positions 1–100 cover 100 units. A coordinate of 1–1 identifies the first base or residue. To place a protein interval on a nucleotide backbone, use a [CDS coordinate model](projection.md).

<a id="features"></a>

## Features on one sequence {#features}

One feature row identifies a coding sequence (CDS), promoter, repeat or another annotation. The installed MT108731.1 GFF3 begins with a whole-genome region, followed by the CDS below. Select coding sequences by type.

<div class="field-definitions">

| Column | When needed | Meaning and unit |
| --- | --- | --- |
| `accver` | Required | Identifier of the reference carrying the feature. |
| `start` | Required | First interval bound in the reference unit. |
| `end` | Required | Last interval bound in the same unit, inclusive. |
| `strand` | For directional shapes | Use `"+"` or `"-"` when known. For a feature drawn without an arrow, use `"."`. |
| `directionality` | Optional in the input to `as_ggchord_feature()` | Use `"forward"`, `"reverse"`, `"bidirectional"` or `"none"`. If `strand` is also present, it must agree. |
| `type` | Optional for a general feature, needed for type-based selection unless another classification field is supplied | Feature class such as CDS. |
| `anno` | Optional | Source annotation or product description. |
| `feature_id` | Needed for joins by ID | Stable feature identifier from the source, useful for linking labels to selected records. |

</div>

| accver | start (bp) | end (bp) | strand | type | GFF3 ID | Product |
| --- | ---: | ---: | --- | --- | --- | --- |
| MT108731.1 | 13,064 | 14,554 | + | CDS | cds-QIQ63803.1 | hypothetical protein |

Use `as_ggchord_feature()` to normalize annotations before plotting. It can map existing column names to the expected fields and records unstranded features as `"."` when strand is missing. For directional shapes, retain the source strand or an explicit direction. The bundled `gene_data_example` keeps the original annotation in `source_strand`.

Select rows by their recorded type before adding the layer. Use `geom_gene()` for coding sequences with a known strand. Use `geom_feature()` for intervals with a generic shape, or map a type column to `shape` when several annotation classes share a layer.

<a id="alignments"></a>

## Alignments between two intervals {#alignments}

A high-scoring segment pair (HSP) is one local alignment. The compact field names reflect ggchord’s origins in plotting BLAST tables. The prefix q means query and s means subject. An accession with its version supplies the name accver. An alignment from another source can use the same columns.

BLAST tabular formats 6 and 7 have 12 fields by default. Format 7 adds comment lines. The original file for the bundled example uses a custom format 7 with 17 fields, including full query and subject lengths and subject strand. Six fields suffice to locate a ribbon: two identifiers and four endpoints. The installed tutorial TSV also includes the same first alignment in the standard 12-column layout.

<div class="field-definitions">

| Column | When needed | Meaning and unit |
| --- | --- | --- |
| `qaccver` | Required | Reference identifier on the query side. |
| `saccver` | Required | Reference identifier on the subject side. |
| `qstart` | Required | First reported query position, in its reference unit. |
| `qend` | Required | Last reported query position, inclusive. |
| `sstart` | Required | First reported subject position, in its reference unit. |
| `send` | Required | Last reported subject position, inclusive. |
| `pident` | Optional, needed to colour by measured identity | Percentage of identical aligned positions, from 0 to 100. |
| `length` | Optional | Alignment length in columns, including gaps. On the reference table this name instead denotes a complete sequence length. |

</div>

The first row of `ribbon_data_example` comes from the first BLAST alignment for the tutorial pair. Both intervals increase in the source. If one interval decreases and the other increases, their relative orientation is reversed. Keep the four endpoints in their reported order. A feature’s biological strand and a backbone’s display orientation are separate from this relationship.

| Query | Query start (bp) | Query end (bp) | Subject | Subject start (bp) | Subject end (bp) | Identity (%) |
| --- | ---: | ---: | --- | ---: | ---: | ---: |
| MT108731.1 | 26,298 | 51,139 | MT118296.1 | 7,121 | 31,959 | 98.612 |

<a id="points"></a>

## Links between two individual positions {#points}

A line can connect a site on one reference to a site on another. If you deliberately represent an interval by one point, record how you selected it. For example, taking the midpoint of each side of the HSP above and rounding down gives 38,718 bp on MT108731.1 and 19,540 bp on MT118296.1. The line conveys the locations. Use a ribbon to retain the full aligned spans.

<div class="field-definitions">

| Column | When needed | Meaning and unit |
| --- | --- | --- |
| `qaccver` | Required | Reference identifier for the first point. |
| `qpos` | Required | Position of the first point, in the reference unit. |
| `saccver` | Required | Reference identifier for the second point. |
| `spos` | Required | Position of the second point, in the reference unit. |
| `sequence_unit` | Optional | When supplied, it must match the backbone unit on both sides. |

</div>

<a id="measurements"></a>

## Measurements, windows and counts {#measurements}

Use a measurement when a value is already known for one site or interval. The installed pUC19 bedGraph reports 60% GC for its first window. The file gives zero-based, half-open bounds 0–100. After import, the same 100 bases run from 1 through 100 inclusive. GC percentage describes the sequence composition in this window.

<div class="field-definitions">

| Column | When needed | Meaning and unit |
| --- | --- | --- |
| `accver` | Required | Reference on which the observation or counting window lies. |
| `start` | Required for an interval or an explicit window | First included base or residue. |
| `end` | Required for an interval or an explicit window | Last included base or residue. |
| `pos` | For a point measurement or an event at a position | Position in the reference unit. You can also map a named column to `x`. |
| `value` | For an already measured signal | The measured quantity and its own unit, such as GC percent. Map the measurement column to height or fill. |
| `window_id` | Optional for explicit windows | A unique label for matching each window to its statistics. |

</div>

For an event count, give one row per event. Two observations at the same position contribute two to the count. For interval coverage, give the source intervals and choose a bin width or supply explicit windows. Covered bases count each position once, even where intervals overlap. Dividing by the actual window length gives the covered fraction. The [track guide](tracks.md) plots measured values, and [window statistics](tracks-statistics.md) follows intervals into coverage estimates.

<a id="prepared"></a>

## Supply these tables to a plot {#prepared}

Start with `ggchord(reference)`. Pass each observation table to the layer's `data` argument. This example bundles reference and CDS tables with `ggchord_data()` and inspects `features$cds` before plotting. Use the same component when supplying the data to a feature layer.

{{example:prepared-tables}}

Continue to the [two-reference chord diagram](chord.md#backbones). If you are starting from files, see [sequences and annotations](import-sequences.md), [alignments and relationships](import-relations.md) or [measured signals](import-signals.md).
