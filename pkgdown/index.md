# ggchord: Multi-Sequence Alignment Chord Diagram Visualization Tool

## Overview
`ggchord` is an R package built on `ggplot2` that visualizes multi-sequence alignment results as intuitive chord diagrams using layered grammar of graphics. Version 0.4.0 moved layout parameters from `ggchord()` into individual `geom_*` layers, aligning closely with `ggplot2`'s design philosophy — each layer manages its own style. The latest v0.5.0 release removes the `ggnewscale`/`RColorBrewer` dependencies, adds ribbon outline customization, and ships updated documentation and example plots.

- Each sequence is presented as an arc with length proportionally mapped.
- Colored ribbons represent alignment regions between sequences, supporting coloring by similarity or source.
- Gene annotations are overlaid as arrows, colorable by strand direction or functional category.
- Customizable axes provide precise annotation of sequence positions and lengths.
- Layout parameters can be distributed across the corresponding `geom_*` layers, or omitted to use sensible defaults.

Suitable for comparative genomics, pan-genome analysis, phage-host sequence relationship studies, and more.

## Key Features
- **Genuine `ggplot2` Layered Style**: `ggchord()` only takes data and global parameters; each `geom_*` receives its own layout parameters, just like `ggplot2`.
- **Deferred Computation Model**: Layout is computed at `print()` time, so parameters specified across `geom`s are all collected and applied during rendering.
- **Multi-sequence Support**: Display alignment relationships of two or more sequences simultaneously.
- **Sequence-level Customization**: Order, orientation, gaps, radius, and curvature — parameters belong in `geom_seq()`.
- **Flexible Ribbon Styles**: 3 coloring schemes, Bézier control points — parameters belong in `geom_ribbon()`.
- **Gene Arrow Layer**: Color by strand or annotation category, with labels — parameters belong in `geom_gene()`.
- **Refined Axes**: Major/minor ticks, label size and orientation — parameters belong in `geom_axis()`.

## Installation
### Dependencies
- R (≥ 3.6.0)
- ggplot2 (≥ 3.3.0)

```r
install.packages("ggplot2")
```

### Installing ggchord

From CRAN:

```r
install.packages("ggchord")
```

From GitHub:

```r
devtools::install_github("DangJem/ggchord")
```

## Quick Start

### Building a Chord Diagram

```r
library(ggchord)

data(seq_data_example)
data(ribbon_data_example)
data(gene_data_example)

# Simplest: all default parameters
p <- ggchord(
  seq_data = seq_data_example,
  ribbon_data = ribbon_data_example,
  gene_data = gene_data_example
) +
  geom_seq() +
  geom_ribbon() +
  geom_gene() +
  geom_axis()

print(p)
```

### Parameters in Geom Layers (v0.4.0)

```r
# Just like ggplot2: parameters go with the layer that uses them
ggchord(seq_data_example, ribbon_data_example, gene_data_example,
        title = "Fine Parameter Control", rotation = 30) +
  geom_seq(
    seq_radius = c(3, 2, 2, 1),
    seq_curvature = c(0, 1, -1, 1.5),
    seq_orientation = c(-1, -1, -1, 1)
  ) +
  geom_ribbon(
    ribbon_color_scheme = "pident",
    ribbon_gap = 0.1
  ) +
  geom_gene(
    gene_offset = list(
      c("+" = 0.2, "-" = -0.2),
      c("+" = 0.2, "-" = -0.2),
      c("+" = 0.2, "-" = 0),
      c("+" = 0.2, "-" = 0.1)
    ),
    gene_width = 0.08,
    gene_label_show = TRUE,
    gene_label_rotation = list(
      c("+" = 45, "-" = -45),
      c("+" = 0.2, "-" = -0.2),
      c("+" = 0.2, "-" = 0),
      c("+" = 0.2, "-" = 0.1)
    )
  ) +
  geom_axis(
    axis_label_orientation = c(0, 45, 80, 130),
    axis_gap = 0,
    axis_tick_major_length = 0.03,
    axis_label_size = 2
  )
```

> All parameters have sensible defaults. You can write `geom_seq()` with no arguments, or fine-tune each detail as shown above.

---

## Usage Instructions

### Data Preparation

Three types of input data are expected:

#### [Required] Sequence Data (`seq_data`)

A data frame describing the sequences to draw. It must contain the following columns:

| Column | Type | Description |
| --- | --- | --- |
| `seq_id` | character | Unique sequence identifier |
| `length` | integer | Sequence length (positive) |

Example:

```r
seq_data <- read.delim("seq_track.tsv", sep = "\t", stringsAsFactors = FALSE)
```

Where `seq_track.tsv` looks like:

| seq_id | length |
| --- | --- |
| MT108731.1 | 64323 |
| MT118296.1 | 32090 |
| OQ646790.1 | 57367 |
| OR222515.1 | 83080 |

Auto-generate from FASTA files:

```bash
seqkit fx2tab -nil *fna | sed '1i seq_id\tlength' > seq_track.tsv
```

#### [Optional] Alignment Data (`ribbon_data`)

A data frame with one row per pairwise alignment block between two sequences. It must contain the following columns (named after common alignment-tool output conventions, e.g., BLAST):

| Column | Type | Description |
| --- | --- | --- |
| `qaccver` | character | Query sequence ID |
| `saccver` | character | Subject sequence ID |
| `length` | integer | Alignment length (bp) |
| `pident` | numeric | Percent identity (0-100) |
| `qstart` | integer | Start position on the query sequence |
| `qend` | integer | End position on the query sequence |
| `sstart` | integer | Start position on the subject sequence |
| `send` | integer | End position on the subject sequence |

Example rows:

| qaccver | saccver | length | pident | qstart | qend | sstart | send |
| --- | --- | --- | --- | --- | --- | --- | --- |
| MT108731.1 | MT118296.1 | 24856 | 98.612 | 26298 | 51139 | 7121 | 31959 |
| MT108731.1 | MT118296.1 | 4412 | 97.031 | 21513 | 25922 | 2365 | 6772 |
| MT108731.1 | MT118296.1 | 464 | 94.181 | 20691 | 21146 | 1032 | 1495 |

For example, the output of a batch BLAST run can be parsed into this table directly:

```bash
seqs=("MT108731.1" "MT118296.1" "OQ646790.1" "OR222515.1")
ext="fna"
for ((i=0; i<${#seqs[@]}-1; i++)); do
  for ((j=i+1; j<${#seqs[@]}; j++)); do
    blastn \
      -outfmt '7 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen sstrand stitle' \
      -query "${seqs[$i]}.${ext}" \
      -subject "${seqs[$j]}.${ext}" \
      -out "${seqs[$i]}__${seqs[$j]}.o7"
  done
done
```

#### [Optional] Gene Data (`gene_data`)

A data frame annotating genes (or other features) on the sequences. It must contain the following columns:

| Column | Type | Description |
| --- | --- | --- |
| `seq_id` | character | Sequence ID the gene belongs to |
| `start` | integer | Gene start position |
| `end` | integer | Gene end position |
| `strand` | character | Strand direction (`+` or `-`) |
| `anno` | character | Gene annotation / functional category |

Example rows:

| seq_id | start | end | strand | anno |
| --- | --- | --- | --- | --- |
| MT108731.1 | 60709 | 63087 | + | hypothetical protein |
| MT118296.1 | 14628 | 16301 | + | virion structural protein |
| OQ646790.1 | 43765 | 46140 | + | integrase |
| OQ646790.1 | 13194 | 15551 | + | tail tape measure protein |

For example, this table can be converted from a GFF3 file:

```r
library(tidyverse)
gff3FilesPath <- list.files(path = ".", pattern = "*.gff3")
gff3Table <- map_df(gff3FilesPath, ~read_tsv(.x, show_col_types = F, comment = "#",
  col_names = F) %>% set_names(c("seq_id", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes")))
geneTrackTable <- gff3Table %>%
  filter(type == "CDS") %>%
  mutate(anno = str_extract(attributes, "(?<=product=)[^;]+(?=;)")) %>%
  select(seq_id, start, end, strand, anno)
write_tsv(geneTrackTable, "gene_track.tsv")
```

---

## Usage Examples

### Reading Data

```r
seq_data <- read.delim("seq_track.tsv", sep = "\t", stringsAsFactors = FALSE)

read_blast <- function(file) {
  df <- read.delim(file, sep = "\t", header = FALSE,
                   stringsAsFactors = FALSE, comment.char = "#")
  colnames(df) <- c("qaccver", "saccver", "pident", "length", "mismatches",
                    "gapopen", "qstart", "qend", "sstart", "send",
                    "evalue", "bitscore", "qcovs", "qlen", "slen",
                    "sstrand", "stitle")
  df
}
blast_files <- list.files(path = ".", pattern = "\\.o7$", full.names = TRUE)
all_blast <- do.call(rbind, lapply(blast_files, read_blast))
ribbon_data <- subset(all_blast, length >= 100)

gene_data <- read.delim("gene_track.tsv", sep = "\t",
  stringsAsFactors = FALSE) |>
  dplyr::slice_max(order_by = end - start, n = 5, by = seq_id)
```

### Sequences Only

```r
# Default: counterclockwise in seq_data order
ggchord(seq_data = seq_data) + geom_seq()
```

![Sequence chord diagram (default)](examples/plots/seq_only_default.png)

```r
# Custom sequence layout
ggchord(seq_data = seq_data) +
  geom_seq(
    seq_order = c("MT118296.1", "OR222515.1", "MT108731.1", "OQ646790.1"),
    seq_orientation = c(1, -1, 1, -1),
    seq_curvature = c(0, 2, -2, 6),
    seq_colors = c("steelblue", "orange", "pink", "yellow")
  )
```

![Sequence chord diagram (custom layout)](examples/plots/seq_only_custom.png)

### Adding Alignment Ribbons

```r
# Default: colored by percent identity
ggchord(seq_data = seq_data, ribbon_data = ribbon_data) +
  geom_seq() + geom_ribbon()
```

![Alignment ribbons colored by percent identity](examples/plots/ribbon_pident.png)

```r
# Color by query sequence
ggchord(seq_data = seq_data, ribbon_data = ribbon_data) +
  geom_seq() +
  geom_ribbon(ribbon_color_scheme = "query")
```

![Alignment ribbons colored by query sequence](examples/plots/ribbon_query.png)

```r
# Single color
ggchord(seq_data = seq_data, ribbon_data = ribbon_data) +
  geom_seq() +
  geom_ribbon(ribbon_color_scheme = "single", ribbon_colors = "orange")
```

![Alignment ribbons in a single color](examples/plots/ribbon_single.png)

### Adding Gene Annotations

```r
# Color by strand
ggchord(seq_data = seq_data, gene_data = gene_data) +
  geom_seq() + geom_gene()
```

![Gene arrows colored by strand](examples/plots/gene_strand.png)

```r
# Color by annotation + show labels
ggchord(seq_data = seq_data, gene_data = gene_data) +
  geom_seq() +
  geom_gene(
    gene_color_scheme = "manual",
    gene_label_show = TRUE,
    gene_label_rotation = 45,
    gene_label_radial_offset = 0.1
  )
```

![Gene arrows colored by annotation with labels](examples/plots/gene_manual_label.png)

### Full Example

```r
# All defaults
ggchord(seq_data, ribbon_data, gene_data) +
  geom_seq() + geom_ribbon() + geom_gene() + geom_axis()
```

![Full chord diagram with all default parameters](examples/plots/combined_default.png)

```r
# v0.4.0: fine-grained control with parameters distributed across geom layers
ggchord(
  seq_data = seq_data,
  ribbon_data = ribbon_data,
  gene_data = gene_data,
  title = "Multi-sequence Chord Diagram with Gene Annotations",
  rotation = 45
) +
  geom_seq(
    seq_radius = c(3, 2, 2, 1),
    seq_orientation = c(-1, -1, -1, 1),
    seq_curvature = c(0, 1, -1, 1.5),
    seq_gap = 0.03
  ) +
  geom_ribbon(
    ribbon_color_scheme = "pident",
    ribbon_gap = 0.1
  ) +
  geom_gene(
    gene_offset = list(
      c("+" = 0.2, "-" = -0.2),
      c("+" = 0.2, "-" = -0.2),
      c("+" = 0.2, "-" = 0),
      c("+" = 0.2, "-" = 0.1)
    ),
    gene_width = 0.08,
    gene_label_show = TRUE,
    gene_label_rotation = list(
      c("+" = 45, "-" = -45),
      c("+" = 0.2, "-" = -0.2),
      c("+" = 0.2, "-" = 0),
      c("+" = 0.2, "-" = 0.1)
    )
  ) +
  geom_axis(
    axis_gap = 0,
    axis_tick_major_length = 0.03,
    axis_label_size = 2,
    axis_label_orientation = c(0, 45, 80, 130)
  )
```

![Full chord diagram with fine-grained control](examples/plots/combined_fine.png)

> Sequence-level parameters (e.g., `seq_radius`, `seq_gap`) support single value, unnamed vector, and named vector formats. Gene-level parameters additionally support per-strand (`+`/`-`) list formats.

---

## Layer Reference

| Layer | Function | Description |
| --- | --- | --- |
| Sequence Arcs | `geom_seq()` | Draws arcs (or lines) for each sequence, with directional arrows |
| Alignment Ribbons | `geom_ribbon()` | Draws colored ribbons from alignment results |
| Gene Arrows | `geom_gene()` | Draws gene annotation arrow polygons and labels |
| Axes | `geom_axis()` | Draws axis lines, major/minor ticks, and tick labels |

---

## Parameter Reference

### ggchord() Parameters

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `seq_data` | data.frame | - | Sequence info; must include `seq_id`, `length` |
| `ribbon_data` | data.frame | NULL | Alignment results (e.g., BLAST output) |
| `gene_data` | data.frame | NULL | Gene annotation data |
| `title` | character | NULL | Plot title |
| `rotation` | numeric | 45 | Global rotation angle (degrees) |
| `panel_margin` | numeric/list | 0 | Panel margins |
| `show_legend` | logical | TRUE | Display legend |
| `debug` | logical | FALSE | Output debug info |

### geom_seq() Parameters

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `seq_order` | character vector | NULL | Drawing order of sequences |
| `seq_labels` | character vector | NULL | Sequence labels |
| `seq_orientation` | numeric (1/-1) | 1 | Sequence direction |
| `seq_gap` | numeric [0, 0.5) | 0.03 | Gap proportion between sequences |
| `seq_radius` | numeric (> 0) | 1.0 | Sequence arc radius |
| `seq_curvature` | numeric | 1.0 | Curvature: 0=straight, 1=standard, >1=more curved |
| `seq_colors` | color vector | Set1 | Sequence arc colors |
| `linewidth` | numeric | 1.2 | Arc line width |
| `show_legend` | logical | TRUE | Show legend |

### geom_ribbon() Parameters

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `ribbon_color_scheme` | character | "pident" | Scheme: `"pident"`, `"query"`, `"single"` |
| `ribbon_colors` | color vector | auto | Ribbon color parameters |
| `ribbon_alpha` | numeric (0-1) | 0.35 | Ribbon transparency |
| `ribbon_ctrl_point` | vector/list | c(0,0) | Bézier control points |
| `ribbon_gap` | numeric/vector | 0.15 | Radial gap between sequences and ribbons |
| `alpha` | numeric | 0.35 | Transparency (overrides ribbon_alpha) |
| `ribbon_outline_color` | character | "black" | Ribbon outline (border) color |
| `ribbon_outline_width` | numeric | 0.05 | Ribbon outline line width |
| `ribbon_outline_linetype` | numeric/character | 1 | Ribbon outline line type (1 = solid) |
| `show_legend` | logical | TRUE | Show legend |

### geom_gene() Parameters

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `gene_offset` | numeric/vector/list | 0.1 | Radial offset of gene arrows |
| `gene_width` | numeric/vector/list | 0.05 | Gene arrow width |
| `gene_color_scheme` | character | "strand" | Scheme: `"strand"` or `"manual"` |
| `gene_colors` | color vector | auto | Gene arrow fill colors |
| `gene_order` | character vector | NULL | Gene display order in legend |
| `gene_label_show` | logical | FALSE | Display gene labels |
| `gene_label_size` | numeric | 2.5 | Label font size |
| `gene_label_rotation` | numeric/vector/list | 0 | Label rotation angle |
| `gene_label_radial_offset` | numeric/vector/list | 0 | Radial offset of labels |
| `gene_label_circum_offset` | numeric/vector/list | 0 | Circumferential offset |
| `gene_label_circum_limit` | logical/vector/list | TRUE | Limit circumferential offset |
| `show_legend` | logical | TRUE | Show legend |
| `show_label` | logical | NULL | Override gene_label_show |
| `label_size` | numeric | NULL | Override gene_label_size |

### geom_axis() Parameters

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `show_axis` | logical | TRUE | Display axes |
| `axis_gap` | numeric/vector | 0.04 | Radial gap to sequences |
| `axis_tick_major_number` | integer/vector | 5 | Number of major ticks |
| `axis_tick_major_length` | numeric/vector | 0.02 | Major tick length proportion |
| `axis_tick_minor_number` | integer/vector | 4 | Number of minor ticks |
| `axis_tick_minor_length` | numeric/vector | 0.01 | Minor tick length proportion |
| `axis_label_size` | numeric/vector | 3 | Tick label font size |
| `axis_label_offset` | numeric/vector | 1.5 | Label offset proportion |
| `axis_label_orientation` | character/numeric/vector | "horizontal" | Label orientation |
| `show_legend` | logical | FALSE | Show legend |

---

## Plot Interpretation
- **Sequence Arcs**: Each colored arc represents a sequence, proportionally mapped by length, with arrows indicating direction.
- **Ribbons**: Colored regions connecting different sequences represent alignment intervals.
- **Gene Arrows**: Arrow polygons annotated on sequences, colored by strand or functional category, with optional labels.
- **Axes**: Ticks and numbers outside each sequence arc, labeling sequence positions.

---

## Version History
### v0.5.0 (Latest)
- **Ribbon outline customization**: `geom_ribbon()` gains `ribbon_outline_color` (default `"black"`), `ribbon_outline_width` (default `0.05`) and `ribbon_outline_linetype` (default `1`, solid).
- **Dependencies removed**: dropped `ggnewscale` (ribbon/gene fill scales are now separated internally) and `RColorBrewer` (Set1 palette built in).
- **Bug fixes**: ribbon fill scale overwritten by the gene fill scale; `ribbon_alpha` rendered at the wrong opacity; `geom_axis(show_axis = FALSE)` error; mixed `axis_label_orientation` vectors; `brewer.pal()` warnings for <3 items; `geom_gene()` added before `geom_ribbon()`; axis-only plots; S3 method registration.
- **Documentation**: comments/messages in English; new man pages; data-preparation tables with example rows; rendered example plots; BLAST de-emphasized; vignette rewritten for the layered API.

### v0.4.0
- **Parameter Redistribution**: Layout parameters moved from `ggchord()` to individual `geom_*` layers. `ggchord()` now only handles data validation and global style (`title`, `rotation`, `panel_margin`, `show_legend`, `debug`).
- **Deferred Computation**: Coordinate layout computed at `print()` time. Parameters from all `geom`s are collected and applied during rendering.
- **Custom `print.ggchord()` Method**: Merge parameters → compute layout → inject data into layers → render.
- **15 unit tests** added.

### v0.3.0
- Layered API refactoring: split from monolithic function to `ggchord() + geom_seq() + geom_ribbon() + geom_gene() + geom_axis()`.
- Custom `+.ggchord` method for automatic list flattening.
- Lightweight `coord_chord()` coordinate system.

### v0.2.0
- Enhanced arc and line mode optimization.
- Precise curvature and gap control.
- Enhanced color customization.

### v0.1.0
- Separate management of sequence, alignment, and gene data.
- Sequence orientation, custom order, gap and radius adjustment.
- Customizable axes. Ribbons support 3 coloring schemes.

### v0.0.2
- Multi-sequence support. Arc/line mode switching.

### v0.0.1
- Initial release. Pairwise BLAST alignment chord diagram visualization.

---

## Contributions & Feedback
Bug reports and feature requests are welcome via GitHub issues. Pull requests are also appreciated.

## License
This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
