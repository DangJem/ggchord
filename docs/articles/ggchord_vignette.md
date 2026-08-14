# ggchord tutorial

``` r

library(ggchord)
library(ggplot2)

data(seq_data_example)
data(ribbon_data_example)
data(gene_data_example)
```

This tutorial walks through the complete `ggchord` workflow: preparing
input data, importing files in R, validating and cleaning data, and
building plots layer by layer.

## 1. Data preparation

Three types of input data are used. All three are plain data frames.

### \[Required\] Sequence data (`seq_data`)

| Column   | Type      | Description                        |
|----------|-----------|------------------------------------|
| `seq_id` | character | Unique sequence identifier         |
| `length` | integer   | Sequence length (must be positive) |

Example:

| seq_id     | length |
|------------|--------|
| MT108731.1 | 64323  |
| MT118296.1 | 32090  |
| OQ646790.1 | 57367  |
| OR222515.1 | 83080  |

A common way to build this table from FASTA files:

``` bash
seqkit fx2tab -nil examples/fasta/*.fna | sed '1i seq_id\tlength' > examples/seq_track.tsv
```

### \[Optional\] Alignment data (`ribbon_data`)

One row per alignment block between two sequences (the column names
follow common alignment-tool conventions):

| Column    | Type      | Description                            |
|-----------|-----------|----------------------------------------|
| `qaccver` | character | Query sequence ID                      |
| `saccver` | character | Subject sequence ID                    |
| `length`  | integer   | Alignment length (bp)                  |
| `pident`  | numeric   | Percent identity (0–100)               |
| `qstart`  | integer   | Start position on the query sequence   |
| `qend`    | integer   | End position on the query sequence     |
| `sstart`  | integer   | Start position on the subject sequence |
| `send`    | integer   | End position on the subject sequence   |

Example rows:

| qaccver    | saccver    | length | pident | qstart | qend  | sstart | send  |
|------------|------------|--------|--------|--------|-------|--------|-------|
| MT108731.1 | MT118296.1 | 24856  | 98.612 | 26298  | 51139 | 7121   | 31959 |
| MT108731.1 | MT118296.1 | 4412   | 97.031 | 21513  | 25922 | 2365   | 6772  |
| MT108731.1 | MT118296.1 | 464    | 94.181 | 20691  | 21146 | 1032   | 1495  |

For example, the standard BLAST `-outfmt 7` output can be parsed into
this table directly:

``` bash
seqs=("MT108731.1" "MT118296.1" "OQ646790.1" "OR222515.1")
ext="fna"
for ((i=0; i<${#seqs[@]}-1; i++)); do
  for ((j=i+1; j<${#seqs[@]}; j++)); do
    blastn \
      -outfmt '7 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen sstrand stitle' \
      -query "examples/fasta/${seqs[$i]}.${ext}" \
      -subject "examples/fasta/${seqs[$j]}.${ext}" \
      -out "examples/blastn/${seqs[$i]}__${seqs[$j]}.o7"
  done
done
```

### \[Optional\] Gene data (`gene_data`)

One row per gene (or feature) on a sequence:

| Column   | Type      | Description                           |
|----------|-----------|---------------------------------------|
| `seq_id` | character | Sequence ID the gene belongs to       |
| `start`  | integer   | Gene start position                   |
| `end`    | integer   | Gene end position                     |
| `strand` | character | Strand direction (`+` or `-`)         |
| `anno`   | character | Gene annotation / functional category |

Example rows:

| seq_id     | start | end   | strand | anno                      |
|------------|-------|-------|--------|---------------------------|
| MT108731.1 | 60709 | 63087 | \+     | hypothetical protein      |
| MT118296.1 | 14628 | 16301 | \+     | virion structural protein |
| OQ646790.1 | 43765 | 46140 | \+     | integrase                 |
| OQ646790.1 | 13194 | 15551 | \+     | tail tape measure protein |

For example, this table can be converted from GFF3 files:

``` r

library(tidyverse)
gff3FilesPath <- list.files(path = "examples/gff3", pattern = "\\.gff3$", full.names = TRUE)
gff3Table <- map_df(gff3FilesPath, ~read_tsv(.x, show_col_types = F, comment = "#",
  col_names = F) %>% set_names(c("seq_id", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes")))
geneTrackTable <- gff3Table %>%
  filter(type == "CDS") %>%
  mutate(anno = str_extract(attributes, "(?<=product=)[^;]+(?=;)")) %>%
  select(seq_id, start, end, strand, anno)
write_tsv(geneTrackTable, "examples/gene_track.tsv")
```

## 2. Importing FASTA / BLAST / GFF3 data in R

External command-line tools (BLAST, seqkit, …) are often used to
*prepare* the data, but the files they produce can also be read directly
in R with the built-in import helpers. This keeps the “prepare data
outside R” and “import + plot in R” steps clearly separated:

``` r

library(ggchord)

# FASTA -> seq_data (read and combine all example FASTA files)
seq_data <- invisible(read_fasta_lengths(files = "examples/fasta/*.fna"))

# BLAST -outfmt 6/7 tabular output -> ribbon_data (12 or 17 columns auto-detected)
ribbon_data <- invisible(read_blast(files = "examples/blastn/*.o7"))

# GFF3 -> gene_data (CDS features by default; anno from product/Name/...)
gene_data <- invisible(read_gff3(files = "examples/gff3/*.gff3"))

ggchord(seq_data, ribbon_data, gene_data) +
  geom_seq() + geom_ribbon() + geom_gene()
```

[`read_blast()`](https://dangjem.github.io/ggchord/reference/read_blast.md)
preserves useful extra columns (`evalue`, `bitscore`, `qcovs`, `qlen`,
`slen`, `sstrand`, `stitle`);
[`read_gff3()`](https://dangjem.github.io/ggchord/reference/read_gff3.md)
preserves `type`, `source`, `score`, `phase` and `attributes`;
[`read_fasta_lengths()`](https://dangjem.github.io/ggchord/reference/read_fasta_lengths.md)
optionally splits NCBI-style headers with `header_delim = "|"`.

## 3. Tutorial: from data to a diagram

This section turns the prepared data into a complete chord diagram. Each
step shows a short, runnable example and its corresponding figure.

### 3.1 Validate and clean input data

``` r

data(seq_data_example)
data(ribbon_data_example)
data(gene_data_example)

validate_ggchord_data(seq_data_example, ribbon_data_example, gene_data_example)

clean_ggchord_data(seq_data_example, ribbon_data_example, gene_data_example,
                   unknown_id = "drop", out_of_range = "clip",
                   reversed_interval = "sort", invalid_pident = "clip")
```

[`ggchord()`](https://dangjem.github.io/ggchord/reference/ggchord.md)
also validates automatically with `validate = "warn"`. Use `"error"` to
stop on severe problems or `"none"` for large inputs.

### 3.2 Filter and merge ribbons

``` r

kept <- filter_ggchord_ribbons(ribbon_data_example, min_pident = 90,
                               drop_self_links = TRUE,
                               sort_by = "pident")
dedup <- deduplicate_ggchord_ribbons(kept$data, by = "exact",
                                     keep = "best_pident")
merged <- merge_ggchord_ribbons(dedup$data, max_gap = 0)
```

The returned `$data` tables can be passed directly to
[`ggchord()`](https://dangjem.github.io/ggchord/reference/ggchord.md).

### 3.3 Start with sequence arcs

``` r

ggchord(seq_data_example) +
  geom_seq()
```

![Sequence arcs.](../reference/figures/seq_only_default.png)

Sequence arcs.

### 3.4 Add alignment ribbons

``` r

ggchord(seq_data_example, ribbon_data_example) +
  geom_seq() + geom_ribbon()
```

![Ribbons coloured by percent
identity.](../reference/figures/ribbon_pident.png)

Ribbons coloured by percent identity.

### 3.5 Add genes and labels

``` r

ggchord(seq_data_example, gene_data = gene_data_example) +
  geom_seq() + geom_gene() + geom_gene_label_repel()
```

![Gene arrows with repelled
labels.](../reference/figures/gene_repel.png)

Gene arrows with repelled labels.

### 3.6 Add axes and sequence labels

``` r

ggchord(seq_data_example) +
  geom_seq() + geom_axis() + geom_seq_label()
```

![Axes and sequence labels.](../reference/figures/axis_seq_label.png)

Axes and sequence labels.

### 3.7 Group sequences

``` r

seq_grouped <- transform(seq_data_example,
                         seq_group = c("host", "host", "phage", "phage"))

ggchord(seq_grouped, ribbon_data_example, gene_data_example) +
  geom_seq(seq_group = "seq_group",
           seq_group_colors = c(host = "#E41A1C", phage = "#377EB8")) +
  geom_ribbon() + geom_gene()
```

![Sequences grouped with an inter-group gap and
labels.](../reference/figures/tutorial_seq_group.png)

Sequences grouped with an inter-group gap and labels.

### 3.8 Map numeric ribbon columns and direction

``` r

rb_scored <- transform(ribbon_data_example,
                       bitscore = seq_len(nrow(ribbon_data_example)) * 10)

ggchord(seq_data_example, rb_scored) +
  geom_seq() +
  geom_ribbon(ribbon_color_by = "bitscore",
              ribbon_alpha_by = "bitscore",
              ribbon_direction = "linetype")
```

![Continuous ribbon fill, alpha and direction
mapping.](../reference/figures/tutorial_ribbon_mapping.png)

Continuous ribbon fill, alpha and direction mapping.

### 3.9 Highlight regions and ribbons

``` r

regions <- data.frame(seq_id = "MT108731.1",
                      start = 1000, end = 4000, color = "orange")

ggchord(seq_data_example, ribbon_data_example) +
  geom_seq() + geom_ribbon() +
  geom_seq_region(regions = regions) +
  geom_ribbon_highlight(ribbon_ids = 1)
```

![Sequence regions and ribbon
highlight.](../reference/figures/tutorial_highlights.png)

Sequence regions and ribbon highlight.

### 3.10 Draw generic features

``` r

features <- data.frame(seq_id = c("MT108731.1", "MT118296.1"),
                       start = c(1000, 500), end = c(4000, 2000),
                       strand = c("+", "-"), type = c("CDS", "tRNA"))

ggchord(seq_data_example, ribbon_data_example) +
  geom_seq() + geom_ribbon() + geom_feature(features)
```

![Generic features drawn with
geom_feature().](../reference/figures/tutorial_features.png)

Generic features drawn with geom_feature().

### 3.11 Apply themes and scales

``` r

ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
  geom_seq() + geom_ribbon() + geom_gene() + geom_axis() +
  scale_color_manual(values = c("MT108731.1" = "#E41A1C",
                                "MT118296.1" = "#377EB8",
                                "OQ646790.1" = "#4DAF4A",
                                "OR222515.1" = "#984EA3")) +
  theme(panel.background = element_rect(fill = "grey95"),
        legend.position = "bottom", legend.box = "horizontal")
```

![A themed plot with all legends
grouped.](../reference/figures/legend_bottom.png)

A themed plot with all legends grouped.

### 3.12 Fine-tuned publication-style plot

``` r

ggchord(seq_data_example, ribbon_data_example, gene_data_example,
        title = "ggchord") +
  geom_seq(seq_radius = c(3.3, 2.5, 1.8, 1.25),
           seq_orientation = c(-1, -1, 1, -1),
           seq_colors = c("MT108731.1" = "#E76F51",
                          "MT118296.1" = "#264653",
                          "OQ646790.1" = "#2A9D8F",
                          "OR222515.1" = "#D9A62E")) +
  geom_ribbon(ribbon_alpha = 0.45) +
  geom_gene() +
  geom_gene_label_repel(gene_label_size = 2, seed = 42) +
  geom_seq_label() +
  geom_axis() +
  theme(plot.background = element_rect(fill = "#FBF9F6", colour = NA),
        panel.background = element_rect(fill = "#FBF9F6", colour = NA))
```

![A complete fine-tuned chord
diagram.](../reference/figures/combined_fine.png)

A complete fine-tuned chord diagram.

## 6. Flexible parameter formats

Sequence-level parameters (`seq_radius`, `seq_gap`, `axis_label_size`,
…) accept **a single value, an unnamed vector, a vector/list named by
sequence ID, a list named by sequence order (`"1"`, `"2"`, …), or an
unnamed list**. Gene-level parameters (`gene_label_rotation`,
`gene_offset`, …) additionally accept per-strand (`+`/`-`) values. All
of the following are valid:

``` r

# 1. One value for everything
gene_label_rotation = 20

# 2. One value per strand, for every sequence
gene_label_rotation = c("+" = -15, "-" = -45)

# 3. By sequence ID
gene_label_rotation = list(
  "MT118296.1" = c("+" = -15, "-" = -45),
  "OR222515.1" = c("+" = 30, "-" = -30),
  "MT108731.1" = c("+" = 15, "-" = -15),
  "OQ646790.1" = c("+" = 0,  "-" = 0)
)

# 4. By sequence order ("1" = first sequence)
gene_label_rotation = list(
  "1" = c("+" = -15, "-" = -45),
  "2" = c("+" = 30, "-" = -30),
  "3" = c("+" = 15, "-" = -15),
  "4" = c("+" = 0,  "-" = 0)
)

# 5. Unnamed list: follows the sequence order (same as #4)
gene_label_rotation = list(
  c("+" = -15, "-" = -45),
  c("+" = 30, "-" = -30),
  c("+" = 15, "-" = -15),
  c("+" = 0,  "-" = 0)
)

# 6. A length-one list recycles to every sequence
gene_label_rotation = list(20)
```

## 7. Layer reference

| Layer | Function | Description |
|----|----|----|
| Sequence arcs | [`geom_seq()`](https://dangjem.github.io/ggchord/reference/geom_seq.md) | Draws an arc (or line) for each sequence, with direction arrows |
| Alignment ribbons | [`geom_ribbon()`](https://dangjem.github.io/ggchord/reference/geom_ribbon.md) | Draws colored ribbons from alignment results |
| Gene arrows | [`geom_gene()`](https://dangjem.github.io/ggchord/reference/geom_gene.md) | Draws gene/feature arrow polygons |
| Gene labels | [`geom_gene_label()`](https://dangjem.github.io/ggchord/reference/geom_gene_label.md) | Draws gene labels at fixed positions |
| Repelled gene labels | [`geom_gene_label_repel()`](https://dangjem.github.io/ggchord/reference/geom_gene_label_repel.md) | ggrepel-style labels with leader lines, wrapping and overlap hiding |
| Axes | [`geom_axis()`](https://dangjem.github.io/ggchord/reference/geom_axis.md) | Draws axis lines, major/minor ticks and tick labels |
| Sequence labels | [`geom_seq_label()`](https://dangjem.github.io/ggchord/reference/geom_seq_label.md) | Places sequence names on/outside the arcs |

## 8. Plot interpretation

- **Sequence arcs** — each colored arc is one sequence, with length
  proportional to the sequence length and arrows showing direction.
- **Ribbons** — colored regions connecting sequences represent aligned/
  homologous intervals; color encodes identity, query or subject by
  default.
- **Gene arrows** — arrow polygons on the sequences; color encodes
  strand or functional category, with optional labels.
- **Axes** — ticks and numbers outside each arc label sequence
  positions.

## Further reading

- [Function
  reference](https://dangjem.github.io/ggchord/reference/index.html)
- [Version history
  (NEWS.md)](https://dangjem.github.io/ggchord/news/index.html)
- [Package homepage](https://dangjem.github.io/ggchord/)
