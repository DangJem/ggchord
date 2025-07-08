# ggchord: Multi-Sequence BLAST Alignment Chord Diagram Visualization Tool

## Overview
`ggchord` is an R function based on `ggplot2` for visualizing BLAST alignment results of multiple sequences as intuitive chord diagrams. It supports extensive style customization, making it easy to display homologous regions and structural relationships between sequences. Version 0.1.0 of `ggchord` represents a breakthrough upgrade from simple multi-sequence chord diagrams to **more feature-rich** multi-sequence chord diagrams, capable of simultaneously showing alignment relationships between multiple sequences:
- Each sequence is presented as an arc or custom track, with length proportionally mapped.
- Colored ribbons represent alignment regions between sequences, supporting coloring by similarity or source.
- Equipped with customizable axes for precise annotation of sequence positions and lengths.
- Supports layout optimizations such as global rotation and sequence orientation adjustment to adapt to different analysis scenarios.

It is suitable for research in comparative genomics, pan-genome analysis, phage-host sequence relationship studies, etc., helping researchers quickly identify homologous patterns between sequences.

## Key Features
- **Multi-sequence Support**: Simultaneously display alignment relationships of 2 or more sequences, no longer limited to pairwise comparisons.
- **Sequence-level Customization**:
  - Customize sequence order, orientation (forward/reverse), gaps, and radii.
  - Automatically or manually specify sequence colors and labels to improve readability.
- **Refined Axes**:
  - Each sequence has independent axes with major/minor ticks, clearly labeling length positions.
  - Adjust tick lengths, label sizes, and offsets to balance aesthetics and information density.
- **Flexible Ribbon Styles**:
  - 3 coloring schemes (single color, by query sequence, gradient by similarity).
  - Adjustable gap between ribbons and sequences; supports customization of Bézier curve control points for smoothness.
- **Layout Optimization**: The entire graph can be rotated to meet different display needs.
- **Debug Mode**: Assists in troubleshooting data issues by displaying counts of valid/invalid alignments.

## Installation
### Dependencies
- R (≥ 3.6.0)
- ggplot2 (≥ 3.3.0)
- ggnewscale (≥ 0.5.0)
- RColorBrewer

```r
install.packages("ggplot2")
install.packages("ggnewscale")
install.packages("RColorBrewer")
```

### How to install ggchord？
Install the stable version of gggenes from CRAN:

`install.packages("ggchord")`

If you want the development version, install it from GitHub:

`devtools::install_github("DangJem/ggchord")`


## Parameter Details
| Parameter Category | Parameter Name | Type | Default Value | Description |
| --- | --- | --- | --- | --- |
| **Core Data** | `seq_data` | data.frame/tibble | - | Data frame containing sequence information, must include columns:<br> - `seq_id`: Unique sequence identifier<br> - `length`: Sequence length |
|  | `ribbon_data` | data.frame/tibble | - | Data frame containing BLAST alignment results, optional. If provided, must include columns:<br> - `qaccver`: Query sequence ID<br> - `saccver`: Subject sequence ID<br> - `length`: Alignment length<br> - `pident`: Percentage of sequence identity<br> - `qstart`: Query sequence start position<br> - `qend`: Query sequence end position<br> - `sstart`: Subject sequence start position<br> - `send`: Subject sequence end position |
|  | `gene_data` | data.frame/tibble | - | Data frame containing gene annotation information, optional. If provided, must include columns:<br> - `seq_id`: Unique sequence identifier<br> - `start`: Gene start position<br> - `end`: Gene end position<br> - `strand`: Strand direction (+ or -)<br> - `anno`: Gene annotation |
| **Basic Style** | `title` | Character | NULL | Main title of the graph |
| **Sequence Layout** | `seq_order` | Character vector | NULL | Specify the drawing order of sequences; if NULL, uses the order in `seq_data` |
|  | `seq_labels` | Character vector or named vector | NULL | Labels for sequences; if NULL, uses `seq_id` |
|  | `seq_orientation` | Numeric vector or single value | 1 | Orientation of each sequence: 1 (forward) or -1 (reverse); default is forward |
|  | `seq_gap` | Numeric or vector | 0.03 | Length consistent with the number of sequences, defining the arc proportion [0,0.5) from the head of one sequence to the tail of the next |
|  | `seq_radius` | Numeric or vector | 1.0 | Radius of sequence arcs, supports single value or vector with length equal to the number of sequences |
|  | `seq_curvature` | Numeric or vector | 1.0 | Curvature of sequence arcs: 1 for standard arc, 0 for straight line, >1 for more curved |
|  | `seq_colors` | Color vector or named vector | NULL | Define colors for each sequence arc; if NULL, automatically generated based on RColorBrewer Set1 |
| **Gene Style** | `gene_offset` | Numeric, vector, or list | 0.03 | Radial offset distance between gene arrows and sequence arcs. Supports:<br> - Single value: same offset for all strands of all sequences<br> - Vector: length consistent with the number of sequences, same offset for all strands of each sequence<br> - List: named list where each element corresponds to a sequence; elements can be a single value (all strands of the sequence) or a named vector containing "+" and "-" (strand-specific) |
|  | `gene_width` | Numeric or vector | 0.1 | Width of gene arrows |
|  | `gene_label_show` | Logical | FALSE | Whether to display gene labels |
|  | `gene_label_rotation` | Numeric, vector, or list | 0 | Rotation angle (degrees) of gene labels, supports the same parameter format as `gene_offset` |
|  | `gene_label_size` | Numeric | 2.5 | Font size of gene annotations |
|  | `gene_label_radial_offset` | Numeric, vector, or list | 0 | Radial offset of gene labels relative to arrows (positive values outward, negative values inward), supports the same parameter format as `gene_offset` |
|  | `gene_label_circum_offset` | Numeric, vector, or list | 0 | Circumferential offset proportion of gene labels along the sequence (relative to gene length), supports the same parameter format as `gene_offset` |
|  | `gene_label_circum_limit` | Logical, vector, or list | TRUE | Whether to limit circumferential offset to no more than half the gene length, supports the same parameter format as `gene_offset` |
|  | `gene_color_scheme` | Character | "strand" | Specify gene color scheme, optional "strand" (by strand direction) or "manual" (manual specification) |
|  | `gene_colors` | Color vector | - | Fill colors for gene arrows, behavior depends on `gene_color_scheme`:<br> - "strand" mode: supports named vectors (only "+"/"-"), unnamed vectors (first "+" then "-"), or single value (same color for both strands); defaults to red for "+" and blue for "-"<br> - "manual" mode: supports named vectors (corresponding to `anno`), unnamed vectors (truncate excess, pad不足); defaults to automatically generated colors |
|  | `gene_order` | Character vector | NULL | Specify the display order of genes in the legend; if NULL, uses the order of genes in the data |
| **Ribbon Style** | `ribbon_color_scheme` | Character | "pident" | Coloring scheme for ribbons, optional "single", "query", or "pident" |
|  | `ribbon_colors` | - | - | Color parameters for ribbons:<br> - single: single color (single value or first element of vector)<br> - query: map colors by query sequence (named/unnamed vector or single value)<br> - pident: gradient color scale vector for generating gradients by similarity percentage, defaults to blue-to-yellow gradient |
|  | `ribbon_alpha` | Numeric | 0.35 | Transparency of ribbons [0,1] |
|  | `ribbon_ctrl_point` | Vector or list | - | Bézier control points for adjusting ribbon shape:<br> - Vector: length 2 (single control point) or 4 (c1x,c1y,c2x,c2y, dual control points)<br> - List: each element is a sublist containing 1-2 control points, defaults to automatic calculation |
|  | `ribbon_gap` | Numeric or vector | 0.15 | Radial distance between sequence arcs and ribbons |
| **Axis Settings** | `axis_gap` | Numeric or vector | 0.04 | Radial distance between axes and sequence arcs, supports negative values |
|  | `axis_tick_major_number` | Integer or vector | 5 | Number of major ticks per sequence |
|  | `axis_tick_major_length` | Numeric or vector | 0.02 | Length proportion of major ticks |
|  | `axis_tick_minor_number` | Integer or vector | 4 | Number of minor ticks between two major ticks |
|  | `axis_tick_minor_length` | Numeric or vector | 0.01 | Length proportion of minor ticks |
|  | `axis_label_size` | Numeric or vector | 3 | Font size of axis tick labels |
|  | `axis_label_offset` | Numeric or vector | 0 | Offset proportion of axis labels relative to ticks |
|  | `axis_label_orientation` | Character, numeric, or vector | "horizontal" | Orientation of axis labels:<br> - "horizontal": horizontal direction<br> - Numeric: rotation angle (degrees)<br> - Vector: length consistent with the number of sequences or named vector |
| **Layout Settings** | `rotation` | Numeric | 45 | Rotation angle (degrees) of the entire graph |
|  | `panel_margin` | List | list(t=0,r=0,b=0,l=0) | Margins around the graph (t=top, r=right, b=bottom, l=left) |
| **Legend & Debug** | `show_legend` | Logical | TRUE | Whether to display the legend |
|  | `show_axis` | Logical | TRUE | Whether to display axes and ticks |
|  | `debug` | Logical | FALSE | Whether to output debug information | 

## Plot Interpretation
- **Sequence Arcs**: Each colored arc represents a sequence, with length proportionally mapped. Arrows indicate direction (forward/reverse).
- **Ribbons**: Colored regions connecting different sequences, representing alignment intervals:
  - When colored by similarity, the color gradient reflects sequence identity (e.g., from blue to red indicates increasing similarity).
  - When colored by query sequence, ribbons of the same color originate from the same query sequence.
- **Axes**: Ticks and numbers outside each sequence arc, labeling sequence positions (units match sequence length) for easy localization of alignment regions.

