🌐 Language Switch: 【[现代汉语（Han）](README-Hans.md) | [英文（English）](README.md)】

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


## Usage Instructions
### Preliminary Data Preparation
Two types of input data need to be prepared:

#### 【Required】Sequence Information Data (`seq_data`)
A TSV (Tab-Separated Values) file containing basic sequence information, must include the following columns:
- `seq_id`: Unique sequence identifier (e.g., gene name, accession number)
- `length`: Sequence length (positive number)

Example:
```r
seq_data <- read.delim("seq_track.tsv", sep = "\t", stringsAsFactors = FALSE)
```
The format of `seq_track.tsv` is as follows (example):
```txt
seq_id	length
MT108731.1	64323
MT118296.1	32090
OQ646790.1	57367
OR222515.1	83080
```
You can automatically generate this table from FASTA files using the following command:
```bash
seqkit fx2tab -nil *fna | sed '1i seq_id\tlength' > seq_track.tsv
```

#### 【Optional】Alignment Data (`ribbon_data`)
A TSV (Tab-Separated Values) file containing BLAST alignment results (convertible from `outfmt6` or `outfmt7` formats), must include the following columns:
- `qaccver`: Query sequence ID (must exist in `seq_data$seq_id`)
- `saccver`: Subject sequence ID (must exist in `seq_data$seq_id`)
- `length`: Alignment length
- `pident`: Sequence similarity (percentage)
- `qstart`/`qend`: Start/end positions of the alignment on the query sequence
- `sstart`/`send`: Start/end positions of the alignment on the subject sequence

You can use the following script to perform BLAST alignments on example sequences and obtain results in outfmt7 format:
```bash
# Script to run BLAST alignments using example FASTA files
seqs=("MT108731.1" "MT118296.1" "OQ646790.1" "OR222515.1")
seqsNum=${#seqs[@]}
ext="fna"
for ((i=0; i<seqsNum-1; i++)); do
  for ((j=i+1; j<seqsNum; j++)); do
    echo -e "Running BLASTN: ${seqs[$i]} vs ${seqs[$j]}"
    blastn \
      -outfmt '7 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen sstrand stitle' \
      -query "${seqs[$i]}.${ext}" \
      -subject "${seqs[$j]}.${ext}" \
      -out "${seqs[$i]}__${seqs[$j]}.o7"
  done
done
```

#### 【Optional】Gene Data (`gene_data`)
A TSV (Tab-Separated Values) file, which must include the following columns:
- `seq_id`: Unique sequence identifier, must correspond to `seq_id` in the sequence information data (`seq_data`), such as gene names, accession numbers, etc.
- `start`: Gene start position
- `end`: Gene end position
- `strand`: Strand direction (usually `+` for forward, `-` for reverse)
- `anno`: Gene annotation, such as functional description of the gene

The format of the example file `gene_track.tsv` is as follows:
```txt
seq_id	start	end	strand	anno
MT108731.1	100	200	+	DNA binding protein
MT118296.1	300	400	-	Transcription factor
```

You can convert GFF3 format files into a gene data table using the `gff2gene_track.R` script. The script content is as follows:
```r
library(tidyverse)

# Get paths of all gff3 files in the current directory
gff3FilesPath <- list.files(path = ".", pattern = "*.gff3")

# Read all gff3 files and merge into a data frame
gff3Table <- map_df(gff3FilesPath,~read_tsv(.x,show_col_types = F,comment = "#",col_names = F) %>% set_names(c("seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes")))

# Filter records of type CDS and extract annotation information
geneTrackTable <- gff3Table %>% filter(type=="CDS") %>% mutate(anno=str_extract(attributes,"(?<=product=)[^;]+(?=;)")) %>% select(seq_id,start,end,strand,anno)

# Save the processed data frame as a TSV file
write_tsv(geneTrackTable,"gene_track.tsv")
```
After running the above script, a `gene_track.tsv` file will be generated in the current directory, which can be used as gene data for subsequent analysis and visualization.


## Usage Examples
### Data Reading
```R
# Read sequence length data
seq_data <- read.delim("seq_track.tsv", sep = "\t", stringsAsFactors = FALSE)

# Read and process BLAST data
read_blast <- function(file) {
  df <- read.delim(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
  colnames(df) <- c("qaccver","saccver","pident","length","mismatches","gapopen",
                    "qstart","qend","sstart","send","evalue","bitscore",
                    "qcovs","qlen","slen","sstrand","stitle")
  df
}
blast_files <- list.files(path = ".", pattern = "*.o7", full.names = TRUE)
all_blast <- do.call(rbind, lapply(blast_files, read_blast))
ribbon_data <- subset(all_blast, length >= 100)

# Read gene annotation data; to make the image more aesthetically pleasing, shorter gene annotations are filtered out here
gene_data <- read.delim("gene_track.tsv", sep = "\t", stringsAsFactors = FALSE) |> dplyr::slice_max(order_by = end-start, n = 5, by = seq_id)
```

### Passing Only Essential `seq_data`
For `ggchord`, sequence data is the most important and indispensable. By default, sequences will be arranged counterclockwise in the order of the input `seq_data`. Of course, these can be modified.
```R
part1_1 <- ggchord(
  seq_data = seq_data,
)
```
![plot](/examples/plots/v1.0.0/part1_1.jpg)


For example, in the following example, you can control the order, orientation, and curvature of sequences using `seq_order`, `seq_orientation`, and `seq_curvature`, and set sequence colors using `seq_colors`.
```R
part1_2 <- ggchord(
  seq_data = seq_data,
  seq_order = c("MT118296.1", "OR222515.1", "MT108731.1", "OQ646790.1"),
  seq_orientation = c(1,-1,1,-1),
  seq_curvature = c(0,2,-2,6),
  seq_colors = c("steelblue", "orange", "pink", "yellow")
)
```
![plot](/examples/plots/v1.0.0/part1_2.jpg)

### Adding Sequence Alignment Data
For gene alignment chord diagrams, sequence alignment is undoubtedly our main focus, so `ribbon_data` is the most important data next to `seq_data`.
By default, the fill color of ribbons is determined by the percentage identity in the BLAST results.
```R
part2_1 <- ggchord(
  seq_data = seq_data,
  ribbon_data = ribbon_data
)
```
![plot](/examples/plots/v1.0.0/part2_1.jpg)

Of course, these can also be modified. For example, you can set the fill color to be based on the query sequence, making it easier for users to identify alignments between different sequences.
```R
part2_2 <- ggchord(
  seq_data = seq_data,
  ribbon_data = ribbon_data,
  ribbon_color_scheme = "query"
)
```
![plot](/examples/plots/v1.0.0/part2_2.jpg)

If you think color is not important, you can also set it to a single color.
```R
part2_3 <- ggchord(
  seq_data = seq_data,
  ribbon_data = ribbon_data,
  ribbon_color_scheme = "single",
  ribbon_colors = "orange"
)
```
![plot](/examples/plots/v1.0.0/part2_3.jpg)


In addition, ribbons will automatically adjust to perfectly match parameters such as sequence orientation, curvature, spacing, and radius (note: the same applies to axes and gene arrows).
> The current version still has some issues; image distortion may occur with certain parameter combinations, which will be fixed in future versions.
```R
part2_4 <- ggchord(
  seq_data = seq_data,
  ribbon_data = ribbon_data,
  seq_orientation = c(1,-1,1,-1),
  seq_curvature = c(0,2,-2,6),
  seq_gap = c(.1,.05,.09,.05),
  seq_radius = c(1,5,1,1)
)
```
![plot](/examples/plots/v1.0.0/part2_4.jpg)


## Adding Gene Annotation Information
With gene annotation information, we can more easily identify consistent regions in the gene regions of sequences, helping to explore evolutionary relationships between different species.
```R
part3_1 <- ggchord(
  seq_data = seq_data,
  gene_data = gene_data
)
```
![plot](/examples/plots/v1.0.0/part3_1.jpg)

By default, the fill color of arrows follows the strand mode, i.e., using the strand direction of the gene sequence. You can also use the manual mode, where colors are filled based on the `anno` categories in `gene_data`.
```R
part3_2 <- ggchord(
  seq_data = seq_data,
  gene_data = gene_data,
  gene_color_scheme = "manual"
)
```
![plot](/examples/plots/v1.0.0/part3_2.jpg)

Of course, gene annotation labels can also be displayed. However, as you can see, adjusting to achieve a perfect effect may take some effort.
```R
part3_3 <- ggchord(
  seq_data = seq_data,
  gene_data = gene_data,
  gene_label_show = T,
  gene_label_rotation = 45,
  gene_label_radial_offset = .1,
  panel_margin = list(l=.2)
)
```
![plot](/examples/plots/v1.0.0/part3_3.jpg)


## Comprehensive Example
We usually plot using `seq_data`, `ribbon_data`, and `gene_data` simultaneously, resulting in a more visually appealing image.
```R
part4_1 <- ggchord(
  seq_data = seq_data,
  ribbon_data = ribbon_data,
  gene_data = gene_data,
)
```
![plot](/examples/plots/v1.0.0/part4_1.jpg)

Of course, `ggchord` also provides rich parameters to control various details of the image.
```R
part4_2 <- ggchord(
  seq_data = seq_data,
  ribbon_data = ribbon_data,
  gene_data = gene_data,
  title = "Multi-sequence Chord Diagram with Gene Annotations",
  seq_gap = .03,
  seq_radius = c(3,2,2,1),
  seq_orientation = c(-1, -1, -1, 1),
  seq_curvature = c(0,1,-1,1.5),
  gene_offset = list(c("+"=.2,"-"=-.2), 
                     c("+"=.2,"-"=-.2), 
                     c("+"=.2,"-"=0), 
                     c("+"=.2,"-"=.1)),
  gene_label_rotation = list(c("+"=45,"-"=-45), 
                             c("+"=.2,"-"=-.2), 
                             c("+"=.2,"-"=0), 
                             c("+"=.2,"-"=.1)),
  gene_label_radial_offset = c(0,0,0,0),
  gene_label_circum_offset = c(1, 0, -2, 0),
  gene_label_circum_limit = c(T,T,T,T),
  gene_width = .08,
  gene_label_show = T,
  gene_color_scheme = "strand",
  ribbon_gap = .1,
  ribbon_color_scheme = "pident",
  ribbon_ctrl_point = c(0,0),
  axis_label_orientation = c(0,45,80,130),
  axis_gap = 0,
  axis_tick_major_number = 5,
  axis_tick_major_length = 0.03,
  axis_tick_minor_number = 5,
  axis_tick_minor_length = 0.01,
  axis_label_size = 2,
  axis_label_offset = 2,
  rotation = 45, 
  show_axis = T,
  panel_margin = list(t=.1),
  debug = TRUE,
)
```
![plot](/examples/plots/v1.0.0/part4_2.jpg)

However, this is still an early version of the software, and there are many imperfect aspects or even bugs. These issues are expected to be resolved in future versions.


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

## Version History
### v0.2.0 (Latest)
- **Advanced Arc and Line Mode Optimization**:
  - Through enhanced curve-fitting algorithms, "arc mode" and "line mode" in pairwise sequence visualization achieve smoother transitions between different alignment regions. For example, when visualizing two closely related sequences with multiple short alignments, arcs or lines can now connect these regions more elegantly, reducing visual clutter.
  - Ensures more accurate representation of alignment data. In previous versions, there might have been slight distortions in visualization, especially for long-distance alignments. Now, the lengths and positions of arcs and lines more accurately match actual alignment coordinates, providing a more realistic view of sequence relationships.
- **Precise Curvature and Gap Control**:
  - Users can now control the curvature of arcs in the chord diagram with finer precision. A new parameter has been introduced to allow step-by-step adjustment of arc curvature, enabling users to highlight different types of alignment patterns. For instance, more pronounced curvature can be used to emphasize highly conserved regions, while flatter curves can be used for less significant alignments.
  - The gap between sequences and ribbons can be adjusted with greater granularity. This is useful for visualizing complex alignment scenarios where different sequences have varying degrees of similarity. Users can now set different gap values for different sequence pairs, ensuring the visualization is both clear and informative.
- **Enhanced Color Customization**:
  - Version 0.2.0 offers a wider range of color palettes for sequences and ribbons. In addition to existing default color schemes, users can now choose from various predefined palettes optimized for different types of data visualization. For example, there are palettes designed to highlight high-contrast regions and those for creating more subtle and aesthetically pleasing visualizations.
  - When using the `pident` color scheme for ribbons, users can now customize gradient colors. This allows for more personalized visualization, better representing the distribution of sequence similarity. For example, users can choose a heatmap-like gradient to clearly show the range of similarity values from low to high.
### v0.1.0
- Supports separate management of sequence, alignment, and gene data via `seq_data`, `ribbon_data`, and `gene_data`.
- Added sequence orientation control, custom order, gap, and radius adjustment.
- Implemented customizable axes (major/minor ticks, label positions).
- Ribbons support 3 coloring schemes (single color, by query sequence, by similarity gradient).
- Added global rotation and debug mode.
### v0.0.2
- Added multi-sequence support (upgraded from pairwise).
- Added "arc mode" and "line mode" switching for pairwise sequences.
- Supported arc curvature adjustment and gap control.
### v0.0.1
- Initial release, supporting chord diagram visualization for pairwise BLAST alignments.

## Contributions & Feedback
Issue reports for bugs or feature requests are welcome, as are contributions via Pull Requests. For usage issues, refer to example code or enable debug mode (`debug = TRUE`) to troubleshoot data problems.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
