🌐 Language Switch: [现代汉语](README-Hans.md)  

# ggchord: Multi-Sequence BLAST Alignment Chord Diagram Visualization Tool  

An R function built on ggplot2 for visualizing BLAST alignment results across multiple sequences as intuitive chord diagrams, with extensive customization options to highlight homologous regions and structural relationships between sequences.  


## Overview  
Version 0.1.0 of `ggchord` marks a breakthrough upgrade from supporting pairwise sequences to **multi-sequence chord diagrams**, enabling simultaneous visualization of alignment relationships across multiple sequences:  
- Each sequence is represented as an arc or customized trajectory, with length proportionally mapped to its actual size.  
- Colored ribbons indicate alignment regions between sequences, with support for coloring by similarity or source.  
- Equipped with customizable axes to precisely label sequence positions and lengths.  
- Supports global rotation, sequence orientation adjustments, and other layout optimizations to suit diverse analytical scenarios.  

Ideal for comparative genomics, pan-genome analysis, phage-host sequence relationship studies, and other research areas, helping researchers quickly identify homologous patterns between sequences.  


## Key Features  
- **Multi-sequence Support**: Visualize alignments across 2+ sequences, no longer limited to pairwise comparisons.  
- **Sequence-level Customization**:  
  - Customize sequence order, orientation (forward/reverse), gaps, and radii.  
  - Automatically generated or manually specified sequence colors and labels for improved readability.  
- **Fine-tuned Axes**:  
  - Each sequence has independent axes with major/minor ticks, clearly labeling length and positions.  
  - Adjust tick lengths, label sizes, and offsets to balance readability and aesthetics.  
  - Radial distance between axes and sequence arcs (`axis_gap`) is adjustable, including support for negative values (inward indentation).  
- **Flexible Ribbon Styling**:  
  - 3 coloring schemes: `single` (uniform color), `query` (by query sequence), and `pident` (gradient by similarity).  
  - Adjustable radial gap between ribbons and sequences; supports custom Bézier curve control points for smoothness.  
- **Layout Optimization**: Global rotation of the entire plot to adapt to different display needs.  
- **Debug Mode**: Assists in troubleshooting data issues by showing counts of valid/invalid alignments.  


## Installation  
### Dependencies  
- R (≥ 3.6.0)  
- ggplot2 (≥ 3.3.0)  
- RColorBrewer (≥ 1.1-3)  
- grDevices (built-in package)  

### Installation Steps  
```r
# Install dependencies
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("RColorBrewer")) install.packages("RColorBrewer")

# Install from GitHub (replace with your repository path)
if (!require("devtools")) install.packages("devtools")
devtools::install_github("your_username/ggchord@v0.1.0")
```


## Usage Instructions  

### Data Preparation  
Two types of input data are required:  

1. **Sequence Information Data (`seq_data`)**  
   A data frame containing basic sequence information, with mandatory columns:  
   - `seq_id`: Unique identifier for the sequence (e.g., gene name, accession number).  
   - `length`: Length of the sequence (positive value).  

   Example:  
   ```r
   seq_data <- data.frame(
     seq_id = c("seqA", "seqB", "seqC"),
     length = c(5000, 8000, 6500)  # Sequence lengths
   )
   ```  

2. **Alignment Data (`ribbon_data`)**  
   A data frame containing BLAST alignment results (convertible from `outfmt6` or `outfmt7` formats), with mandatory columns:  
   - `qaccver`: Query sequence ID (must exist in `seq_data$seq_id`).  
   - `saccver`: Subject sequence ID (must exist in `seq_data$seq_id`).  
   - `length`: Length of the alignment.  
   - `pident`: Percentage sequence identity.  
   - `qstart`/`qend`: Start/end positions on the query sequence.  
   - `sstart`/`send`: Start/end positions on the subject sequence.  

   Example: Generating from BLAST result files:  
   ```r
   # Read a single BLAST result
   read_blast <- function(file) {
     df <- read.delim(file, sep = "\t", header = FALSE, comment.char = "#")
     colnames(df) <- c("qaccver","saccver","pident","length","mismatches","gapopen",
                      "qstart","qend","sstart","send","evalue","bitscore",
                      "qcovs","qlen","slen","sstrand","stitle")
     df
   }
   
   # Batch read and merge multiple BLAST results
   blast_files <- list.files(pattern = "*.o7")  # Replace with your BLAST file paths
   all_blast <- do.call(rbind, lapply(blast_files, read_blast))
   ribbon_data <- subset(all_blast, length >= 100)  # Filter short alignments (optional)
   ```  


### Basic Usage Example  
```r
# Load packages
library(ggchord)
library(ggplot2)

# Generate multi-sequence chord diagram with ggchord
p <- ggchord(
  seq_data = seq_data,                # Sequence information
  ribbon_data = ribbon_data,          # Alignment data
  title = "Multi-sequence Alignment Chord Diagram",  # Plot title
  seq_order = c("seqA", "seqB", "seqC"),  # Custom sequence order
  seq_gap = 0.03,                     # Gap proportion between sequences
  seq_orientation = c(1, -1, 1),      # Sequence orientation (1=forward, -1=reverse)
  ribbon_color_scheme = "pident",     # Color by sequence identity
  ribbon_alpha = 0.7,                 # Ribbon transparency
  axis_gap = 0.1,                     # Distance between axes and sequences
  rotation = 45,                      # Global rotation (45 degrees)
  show_legend = TRUE                  # Show legend
)

# Display the plot
print(p)

# Save as high-resolution image
ggsave("multi_sequence_chord.png", plot = p, width = 10, height = 10, dpi = 300)
```  


## Parameter Details  

| Parameter Category   | Parameter Name               | Type          | Default Value                  | Description                                                                                     |
|----------------------|------------------------------|---------------|--------------------------------|-------------------------------------------------------------------------------------------------|
| **Core Data**        | `seq_data`                   | data.frame    | -                              | Data frame with sequence info, containing `seq_id` (sequence ID) and `length` (sequence length).|
|                      | `ribbon_data`                | data.frame    | -                              | Data frame with alignment info (must include `qaccver`, `saccver`, etc.).                       |
| **Basic Styling**    | `title`                      | character     | "Multi-sequence Chord Diagram" | Plot title.                                                                                     |
|                      | `show_legend`                | logical       | TRUE                           | Whether to display the legend.                                                                  |
|                      | `rotation`                   | numeric       | 45                             | Global rotation angle (degrees, positive = counterclockwise).                                  |
| **Sequence Layout**  | `seq_order`                  | character     | NULL                           | Order of sequences (defaults to order in `seq_data`).                                           |
|                      | `seq_labels`                 | character     | NULL                           | Labels for sequences (defaults to `seq_id`).                                                   |
|                      | `seq_orientation`            | numeric       | 1                             | Sequence orientation (1=forward, -1=reverse); supports single value or per-sequence settings.   |
|                      | `seq_gap`                    | numeric/vector| 0.05                          | Gap proportion between sequences ([0,0.5)), controlling arc spacing.                           |
|                      | `seq_radius`                 | numeric/vector| 1.0                           | Radius of sequence arcs; supports single value or per-sequence settings.                        |
|                      | `seq_colors`                 | character     | NULL                           | Colors for sequence arcs (auto-generated via RColorBrewer's Set1 by default).                   |
| **Ribbon Styling**   | `ribbon_color_scheme`        | character     | "single"                       | Ribbon coloring scheme: `single` (uniform), `query` (by query sequence), `pident` (gradient by identity). |
|                      | `ribbon_colors`              | character     | NULL                           | Parameters for coloring (single value/vector/gradient colors).                                  |
|                      | `ribbon_alpha`               | numeric       | 0.6                           | Transparency of ribbons ([0,1]).                                                                |
|                      | `ribbon_gap`                 | numeric/vector| 0.1                           | Radial gap between ribbons and sequence arcs.                                                  |
|                      | `ribbon_ctrl_point`          | vector/list   | NULL                           | Bézier curve control points (defaults to origin) for adjusting ribbon smoothness.               |
| **Axis Settings**    | `axis_gap`                   | numeric/vector| 0.05                          | Radial distance between axes and sequences (supports negative values for inward placement).     |
|                      | `axis_tick_major`            | integer/vector| 5                             | Number of major ticks.                                                                        |
|                      | `axis_tick_major_length`     | numeric/vector| 0.02                          | Length proportion of major ticks.                                                              |
|                      | `axis_tick_minor`            | integer/vector| 4                             | Number of minor ticks between major ticks.                                                     |
|                      | `axis_tick_minor_length`     | numeric/vector| 0.01                          | Length proportion of minor ticks.                                                              |
|                      | `axis_label_size`            | numeric/vector| 3                             | Size of axis label text.                                                                       |
|                      | `axis_label_offset`          | numeric/vector| 0                             | Offset of axis labels (positive = outward, negative = inward).                                  |
| **Debug**            | `debug`                      | logical       | FALSE                          | Whether to print debug info (counts of valid/invalid alignments).                               |


## Plot Interpretation  
- **Sequence Arcs**: Colored arcs representing each sequence, with length proportionally mapped. Arrows indicate orientation (forward/reverse).  
- **Ribbons**: Colored polygons connecting sequences, representing alignment regions:  
  - When colored by `pident`, the gradient reflects sequence similarity (e.g., from low to high identity).  
  - When colored by `query`, ribbons from the same query sequence share the same color.  
- **Axes**: Ticks and numbers outside each sequence arc, labeling positions (units match sequence length) for easy alignment localization.  


## Dependencies  
- [ggplot2](https://ggplot2.tidyverse.org/): Core plotting engine  
- [RColorBrewer](https://cran.r-project.org/package=RColorBrewer): Color scheme support  
- [grDevices](https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/00Index.html): Color handling tools  


## License  
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.  


## Version History  
### v0.1.0 (Latest)  
- Added multi-sequence support (upgraded from pairwise), with separate management of sequence and alignment data via `seq_data` and `ribbon_data`.  
- Added sequence orientation control, custom order, gap, and radius adjustments.  
- Implemented customizable axes (major/minor ticks, label positions).  
- Ribbons support 3 coloring schemes (uniform, by query, by identity gradient).  
- Added global rotation and debug mode.  

### v0.0.2  
- Added "arc mode" and "line mode" for pairwise sequences.  
- Supported arc curvature adjustment and gap control.  

### v0.0.1  
- Initial release, supporting chord diagram visualization for pairwise BLAST alignments.  


## Contributions & Feedback  
Issue reports for bugs or feature requests are welcome, as are contributions via Pull Requests. For usage issues, refer to example code or enable debug mode (`debug=TRUE`) to troubleshoot data problems.
