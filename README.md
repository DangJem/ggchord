# ggchord: Chord Diagram for BLAST Alignment Visualization  

A ggplot2-based R function to visualize pairwise sequence alignment results from BLAST as intuitive chord diagrams, supporting both **arc** and **line** modes for flexible representation of homologous regions between query and subject sequences.  


## Overview  
`ggchord` transforms BLASTN outfmt6 format results into clear chord diagrams, with two visualization modes:  
- **Arc Mode** (default): Sequences mapped to semi-circular arcs, ideal for circular genomes (e.g., plasmids, phages).  
- **Line Mode**: Sequences displayed as horizontal lines, suitable for linear genomes or simplified comparisons.  

This tool is designed for researchers in bioinformatics, genomics, or molecular biology to quickly visualize sequence homology and structural relationships.  


## Features  
- **Dual Visual Modes**: Switch between arc and line modes via `arc_mode`.  
- **Adjustable Curvature**: Control arc smoothness with `curvature` (0–1).  
- **Precision Ribbon Alignment**: Ensures ribbons connect exact sequence coordinates in arc mode.  
- **Customizable Gaps**: Adjust horizontal gaps in arc mode (`gap_frac_arc`) or vertical gaps in line mode (`line_gap_frac`).  
- **Filtering**: Focus on meaningful alignments with `min_len` parameter.  
- **Aesthetic Customization**: Modify ribbon color, arc radii, and plot title.  
- **ggplot2 Compatibility**: Outputs a ggplot2 object for further styling or export.  


## Installation  
### Prerequisites  
- R (≥ 3.6.0)  
- `ggplot2` (≥ 3.3.0)  

### Install from GitHub  
```r
# Install dependencies
if (!require("ggplot2")) install.packages("ggplot2")

# Install ggchord (using devtools)
if (!require("devtools")) install.packages("devtools")
devtools::install_github("your_username/ggchord")  # Replace with your GitHub repo path
```


## Usage  

### Data Preparation  
Input data must be a BLASTN result in `outfmt6` format (or compatible), with the following columns:  
- `qstart`, `qend`: Alignment start/end positions on the query sequence  
- `sstart`, `send`: Alignment start/end positions on the subject sequence  
- `qlen`, `slen`: Total lengths of query and subject sequences (**must be constant across all rows**)  
- `length`: Length of the alignment  


### Basic Example  
```r
# Load the function
library(ggchord)

# Read BLAST results (adjust file path and column names as needed)
blast_df <- read.table(
  "blast_results.outfmt6",  # Your BLAST output file
  header = FALSE, 
  comment.char = "#", 
  stringsAsFactors = FALSE
)

# Rename columns to match required fields (adjust based on your outfmt6 format)
colnames(blast_df) <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore",
  "qlen", "slen"  # Ensure these columns exist
)

# Generate chord diagram in arc mode (default)
p_arc <- ggchord(
  blast_df = blast_df,
  min_len = 500,          # Filter alignments shorter than 500 bp
  title = "Arc Mode: Query vs Subject",
  arc_mode = TRUE,        # Arc mode
  curvature = 1.0,        # Full circular arc
  gap_frac_arc = 0.02     # Small gap between arcs
)

# Generate chord diagram in line mode
p_line <- ggchord(
  blast_df = blast_df,
  min_len = 500,          # Filter alignments shorter than 500 bp
  title = "Line Mode: Query vs Subject",
  arc_mode = FALSE,       # Line mode
  line_gap_frac = 0.5     # Vertical gap between lines
)

# Display the plots
print(p_arc)
print(p_line)

# Save to file
ggsave("blast_chord_arc.png", plot = p_arc, width = 8, height = 8, dpi = 300)
ggsave("blast_chord_line.png", plot = p_line, width = 8, height = 8, dpi = 300)
```


## Generating BLASTN Results for `ggchord`  
To use `ggchord`, you first need pairwise alignment results from BLASTN in `outfmt7` format (with comments) or `outfmt6` (without comments). Below is a bash script to automate pairwise BLASTN alignments and generate results compatible with `ggchord`.  

### Step 1: Prepare Your Sequence Files  
Ensure your sequences are in FASTA format (e.g., `.fna`), with filenames matching the sequence IDs (e.g., `vB_AbaM_CP14.fna` for sequence `vB_AbaM_CP14`).  

### Step 2: Run the BLASTN Script  
Use this bash script to perform pairwise alignments between all sequences in your list. It will generate output files formatted for `ggchord`:  

```bash
# Define your sequence IDs (match FASTA filenames without .fna)
seqs=("vB_AbaM_CP14" "PQ859668.1" "your_other_seq" "another_seq")  # Add/remove sequences as needed

# Get the number of sequences
seqsNum=${#seqs[@]}

# Run pairwise BLASTN alignments (i vs j, where i < j to avoid duplicates)
for ((i=0; i<seqsNum-1; i++)); do
  for ((j=i+1; j<seqsNum; j++)); do
    echo -e "Running BLASTN: ${seqs[$i]} vs ${seqs[$j]}"
    
    # Run BLASTN with output format compatible with ggchord
    blastn \
      -outfmt '7 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen sstrand stitle' \
      -query "${seqs[$i]}.fna" \  # Query sequence file (FASTA format)
      -subject "${seqs[$j]}.fna" \  # Subject sequence file (FASTA format)
      -out "${seqs[$i]}__${seqs[$j]}.o7"  # Output file (compatible with ggchord)
  done
done
```  

### What the Script Does:  
- **Sequence List**: The `seqs` array defines your sequence IDs (e.g., `vB_AbaM_CP14`, `PQ859668.1`). Replace these with your actual sequence names.  
- **Pairwise Alignments**: The nested loops run BLASTN for every unique pair of sequences (e.g., `seq1 vs seq2`, `seq1 vs seq3`, etc.) to avoid redundant comparisons.  
- **Output Format**: The `-outfmt '7 ...'` flag specifies a BLAST format with:  
  - Comments (via `7`) to make the file human-readable  
  - All fields required by `ggchord`: `qstart`, `qend`, `sstart`, `send`, `qlen`, `slen`, and `length` (alignment length).  
- **Output Files**: Results are saved as `[query]__[subject].o7` (e.g., `vB_AbaM_CP14__PQ859668.1.o7`), matching the example file included in this repo.  

### Requirements:  
- Install [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (ensure `blastn` is in your system PATH).  
- Your sequence files must exist in the working directory with names like `[seqID].fna` (e.g., `vB_AbaM_CP14.fna`).  


## Parameter Details  

| Parameter      | Type    | Default         | Description                                                                 |
|----------------|---------|-----------------|-----------------------------------------------------------------------------|
| `blast_df`     | data.frame | -               | Input data frame with BLAST results (must contain required fields)          |
| `min_len`      | numeric | 100             | Minimum alignment length to include (shorter alignments are filtered out)   |
| `title`        | string  | "BLAST Chord Diagram" | Plot title                                                                 |
| `ribbon_col`   | string  | "steelblue"     | Fill color for alignment ribbons                                             |
| `arc_mode`     | logical | TRUE            | TRUE for arc mode, FALSE for line mode                                      |
| `curvature`    | numeric | 1.0             | Curvature of arcs in arc mode [0,1]: 0=polyline, 1=full arc                |
| `gap_frac_arc` | numeric | 0.02            | Horizontal gap fraction between sequence arcs in arc mode                    |
| `line_gap_frac`| numeric | 0.5             | Vertical gap fraction between sequence lines in line mode                    |
| `r_query`      | numeric | 1.0             | Radius for query arc (arc mode) or ignored in line mode                     |
| `r_subject`    | numeric | 0.8             | Radius for subject arc (arc mode) or ignored in line mode                  |


## Plot Interpretation  

### Arc Mode  
- **Arcs**: The upper arc (0~π) represents the query sequence; the lower arc (π~2π) represents the subject sequence.  
- **Ribbons**: Translucent polygons connecting arcs represent alignment intervals, indicating homologous regions.  
- **Curvature**: Adjusts arc smoothness (`curvature = 1` for full arcs, `0` for segmented lines).  

### Line Mode  
- **Lines**: Horizontal lines represent query (top) and subject (bottom) sequences.  
- **Ribbons**: Rectangles connecting lines show alignment intervals.  
- **Line Gap**: Controlled by `line_gap_frac` for vertical spacing.  


## Example File  
This repository includes a sample BLAST result file `vB_AbaM_CP14__PQ859668.1.o7` for testing purposes. You can use it directly to try out `ggchord`:  

```r
# Load the sample data
blast_df <- read.table(
  "vB_AbaM_CP14__PQ859668.1.o7",  # Path to the sample file
  header = FALSE, 
  comment.char = "#", 
  stringsAsFactors = FALSE
)

# Rename columns and generate plots
colnames(blast_df) <- c("qaccver","saccver","pident","length","mismatch","gapopen",
                        "qstart","qend","sstart","send","evalue","bitscore",
                        "qcovs","qlen","slen","sstrand","stitle")

# Arc mode with full curvature
p_arc <- ggchord(blast_df, min_len=5000, title="Arc Mode", arc_mode=TRUE, curvature=1)
print(p_arc)

# Line mode with custom gap
p_line <- ggchord(blast_df, min_len=5000, title="Line Mode", arc_mode=FALSE, line_gap_frac=0.3)
print(p_line)
```  


## Notes  
- Ensure `qlen` and `slen` are constant across all rows of `blast_df` (the function will throw an error if not).  
- If no alignments meet `min_len`, the function will stop with a warning.  
- The output is a ggplot2 object, so you can customize it further (e.g., add themes, adjust labels):  
  ```r
  p_arc + theme(plot.title = element_text(hjust = 0.5))  # Center the title
  ```  


## Dependencies  
- [ggplot2](https://ggplot2.tidyverse.org/): For generating the visualization  


## License  
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.  


## Contributing  
Issues, bug reports, and pull requests are welcome! Feel free to open an issue if you encounter problems or have suggestions for improvement.  


## Release Notes  
### v0.0.2 (2025-07-01)  
- Added dual visualization modes (arc and line) for flexible sequence representation.  
- Implemented adjustable curvature for arcs in arc mode.  
- Enhanced ribbon alignment precision in arc mode.  
- Added gap control parameters for both visualization modes.  
- Improved computational efficiency.  

### v0.0.1 (Initial Release)  
- Core functionality for generating BLAST chord diagrams.  
- Support for pairwise sequence alignments.  
- Customizable plot aesthetics.  
- Example dataset (`vB_AbaM_CP14__PQ859668.1.o7`) for quick testing.  
- BLASTN script for generating compatible input files.  

