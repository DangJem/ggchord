# ggchord: Chord Diagram for BLAST Alignment Visualization  

A ggplot2-based R function to visualize pairwise sequence alignment results from BLAST as intuitive chord diagrams, highlighting homologous regions between query and subject sequences.  


## Overview  
`ggchord` transforms BLASTN outfmt6 format results into clear chord diagrams, where:  
- The upper semicircle represents the query sequence  
- The lower semicircle represents the subject sequence  
- Colored ribbons connect aligned intervals, making it easy to identify homologous regions  

This tool is designed for researchers in bioinformatics, genomics, or molecular biology to quickly visualize sequence homology and structural relationships.  


## Features  
- Parses BLASTN outfmt6 results directly, requiring minimal data preprocessing  
- Maps sequence lengths proportionally to circular arcs for accurate scaling  
- Filters alignments by minimum length to focus on meaningful regions  
- Customizable aesthetics (ribbon color, radii, gap size, title)  
- Outputs a ggplot2 object for further styling or export  


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

# Generate chord diagram
p <- ggchord(
  blast_df = blast_df,
  min_len = 500,          # Filter alignments shorter than 500 bp
  title = "Query vs Subject Alignment",
  ribbon_col = "darkorange",  # Custom ribbon color
  r_query = 1.2,          # Radius for query arc
  r_subject = 0.9         # Radius for subject arc
)

# Display the plot
print(p)

# Save to file
ggsave("blast_chord.png", plot = p, width = 8, height = 8, dpi = 300)
```


## Generating BLASTN Results for `ggchord`  
To use `ggchord`, you first need pairwise alignment results from BLASTN in `outfmt7` format (with comments) or `outfmt6` (without comments). Below is a bash script to automate pairwise BLASTN alignments and generate results compatible with `ggchord`.  

### Step 1: Prepare Your Sequence Files  
Ensure your sequences are in FASTA format (e.g., `.fna`), with filenames matching the sequence IDs (e.g., `vB_AbaM_CP14.fna` for sequence `vB_AbaM_CP14`).  

### Step 2: Run the BLASTN Script  
Use this bash script to perform pairwise alignments between all sequences in your list. It will generate output files formatted for `ggchord`:  

```bash
# Define your sequence IDs (match FASTA filenames without .fna)
seqs=("vB_AbaM_CP14" "PQ859668.1" "another_seq")  # Add/remove sequences as needed

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
| `gap_frac`     | numeric | 0.02            | Fraction of the circle allocated to gaps between query and subject arcs     |
| `r_query`      | numeric | 1.0             | Radius of the query sequence arc                                            |
| `r_subject`    | numeric | 0.8             | Radius of the subject sequence arc                                         |


## Plot Interpretation  
- **Arcs**: The upper arc (0~π) represents the query sequence; the lower arc (π~2π) represents the subject sequence. Arc length is proportional to sequence length.  
- **Ribbons**: Translucent polygons connecting arcs represent alignment intervals, indicating homologous regions.  
- **Gaps**: The blank space between arcs (controlled by `gap_frac`) prevents overlap.  


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

# Rename columns and generate plot (as shown in the basic example)
colnames(blast_df) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                        "qstart","qend","sstart","send","evalue","bitscore",
                        "%qcov_per_subject","qlen","slen","strand","title")
p <- ggchord(blast_df, min_len=500, title="Demo: vB_AbaM_CP14 vs PQ859668")
print(p)
```  


## Notes  
- Ensure `qlen` and `slen` are constant across all rows of `blast_df` (the function will throw an error if not).  
- If no alignments meet `min_len`, the function will stop with a warning.  
- The output is a ggplot2 object, so you can customize it further (e.g., add themes, adjust labels):  
  ```r
  p + theme(plot.title = element_text(hjust = 0.5))  # Center the title
  ```  


## Dependencies  
- [ggplot2](https://ggplot2.tidyverse.org/): For generating the visualization  


## License  
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.  


## Contributing  
Issues, bug reports, and pull requests are welcome! Feel free to open an issue if you encounter problems or have suggestions for improvement.  


## Release Notes  
### v0.0.1 (Initial Release)  
- Core functionality for generating BLAST chord diagrams  
- Support for pairwise sequence alignments  
- Customizable plot aesthetics  
- Example dataset (`vB_AbaM_CP14__PQ859668.1.o7`) for quick testing  
- BLASTN script for generating compatible input files  
