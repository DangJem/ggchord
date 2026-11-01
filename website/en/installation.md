<div class="lesson-intro">
Build a sequence comparison from bundled data and save a figure. You do not need to run BLAST or prepare your own files first.
</div>

<div class="lesson-goal">
By the end of this page, you will have a collinear plot containing sequence backbones, gene arrows and observed alignment intervals, saved as a PNG at an explicit size. The real-data tutorials use the same layer composition.
</div>

<div class="lesson-path">Install and load → inspect three tables → compose layers → read the result → save</div>

## 1. Install and load {#requirements}
Use R ≥4.1.0 and ggplot2 ≥4.0.0. Open R or RStudio and run the installation commands in its console.
Skip installation if ggchord is already installed and loads successfully.

```r
install.packages("remotes")
remotes::install_github("DangJem/ggchord@v0.16.5")
library(ggchord)
packageVersion("ggchord")
```

### Which version are you using? {#source}
This tutorial documents **ggchord 0.16.5**. The installation above pins its formal tag; see [GitHub Releases](https://github.com/DangJem/ggchord/releases/tag/v0.16.5) for release assets.
For a source checkout, run from its repository root:

```r
install.packages("pkgload")
pkgload::load_all(".")
packageVersion("ggchord")
```

Use the loading route appropriate to your source; you do not need both installation paths.
Record the version when reproducing a plot. This lesson needs neither optional Bioconductor packages nor Biopython.

## 2. Inspect the three input tables {#input}
Select the bundled sequences `MT108731.1` and `MT118296.1`. Keep the complete source objects;
the code below selects only the columns used here. Run the complete block first:

{{code:firstdata}}

| Object | Required columns | One row represents |
| --- | --- | --- |
| `first_sequences` | `accver`, `length` | A sequence and its full length |
| `first_genes` | `accver`, `start`, `end`, `strand` | A gene annotation |
| `first_links` | Two IDs and four query/subject endpoints | An observed alignment interval between sequences |

`accver` joins the tables and must agree; it need not be a GenBank accession. Coordinates are one-based inclusive bp.
The bundled gene fixture separately retains `source_strand`; restoring it uses biological direction rather than alternating strands for display.

Check column names, then inspect the first annotation rows:

```r
names(first_sequences)
names(first_genes)
names(first_links)
head(first_genes)
```

<div class="lesson-check">
Checkpoint: the sequence table has 2 columns, genes 4 and links 6. Neither pident nor qlen is required for this basic plot;
identity colouring or coverage calculations require their corresponding fields later.
</div>

## 3. Compose your first plot {#first-plot}
Start with a linear view so genomic distances are easy to read.
`ggchord()` supplies the data, three geometry layers draw links, backbones and genes, and the coordinate selects a linear arrangement.
Named arguments avoid accidentally passing genes into the second positional argument, `ribbon_data`.

{{code:firstplot}}
{{figure:first}}

The two colours identify the sequence containing each gene; **they do not assert homology-cluster membership**.
Neutral alignment bands keep attention on genes and positions.

## 4. Read the result {#read}
Inspect the plot in this order:

1. Names on the left identify the sequences. Backbones and ticks show full lengths and bp positions.
2. Arrows show selected bundled gene annotations and their source strands. An unshown gene is not evidence of absence.
3. Bands between rows show supplied alignment intervals. They are neither newly computed alignments nor base-by-base alignment displays.

This view teaches positions and links. Use a chord view for overall relationships among several sequences,
or a single-sequence circular view for plasmid annotations. Choose the question before the geometry;
see [which view to use](modes.html).

## 5. Save at an explicit size {#save}
`p_first` is an ordinary ggplot object. Save it to the current working directory:

{{code:firstsave}}

```r
file.exists("first-comparison.png")
```

`TRUE` confirms that the file was created. RStudio's Plots pane can change the layout when resized.
To inspect the intended export, preview at the same dimensions:

```r
view_ggchord(p_first, width = 10, height = 5, units = "in")
```

## 6. Troubleshoot a first run {#troubleshooting}
| Symptom | Check first |
| --- | --- |
| `first_sequences` or `p_first` is missing | Run the blocks in order from step 2, rather than copying only the last line |
| A function is unavailable | Check that the loaded version matches this tutorial with `packageVersion()` |
| Links are missing | Check that both IDs occur in the sequence table and that the selected link table is nonempty |
| The plot is small or labels differ | Preview at the explicit 10 × 5 in size; see [output and sizing](output.html) |
| You cannot find the saved file | Check the working directory with `getwd()` |

## Choose your next task {#first-steps}
You have prepared data, composed layers and saved a figure. Continue according to your next task:

- **Use your own records**: [import and adaptation](formats.html); BLAST users can continue with [standard tables and protein coordinates](alignments.html).
- **Compare complete real records**: [the FASTA, GFF3 and alignment workflow](real-data.html).
- **Make a plasmid map**: [pUC19 genes, primers and restriction annotations](plasmids.html).
- **Control the appearance**: [colours and guides](styling.html), [label placement](labels.html).
- **Look up a parameter**: [complete function reference](reference.html).
