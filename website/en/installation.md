## Requirements {#requirements}
Use R ≥ 4.1.0 and ggplot2 ≥ 4.0.0. ggchord itself needs no compiled code.
Optional interval adapters use Bioconductor packages; spreadsheet imports need
readxl. Install these only for the input formats you actually use.

## Choose a source {#source}
The GitHub repository is the installation source. A development checkout can
contain changes newer than a formal release; record its commit for reproducibility.
A CRAN release is not claimed by this guide.

```r
install.packages("remotes")
remotes::install_github("DangJem/ggchord")
library(ggchord)
packageVersion("ggchord")
```

For the checked-out 0.16.5 source, run from the repository root:

```r
install.packages("pkgload")
pkgload::load_all(".")
```

Repository paths used in the [file-import example](data.html#import) are not
installed package data. The introductory examples select minimal columns from bundled data and work after installation.

## First steps {#first-steps}
Read the [data contract](data.html), copy the [shared setup](index.html#quick-start),
then choose [a coordinate system](modes.html). Save with explicit width and height;
RStudio's Plots pane can rearrange labels when resized.
