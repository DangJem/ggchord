<div class="lesson-intro">Install and confirm the environment, then choose one progressive module. Reading the entire reference is not a prerequisite.</div>

## Install and inspect the environment {#requirements}
R ≥4.1.0 and ggplot2 ≥4.0.0 are required. The formal installation is pinned to v0.16.5; current source examples include unreleased v0.16.6 work. Use the matching checkout to reproduce development figures.

```r
install.packages(c("remotes", "dplyr", "purrr", "tidyr"))
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

### Development checkout and support evidence {#support-evidence}
The formal release remains **v0.16.5**. This source checkout is the unreleased
**0.16.6.9000** RC development version. Its finite local checks used R 4.6.0 and
ggplot2 4.0.3 on arm64 macOS. The declared lower bounds (R 4.1.0, ggplot2 4.0.0)
have not been exercised by these checks. Earlier v0.16.5 Windows/macOS/Ubuntu
release and Ubuntu devel CI is historical evidence, not a v0.16.6 result.
Optional readers/adapters need their documented dependencies; full release tests
and other device/platform combinations remain unverified for this checkout.


## Choose a learning path {#first-steps}
1. [Module I · Chord: learn the three tables](chord-data.html), then import, backbones, layers and geometry.
2. [Module II · Collinear: reuse the chord inputs](collinear-data.html), then switch coordinates, route links and focus locally.
3. [Module III · Plasmid: start with one sequence and DNA](plasmid-data.html), then add layers and explore 13 reference cases.

Starting with chord follows ggchord's development order. See the [learning map](articles.html). The two-sequence exercise below is an optional environment check rather than a substitute for the modules.

## Which columns come first? {#input}
Sequences need accver/length; links need two IDs and four endpoints; genes need accver/start/end/strand. [Module I's first lesson](chord-data.html) explains real tables, optional fields and BLAST naming.

## Optional: check one complete plot {#first-plot}
Run the whole folded block, then inspect the figure. Bundled source-backed records need no external files or BLAST run.

{{prepare:firstdata,firstplot}}
{{figure:first}}

## Read the quick result {#read}
Arrows retain biological strands; colours identify sequences only. Ribbons are existing comparison intervals, not new homology results. This small linear example makes positions easy to read; the full chord module starts with backbones.

## Optional: save the quick result {#save}
Run this only after the block above has created p_first:

{{code:firstsave}}

Use the same 10×5 in dimensions for preview and export; see [output](output.html).

## Common execution problems {#troubleshooting}
For missing objects, run the preparation block. For missing functions, check the package version. getwd() determines the save directory. Optional readers need dependencies only in their own topics; basic modules do not require Bioconductor or Biopython.
