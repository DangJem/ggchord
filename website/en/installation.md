ggchord uses sequence coordinates to show where features lie and how sequences match. Suppose you have two phage genomes and a local alignment between them. You can start with their lengths, add the aligned intervals and then see which coding sequences occupy those regions. You can also lay the same genomes out in parallel to compare their positions from left to right.

Map an annotation category to fill or an alignment score to colour, and the scale supplies a legend for reading those values. A [quantitative track](tracks.md) can show GC content along a backbone. A [Fusion composition](fusion.md) can connect one interval in the overview to its own detail panel.

<a id="requirements"></a>

## Install ggchord {#requirements}

For the CRAN release, run the following in R. CRAN currently supplies ggchord 0.28.2. This manual uses the 0.50.0 development interface.

```{.r}
install.packages("ggchord")
```

To run the examples in this manual, [download the matching package source](../../) and unzip it. Open R in the extracted ggchord folder that contains `DESCRIPTION`. With R 4.1.0 or later, install the source and the dplyr package used to select example records. The source installer resolves the required dependencies, including ggplot2 4.0.0 or later.

```{.r}
install.packages(c("remotes", "dplyr"))
remotes::install_local(".", dependencies = NA, upgrade = "never", build = FALSE)
```

Restart R and run `packageVersion("ggchord")`. For this manual the result should be 0.50.0.

<a id="layers"></a>

## Read a sequence comparison in three steps {#layers}

The first picture shows two complete genomes. The longer arc is MT108731.1 (64,323 bp) and the shorter arc is MT118296.1 (32,090 bp). In the second picture, one ribbon joins the intervals reported by their first local alignment. In the third, green arrows mark coding sequences at their annotated positions. Follow an arrow along a backbone to read its strand, or look across the ribbon to see what the matching interval contains.

{{preview:first-backbones}}

{{preview:first-chord}}

{{preview:first-features}}

<a id="layouts"></a>

## Read the same alignment on arcs and on parallel lines {#layouts}

In a chord diagram you can follow either end of the ribbon around the circle. In a collinear plot you can compare the interval with the start and end of each genome along a straight backbone. The two pictures use the same genomes, alignment and endpoints.

{{preview:views-chord}}

{{preview:views-collinear}}

[Draw your first chord diagram with the bundled data](chord.md#backbones). To learn what to put in your own tables, read the [data guide](data.md). To choose another view, open the [plotting guide](modes.md).

<a id="first-plot"></a>

<a id="example-first-plot"></a>

To run the example with one alignment directly, go to [Join one pair of aligned intervals](chord.md#draw).
