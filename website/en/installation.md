<a id="requirements"></a>

## Install ggchord {#requirements}

Already installed ggchord? [Try the bundled example](#first-plot).

ggchord requires R ≥ 4.1.0 and ggplot2 ≥ 4.0.0. For a development preview, use the checkout used to build it. For a released version, select the version matching the documentation from [GitHub Releases](https://github.com/DangJem/ggchord/releases), then download and extract the source.

Set your R working directory to the ggchord source directory containing `DESCRIPTION`, then run:

```{.r}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_local(".", dependencies = NA, upgrade = "never", build = FALSE)
```

For scripts from an earlier version, see the [API migration guide](migration.md).

<a id="first-plot"></a>

## Draw a chord diagram with bundled data {#first-plot}

Copy the whole block into R and run it. The data are included with the package.

```r
library(ggchord)

p <- ggchord(data = seq_data_example) +
  geom_seq() +
  geom_seq_label() +
  coord_chord()
p
```

You will see four labelled arcs. Each arc represents a sequence, and its span reflects the full sequence length. Use `geom_seq()` to draw the backbones, `geom_seq_label()` to label them, and `coord_chord()` to select the circular arrangement.

<a id="first-steps"></a>

## Continue with your own data {#first-steps}

- [Add links and annotations](multi-sequence-start.md#links)
- [Choose another plot type](modes.md#modules)
- [Work with your own files](formats.md#choose-format)
- [Preview and save a plot](output.md#size)

<a id="environment"></a>

## Check the installed version {#environment}

```r
getRversion()
packageVersion("ggplot2")
packageVersion("ggchord")
find.package("ggchord")
```

<a id="support-evidence"></a>

## Which additional dependencies do examples use? {#support-evidence}

Each example lists any extra R packages or file decoders it needs. To read your own files, check the [supported formats](formats.md#choose-format).

<a id="troubleshooting"></a>

## Installation or plotting problems {#troubleshooting}

Check that the dependency named in the error is available and that R loads ggchord from the intended library. If a source install fails, keep the first compiler error. Use [Diagnostics and troubleshooting](diagnostics.md) to investigate plotting warnings and errors.

<span id="input"></span><span id="read"></span><span id="save"></span><span id="source"></span><span id="stable"></span><span id="verified-versions"></span>
