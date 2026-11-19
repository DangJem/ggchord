<a id="first-plot"></a>

## Draw your first chord diagram {#first-plot}

The sequence table supplies the lengths of four phage genome sequences. The first four rows of the alignment table supply the interval pairs joined by ribbons. With ggchord installed, run these five lines in R.

{{example:first-plot}}

<a id="requirements"></a>

## Install the package from source {#requirements}

The examples on this site use the source version supplied with this manual. It requires R 4.1.0 or later and ggplot2 4.0.0 or later. Open R in the ggchord source folder containing `DESCRIPTION`, install the remotes package if needed, and run the local installer.

```{.r}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_local(".", dependencies = NA, upgrade = "never", build = FALSE)
```

Restart R after installation. Check the loaded version with `packageVersion("ggchord")`. To use your own sequence lengths and interval pairs, see the [chord guide](chord.md).
