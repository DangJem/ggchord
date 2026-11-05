## Download and install {#requirements}

Use R ≥ 4.1.0 and ggplot2 ≥ 4.0.0. Run these commands in R. CRAN supplies ggchord and its required dependencies.

```r
install.packages("ggchord")
library(ggchord)
packageVersion("ggchord")
```

The last line reports your installed version. **This site targets development version 0.16.7.9000**; CRAN provides 0.16.6 (checked 2026-11-05). They differ in bundled protein objects, dedicated cross-row border width, the circular zero-rotation reference and some default appearances.

The development changes are not yet pushed or released. Install from the local source checkout matching this site; do not assume remote main includes them. From the repository root run:

```r
install.packages(".", repos = NULL, type = "source")
library(ggchord)
stopifnot(as.character(packageVersion("ggchord")) == "0.16.7.9000")
```

Install dependencies first, then use a fresh R session. A version number alone cannot distinguish commits; the site `build-manifest.json` and audit record retain source SHA. Readers using CRAN should use the [0.16.6 offline tutorials and reference](https://stat.ethz.ch/CRAN/web/packages/ggchord/index.html), and avoid development protein examples and `nonadjacent_linewidth`.

### Verified differences {#verified-versions}

The three introductory plotting sequences were actually rendered in six fresh R processes: each lesson on CRAN 0.16.6 and on this development checkout. This verifies the basic plotting calls, not identical appearances or all later lessons. The following probes were also run against both namespaces.

| Capability | 0.16.6 | 0.16.7.9000 |
| --- | --- | --- |
| `geom_link_ribbon(nonadjacent_linewidth=...)` | Unavailable | Available |
| Explicit `read_fasta_lengths(sequence_unit=...)` | Unavailable | Available |
| Bundled native protein data | Unavailable | Available |

The stable source archive SHA256 and exact process results are retained in the [version audit](../assets/source/tools/teaching-audit/version-check.json).

### Source downloads {#source}

Source packages are available from [CRAN](https://cran.r-project.org/web/packages/ggchord/index.html) and [GitHub v0.16.6](https://github.com/DangJem/ggchord/releases/tag/v0.16.6).

### Optional dependencies {#support-evidence}

Getting started and the basic chord lessons require only ggchord. The case study lists its extra dependencies.

## Confirm the installation {#troubleshooting}

Load the package in a new R session with `library(ggchord)`. If installation or loading fails, check the R or dependency version named in the error.

## Next steps {#first-steps}

- [Basic input tables](chord-data.html) shows real sequence, link and feature rows and explains each column.
- [Plotting concepts](modes.html) introduces backbones, layers and coordinate systems.
- <span id="first-plot"></span><span id="input"></span><span id="read"></span>[Draw sequence backbones](chord-start.html) opens the chord chapter.
- <span id="save"></span>[Save plots](output.html) explains dimensions and file formats.
