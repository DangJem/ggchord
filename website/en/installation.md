## Download and install {#requirements}

Use R ≥ 4.1.0 and ggplot2 ≥ 4.0.0. Run these commands in R. CRAN supplies ggchord and its required dependencies.

```r
install.packages("ggchord")
library(ggchord)
packageVersion("ggchord")
```

The last line reports the installed version. This tutorial uses the public interface from 0.16.6.

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
