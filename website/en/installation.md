## Choose a version and dependencies {#requirements}

These lessons use the **0.27.0.9000 development API**. Choose a matching source; the stable 0.16.6 package has a different teaching scope.

R ≥ 4.1.0; ggplot2 ≥ 4.0.0.

```r
install.packages("remotes")
```

### Stable release {#stable}

The latest recorded formal release is 0.16.6. Its source and versioned help are available at [GitHub v0.16.6](https://github.com/DangJem/ggchord/releases/tag/v0.16.6). Use its bundled help after installing the fixed tag; these lessons include newer APIs.

```r
remotes::install_github("DangJem/ggchord@v0.16.6", dependencies = NA)
```

### Development snapshot {#source}

The development runtime is pinned to public audit commit `0f5dd96da57876da74cd1234c7de9cf5579c80b8`. This local documentation refactor does not change runtime behavior. The separate documentation revision, build time and input hashes are recorded in the [build manifest](../build-manifest.json). Installation downloads source and dependencies.

```r
runtime_commit <- "0f5dd96da57876da74cd1234c7de9cf5579c80b8"
remotes::install_github(
  paste0("DangJem/ggchord@", runtime_commit),
  dependencies = NA, upgrade = "never"
)
```

<span id="verified-versions"></span>

## Confirm in a fresh session and draw {#troubleshooting}

```r
library(ggchord)
stopifnot(as.character(packageVersion("ggchord")) == "0.27.0.9000")
stopifnot(packageDescription("ggchord")$RemoteSha ==
  "0f5dd96da57876da74cd1234c7de9cf5579c80b8")
p <- ggchord(seq_data_example) +
  geom_seq() +
  geom_seq_label() +
  coord_chord()
p
```

### Optional dependencies {#support-evidence}

The first plot needs only ggchord. Data lessons explicitly use dplyr; import recipes name their optional decoders.

## Continue learning {#first-steps}

<span id="first-plot"></span><span id="input"></span><span id="read"></span>

- [Input tables](chord-data.html)
- [First multi-sequence plot](chord-start.html)
- <span id="save"></span>[Save at target size](output.html)
