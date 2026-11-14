## Install from GitHub {#requirements}

Install ggchord and its required dependencies from the project repository. Source installation requires a working R development environment.

```r
install.packages("remotes")
remotes::install_github("DangJem/ggchord", dependencies = NA)
```

## Environment and installation checks {#environment}

Requires R ≥ 4.1.0 and ggplot2 ≥ 4.0.0. Inspect the packages actually used by the current session and include this output when reporting a problem.

```r
getRversion()
packageVersion("ggplot2")
packageVersion("ggchord")
find.package("ggchord")
```

## Load and draw {#first-plot}

```r
library(ggchord)

p <- ggchord(seq_data_example) +
  geom_seq() +
  geom_seq_label() +
  coord_chord()
p
```

Success produces four arcs with different lengths and sequence IDs, without a build error. For missing functions, inspect the library path and loaded package before checking the matching [function reference](reference/index.html).

## Optional packages and input files {#support-evidence}

Examples that use dplyr, Biostrings or a format decoder declare that dependency beside the code. Install it only for the corresponding workflow. Built-in data works offline; file examples use the [input archive](formats.html).

## If installation fails {#troubleshooting}

| Symptom | First check | Action |
| --- | --- | --- |
| Repository or index cannot download | First network error, proxy and GitHub access | Restore access and retry without changing data |
| Missing or outdated dependency | Dependency named in the error and minimum R | Install/update that dependency, then retry |
| Source compilation fails | First compiler error and platform build tools | Install the appropriate R development tools and retain the log |
| Installation succeeds but functions are missing | `find.package()` and the loaded library | Restart R, load the correct library and retry the minimal plot |

For input and plotting problems, see [diagnostics](diagnostics.html).

## Related sections {#first-steps}

[Plotting concepts](modes.html) · [Input tables](data.html) · [Save a figure](output.html)


<span id="input"></span><span id="read"></span><span id="save"></span><span id="source"></span><span id="stable"></span><span id="verified-versions"></span>
