## Install from GitHub {#requirements}

Install ggchord and its required dependencies from the project repository. Source installation requires a working R development environment.

```r
install.packages("remotes")
remotes::install_github("DangJem/ggchord", dependencies = NA)
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

## Optional packages and input files {#support-evidence}

Examples that use dplyr, Biostrings or a format decoder declare that dependency beside the code. Install it only for the corresponding workflow. Built-in data works offline; file examples use the [input archive](formats.html).

## If installation fails {#troubleshooting}

Read the first dependency or compiler error, restart R after installation, then use `help(package = "ggchord")` to inspect available functions. See [diagnostics](diagnostics.html) for plotting errors.

## Related sections {#first-steps}

[Plotting concepts](modes.html) · [Input tables](data.html) · [Save a figure](output.html)


<span id="input"></span><span id="read"></span><span id="save"></span><span id="source"></span><span id="stable"></span><span id="verified-versions"></span>
