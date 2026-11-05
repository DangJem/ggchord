## Design references and acknowledgements

ggchord draws on these open-source projects and public tools for graphical grammar, genomic representation and figure preview:

- [ggplot2](https://ggplot2.tidyverse.org/): composing layers while separating mappings, scales, geometry, coordinates and themes.
- [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/): query/subject table naming that helps trace alignment endpoints and sequence identifiers.
- [Circos](https://circos.ca/) and [circlize](https://github.com/jokergoo/circlize): circular information hierarchy, ribbons and readable ordering of multiple tracks.
- [clinker](https://github.com/gamcil/clinker), [gggenes](https://github.com/wilkox/gggenes) and [genoPlotR](https://genoplotr.r-forge.r-project.org/): readable gene arrows, gene clusters and linear synteny comparisons.
- [DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) and [gggenomes](https://thackl.github.io/gggenomes/): annotation placement, multi-table input, feature stacking and local comparison workflows.
- [ggview](https://github.com/idmn/ggview): previewing figures at their intended export dimensions; ggchord previews with `view_ggchord()` and saves with standard `ggsave()`.

These projects inform design and comparison. Consult ggchord's function documentation and `DESCRIPTION` for its supported features, parameters and runtime dependencies.
