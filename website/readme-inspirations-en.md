## Design References

ggchord draws on established open-source projects and scientific visualization approaches. These references have informed its graphical grammar, genomic layouts, annotation strategies, and figure composition. Each project is acknowledged for its specific contributions to the design process; this does not imply direct reuse of its algorithms or implementations.

### Graphical grammar and API design

- [ggplot2](https://ggplot2.tidyverse.org/) — Provides the foundational layered grammar, separating data mappings, geometries, statistics, positions, scales, guides, coordinates, and themes. ggchord builds on this architecture to support composable genomic visualizations and consistent user-defined styling.
- [gggenomes](https://thackl.github.io/gggenomes/) — Offers important references for organizing genomic sequences and annotations, strand-aware layouts, feature stacking, sequence orientation, anchoring, and focus-based workflows. ggchord follows its own reference-coordinate and layer-data contracts while drawing on these domain-specific design concepts.
- [gggenes](https://wilkox.org/gggenes/) — Informs directional gene-arrow geometry, strand-aware positioning, and the distinction between genomic interval lengths and physical drawing dimensions.

### Circular and comparative genomic visualization

- [Circos](https://circos.ca/) — Inspires the visual organization of circular genomic information, including sequence sectors, layered annotations, and ribbons representing relationships between genomic intervals.
- [circlize](https://github.com/jokergoo/circlize) — Serves as a reference for flexible circular track organization, ring arrangement, and multi-layer genomic visualization.
- [clinker](https://github.com/gamcil/clinker) — Provides a visual reference for gene-cluster comparisons, directional gene arrows, and the readable presentation of genomic similarity and synteny.
- [genoPlotR](https://genoplotr.r-forge.r-project.org/) — Offers design references for linear comparative genomics, aligned sequence regions, genomic coordinate presentation, and relationships across multiple genomic segments.

### Feature annotations and multi-scale visualization

- [DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) — Informs the presentation of genomic features, annotation geometry, label placement, and collision management in densely annotated regions.
- [ggmagnify](https://github.com/hughjonesd/ggmagnify) — Provides a reference for linking overview figures to magnified details. ggchord develops this concept through explicit source intervals and connections between independently constructed views.
- [ggmsa](https://bioconductor.org/packages/ggmsa/) — Provides references for aligned nucleotide and amino-acid visualization, compact residue cells, and chemistry-based colouring. Selected residue palette values are adapted with attribution in [NOTICE](../inst/NOTICE); sequence registration and rendering are independently implemented.
- [ggview](https://github.com/idmn/ggview) — Informs the workflow of previewing figures at their intended output dimensions. ggchord provides size-aware previews through `view_ggchord()` while retaining standard `ggsave()` output.

### Data conventions and related research

- [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) — Provides established query/subject terminology reflected in ggchord's alignment data conventions. ggchord visualizes existing alignment results rather than performing sequence alignment.
- [trackViewer](https://bioconductor.org/packages/trackViewer/) — Serves as a research reference for genomic signal tracks, interval-based measurements, and quantitative annotation displays.
- [ggseqalign](https://bioconductor.org/packages/ggseqalign/), [ggseqlogo](https://omarwagih.github.io/ggseqlogo/), [Jalview](https://www.jalview.org/), and [ESPript](https://espript.ibcp.fr/ESPript/ESPript/) — Have been examined as related approaches to alignment visualization, sequence logos, residue annotations, and publication-oriented sequence figures. Their inclusion indicates research and comparison, not implemented integration or adoption of their rendering algorithms.
