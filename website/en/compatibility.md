Check units, layers and parameters before changing coordinates. This table combines current Rd/source contracts with actual renders. Support describes an implementation, not every parameter combination. Layouts can use native bp/aa units; genes, primers and restriction layers retain nucleotide constraints.

## Layers across three coordinate systems {#layers}

| Layer/reference | Chord | Collinear | Single circular | Conditions and boundaries |
| --- | --- | --- | --- | --- |
| [`geom_seq()`](reference/geom_seq.html) | Supported | Conditional | Supported | Linear mode rejects explicit radius/curvature/ring |
| [`geom_seq_label()`](reference/geom_seq_label.html) | Supported | Conditional | Supported | No explicit radial label radius in linear mode |
| [`geom_gene()`](reference/geom_gene.html) | Conditional | Conditional | Conditional | Nucleotide units and source strands; use generic features for aa |
| [`geom_feature()`](reference/geom_feature.html) | Supported | Supported | Supported | Valid IDs, intervals and matching native units |
| [`geom_feature_plasmid()`](reference/geom_feature_plasmid.html) | Supported | Supported | Supported | Generic feature preset; circular Position defaults need review when transferred |
| [`geom_gene_label()`](reference/geom_gene_label.html) | Supported | Conditional | Supported | Share the gene Position; radial/tangent become vertical/horizontal in linear mode |
| [`geom_feature_label()`](reference/geom_feature_label.html) | Supported | Conditional | Supported | Share the feature Position; radial/tangent become vertical/horizontal in linear mode |
| [`geom_gene_label_repel()`](reference/geom_gene_label_repel.html) | Supported | Conditional | Supported | Linear rows reject external=FALSE; inspect hiding at target size |
| [`geom_feature_label_repel()`](reference/geom_feature_label_repel.html) | Supported | Conditional | Supported | Linear rows reject external=FALSE; inspect fitting and hiding |
| [`geom_link_ribbon()`](reference/geom_link_ribbon.html) | Supported | Conditional | Unsupported | Linear adjacency follows display rows; segmented same-row ribbons are omitted and reported |
| [`geom_link_line()`](reference/geom_link_line.html) | Supported | Conditional | Unsupported | Point coordinates required; explicit segmented same-row links unsupported |
| [`geom_seq_region()`](reference/geom_seq_region.html) | Supported | Supported | Supported | Valid source intervals, IDs and local side |
| [`geom_seq_scale()`](reference/geom_seq_scale.html) | Unsupported | Supported | Unsupported | Straight scale bars require collinear mode |
| [`geom_seq_break()`](reference/geom_seq_break.html) | Unsupported | Conditional | Unsupported | Collinear window boundaries required; does not edit source sequence |
| [`geom_seq_center_label()`](reference/geom_seq_center_label.html) | Conditional | Unsupported | Conditional | One displayed sequence; length uses its native bp/aa unit |
| [`geom_primer()`](reference/geom_primer.html) | Conditional | Unsupported | Conditional | One displayed DNA sequence and valid binding coordinates |
| [`geom_primer_label_repel()`](reference/geom_primer_label_repel.html) | Conditional | Unsupported | Conditional | One displayed DNA sequence; share primer placement |
| [`geom_restriction_site()`](reference/geom_restriction_site.html) | Conditional | Unsupported | Conditional | One displayed DNA sequence and actual computed cleavage sites |
| [`geom_track_bar()`](reference/geom_track_tile.html) | Supported | Supported | Supported | Explicit position_track; numeric heights and declared domains |
| [`geom_track_line()`](reference/geom_track_tile.html) | Supported | Supported | Supported | Explicit position_track; closure must be requested and missing samples respected |
| [`geom_track_tile()`](reference/geom_track_tile.html) | Supported | Supported | Supported | Explicit position_track; no numeric height required |

All 21 geoms were rendered in all three modes using 63 explicitly synthetic geometry fixtures: 53 succeeded and 10 were rejected by contract, including automatic primer labels. Centre labels support bp/aa in single-sequence chord/circular views; primers and restriction layers require DNA. Help now clarifies these boundaries. This verifies mode boundaries, not all parameter combinations, biological recognition or devices. See the [complete runtime record](../assets/source/tools/teaching-audit/layer-compatibility.csv).

## Understand conditional support {#conditions}

`ggchord_validate_coord_contract()` in `R/coord-collinear.R` explicitly checks circular parameters and conflicts involving automatic orientation or segments. Runtime validation also checks general data contracts and nucleotide-only restrictions. Current source files accompany the site under `assets/source/`; reference pages are generated from current roxygen/Rd.

A Position controls local placement; it cannot convert an unsupported layer into a different coordinate layer. Outside/inside means above/below in collinear mode and follows the outward/inward local normal in circular mode. Re-read the result when transferring settings. See [four direction concepts](collinear-layout.html#orientation) for display direction versus biological strand.

## Versions and reproducibility {#versions}

This table targets development version 0.16.7.9000; [installation](installation.html) explains local source installation and stable documentation. Single circular mode requires one displayed sequence. Circular protein layout is a display choice, not a claim of molecular circularity. Project mixed bp/aa inputs first or plot them separately.
