Check units, layers and parameter boundaries before changing coordinates. This table follows current Rd and source contracts. Support means an implementation exists, not that every parameter combination has been tested. Runtime evidence is recorded in the version audit; static documentation is not device acceptance.

## Layers across three coordinate systems {#layers}

| Layer | Chord | Collinear | Single circular | Evidence and conditions |
| --- | --- | --- | --- | --- |
| `geom_seq()` | Supported | Conditional | Supported | [Reference](reference/geom_seq.html); linear mode rejects explicit radius, curvature and ring |
| `geom_seq_label()` | Supported | Conditional | Supported | [Reference](reference/geom_seq_label.html); no explicit radial label radius in linear mode |
| `geom_gene()` | Conditional | Conditional | Conditional | [Reference](reference/geom_gene.html); bp and known strands required; use generic features for aa |
| `geom_feature()` | Supported | Supported | Supported | [Reference](reference/geom_feature.html); IDs, coordinates and units must match |
| Fixed gene/feature labels | Supported | Conditional | Supported | [Feature labels](reference/geom_feature_label.html); radial/tangent become vertical/horizontal in linear mode |
| Automatic gene/feature labels | Supported | Conditional | Supported | [Automatic labels](reference/geom_feature_label_repel.html); linear rows reject external=FALSE |
| `geom_link_ribbon()` | Supported | Conditional | Unsupported | [Reference](reference/geom_link_ribbon.html); valid links required; linear adjacency follows display rows, same-row segmented ribbons are omitted and reported |
| `geom_link_line()` | Supported | Conditional | Unsupported | [Reference](reference/geom_link_line.html); explicit same-row segmented point links unsupported; single circular mode rejects ribbon/link layers |
| `geom_seq_region()` | Supported | Supported | Supported | [Reference](reference/geom_seq_region.html); valid intervals, IDs and side |
| `geom_seq_scale()` | Unsupported | Supported | Unsupported | [Reference](reference/geom_seq_scale.html); straight scale bars require collinear mode; use local axes elsewhere |
| `geom_track_bar/line/tile()` | Supported | Supported | Supported | [Track reference](reference/geom_track_tile.html); explicit position_track required, declarations sharing an ID must agree |
| `geom_seq_break()` | Unsupported | Conditional | Unsupported | [Reference](reference/geom_seq_break.html); requires identifiable window boundaries; does not edit source sequences |
| `geom_seq_center_label()` | Conditional¹ | Unsupported | Supported | [Reference](reference/geom_seq_center_label.html); single sequence; tutorials use circular mode |
| `geom_primer()` / automatic primer labels | Conditional¹ | Unsupported | Conditional | [Reference](reference/geom_primer.html); circular and nucleotide units |
| `geom_restriction_site()` | Conditional¹ | Unsupported | Conditional | [Reference](reference/geom_restriction_site.html); circular, nucleotide units and valid cut sites |

`geom_feature_plasmid()` is a [generic feature preset](reference/geom_feature_plasmid.html), not a separate geometry family. It primarily serves circular plasmid lessons. Start with generic features in other coordinates rather than assuming all circular layout settings transfer unchanged.

¹ Targeted builds of single-sequence pUC19 under `coord_chord()` accept centre labels, primers and restriction layers, although existing collinear documentation calls them circular-only. This records a contract/implementation discrepancy, not visual acceptance or an official cross-mode promise. Automatic primer labels were not separately run. Tutorials continue to use circular mode; runtime code is unchanged.

## Understand conditional support {#conditions}

`ggchord_validate_coord_contract()` in `R/coord-collinear.R` explicitly checks circular parameters and conflicts involving automatic orientation or segments. Runtime validation also checks general data contracts and nucleotide-only restrictions. Current source files accompany the site under `assets/source/`; reference pages are generated from current roxygen/Rd.

A Position controls local placement; it cannot convert an unsupported layer into a different coordinate layer. Outside/inside means above/below in collinear mode and follows the outward/inward local normal in circular mode. Re-read the result when transferring settings. See [four direction concepts](collinear-layout.html#orientation) for display direction versus biological strand.

## Versions and reproducibility {#versions}

This table targets development version 0.16.7.9000; [installation](installation.html) explains local source installation and stable documentation. Single circular mode requires one displayed sequence. Circular protein layout is a display choice, not a claim of molecular circularity. Project mixed bp/aa inputs first or plot them separately.
