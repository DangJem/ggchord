Check units, layers and parameters before changing coordinates. This table combines current Rd/source contracts with actual renders. Support describes an implementation, not every parameter combination. Layouts can use native bp/aa units; genes, primers and restriction layers retain nucleotide constraints.

<a id="layers"></a>

## Layers across six coordinate systems {#layers}

| Layer/reference | Chord | Collinear | Single circular | Radial | Sankey | Dotplot | Conditions and boundaries |
| --- | --- | --- | --- | --- | --- | --- | --- |
| [`geom_seq()`](../../man/geom_seq.Rd) | Supported | Conditional | Supported | Conditional | Conditional | Unsupported | Linear mode rejects explicit radius/curvature/ring |
| [`geom_seq_label()`](../../man/geom_seq_label.Rd) | Supported | Conditional | Supported | Conditional | Conditional | Unsupported | No explicit radial label radius in linear mode |
| [`geom_feature(preset = 'gene')`](../../man/geom_feature.Rd) | Conditional | Conditional | Conditional | Conditional | Conditional | Unsupported | Nucleotide units and source strands; use generic features for aa |
| [`geom_feature()`](../../man/geom_feature.Rd) | Supported | Supported | Supported | Conditional | Conditional | Unsupported | Valid IDs, intervals and matching native units |
| [`geom_feature(preset = 'plasmid')`](../../man/geom_feature.Rd) | Supported | Supported | Supported | Conditional | Conditional | Unsupported | Generic feature preset; circular Position defaults need review when transferred |
| [`geom_gene_label()`](../../man/geom_gene_label.Rd) | Supported | Conditional | Supported | Conditional | Conditional | Unsupported | Share the gene Position; radial/tangent become vertical/horizontal in linear mode |
| [`geom_feature_label()`](../../man/geom_feature_label.Rd) | Supported | Conditional | Supported | Conditional | Conditional | Unsupported | Share the feature Position; radial/tangent become vertical/horizontal in linear mode |
| [`geom_gene_label_repel()`](../../man/geom_gene_label_repel.Rd) | Supported | Conditional | Supported | Conditional | Conditional | Unsupported | Linear rows reject external=FALSE; inspect hiding at target size |
| [`geom_feature_label(placement = 'repel')`](../../man/geom_feature_label.Rd) | Supported | Conditional | Supported | Conditional | Conditional | Unsupported | Linear rows reject external=FALSE; inspect fitting and hiding |
| [`geom_link_ribbon()`](../../man/geom_link_ribbon.Rd) | Supported | Conditional | Unsupported | Conditional | Conditional | Unsupported | Linear adjacency follows display rows; segmented same-row ribbons are omitted and reported |
| [`geom_link_line()`](../../man/geom_link_line.Rd) | Supported | Conditional | Unsupported | Conditional | Conditional | Unsupported | Point coordinates required; explicit segmented same-row links unsupported |
| [`geom_seq_region()`](../../man/geom_seq_region.Rd) | Supported | Supported | Supported | Conditional | Conditional | Unsupported | Valid source intervals, IDs and local side |
| [`geom_seq_scale()`](../../man/geom_seq_scale.Rd) | Unsupported | Supported | Unsupported | Unsupported | Conditional | Unsupported | Straight scale bars require collinear mode |
| [`geom_seq_break()`](../../man/geom_seq_break.Rd) | Unsupported | Conditional | Unsupported | Conditional | Conditional | Unsupported | Collinear window boundaries required; does not edit source sequence |
| [`geom_seq_center_label()`](../../man/geom_seq_center_label.Rd) | Conditional | Unsupported | Conditional | Unsupported | Unsupported | Unsupported | One displayed sequence; length uses its native bp/aa unit |
| [`geom_primer()`](../../man/geom_primer.Rd) | Conditional | Unsupported | Conditional | Unsupported | Unsupported | Unsupported | One displayed DNA sequence and valid binding coordinates |
| [`geom_primer_label_repel()`](../../man/geom_primer_label_repel.Rd) | Conditional | Unsupported | Conditional | Unsupported | Unsupported | Unsupported | One displayed DNA sequence; share primer placement |
| [`geom_restriction_site()`](../../man/geom_restriction_site.Rd) | Conditional | Unsupported | Conditional | Unsupported | Unsupported | Unsupported | One displayed DNA sequence and actual computed cleavage sites |
| [`geom_track()`](../../man/geom_track.Rd) | Supported | Supported | Supported | Conditional | Conditional | Unsupported | Source-space Stat/Position; [base matrix](tracks-statistics.md#tracks) and [estimation combinations](tracks-statistics.md#estimates) |
| [`geom_track_bar()`](../../man/geom_track.Rd) | Supported | Supported | Supported | Conditional | Conditional | Unsupported | Count observations; use col for supplied values |
| [`geom_track_col()`](../../man/geom_track.Rd) | Supported | Supported | Supported | Conditional | Conditional | Unsupported | Default or named band; numeric heights and declared domains |
| [`geom_track_line()`](../../man/geom_track.Rd) | Supported | Supported | Supported | Conditional | Conditional | Unsupported | Default or named band; closure must be requested and missing samples respected |
| [`geom_track_tile()`](../../man/geom_track.Rd) | Supported | Supported | Supported | Conditional | Conditional | Unsupported | Default or named band; no numeric height required |
| [`geom_alignment_segment()`](../../man/geom_alignment_segment.Rd) | Unsupported | Unsupported | Unsupported | Unsupported | Unsupported | Identity HSP segments | Dedicated two-axis geometry |


23 representative layers were rendered in the three main modes using 69 explicitly synthetic geometry fixtures: 59 succeeded and 10 were rejected by contract, including automatic primer labels. Centre labels support bp/aa in single-sequence chord/circular views; primers and restriction layers require DNA. Help now clarifies these boundaries. This verifies mode boundaries, not all parameter combinations, biological recognition or devices. See the [complete runtime record](../../tools/teaching-audit/layer-compatibility.csv).

<a id="conditions"></a>

## Understand conditional support {#conditions}

`ggchord_validate_coord_contract()` in `R/coord-collinear.R` explicitly checks circular parameters and conflicts involving automatic orientation or segments. Runtime validation also checks general data contracts and nucleotide-only restrictions. Current source files accompany the site under `assets/source/`; reference pages are generated from current roxygen/Rd.

A Position controls local placement; it cannot convert an unsupported layer into a different coordinate layer. Outside/inside means above/below in collinear mode and follows the outward/inward local normal in circular mode. Re-read the result when transferring settings. See [four direction concepts](collinear-layout.md#orientation) for display direction versus biological strand.

<a id="versions"></a>

## Input units and script reuse {#versions}

See [installation](installation.md) for setup instructions and [migration](../../MIGRATION.md) when updating existing scripts. Single circular mode requires one displayed sequence. Circular protein layout is a display choice, not a claim of molecular circularity. Project mixed bp/aa inputs first or plot them separately.


<a id="new-views"></a>

## New coordinate systems {#new-views}

The [Dotplot/Radial support matrix](dotplot-start.md#support) records the new contracts. Dotplot uses a dedicated identity HSP mark. Radial adapts the existing feature/Track frames with cyclic adjacency; circular-only layers, shared straight scales and branches are explicitly rejected. The current layer probes are recorded separately from the earlier three-view evidence.


| Layer/reference | Dotplot | Radial |
| --- | --- | --- |
| [`geom_alignment_segment()`](../../man/geom_alignment_segment.Rd) | Supported identity HSP segment | Unsupported |


<a id="fusion"></a>

## Static cross-panel composition {#fusion}

[Fusion](fusion.md) accepts registered Circular/Chord/Collinear/Sankey/Radial plots, nested Fusions and explicitly registered Cartesian character details. Eight named sides and explicit angle/box placements share the source-frame resolver. Each panel keeps its existing layer contract. Cross-panel interval
fronts use one-based inclusive boundaries with continuous curve mapping;
fully visible same-unit source intervals and unambiguous panel/instance routes
are required. Generic Fusion scales affect cross-panel layers only. Dotplot, facets and implicit third-party coordinate adapters remain unregistered. Physical panels are remeasured at save time; this does not imply
global collision avoidance or general font scaling.


<a id="aligned-cells"></a>

## Local aligned-character cells {#aligned-cells}

| Layer/reference | Chord / Collinear / Circular / Dotplot / Radial / Sankey | Ordinary Cartesian ggplot |
| --- | --- | --- |
| [`geom_residue()`](../../man/geom_residue.Rd) | Unsupported biological backbone projection | Supported local character cells with identity/reverse position scales |

The [aligned-character workflow](aligned-sequences.md#aligned-characters) uses ordinary ggplot x/y coordinates. Source positions and alignment columns are separate fields. Final physical glyph measurement selects readable text or explicit tile/error behavior; it does not change source data or shrink fonts. Explicit source registration supports Fusion/magnify details. Optional source-backed CDS translation uses Biostrings and the existing CDS map; no alignment or missing phase is inferred. See [CDS translation](aligned-sequences.md#cds-translation) and [character details](fusion-details.md#character-detail).

<a id="adaptive"></a>

## Optional target-size policies {#adaptive}

[`adapt_ggchord()`](styling.md#adaptive) supports individual Chord, Circular,
Collinear, Sankey, Radial and Dotplot plots, and opted-in Fusion panels.
It scales default theme styles and arranges native guides under explicit user
constraints; no policy is applied by default. Facets, arbitrary custom grobs,
outer Fusion allocations and a global label solver are outside this scope.

[Magnification](fusion-details.md#multiple-details) supports independent registered details with explicit per-row targets, eight named sides, nested leaf paths and free boxes. Deliberate overlap is opt-in; global avoidance is separate.

Fusion also supports Sankey/Radial, layered compositions and eight named sides or explicit boxes. Leaf-panel paths preserve independent sources, and deliberate overlaps export their paint order. Dotplot still requires independent axis registration. See [Advanced boxes and nesting](fusion-layout.md#advanced-placement).
