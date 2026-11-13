{{code:courseblockbw}}

{{prepare:fusiondata,fusionpanels,fusionlinear,fusionview}}

## Nested source regions retain their own namespaces {#nested-detail}

A detail may itself be a Fusion. This example uses two explicit windows on the unchanged pUC19 source. The inner source is the larger window; its child is the smaller window. The outer selection targets the inner source through `detail/center`, not an arbitrary rectangle in an image.

{{code:courseblockbx}}

![Nested details with explicit leaf paths](../assets/courseouter.png)

`composition_id` plus `region_row` connects derived regions/fronts/leaders to each original selection in `focus_region_data_by_composition`. Top-level `focus_region_data` remains the original outer input. Nested leaf names use paths such as `detail/center` and `detail/detail`; a composition container is never treated as one genomic axis. If an accession occurs in multiple leaves, select the full panel path and exact instance. `source_panel` can select a source leaf under `center`. Layout exports include actual mm boundaries and requested boxes/angles.

Each plot keeps its own layers, scales and statistical population; the compositor does not silently copy layers or recompute statistics. Origin wraps use the explicit source-length declaration below. The complete negative-strand CDS-to-genome-focus-to-real-character example is in `examples/m3-nested.R`. Registered characters retain their physical tile/error policy; a small nested allocation cannot make 60 letters readable merely by increasing DPI. Full multi-device publication acceptance and broader layout optimisation remain separate concerns.

## Explicit origin-crossing selection {#origin-detail}

Circular layout is not evidence of a circular source. Set `wrap = TRUE` on a
selection with `start > end`. Both selected leaf panels must retain the same
full source length. Character adapters additionally require `circular = TRUE`
and `source_data` with a full sequence or an exact `source_length`; without the
full string, characters remain unverified. At most one source tour is allowed,
so inverse source-to-column mapping stays unique. Negative strands wrap in
their own direction and retain gaps. Different focused instances are never
joined implicitly.

Run the complete real pUC19 workflow in `examples/origin-detail.R`. It extracts
82 observed bases (2640–2686, then 1–35) from the unchanged bundled FASTA.
`focus_region_data` retains the original single row. Derived tables identify
two pieces with `region_piece`, sharing one `region_row`; only the two external
boundaries receive leaders. The seam is not another selection or relationship.

{{code:courseblockby}}

![An explicit origin-crossing selection](../assets/courseorigindetail.png)
