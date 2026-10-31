#' Layout results, units and ownership
#'
#' [export_ggchord_layout()] is the recommended interface for consuming
#' layer-aware geometry tables. [get_chord_layout()] retains the complete
#' layout list and its documented geometry, bounds, colour and label members.
#' Neither entrance is deprecated. [view_ggchord()] returns a temporary image
#' path; [diagnose_ggchord()] returns evidence for one requested render.
#'
#' @section Geometry and provenance:
#' Selected export components are data frames; absent components may be empty
#' data frames without columns. Nonempty geometry has `layer_id`, `component`
#' and `source_row`. An ID identifies one layer occurrence in that plot;
#' editing/reordering layers can change IDs. A source row indexes the resolved
#' input for that layer/component, not necessarily a physical file line.
#' Generated decorations use missing source rows. Aggregates may carry
#' `source_rows` membership instead of identifying one observation.
#' `original_data = TRUE` includes resolved inputs for selected nonempty
#' components, after role mapping and input adaptation. It is not a raw-file
#' archive. File line/record provenance and conversion attributes retain their
#' separate meanings from [ggchord_data_contract].
#'
#' Cartesian geometry and bounding boxes use native layout data units.
#' Genomic source coordinates remain one-based inclusive bp; they are not
#' converted to Cartesian distances. Angles use degrees except explicitly
#' internal sequence-reference angles. Text sizes use mm. Metadata retains
#' `coordinate_space`, `units`, `coordinate`, `rotation`, `rotation_applied`,
#' `coord_transform_applied`, `ratio`, `fit`, `xlim`, `ylim` and `clip`.
#' Exported coordinates already include layout rotation, before the ggplot2
#' coord transform. Native data units are not millimetres or genomic bp.
#'
#' @section Coordinate-specific tables:
#' Circular `annotation_registry` retains its documented identity, owner,
#' component/source row, inner track/radius/width/boundaries, outer region,
#' sector, band, dominant band, slot, spill reason, crossing count, bounding box
#' and leader corridor fields. Radii, widths and boxes use layout data units;
#' tracks/slots are build-local resource indices. Missing values denote fields
#' that do not apply. These resources are distinct from named quantitative
#' track IDs and Position lanes.
#'
#' Collinear `sequence_layout` retains reversible window/source origins,
#' horizontal range, backbone y, display orientation, units per bp, x offset
#' and anchor. Segmented rows also retain row/window order, breaks and display
#' `gap_bp`; omitted gaps have no source mapping. `orientation_diagnostics`,
#' `connection_report` and `label_diagnostics` retain support, selection and
#' placed/hidden-label evidence, including layer/row provenance where applicable.
#'
#' `track_registry` retains track ID, accession, placement, source layers,
#' numerical `limits` in original value units, optional `unit`, domain status
#' and measurement status. Native offsets/widths use layout units; `_mm`
#' fields are physical mm contributions. `normal_min`/`normal_max` use layout
#' units. `track_diagnostics` records input rows, retention/reason, coordinate
#' clipping and the original bar baseline. `axis_selection` records selected
#' sequence axes. `label_fit_diagnostics` retains original text, requested,
#' minimum/final sizes in mm and fallback/hiding reasons. Empty diagnostic
#' tables need not have a fixed column prototype.
#'
#' @section Measurement and cache scope:
#' `render_context` and diagnostic `measurement` distinguish `estimated`,
#' `unavailable` and `measured` panel status. `panel_mm`, `canvas_mm` and,
#' when measured, `panel_bounds_mm` use physical mm. `units_per_mm` and
#' `axis_units_per_inch` convert layout units along x/y; `units_per_inch` is
#' a conservative scalar text metric. `geometry_metric_status` and `converged`
#' describe refinement accuracy separately from panel measurement. Estimated
#' contexts can have missing dimensions and convergence. Device/basis/font
#' conditions state the source and limits of the evidence. Refinement counts
#' and private solver/candidate/cache structures are implementation details.
#'
#' A plot composition owns a fresh layout cache; a build/getter may populate
#' it. Device/context changes invalidate estimates, even for ordinary labels.
#' `get_chord_layout(build = FALSE)` requires a current cache. Previews and
#' diagnostics render independent copies and preserve the source cache.
#' Reference-path reuse lasts only for one geometry calculation. Temporary
#' render contexts restore their predecessor on normal or exceptional exit.
#' Ordinary R assignment retains ggplot2 reference semantics for ggproto
#' objects; arbitrary direct mutation of shared layer environments is not
#' tracked as a new plot composition.
#'
#' @section Stability and diagnostic limits:
#' Existing documented members above remain supported. Additional columns may
#' be attached for input metadata or geometry-specific needs; consumers should
#' select named fields, not fixed column positions/counts. Numeric layouts can
#' change with size, device, font and legitimate collision handling. Resource
#' numbers and diagnostic reason details are not persistent biological IDs.
#' No candidate layout, private cache or iteration count is a stable result.
#'
#' Diagnostics retain `measurement`, `checks`, `text`, `guides`,
#' `colour_suggestions`, `label_diagnostics`, `annotation_registry`,
#' `track_registry`, `track_diagnostics` and `warnings`. Text/guide bounds and
#' line/text sizes use mm; colour distances use the named model in each row.
#' A measured panel does not certify collision-free geometry, font embedding,
#' fallback, arbitrary grobs or cross-application portability. Preview paths
#' live in a temporary directory and may be removed by cleanup; use ordinary
#' [ggplot2::ggsave()] to keep a final file.
#'
#' @md
#' @name ggchord_layout_contract
#' @family layout results
#' @seealso [ggchord_data_contract], [export_ggchord_layout()],
#'   [get_chord_layout()], [view_ggchord()], [diagnose_ggchord()]
#' @examples
#' p <- ggchord(data.frame(accver = "A", length = 1000)) + geom_seq()
#' invisible(ggplot2::ggplot_build(p))
#' geometry <- export_ggchord_layout(p, include = "seq")
#' unique(geometry$seq[c("layer_id", "component", "source_row")])
#' geometry$metadata$coordinate
NULL
