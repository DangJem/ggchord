Browse complete current parameters by purpose. Defaults come from actual formals and descriptions from matching Rd. NULL can resolve during layout. A choice vector may mean match.arg alternatives or multiple enabled options; read its description. Upstream reexports follow their own help.

Start with the [task index](task-index.html). Lessons compare important parameters; references define full signatures and forwarding. Expand functions as needed; upstream dots are not presented as separate ggchord APIs.


## Coordinates and global views {#coordinates}

<span id="api-coord_chord"></span>
<details><summary><code>coord_chord()</code></summary>


[Full reference](reference/coord_chord.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `rotation` | `0` | Numeric/integer; range in description | Global counterclockwise layout rotation in degrees, default 0, measured from the positive x axis. This rotates the complete layout; it does not change the biological strand or sequence orientation. |
| `ratio` | `1` | Numeric/integer; range in description | Fixed y/x aspect ratio, default 1. |
| `xlim` | `NULL` | Nullable R object; type/inheritance below | Optional user limits. Explicit limits take priority over automatically fitted limits. |
| `ylim` | `NULL` | Nullable R object; type/inheritance below | Optional user limits. Explicit limits take priority over automatically fitted limits. |
| `expand` | `FALSE` | Logical | Logical. Expand coordinate limits, default FALSE. Automatic fitting already includes a small safety margin; set TRUE to request the additional ggplot2 coordinate expansion. |
| `clip` | `"off"` | Character/enum; choices below | Whether drawing is clipped to the panel, default "off". |
| `fit` | `c("labels", "geometry", "manual")` | Vector; choices versus multiple selection below | Fitting strategy: "labels" includes measured label boxes, "geometry" fits geometric elements only, and "manual" requires explicit xlim and ylim. |

</details>

<span id="api-coord_circular"></span>
<details><summary><code>coord_circular()</code></summary>


[Full reference](reference/coord_circular.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `rotation` | `0` | Numeric/integer; range in description | Counterclockwise rotation in degrees from a 12-o'clock genomic origin. The default 0 keeps the origin at 12 o'clock; 90 places it at 9 o'clock, and -90 at 3 o'clock. |
| `direction` | `c("clockwise", "counterclockwise")` | Vector; choices versus multiple selection below | Direction of increasing genomic coordinates. |
| `gap` | `0` | Numeric/integer; range in description | Opening in degrees. Zero draws a closed circle. |
| `ratio` | `1` | Numeric/integer; range in description | Fixed y/x aspect ratio, default 1. |
| `xlim` | `NULL` | Nullable R object; type/inheritance below | Optional user limits. Explicit limits take priority over automatically fitted limits. |
| `ylim` | `NULL` | Nullable R object; type/inheritance below | Optional user limits. Explicit limits take priority over automatically fitted limits. |
| `expand` | `FALSE` | Logical | Logical. Expand coordinate limits, default FALSE. Automatic fitting already includes a small safety margin; set TRUE to request the additional ggplot2 coordinate expansion. |
| `clip` | `"off"` | Character/enum; choices below | Whether drawing is clipped to the panel, default "off". |
| `fit` | `c("labels", "geometry", "manual")` | Vector; choices versus multiple selection below | Fitting strategy: "labels" includes measured label boxes, "geometry" fits geometric elements only, and "manual" requires explicit xlim and ylim. |

</details>

<span id="api-coord_collinear"></span>
<details><summary><code>coord_collinear()</code></summary>


[Full reference](reference/coord_collinear.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `ratio` | `1` | Numeric/integer; range in description | Fixed y/x aspect ratio, default 1. |
| `xlim` | `NULL` | Nullable R object; type/inheritance below | Optional user limits. Explicit limits take priority over automatically fitted limits. |
| `ylim` | `NULL` | Nullable R object; type/inheritance below | Optional user limits. Explicit limits take priority over automatically fitted limits. |
| `expand` | `FALSE` | Logical | Logical. Expand coordinate limits, default FALSE. Automatic fitting already includes a small safety margin; set TRUE to request the additional ggplot2 coordinate expansion. |
| `clip` | `"off"` | Character/enum; choices below | Whether drawing is clipped to the panel, default "off". |
| `fit` | `c("labels", "geometry", "manual")` | Vector; choices versus multiple selection below | Fitting strategy: "labels" includes measured label boxes, "geometry" fits geometric elements only, and "manual" requires explicit xlim and ylim. |
| `align` | `c("left", "center", "right")` | Vector; choices versus multiple selection below | Display alignment: "left" (default), "center", or "right". Applied to the full sequence or the interval supplied by anchor. |
| `anchor` | `NULL` | Explicit window/anchor table or documented equivalent | Optional anchor positions as a numeric vector named by accver, or a data frame with accver and either position or start/end. Every displayed sequence must have exactly one anchor. Missing and duplicate anchors are errors; multiple-copy targets must be selected explicitly before constructing the coordinate system. |
| `orientation` | `c("manual", "auto")` | Vector; choices versus multiple selection below | "manual" (default) uses sequence directions supplied by geom_seq(). "auto" derives display directions from interval-link evidence and records support and unresolved cases. |
| `min_support` | `0.1` | Numeric/integer; range in description | Minimum directional imbalance in [0, 1] for an automatic change. It is the absolute difference between forward and reverse support divided by their sum, default 0.1. A tie is always unresolved. |
| `segments` | `NULL` | Explicit window/anchor table or documented equivalent | Optional data frame with unique window accver, row_id and orientation (1 or -1). Every displayed accession must be declared. Each window has its own local coordinates and reversible source metadata. |
| `segment_gap` | `100` | Numeric/integer; range in description | Positive display gap in native sequence units when segments is supplied, default 100. It never implies genomic continuity. |
| `seq_axis` | `c("horizontal", "vertical")` | Vector; choices versus multiple selection below | "horizontal" (default) or "vertical". Vertical mode projects sequence-local geometry onto upright columns and lays out text in its own readable orientation; it does not change source strands. |

</details>

<span id="api-coord_dotplot"></span>
<details><summary><code>coord_dotplot()</code></summary>


[Full reference](reference/coord_dotplot.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `x_seqs` | `NULL` | Nullable R object; type/inheritance below | Optional unique accession vectors in axis order. Defaults to the query/subject accessions occurring in alignment layers, or all sequences when there are no alignments. Rows outside a selected axis set are omitted and recorded; self and repeated HSP rows are retained. |
| `y_seqs` | `NULL` | Nullable R object; type/inheritance below | Optional unique accession vectors in axis order. Defaults to the query/subject accessions occurring in alignment layers, or all sequences when there are no alignments. Rows outside a selected axis set are omitted and recorded; self and repeated HSP rows are retained. |
| `x_orientation` | `1` | Numeric/integer; range in description | Display direction: 1 or -1, either a scalar or a vector named by every selected accession. Input is unchanged. |
| `y_orientation` | `1` | Numeric/integer; range in description | Display direction: 1 or -1, either a scalar or a vector named by every selected accession. Input is unchanged. |
| `gap` | `NULL` | Nullable R object; type/inheritance below | Nonnegative spacing in native sequence units between axis segments. NULL uses five percent of the longest sequence. |
| `ratio` | `1` | Numeric/integer; range in description | Ratio of physical y units to physical x units, default 1. |
| `xlim` | `NULL` | Nullable R object; type/inheritance below | Optional global concatenated position limits. These zoom the panel without changing exported HSP endpoints. |
| `ylim` | `NULL` | Nullable R object; type/inheritance below | Optional global concatenated position limits. These zoom the panel without changing exported HSP endpoints. |
| `expand` | `FALSE` | Logical | Logical. Expand panel limits, default FALSE. |
| `clip` | `"on"` | Character/enum; choices below | Panel clipping, "on" (default) or "off". |

</details>

<span id="api-coord_radial_synteny"></span>
<details><summary><code>coord_radial_synteny()</code></summary>


[Full reference](reference/coord_radial_synteny.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `start_angle` | `90` | Numeric/integer; range in description | Direction of the first outward spoke in degrees, counterclockwise from the positive x axis. Default 90 starts at noon. |
| `clockwise` | `TRUE` | Logical | Logical, order subsequent spokes clockwise, default TRUE. |
| `inner_radius` | `0.6` | Numeric/integer; range in description | Positive inner radius in layout units, default 0.6. The longest proportional spoke has length 2. The central hole is display spacing and does not remove any source positions. |
| `ratio` | `1` | Numeric/integer; range in description | Fixed y/x aspect ratio, default 1. |
| `xlim` | `NULL` | Nullable R object; type/inheritance below | Optional user limits. Explicit limits take priority over automatically fitted limits. |
| `ylim` | `NULL` | Nullable R object; type/inheritance below | Optional user limits. Explicit limits take priority over automatically fitted limits. |
| `expand` | `FALSE` | Logical | Logical. Expand coordinate limits, default FALSE. Automatic fitting already includes a small safety margin; set TRUE to request the additional ggplot2 coordinate expansion. |
| `clip` | `"off"` | Character/enum; choices below | Whether drawing is clipped to the panel, default "off". |
| `fit` | `c("labels", "geometry", "manual")` | Vector; choices versus multiple selection below | Fitting strategy: "labels" includes measured label boxes, "geometry" fits geometric elements only, and "manual" requires explicit xlim and ylim. |

</details>

<span id="api-coord_sankey"></span>
<details><summary><code>coord_sankey()</code></summary>


[Full reference](reference/coord_sankey.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `column_gap` | `3` | Numeric/integer; range in description | Positive separation of the two backbones in layout units. |
| `self_links` | `FALSE` | Logical | Show relationships whose query and subject entity are equal. |
| `ratio` | `1` | Numeric/integer; range in description | Fixed y/x aspect ratio, default 1. |
| `xlim` | `NULL` | Nullable R object; type/inheritance below | Optional user limits. Explicit limits take priority over automatically fitted limits. |
| `ylim` | `NULL` | Nullable R object; type/inheritance below | Optional user limits. Explicit limits take priority over automatically fitted limits. |
| `expand` | `FALSE` | Logical | Logical. Expand coordinate limits, default FALSE. Automatic fitting already includes a small safety margin; set TRUE to request the additional ggplot2 coordinate expansion. |
| `clip` | `"off"` | Character/enum; choices below | Whether drawing is clipped to the panel, default "off". |
| `fit` | `c("labels", "geometry", "manual")` | Vector; choices versus multiple selection below | Fitting strategy: "labels" includes measured label boxes, "geometry" fits geometric elements only, and "manual" requires explicit xlim and ylim. |

</details>


## Geometry and labels {#layers}

<span id="api-geom_alignment_segment"></span>
<details><summary><code>geom_alignment_segment()</code></summary>


[Full reference](reference/geom_alignment_segment.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard ggplot2 mapping and layer data. Data defaults to the constructor's ribbon_data. Geometry roles are qaccver, saccver, qstart, qend, sstart, and send. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard ggplot2 mapping and layer data. Data defaults to the constructor's ribbon_data. Geometry roles are qaccver, saccver, qstart, qend, sstart, and send. |
| `colour` | `"#59636D"` | Character/enum; choices below | Fixed styles. Mapped styles use link_colour, link_alpha, link_linewidth and link_linetype, with the corresponding scale_link functions. Explicit fixed styles win. |
| `alpha` | `0.6` | Numeric/integer; range in description | Fixed styles. Mapped styles use link_colour, link_alpha, link_linewidth and link_linetype, with the corresponding scale_link functions. Explicit fixed styles win. |
| `linewidth` | `0.4` | Numeric/integer; range in description | Fixed styles. Mapped styles use link_colour, link_alpha, link_linewidth and link_linetype, with the corresponding scale_link functions. Explicit fixed styles win. |
| `linetype` | `1` | Numeric/integer; range in description | Fixed styles. Mapped styles use link_colour, link_alpha, link_linewidth and link_linetype, with the corresponding scale_link functions. Explicit fixed styles win. |
| `position` | `"identity"` | Character/enum; choices below | Only "identity" is supported: data adjustment would change the alignment endpoints. |
| `show.legend` | `NA` | R object; constraints below | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `allow_repeated` | `FALSE` | Logical | Allow explicitly identified source_hsp_id records to be displayed more than once. All participating layers must opt in. Rows without declared IDs are never deduplicated by their values. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional arguments passed to ggplot2's segment geometry. |

</details>

<span id="api-geom_feature"></span>
<details><summary><code>geom_feature()</code></summary>


[Full reference](reference/geom_feature.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Optional aesthetic mapping. Role aesthetics such as accver, start, end and strand may rename input columns; ordinary visual mappings are evaluated after geometry is generated. |
| `data` | `NULL` | Table/supported object; data may be a function | data.frame with accver, start, end and strand; optional type, category and label. Factor start/end columns use their numeric labels, not their internal level codes. Coordinates must resolve to finite numbers. display_priority or logical prioritized_display. As an alternative to strand, directionality or direction may contain "forward", "reverse", "bidirectional", or "none". When omitted, uses the full constructor feature_data (or compatible gene_data). Other layers never remove records from this default input. |
| `shape` | `"arrow"` | Character/enum; choices below | Fixed feature geometry. An explicitly supplied shape takes precedence over aes(feature_shape = ...): "arrow", "compact_arrow", "promoter_arrow", "primer_arrow", "marker", "block", "capped_line", "primer_arc", "chevron", or "lollipop". The default is "arrow". Use aes(feature_shape = type) together with scale_feature_shape_manual() to map categories to geometry. |
| `width` | `NULL` | Geometry number/vector/documented formats | Optional numeric/vector/list controlling width in the shared feature geometry engine. In a circular feature-label composition, the base width responds to final label size within bounded limits. |
| `arrow_head_length` | `0.04` | Numeric/integer; range in description | Arrow-head dimensions in the sequence-local frame, shared with geom_gene(). Values above one for arrow_head_width create a protruding block-arrow head. |
| `arrow_head_width` | `1` | Numeric/integer; range in description | Arrow-head dimensions in the sequence-local frame, shared with geom_gene(). Values above one for arrow_head_width create a protruding block-arrow head. |
| `arrow_head_style` | `c("shouldered", "flush", "triangle")` | Vector; choices versus multiple selection below | Arrowhead construction: "shouldered" keeps the requested protruding head, "flush" keeps it level with the body, and "triangle" draws a full triangular wedge. |
| `short_feature` | `c("auto", "wedge", "block")` | Vector; choices versus multiple selection below | Fallback for arrows too short to hold their requested head: an automatically compressed directional glyph, a wedge, or a block. |
| `segment_boundaries` | `TRUE` | Logical | Draw explicit internal marks carried by a segmented feature. Contiguous source-segment joins use dotted dividers; a non-solid segment-level line_style may override them, while a cleavage_arrows coordinate adds or merges a biological cut mark. These marks never create additional arrowheads. |
| `boundary_colour` | `"#666666"` | Character/enum; choices below | Appearance of internal segment joins. |
| `boundary_linewidth` | `0.25` | Numeric/integer; range in description | Appearance of internal segment joins. |
| `boundary_linetype` | `"22"` | Character/enum; choices below | Appearance of internal segment joins. |
| `position` | `"identity"` | Character/enum; choices below | Feature placement. Use "identity", "strand", "plasmid", position_strand(), position_plasmid(), or position_feature_stack(). |
| `show.legend` | `TRUE` | Logical | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional fixed polygon aesthetics. |

</details>

<span id="api-geom_feature_label"></span>
<details><summary><code>geom_feature_label()</code></summary>


[Full reference](reference/geom_feature_label.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard layer inputs. Map feature text with label. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard layer inputs. Map feature text with label. |
| `label_orientation` | `c("feature", "horizontal", "radial", "tangent")` | Vector; choices versus multiple selection below | One of "feature", "horizontal", "radial", or "tangent". "feature" keeps both internal and adjacent labels tangent to the feature direction; compact labels move beside the feature. |
| `label_side` | `c("outside", "inside", "auto")` | Vector; choices versus multiple selection below | One of "inside", "outside", or "auto". |
| `label_overlap` | `c("hide", "nudge", "allow")` | Vector; choices versus multiple selection below | Fixed-label collision policy. |
| `label_wrap` | `NULL` | Nullable R object; type/inheritance below | Optional character wrapping width. |
| `label_rotation` | `NULL` | Geometry number/vector/documented formats | Optional numeric/vector/list. Label rotation angle, default 0 |
| `label_radial_offset` | `0.04` | Numeric/integer; range in description | Optional numeric/vector/list. Radial offset of labels, default 0.04 |
| `label_circum_offset` | `NULL` | Geometry number/vector/documented formats | Optional numeric/vector/list. Circumferential offset of labels, default 0 |
| `label_circum_limit` | `NULL` | Nullable R object; type/inheritance below | Optional logical/vector/list. Whether to limit circumferential offset, default TRUE |
| `position` | `"identity"` | Character/enum; choices below | Feature Position used for the label anchor. In coord_collinear(), an omitted Position shares the matching feature layer's lanes; an explicit Position changes only label anchors, leaving feature polygons in their own Position. If both inputs contain unique feature_id values, label inputs may be subsets or reordered; their accession, interval and strand must still agree. Without IDs, inheritance requires the same interval rows in order. Multiple matching geometry layers require an explicit label Position or unambiguous feature IDs. |
| `show.legend` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional text parameters. |

</details>

<span id="api-geom_feature_label_repel"></span>
<details><summary><code>geom_feature_label_repel()</code></summary>


[Full reference](reference/geom_feature_label_repel.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard layer inputs. Map feature text with label. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard layer inputs. Map feature text with label. |
| `label_layout` | `c("feature", "radial", "auto", "callout")` | Vector; choices versus multiple selection below | One of "feature", "radial", "auto", or "callout". The default "feature" follows plasmid-map convention: text is tried inside first, then assigned to a fixed-radius internal track. All non-external text is laid out character by character along its selected track radius rather than rotated as one straight string. A light leader appears when one or more intervening tracks separate the feature and its label. Choose another mode for unrestricted external callouts. |
| `label_wrap` | `NULL` | Nullable R object; type/inheritance below | Text fitting controls. label_fit = "shrink" explicitly uses internal-first placement without truncating the text. Text size is bounded by label_min_size; unresolved labels keep their requested font and use adjacent/external placement or hiding. Other fitting modes retain their existing behavior. |
| `label_fit` | `c("wrap", "none", "ellipsis", "auto", "shrink")` | Vector; choices versus multiple selection below | Text fitting controls. label_fit = "shrink" explicitly uses internal-first placement without truncating the text. Text size is bounded by label_min_size; unresolved labels keep their requested font and use adjacent/external placement or hiding. Other fitting modes retain their existing behavior. |
| `label_max_lines` | `2L` | Numeric/integer; range in description | Text fitting controls. label_fit = "shrink" explicitly uses internal-first placement without truncating the text. Text size is bounded by label_min_size; unresolved labels keep their requested font and use adjacent/external placement or hiding. Other fitting modes retain their existing behavior. |
| `label_side` | `c("outside", "inside", "auto")` | Vector; choices versus multiple selection below | One of "inside", "outside", or "auto". |
| `max_overlaps` | `Inf` | R object; constraints below | Maximum unresolved overlaps. |
| `external` | `TRUE` | Logical | Whether the staged feature layout may place a label outside the circular backbone after inside and adjacent placement fail. When FALSE, unresolved labels are hidden instead of pushing inward forever. External labels are rendered as rounded callouts using a lightened form of the feature's resolved fill. When restriction-site labels coexist on a circular map, both retain their own layout semantics but register final boxes in one coordinate-owned exterior track region so they cannot occupy the same space. |
| `position` | `"identity"` | Character/enum; choices below | Feature Position used for the label anchor. In coord_collinear(), an omitted Position shares the matching feature layer's lanes; an explicit Position changes only label anchors, leaving feature polygons in their own Position. If both inputs contain unique feature_id values, label inputs may be subsets or reordered; their accession, interval and strand must still agree. Without IDs, inheritance requires the same interval rows in order. Multiple matching geometry layers require an explicit label Position or unambiguous feature IDs. |
| `show.legend` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `label_min_size` | `1.5` | Numeric/integer; range in description | Minimum text size in millimetres for explicit shrink, default 1.5. Must not exceed the requested text size. Original text, requested/final sizes and fitting reason are exported with the label. |
| `label_segment_overlap` | `"fade"` | Character/enum; choices below | Character, default "fade". How portions of a leader line covered by another label are drawn: "fade" keeps them at reduced opacity, "clip" removes them, and "show" leaves the complete leader visible. This does not affect text wrapping, ellipsis or label hiding. |
| `label_segment_overlap_alpha` | `0.18` | Numeric/integer; range in description | Numeric in [0, 1], default 0.18. Relative opacity applied only to covered leader-line portions when label_segment_overlap = "fade". |
| `label_segment_linetype` | `"auto"` | Character/enum; choices below | Character or numeric, default "auto". Leader-line linetype. "auto" draws solid lines, except for labels that were moved to the other side of their arc, which are drawn dashed. Any other valid ggplot2 linetype (e.g. "solid", "dashed", "dotted", or a numeric dash pattern) is used for all leader lines. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional text parameters. |

</details>

<span id="api-geom_feature_plasmid"></span>
<details><summary><code>geom_feature_plasmid()</code></summary>


[Full reference](reference/geom_feature_plasmid.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Optional aesthetic mapping. Role aesthetics such as accver, start, end and strand may rename input columns; ordinary visual mappings are evaluated after geometry is generated. |
| `data` | `NULL` | Table/supported object; data may be a function | data.frame with accver, start, end and strand; optional type, category and label. Factor start/end columns use their numeric labels, not their internal level codes. Coordinates must resolve to finite numbers. display_priority or logical prioritized_display. As an alternative to strand, directionality or direction may contain "forward", "reverse", "bidirectional", or "none". When omitted, uses the full constructor feature_data (or compatible gene_data). Other layers never remove records from this default input. |
| `shape` | `"arrow"` | Character/enum; choices below | Fixed feature geometry. An explicitly supplied shape takes precedence over aes(feature_shape = ...): "arrow", "compact_arrow", "promoter_arrow", "primer_arrow", "marker", "block", "capped_line", "primer_arc", "chevron", or "lollipop". The default is "arrow". Use aes(feature_shape = type) together with scale_feature_shape_manual() to map categories to geometry. |
| `width` | `0.058` | Numeric/integer; range in description | Track width; the default matches the compact feature bands used by desktop plasmid editors. |
| `arrow_head_length` | `0.068` | Numeric/integer; range in description | Arrow-head dimensions in the sequence-local frame, shared with geom_gene(). Values above one for arrow_head_width create a protruding block-arrow head. |
| `arrow_head_width` | `1.42` | Numeric/integer; range in description | Arrow-head dimensions in the sequence-local frame, shared with geom_gene(). Values above one for arrow_head_width create a protruding block-arrow head. |
| `arrow_head_style` | `"shouldered"` | Character/enum; choices below | Arrowhead construction: "shouldered" keeps the requested protruding head, "flush" keeps it level with the body, and "triangle" draws a full triangular wedge. |
| `short_feature` | `"auto"` | Character/enum; choices below | Fallback for arrows too short to hold their requested head: an automatically compressed directional glyph, a wedge, or a block. |
| `segment_boundaries` | `TRUE` | Logical | Draw explicit internal marks carried by a segmented feature. Contiguous source-segment joins use dotted dividers; a non-solid segment-level line_style may override them, while a cleavage_arrows coordinate adds or merges a biological cut mark. These marks never create additional arrowheads. |
| `boundary_colour` | `"#666666"` | Character/enum; choices below | Appearance of internal segment joins. |
| `boundary_linewidth` | `0.25` | Numeric/integer; range in description | Appearance of internal segment joins. |
| `boundary_linetype` | `"22"` | Character/enum; choices below | Appearance of internal segment joins. |
| `position` | `NULL` | Position object or documented string | Feature position. NULL (the default) automatically stacks overlapping features on compact lanes inside the plasmid backbone. Reuse the same explicit position_feature_stack() object for a separate label layer so polygons and labels share their lanes. Optional numeric display_priority (larger first) or logical prioritized_display can override the default longest-first ordering for overlapping features. |
| `show.legend` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional fixed polygon aesthetics. |

</details>

<span id="api-geom_gene"></span>
<details><summary><code>geom_gene()</code></summary>


[Full reference](reference/geom_gene.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Default NULL (uses pre-computed data) |
| `data` | `NULL` | Table/supported object; data may be a function | Optional layer table or function, overriding constructor data. With unified feature_data, the default selects explicitly classified gene/CDS records; see ggchord(). Legacy gene_data retains all rows. |
| `width` | `NULL` | Geometry number/vector/documented formats | Optional numeric/vector/list. Width of gene arrows, default 0.05 |
| `arrow_head_length` | `0.04` | Numeric/integer; range in description | Shared arrow-head dimensions in the sequence-local frame. |
| `arrow_head_width` | `1` | Numeric/integer; range in description | Shared arrow-head dimensions in the sequence-local frame. |
| `arrow_head_style` | `c("shouldered", "flush", "triangle")` | Vector; choices versus multiple selection below | Arrowhead construction: shouldered, flush, or a full triangular wedge. |
| `short_feature` | `c("auto", "wedge", "block")` | Vector; choices versus multiple selection below | Fallback for arrows too short to hold their requested head: automatic wedge/block selection, a wedge, or a block. |
| `position` | `"identity"` | Character/enum; choices below | Feature placement. Use "identity", "strand", "plasmid", position_strand(), position_plasmid(), or position_feature_stack(). |
| `show.legend` | `TRUE` | Logical | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional arguments passed to geom_polygon() |

</details>

<span id="api-geom_gene_label"></span>
<details><summary><code>geom_gene_label()</code></summary>


[Full reference](reference/geom_gene_label.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Default NULL (uses pre-computed data) |
| `data` | `NULL` | Table/supported object; data may be a function | Optional layer table or function. Default data follows the gene/CDS selection described in ggchord(). |
| `label_orientation` | `"horizontal"` | Character/enum; choices below | Character, default "horizontal". Text orientation relative to each sequence path: "radial", "tangent", or "horizontal". |
| `label_side` | `"outside"` | Character/enum; choices below | Character, default "outside". Place labels outside the chord, use the strand-based "auto" side, or force "inside" placement. |
| `label_overlap` | `"hide"` | Character/enum; choices below | Character, default "hide". Fixed-label collision policy: omit later conflicting labels ("hide"), apply the small position adjustment ("nudge"), or keep all requested positions ("allow"). Input row order therefore provides a simple way to prioritize labels in "hide" mode. |
| `label_rotation` | `NULL` | Geometry number/vector/documented formats | Optional numeric/vector/list. Label rotation angle, default 0 |
| `label_radial_offset` | `0.04` | Numeric/integer; range in description | Optional numeric/vector/list. Radial offset of labels, default 0.04 |
| `label_circum_offset` | `NULL` | Geometry number/vector/documented formats | Optional numeric/vector/list. Circumferential offset of labels, default 0 |
| `label_circum_limit` | `NULL` | Nullable R object; type/inheritance below | Optional logical/vector/list. Whether to limit circumferential offset, default TRUE |
| `label_wrap` | `NULL` | Nullable R object; type/inheritance below | Numeric or NULL, default NULL. When set, long gene annotations are wrapped at this many characters (e.g. 15), which makes the labels narrower and less prone to overlap. |
| `position` | `"identity"` | Character/enum; choices below | Position adjustment. Supply the same position_feature_stack() used by a gene/feature layer when labels are drawn from separate data or without a geometry layer. In coord_collinear(), labels without an explicit Position share a matching geometry layer's lanes; explicit label Positions do not move the geometry layer. Unique feature_id columns allow label subsets and reordered labels when their accession, interval and strand agree with the geometry input. Multiple matching geometry layers require an explicit label Position or unambiguous feature IDs. |
| `show.legend` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional arguments passed to geom_text() |

</details>

<span id="api-geom_gene_label_repel"></span>
<details><summary><code>geom_gene_label_repel()</code></summary>


[Full reference](reference/geom_gene_label_repel.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Default NULL (uses pre-computed data) |
| `data` | `NULL` | Table/supported object; data may be a function | Optional layer table or function. Default data follows the gene/CDS selection described in ggchord(). |
| `label_layout` | `"radial"` | Character/enum; choices below | Character, default "radial". Label layout: "radial" arranges horizontal external labels along sequence-relative contours, with straight normal leaders or a short normal departure followed by a long connector when displaced; "auto" uses shared vertical columns on the left/right and radial placement elsewhere; "arc" rotates text along the sequence tangent. The former "aligned" mode has been removed; use "auto" for side columns. |
| `label_wrap` | `NULL` | Nullable R object; type/inheritance below | Numeric or NULL, default NULL. When set, long gene annotations are wrapped at this many characters (e.g. 15), overriding label_fit. |
| `label_fit` | `"wrap"` | Character/enum; choices below | Character, default "wrap". How labels whose initial boxes are crowded or unusually wide are fitted before deterministic layout: "wrap" preserves the complete annotation over multiple lines, "ellipsis" uses a single shortened line, "auto" wraps and adds an ellipsis only if text exceeds label_max_lines, and "none" disables adaptive fitting. Explicit "shrink" uses internal-first feature placement and reduces the complete text only as far as label_min_size. Labels that still cannot fit use the existing adjacent/external placement or hiding policy, at their requested size. This choice uses feature placement rather than the external layout selected by label_layout; other fitting modes are unchanged. |
| `label_max_lines` | `2L` | Numeric/integer; range in description | Positive integer, default 2. Maximum number of lines used by adaptive "wrap" or "auto" fitting. |
| `label_side` | `"outside"` | Character/enum; choices below | Character, default "outside". Which side of the arc the labels sit on. "auto" keeps the strand-based placement (same as before); "outside" moves labels that would be inside the chord (where they can overlap the ribbons) to the outside of their arc; "inside" does the opposite. Labels moved to the other side are connected with a dashed leader line (see label_segment_linetype). |
| `max_overlaps` | `Inf` | R object; constraints below | Numeric, default Inf. Hide labels that still overlap more than this many other labels after repulsion (ggrepel-style decluttering). Use a finite value to clean up crowded plots. |
| `label_segment_overlap` | `"fade"` | Character/enum; choices below | Character, default "fade". How portions of a leader line covered by another label are drawn: "fade" keeps them at reduced opacity, "clip" removes them, and "show" leaves the complete leader visible. This does not affect text wrapping, ellipsis or label hiding. |
| `label_segment_overlap_alpha` | `0.18` | Numeric/integer; range in description | Numeric in [0, 1], default 0.18. Relative opacity applied only to covered leader-line portions when label_segment_overlap = "fade". |
| `label_segment_linetype` | `"auto"` | Character/enum; choices below | Character or numeric, default "auto". Leader-line linetype. "auto" draws solid lines, except for labels that were moved to the other side of their arc, which are drawn dashed. Any other valid ggplot2 linetype (e.g. "solid", "dashed", "dotted", or a numeric dash pattern) is used for all leader lines. |
| `position` | `"identity"` | Character/enum; choices below | Position adjustment. Supply the same position_feature_stack() used by a gene/feature layer when labels are drawn from separate data or without a geometry layer. |
| `show.legend` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `label_min_size` | `1.5` | Numeric/integer; range in description | Minimum text size in millimetres for explicit label_fit = "shrink", default 1.5. It must not exceed the requested text size. Source text, requested/final size and fitting reason are retained in exported label geometry. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional text styles and drawing options. Automatic layouts own rotation, offsets, orientation and overlap policy; use geom_gene_label() for manual label_rotation, label_radial_offset, label_circum_offset, label_circum_limit, label_orientation or label_overlap. These manual controls are rejected here. |

</details>

<span id="api-geom_link_line"></span>
<details><summary><code>geom_link_line()</code></summary>


[Full reference](reference/geom_link_line.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard ggplot2 mapping and layer data, including data functions. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard ggplot2 mapping and layer data, including data functions. |
| `link_type` | `c("curve", "straight")` | Vector; choices versus multiple selection below | Either "curve" (default) or "straight". |
| `link_gap` | `NULL` | Nullable R object; type/inheritance below | Explicit endpoint spacing. With the default link_avoid = "none", NULL uses a fixed spacing of 0.035. Set link_avoid to "smooth" or "uniform" to opt in to entity-aware endpoint spacing. |
| `link_ctrl_point` | `NULL` | Nullable R object; type/inheritance below | Bezier control points: two numbers, four numbers, or a list per input row. |
| `link_avoid` | `c("none", "smooth", "uniform")` | Vector; choices versus multiple selection below | One of "none" (default), "smooth", or "uniform". Explicit gap takes precedence. |
| `link_branch` | `c("none", "query", "subject")` | Vector; choices versus multiple selection below | Share equal query or subject endpoints: "none", "query", "subject". |
| `link_branch_fraction` | `0.15` | Numeric/integer; range in description | Shared length relative to the shortest member, in (0, 0.5]. |
| `arrow` | `NULL` | Nullable R object; type/inheritance below | A grid arrow specification, or NULL. |
| `colour` | `"#59636D"` | Character/enum; choices below | Fixed styles; mappings use link_* aesthetics. An explicitly supplied fixed value takes precedence over the matching mapping. Omit the fixed argument to retain a mapped style, including an inherited mapping when inherit.aes = TRUE. |
| `alpha` | `0.6` | Numeric/integer; range in description | Fixed styles; mappings use link_* aesthetics. An explicitly supplied fixed value takes precedence over the matching mapping. Omit the fixed argument to retain a mapped style, including an inherited mapping when inherit.aes = TRUE. |
| `linewidth` | `0.4` | Numeric/integer; range in description | Fixed styles; mappings use link_* aesthetics. An explicitly supplied fixed value takes precedence over the matching mapping. Omit the fixed argument to retain a mapped style, including an inherited mapping when inherit.aes = TRUE. |
| `linetype` | `1` | Numeric/integer; range in description | Fixed styles; mappings use link_* aesthetics. An explicitly supplied fixed value takes precedence over the matching mapping. Omit the fixed argument to retain a mapped style, including an inherited mapping when inherit.aes = TRUE. |
| `position` | `"identity"` | Character/enum; choices below | Standard ggplot2 layer arguments. |
| `show.legend` | `NA` | R object; constraints below | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `between` | `NULL` | Nullable R object; type/inheritance below | Two display-instance names or declared-order indices. |
| `allow_repeated` | `FALSE` | Logical | Allow repeated display of an explicitly declared source_hsp_id. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional path parameters. |

</details>

<span id="api-geom_link_ribbon"></span>
<details><summary><code>geom_link_ribbon()</code></summary>


[Full reference](reference/geom_link_ribbon.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Default NULL (uses pre-computed data) |
| `data` | `NULL` | Table/supported object; data may be a function | Default NULL (retrieved automatically from the layout) |
| `link_ctrl_point` | `NULL` | Nullable R object; type/inheritance below | Optional vector/list of Bezier edge controls. Two numbers specify one shared control; four numbers specify controls for the two ribbon edges. A per-row list may contain pairs of coordinate vectors. NULL uses the origin c(0, 0) in chord layouts and derives controls from the sequence endpoints in collinear layouts. With shared branches, controls shape the far ends while the junction follows the shared trunk. |
| `link_gap` | `NULL` | Nullable R object; type/inheritance below | Optional numeric/vector controlling spacing between sequences and ribbon endpoints. With the default link_avoid = "uniform", NULL uses one obstacle-aware clearance across each endpoint. Use "smooth" for clearance that changes only next to geom_gene() and geom_feature() polygons, or "none" for fixed spacing 0.035. Text and leader lines are ignored. Supplying a number disables the automatic rule and uses that exact spacing. Linear layouts reduce automatic baseline clearance to fit narrow channels; an explicit gap or feature clearance that cannot fit raises an error. |
| `link_avoid` | `c("uniform", "none", "smooth")` | Vector; choices versus multiple selection below | One of "uniform" (default), "none", or "smooth". Smooth avoidance changes clearance locally; uniform avoidance uses one clearance across an endpoint. |
| `link_branch` | `c("none", "query", "subject")` | Vector; choices versus multiple selection below | Share equal query or subject endpoints: "none", "query", "subject". |
| `link_branch_fraction` | `0.15` | Numeric/integer; range in description | Shared length relative to the shortest member, in (0, 0.5]. |
| `fill` | `NULL` | Nullable R object; type/inheritance below | Optional fixed ribbon fill. NULL uses the default Identity fill. |
| `alpha` | `0.42` | Numeric/integer; range in description | Standard fixed ribbon styles. Explicit values take precedence over matching mappings; omitted defaults do not suppress user mappings. For example, omit alpha when mapping ribbon_alpha. |
| `colour` | `"#59636D"` | Character/enum; choices below | Standard fixed ribbon styles. Explicit values take precedence over matching mappings; omitted defaults do not suppress user mappings. For example, omit alpha when mapping ribbon_alpha. |
| `linewidth` | `0.08` | Numeric/integer; range in description | Standard fixed ribbon styles. Explicit values take precedence over matching mappings; omitted defaults do not suppress user mappings. For example, omit alpha when mapping ribbon_alpha. |
| `linetype` | `NULL` | Nullable R object; type/inheritance below | Fixed outline line type. The default NULL uses solid outlines for adjacent displayed sequences and dashed outlines for non-adjacent sequences. Collinear outline widths use nonadjacent_linewidth. Chord adjacency includes the first/last pair. This follows the final display order, independently of alignment direction. Explicit fixed styles or aes(ribbon_linetype = ...) (also aes(linetype = ...)) replace the automatic styles. Use a matching scale for data mappings; automatic styles do not train a linetype scale. Global line-type mappings are used when inherit.aes = TRUE. |
| `adjacent_only` | `NULL` | Nullable R object; type/inheritance below | Whether only neighbouring sequences are connected. NULL uses TRUE in coord_collinear() and FALSE in coord_chord(). Set FALSE to retain cross-track ribbons. Chord neighbours include the first/last sequence pair. |
| `pairs` | `NULL` | Nullable R object; type/inheritance below | Optional data frame with qaccver and saccver selecting unordered comparison pairs. With adjacent_only = NULL, explicit pairs take precedence over the coordinate default. Omitted rows are reported by export_ggchord_layout(). |
| `position` | `"identity"` | Character/enum; choices below | Standard ggplot2 layer arguments. |
| `show.legend` | `TRUE` | Logical | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `nonadjacent_linewidth` | `0.4` | Numeric/integer; range in description | Fixed outline width for non-adjacent displayed rows in coord_collinear(), default 0.4 mm. These ribbons use dashed outlines by default. Zero uses the layer's ordinary linewidth for every ribbon. This changes width independently of line type. An explicit aes(linewidth = ...) controls widths instead of this fixed default. |
| `between` | `NULL` | Nullable R object; type/inheritance below | Two display-instance names or declared-order indices for a unique route. |
| `allow_repeated` | `FALSE` | Logical | Explicitly allow a declared source_hsp_id in more than one displayed relation. |
| `tip` | `c("none", "subject")` | Vector; choices versus multiple selection below | Full interval fronts (none, default), or an optional subject taper. The shoulder retains the translated target endpoint span at 88 percent of the channel. The tip is a decorative front-midpoint marker, not a single matched base or biological transfer direction. Degenerate fronts fall back to the full front with a warning. Source intervals are unchanged. |
| `gradient` | `FALSE` | Logical | Logical. Opt into endpoint sequence-colour interpolation; default FALSE. Explicit fill, mapping and scales still take priority. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional arguments passed to geom_polygon() |

</details>

<span id="api-geom_primer"></span>
<details><summary><code>geom_primer()</code></summary>


[Full reference](reference/geom_primer.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard layer inputs. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard layer inputs. |
| `side` | `c("between", "outside", "inside")` | Vector; choices versus multiple selection below | Backbone lane used by the primer arc. "between" uses the backbone centre line; "outside" and "inside" use arc_offset. |
| `arc_offset` | `0.045` | Numeric/integer; range in description | Distance from the backbone centre when side is "outside" or "inside". |
| `arc_linewidth` | `0.026` | Numeric/integer; range in description | Radial width of the primer arc. |
| `colour` | `"#A020F0"` | Character/enum; choices below | Primer arc colour. |
| `position` | `NULL` | Position object or documented string | Standard layer arguments. |
| `show.legend` | `FALSE` | Logical | Standard layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional parameters passed to geom_feature(). |

</details>

<span id="api-geom_primer_label_repel"></span>
<details><summary><code>geom_primer_label_repel()</code></summary>


[Full reference](reference/geom_primer_label_repel.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard layer inputs. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard layer inputs. |
| `side` | `c("outside", "inside")` | Vector; choices versus multiple selection below | Side of the backbone on which the primer callout is arranged. "inside" retains the inner label layout, including when exterior restriction or feature labels share the map. The default is "outside". |
| `show_location` | `TRUE` | Logical | Include the binding range in the label. |
| `colour` | `"#A020F0"` | Character/enum; choices below | Primer arc colour. |
| `max_overlaps` | `Inf` | R object; constraints below | Maximum unresolved overlaps. |
| `position` | `NULL` | Position object or documented string | Standard layer arguments. |
| `show.legend` | `FALSE` | Logical | Standard layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional parameters passed to geom_feature(). |

</details>

<span id="api-geom_residue"></span>
<details><summary><code>geom_residue()</code></summary>


[Full reference](reference/geom_residue.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard ggplot2 layer mapping and data. Required aesthetics are x, y and label. Ordinary colour, fill, alpha, size, family and fontface mappings/scales are supported. size is in mm. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard ggplot2 layer mapping and data. Required aesthetics are x, y and label. Ordinary colour, fill, alpha, size, family and fontface mappings/scales are supported. size is in mm. |
| `overflow` | `c("tile", "error")` | Vector; choices versus multiple selection below | "tile" (default) or "error" when letters do not fit. |
| `position` | `"identity"` | Character/enum; choices below | Standard ggplot2 layer arguments. |
| `show.legend` | `NA` | R object; constraints below | Standard ggplot2 layer arguments. |
| `inherit.aes` | `TRUE` | Logical | Standard ggplot2 layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Fixed text/cell aesthetics and na.rm, passed to the layer. |

</details>

<span id="api-geom_restriction_site"></span>
<details><summary><code>geom_restriction_site()</code></summary>


[Full reference](reference/geom_restriction_site.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard layer inputs. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard layer inputs. |
| `label` | `TRUE` | Logical | Draw labels. |
| `label_style` | `c("enzyme_position", "enzyme", "position")` | Vector; choices versus multiple selection below | Label enzyme and position, enzyme only, or position only. |
| `label_order` | `c("auto", "enzyme_position", "position_enzyme")` | Vector; choices versus multiple selection below | Automatic mirrored label order, enzyme first, or position first. |
| `label_side` | `c("outside", "inside", "auto")` | Vector; choices versus multiple selection below | Draw labels outside, inside, or automatically. |
| `leader` | `c("radial", "elbow", "straight")` | Vector; choices versus multiple selection below | Leader routing style. The default "radial" measures actual text boxes, separates site anchors from final radial label positions, and routes either a direct connector or an independent radial stub plus fan segment. All enzyme-facing text edges follow one close circular contour in genomic order. Dense lateral clusters use evenly spaced contour fans, while sparse clusters retain their natural angular positions. The former "trunk" alias was removed in v0.16.1; use "radial". "straight" connects the tick tip directly to the final label edge, without detours around the sequence or other leaders. "elbow" uses one bend, starting with a local-normal departure where possible. Label positions and grouping are shared by all routing styles. |
| `min_label_gap` | `NULL` | Nullable R object; type/inheritance below | Minimum genomic fraction between label slots. NULL derives a compact device-aware default from the rendered labels. |
| `tick_length` | `0.04` | Numeric/integer; range in description | Local-normal distances. |
| `label_offset` | `0.18` | Numeric/integer; range in description | Local-normal distances. |
| `colour` | `"#202020"` | Character/enum; choices below | Fixed appearance. With the default fontface = NULL, enzyme names are bold only when that enzyme cuts once in the supplied sequence; coordinates remain regular. An explicit fontface disables this automatic distinction. |
| `linewidth` | `0.28` | Numeric/integer; range in description | Fixed appearance. With the default fontface = NULL, enzyme names are bold only when that enzyme cuts once in the supplied sequence; coordinates remain regular. An explicit fontface disables this automatic distinction. |
| `label_size` | `3.2` | Numeric/integer; range in description | Fixed appearance. With the default fontface = NULL, enzyme names are bold only when that enzyme cuts once in the supplied sequence; coordinates remain regular. An explicit fontface disables this automatic distinction. |
| `fontface` | `NULL` | Nullable R object; type/inheritance below | Fixed appearance. With the default fontface = NULL, enzyme names are bold only when that enzyme cuts once in the supplied sequence; coordinates remain regular. An explicit fontface disables this automatic distinction. |
| `family` | `NULL` | Nullable R object; type/inheritance below | Fixed appearance. With the default fontface = NULL, enzyme names are bold only when that enzyme cuts once in the supplied sequence; coordinates remain regular. An explicit fontface disables this automatic distinction. |
| `position` | `"identity"` | Character/enum; choices below | Standard layer arguments. |
| `show.legend` | `FALSE` | Logical | Standard layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional geom parameters. |

</details>

<span id="api-geom_seq"></span>
<details><summary><code>geom_seq()</code></summary>


[Full reference](reference/geom_seq.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Default NULL (uses pre-computed data) |
| `data` | `NULL` | Table/supported object; data may be a function | Default NULL (retrieved automatically from the layout) |
| `seq_order` | `NULL` | Geometry number/vector/documented formats | Source IDs in display order. Repeats create display instances; use names such as c(A_high='A', B='B', A_low='A') for stable instance IDs. Repeated-entity HSPs require a unique between route on each link layer. |
| `seq_orientation` | `NULL` | Geometry number/vector/documented formats | Optional numeric (1 or -1). Sequence orientation, default 1 |
| `seq_gap` | `NULL` | Geometry number/vector/documented formats | Optional numeric. Gap proportion between sequences, default 0.03 in chord coordinates. A named vector can set the gap after each sequence; omitted IDs use the coordinate default. In coord_collinear(), a positive distance between consecutive backbones as a proportion of the longest backbone's width, default 0.3. The final sequence's gap is unused in linear layouts. Segmented linear displays measure it against the longest display row; windows in one row must request the same gap, and only gaps between rows affect vertical placement. |
| `seq_radius` | `NULL` | Geometry number/vector/documented formats | Optional numeric (> 0). Sequence arc radius, default 1.0 |
| `seq_curvature` | `NULL` | Geometry number/vector/documented formats | Optional numeric. Signed arc bow (0=straight, 1=standard arc, negative=opposite bow), default 1.0. Finite positive and negative values are accepted without clipping; magnitudes above 1 amplify the bow. |
| `seq_style` | `c("auto", "single", "double", "band")` | Vector; choices versus multiple selection below | Backbone geometry. "auto" uses a double backbone in coord_circular() and a single backbone in coord_chord(). |
| `seq_backbone_gap` | `0.025` | Numeric/integer; range in description | Gap between double backbone lines. |
| `seq_backbone_width` | `0.035` | Numeric/integer; range in description | Width of a band backbone. |
| `linewidth` | `0.9` | Numeric/integer; range in description | Arc line width, default 0.9 |
| `position` | `"identity"` | Character/enum; choices below | Standard ggplot2 layer arguments. |
| `show.legend` | `TRUE` | Logical | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `length_mode` | `c("proportional", "equal")` | Vector; choices versus multiple selection below | "proportional" (default) retains source-length ratios; "equal" gives each displayed sequence the same angular allocation or linear span without modifying source lengths. Shared absolute scale bars are unavailable in multi-sequence equal mode; local axes retain source units. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional fixed path aesthetics. colour (or color) sets seq_colour. Explicit fixed values take precedence over matching mappings; omit fixed arguments to retain mapped styles. |

</details>

<span id="api-geom_seq_break"></span>
<details><summary><code>geom_seq_break()</code></summary>


[Full reference](reference/geom_seq_break.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard layer mapping and sequence data. By default, use the plot's sequence table. Truncation is read from .truncated_left/.truncated_right, with source interval and length metadata as fallback when available. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard layer mapping and sequence data. By default, use the plot's sequence table. Truncation is read from .truncated_left/.truncated_right, with source interval and length metadata as fallback when available. |
| `break_width` | `0.02` | Numeric/integer; range in description | Width and height of the complete break mark in linear layout units. These are layout distances, not millimetres. |
| `break_height` | `0.04` | Numeric/integer; range in description | Width and height of the complete break mark in linear layout units. These are layout distances, not millimetres. |
| `linewidth` | `0.4` | Numeric/integer; range in description | Fixed line appearance. Ordinary visual aesthetics can also be mapped and controlled by their ggplot2 scales. |
| `colour` | `"#3D4348"` | Character/enum; choices below | Fixed line appearance. Ordinary visual aesthetics can also be mapped and controlled by their ggplot2 scales. |
| `position` | `"identity"` | Character/enum; choices below | Standard ggplot2 layer arguments. |
| `show.legend` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional parameters passed to ggplot2::geom_path(). |

</details>

<span id="api-geom_seq_center_label"></span>
<details><summary><code>geom_seq_center_label()</code></summary>


[Full reference](reference/geom_seq_center_label.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard layer inputs. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard layer inputs. |
| `name` | `TRUE` | Logical | Show the sequence name. A character value overrides the name. |
| `length` | `TRUE` | Logical | Show the formatted sequence length. |
| `separator` | `"\n"` | Character/enum; choices below | Separator between name and length. |
| `name_style` | `list()` | Named list | Named lists of fixed text-style overrides for the name and length grobs. Supported fields are colour, size, family, fontface, and lineheight. |
| `length_style` | `list()` | Named list | Named lists of fixed text-style overrides for the name and length grobs. Supported fields are colour, size, family, fontface, and lineheight. |
| `show.legend` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional arguments passed to geom_text(). |

</details>

<span id="api-geom_seq_label"></span>
<details><summary><code>geom_seq_label()</code></summary>


[Full reference](reference/geom_seq_label.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Default NULL (uses pre-computed data) |
| `data` | `NULL` | Table/supported object; data may be a function | Default NULL (retrieved automatically from the layout) |
| `seq_label_radius` | `1` | Numeric/integer; range in description | Optional numeric/vector. Radial position of the labels as a multiplier of the sequence arc radius: 1 is on the arc, > 1 places the label outside (away from the chord center) and < 1 places it inside, default 1 Explicit values are incompatible with coord_collinear(), which places sequence labels beside the left endpoint of each horizontal track. |
| `seq_label_rotation` | `NULL` | Nullable R object; type/inheritance below | Optional numeric/vector. Additional label rotation (degrees) on top of the arc-aligned orientation, default NULL (0). Ignored when seq_label_orientation = "horizontal". |
| `labels` | `NULL` | Nullable R object; type/inheritance below | Optional character vector. Override the label texts (defaults to the sequence labels from geom_seq() or the sequence IDs) |
| `seq_label_orientation` | `c("arc", "horizontal")` | Vector; choices versus multiple selection below | Character, default "arc". Label text orientation: "arc" rotates the text along the sequence arc (and keeps it readable), "horizontal" draws every label horizontally, extending away from the chord center. |
| `seq_label_hjust` | `NULL` | Nullable R object; type/inheritance below | Optional numeric/vector. Horizontal justification of the labels, default NULL. The default arc orientation uses -0.2 so the text sits just inside the sequence; with seq_label_orientation = "horizontal" the justification is chosen automatically so the text extends away from the chord center. |
| `seq_label_vjust` | `NULL` | Nullable R object; type/inheritance below | Optional numeric/vector. Vertical justification of the labels, default NULL (0.5). |
| `check_overlap` | `FALSE` | Logical | Logical, default FALSE. When TRUE, labels that would overlap a previously drawn label are skipped (ggplot2's geom_text() option). |
| `position` | `"identity"` | Character/enum; choices below | Standard ggplot2 layer arguments. |
| `show.legend` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional arguments passed to geom_text() |

</details>

<span id="api-geom_seq_region"></span>
<details><summary><code>geom_seq_region()</code></summary>


[Full reference](reference/geom_seq_region.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Default NULL (uses pre-computed data) |
| `data` | `NULL` | Table/supported object; data may be a function | data.frame with at least accver, start, end; optional label, category and color. Factor start/end columns use their numeric labels, not level codes. |
| `fill` | `"#F59E0B"` | Character/enum; choices below | Character. Default fill colour for regions, default "#F59E0B". |
| `colour` | `"#B45309"` | Character/enum; choices below | Outline colour, default "#B45309". |
| `alpha` | `0.25` | Numeric/integer; range in description | Numeric (0-1). Region alpha, default 0.25. |
| `region_width` | `0.08` | Numeric/integer; range in description | Numeric. Band width in chord radius units, default 0.08. |
| `region_offset` | `0` | Numeric/integer; range in description | Numeric. Radial offset from the sequence arc, default 0. |
| `region_side` | `c("inside", "outside", "auto")` | Vector; choices versus multiple selection below | Character. "inside", "outside" or "auto"; "auto" places regions inside the chord when possible. |
| `position` | `"identity"` | Character/enum; choices below | Standard ggplot2 layer arguments. |
| `show.legend` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional arguments passed to geom_polygon() |

</details>

<span id="api-geom_seq_scale"></span>
<details><summary><code>geom_seq_scale()</code></summary>


[Full reference](reference/geom_seq_scale.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard ggplot2 layer inputs. Optional scale_length maps one positive integer length per bar, in the sequence unit; label supplies its text. Multiple bars require explicit x and y layout positions in data (or their role mappings), to avoid accidental overprinting. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard ggplot2 layer inputs. Optional scale_length maps one positive integer length per bar, in the sequence unit; label supplies its text. Multiple bars require explicit x and y layout positions in data (or their role mappings), to avoid accidental overprinting. |
| `length` | `NULL` | Nullable R object; type/inheritance below | Optional fixed positive integer length in the sequence unit. Cannot be combined with a scale_length mapping. NULL uses input lengths or an automatic reference length. |
| `label` | `NULL` | Nullable R object; type/inheritance below | Optional fixed character label, or a function receiving the chosen lengths. NULL uses input labels or formatted bp/aa text. |
| `location` | `c("bottomleft", "bottomright", "topleft", "topright")` | Vector; choices versus multiple selection below | Corner of the sequence envelope used for automatic placement. |
| `offset` | `grid::unit(c(0, 6), "mm")` | grid physical-unit object | Two non-negative distances, horizontal then vertical. Absolute grid units are supported; plain numeric values are layout-unit distances. The default is 0 mm horizontally and 6 mm vertically. |
| `tick_height` | `grid::unit(1.5, "mm")` | grid physical-unit object | Non-negative absolute grid units or layout-unit numbers for the end ticks and space between the bar and caption. Physical conversion uses the current render context when available; device-free builds retain an explicitly estimated measurement status. |
| `label_gap` | `grid::unit(1.5, "mm")` | grid physical-unit object | Non-negative absolute grid units or layout-unit numbers for the end ticks and space between the bar and caption. Physical conversion uses the current render context when available; device-free builds retain an explicitly estimated measurement status. |
| `linewidth` | `0.4` | Numeric/integer; range in description | Fixed line colour/width and text size. Text size follows ordinary ggplot2 millimetres. These may also use normal mappings and scales when not supplied as fixed arguments. |
| `colour` | `"#3D4348"` | Character/enum; choices below | Fixed line colour/width and text size. Text size follows ordinary ggplot2 millimetres. These may also use normal mappings and scales when not supplied as fixed arguments. |
| `size` | `3` | Numeric/integer; range in description | Fixed line colour/width and text size. Text size follows ordinary ggplot2 millimetres. These may also use normal mappings and scales when not supplied as fixed arguments. |
| `position` | `"identity"` | Character/enum; choices below | Standard ggplot2 layer arguments. |
| `show.legend` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional fixed aesthetics passed to ggplot2::layer(). |

</details>

<span id="api-geom_track"></span>
<details><summary><code>geom_track()</code></summary>


[Full reference](reference/geom_track.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard layer input. Supply accver and source positions through x/pos, or interval start/end. Map y/track_value for identity values. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard layer input. Supply accver and source positions through x/pos, or interval start/end. Map y/track_value for identity values. |
| `geom` | `"point"` | Character/enum; choices below | Native geometry/statistic names. bar defaults to count; col and other geoms default to identity. |
| `stat` | `NULL` | Nullable R object; type/inheritance below | Native geometry/statistic names. bar defaults to count; col and other geoms default to identity. |
| `track` | `NULL` | Nullable R object; type/inheritance below | A physical position_track() declaration. |
| `position` | `NULL` | Position object or documented string | Data position such as identity, stack, dodge. The legacy position=position_track(...) spelling still declares the physical band. |
| `bounds` | `c("error", "clip")` | Vector; choices versus multiple selection below | Coordinate error/clip policy. |
| `baseline` | `0` | Numeric/integer; range in description | Numeric identity-column baseline. Nonzero baselines require position='identity' so stacking does not change their meaning. |
| `show.legend` | `NA` | R object; constraints below | Standard layer controls. |
| `inherit.aes` | `FALSE` | Logical | Standard layer controls. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Native geom/stat parameters and fixed styles. |

</details>

<span id="api-geom_track_bar"></span>
<details><summary><code>geom_track_bar()</code></summary>


[Full reference](reference/geom_track.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard layer input. Supply accver and source positions through x/pos, or interval start/end. Map y/track_value for identity values. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard layer input. Supply accver and source positions through x/pos, or interval start/end. Map y/track_value for identity values. |
| `position` | `NULL` | Position object or documented string | Data position such as identity, stack, dodge. The legacy position=position_track(...) spelling still declares the physical band. |
| `baseline` | `0` | Numeric/integer; range in description | Numeric identity-column baseline. Nonzero baselines require position='identity' so stacking does not change their meaning. |
| `fill` | `"#86AFB9"` | Character/enum; choices below | Fixed styles for the count shortcut. |
| `colour` | `NA` | R object; constraints below | Fixed styles for the count shortcut. |
| `alpha` | `1` | Numeric/integer; range in description | Fixed styles for the count shortcut. |
| `linewidth` | `0.25` | Numeric/integer; range in description | Fixed styles for the count shortcut. |
| `linetype` | `1` | Numeric/integer; range in description | Fixed styles for the count shortcut. |
| `bounds` | `c("error", "clip")` | Vector; choices versus multiple selection below | Coordinate error/clip policy. |
| `show.legend` | `NA` | R object; constraints below | Standard layer controls. |
| `inherit.aes` | `FALSE` | Logical | Standard layer controls. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Native geom/stat parameters and fixed styles. |

</details>

<span id="api-geom_track_col"></span>
<details><summary><code>geom_track_col()</code></summary>


[Full reference](reference/geom_track.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard layer input. Supply accver and source positions through x/pos, or interval start/end. Map y/track_value for identity values. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard layer input. Supply accver and source positions through x/pos, or interval start/end. Map y/track_value for identity values. |
| `position` | `NULL` | Position object or documented string | Data position such as identity, stack, dodge. The legacy position=position_track(...) spelling still declares the physical band. |
| `track` | `NULL` | Nullable R object; type/inheritance below | A physical position_track() declaration. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Native geom/stat parameters and fixed styles. |
| `show.legend` | `NA` | R object; constraints below | Standard layer controls. |
| `inherit.aes` | `FALSE` | Logical | Standard layer controls. |

</details>

<span id="api-geom_track_line"></span>
<details><summary><code>geom_track_line()</code></summary>


[Full reference](reference/geom_track_tile.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard ggplot2 mapping and layer input. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard ggplot2 mapping and layer input. |
| `position` | `Required/alternative entry; see description` | Position object or documented string | A position_track() declaration. Every layer requires one. |
| `colour` | `"#387D88"` | Character/enum; choices below | Fixed styles; fill and colour mappings use the dedicated track roles and scoped scales. |
| `alpha` | `1` | Numeric/integer; range in description | Fixed styles; fill and colour mappings use the dedicated track roles and scoped scales. |
| `linewidth` | `0.4` | Numeric/integer; range in description | Fixed styles; fill and colour mappings use the dedicated track roles and scoped scales. |
| `linetype` | `1` | Numeric/integer; range in description | Fixed styles; fill and colour mappings use the dedicated track roles and scoped scales. |
| `closed` | `FALSE` | Logical | Explicitly close a complete circular line series. |
| `bounds` | `c("error", "clip")` | Vector; choices versus multiple selection below | "error" (default) or "clip" for genomic coordinates. Clipping is reported and never changes source values. |
| `show.legend` | `NA` | R object; constraints below | Standard layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Extra geom parameters. |

</details>

<span id="api-geom_track_tile"></span>
<details><summary><code>geom_track_tile()</code></summary>


[Full reference](reference/geom_track_tile.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard ggplot2 mapping and layer input. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard ggplot2 mapping and layer input. |
| `position` | `Required/alternative entry; see description` | Position object or documented string | A position_track() declaration. Every layer requires one. |
| `fill` | `"#86AFB9"` | Character/enum; choices below | Fixed styles; fill and colour mappings use the dedicated track roles and scoped scales. |
| `colour` | `NA` | R object; constraints below | Fixed styles; fill and colour mappings use the dedicated track roles and scoped scales. |
| `alpha` | `1` | Numeric/integer; range in description | Fixed styles; fill and colour mappings use the dedicated track roles and scoped scales. |
| `linewidth` | `0.25` | Numeric/integer; range in description | Fixed styles; fill and colour mappings use the dedicated track roles and scoped scales. |
| `linetype` | `1` | Numeric/integer; range in description | Fixed styles; fill and colour mappings use the dedicated track roles and scoped scales. |
| `bounds` | `c("error", "clip")` | Vector; choices versus multiple selection below | "error" (default) or "clip" for genomic coordinates. Clipping is reported and never changes source values. |
| `show.legend` | `NA` | R object; constraints below | Standard layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Extra geom parameters. |

</details>


## Local placement {#positions}

<span id="api-position_feature_stack"></span>
<details><summary><code>position_feature_stack()</code></summary>


[Full reference](reference/position_feature_stack.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `spacing` | `0.1` | Numeric/integer; range in description | Positive minimum radial distance between adjacent feature lanes. Circular plasmid layout enlarges the first transition into a main label corridor and reserves at least one label track between deeper overlapping feature lanes. |
| `side` | `c("strand", "outside", "inside")` | Vector; choices versus multiple selection below | Legacy track-side selector retained for compatibility when base_position is NULL. Its historical names retain their existing geometry. |
| `base_position` | `NULL` | Position object or documented string | Optional identity, strand, or plasmid Position applied before interval-overlap detection and lane allocation. |

</details>

<span id="api-position_plasmid"></span>
<details><summary><code>position_plasmid()</code></summary>


[Full reference](reference/position_plasmid.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `offset` | `-0.1` | Numeric/integer; range in description | Signed local-normal offset. Negative values place the shared band inside a circular sequence and positive values place it outside. |

</details>

<span id="api-position_strand"></span>
<details><summary><code>position_strand()</code></summary>


[Full reference](reference/position_strand.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `offset` | `0.1` | Numeric/integer; range in description | A finite scalar or a flexible per-sequence/per-strand specification. |

</details>

<span id="api-position_track"></span>
<details><summary><code>position_track()</code></summary>


[Full reference](reference/position_track.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `track` | `Required/alternative entry; see description` | R object; constraints below | One stable nonempty character ID. |
| `side` | `c("outside", "inside")` | Vector; choices versus multiple selection below | "outside" or "inside". |
| `order` | `1` | Numeric/integer; range in description | Finite ordering priority, nearest band first. Ties use ID order. |
| `width` | `grid::unit(4, "mm")` | grid physical-unit object | Positive band width and nonnegative preceding gap. Absolute grid units (mm, cm, inches or points) are resolved against the render context. Numeric values retain native layout units. Numeric per-sequence vectors/lists follow the usual sequence parameter formats. |
| `gap` | `grid::unit(1, "mm")` | grid physical-unit object | Positive band width and nonnegative preceding gap. Absolute grid units (mm, cm, inches or points) are resolved against the render context. Numeric values retain native layout units. Numeric per-sequence vectors/lists follow the usual sequence parameter formats. |

</details>


## Link statistics {#statistics}

<span id="api-stat_ribbon_bundle"></span>
<details><summary><code>stat_ribbon_bundle()</code></summary>


[Full reference](reference/stat_ribbon_bundle.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard layer mapping and data arguments. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard layer mapping and data arguments. |
| `bins` | `80L` | Numeric/integer; range in description | Passed to bundle_ggchord_ribbons(). |
| `min_bundle` | `2L` | Numeric/integer; range in description | Passed to bundle_ggchord_ribbons(). |
| `weight` | `c("length", "count", "pident")` | Vector; choices versus multiple selection below | Passed to bundle_ggchord_ribbons(). |
| `group_by` | `NULL` | Nullable R object; type/inheritance below | Passed to bundle_ggchord_ribbons(). |
| `position` | `"identity"` | Character/enum; choices below | Standard ggplot2 layer arguments. |
| `show.legend` | `TRUE` | Logical | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Ribbon geometry and appearance arguments passed to geom_link_ribbon(). |

</details>

<span id="api-stat_ribbon_density"></span>
<details><summary><code>stat_ribbon_density()</code></summary>


[Full reference](reference/stat_ribbon_density.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes mapping | Standard layer mapping and data arguments. |
| `data` | `NULL` | Table/supported object; data may be a function | Standard layer mapping and data arguments. |
| `bins` | `80L` | Numeric/integer; range in description | Positive integer number of normalized midpoint bins per sequence. |
| `weight` | `c("length", "count", "pident")` | Vector; choices versus multiple selection below | Density weight: "length", "count", or "pident". |
| `group_by` | `NULL` | Nullable R object; type/inheritance below | Optional columns that define independent density groups. |
| `position` | `"identity"` | Character/enum; choices below | Standard ggplot2 layer arguments. |
| `show.legend` | `TRUE` | Logical | Standard ggplot2 layer arguments. |
| `inherit.aes` | `FALSE` | Logical | Standard ggplot2 layer arguments. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Ribbon geometry and appearance arguments passed to geom_link_ribbon(). |

</details>


## Role scales {#scales}

<span id="api-scale_feature_fill_manual"></span>
<details><summary><code>scale_feature_fill_manual()</code></summary>


[Full reference](reference/scale_gene_fill_manual.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to ggplot2::scale_fill_manual(). |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | A set of fill values to map data values to. |
| `limits` | `NULL` | Vector/scale rule; interpreted by scale type | Optional display order. Named feature values use their name order by default. |

</details>

<span id="api-scale_feature_fill_plasmid"></span>
<details><summary><code>scale_feature_fill_plasmid()</code></summary>


[Full reference](reference/scale_feature_fill_plasmid.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional scale arguments. |
| `limits` | `NULL` | Vector/scale rule; interpreted by scale type | Optional feature-type order. |

</details>

<span id="api-scale_feature_shape_manual"></span>
<details><summary><code>scale_feature_shape_manual()</code></summary>


[Full reference](reference/scale_feature_shape_manual.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to ggplot2::discrete_scale(). |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | A named or unnamed character vector of feature geometry names. |
| `name` | `"Feature"` | Character/enum; choices below | Scale and guide title, default "Feature". |
| `limits` | `NULL` | Vector/scale rule; interpreted by scale type | Optional category order. Named values use their name order by default so shape and fill guides can merge cleanly. |
| `guide` | `"none"` | Character/enum; choices below | Guide specification. The default is "none" because geom_feature() combines shape silhouettes into its fill legend when both aesthetics describe the same feature categories. Supply guide_ggchord_legend() for a separate shape guide. |

</details>

<span id="api-scale_feature_shape_plasmid"></span>
<details><summary><code>scale_feature_shape_plasmid()</code></summary>


[Full reference](reference/scale_feature_fill_plasmid.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional scale arguments. |
| `limits` | `NULL` | Vector/scale rule; interpreted by scale type | Optional feature-type order. |
| `guide` | `"none"` | Character/enum; choices below | Guide used by the shape preset. |

</details>

<span id="api-scale_gene_fill_manual"></span>
<details><summary><code>scale_gene_fill_manual()</code></summary>


[Full reference](reference/scale_gene_fill_manual.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to ggplot2::scale_fill_manual(). |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | A set of fill values to map data values to. |

</details>

<span id="api-scale_link_alpha_continuous"></span>
<details><summary><code>scale_link_alpha_continuous()</code></summary>


[Full reference](reference/scale_link.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |

</details>

<span id="api-scale_link_alpha_manual"></span>
<details><summary><code>scale_link_alpha_manual()</code></summary>


[Full reference](reference/scale_link.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Manual aesthetic values. |

</details>

<span id="api-scale_link_color_gradientn"></span>
<details><summary><code>scale_link_color_gradientn()</code></summary>


[Full reference](reference/scale_link.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |

</details>

<span id="api-scale_link_color_manual"></span>
<details><summary><code>scale_link_color_manual()</code></summary>


[Full reference](reference/scale_link.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Manual aesthetic values. |

</details>

<span id="api-scale_link_colour_gradientn"></span>
<details><summary><code>scale_link_colour_gradientn()</code></summary>


[Full reference](reference/scale_link.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |

</details>

<span id="api-scale_link_colour_manual"></span>
<details><summary><code>scale_link_colour_manual()</code></summary>


[Full reference](reference/scale_link.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Manual aesthetic values. |

</details>

<span id="api-scale_link_linetype_manual"></span>
<details><summary><code>scale_link_linetype_manual()</code></summary>


[Full reference](reference/scale_link.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Manual aesthetic values. |

</details>

<span id="api-scale_link_linewidth_continuous"></span>
<details><summary><code>scale_link_linewidth_continuous()</code></summary>


[Full reference](reference/scale_link.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |

</details>

<span id="api-scale_link_linewidth_manual"></span>
<details><summary><code>scale_link_linewidth_manual()</code></summary>


[Full reference](reference/scale_link.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Manual aesthetic values. |

</details>

<span id="api-scale_region_fill_manual"></span>
<details><summary><code>scale_region_fill_manual()</code></summary>


[Full reference](reference/scale_gene_fill_manual.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to ggplot2::scale_fill_manual(). |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | A set of fill values to map data values to. |

</details>

<span id="api-scale_ribbon_alpha_continuous"></span>
<details><summary><code>scale_ribbon_alpha_continuous()</code></summary>


[Full reference](reference/scale_ribbon_alpha_continuous.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |
| `range` | `NULL` | Vector/scale rule; interpreted by scale type | Output alpha range for a continuous scale. NULL uses the corresponding ggplot2 scale's default range. |

</details>

<span id="api-scale_ribbon_alpha_manual"></span>
<details><summary><code>scale_ribbon_alpha_manual()</code></summary>


[Full reference](reference/scale_ribbon_alpha_continuous.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Manual aesthetic values. |

</details>

<span id="api-scale_ribbon_color_manual"></span>
<details><summary><code>scale_ribbon_color_manual()</code></summary>


[Full reference](reference/scale_ribbon_alpha_continuous.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Manual aesthetic values. |

</details>

<span id="api-scale_ribbon_colour_manual"></span>
<details><summary><code>scale_ribbon_colour_manual()</code></summary>


[Full reference](reference/scale_ribbon_alpha_continuous.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Manual aesthetic values. |

</details>

<span id="api-scale_ribbon_fill_gradientn"></span>
<details><summary><code>scale_ribbon_fill_gradientn()</code></summary>


[Full reference](reference/scale_ribbon_fill_stepsn.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |
| `colours` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Gradient colours. |
| `values` | `NULL` | Vector/scale rule; interpreted by scale type | Manual values, or positions for gradient colours where supported by the corresponding ggplot2 scale. |
| `colors` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Gradient colours. |
| `guide` | `ggplot2::guide_colourbar(available_aes = "ribbon_fill")` | Guide object or name | A guide object; defaults accept the ribbon_fill aesthetic. |

</details>

<span id="api-scale_ribbon_fill_identity"></span>
<details><summary><code>scale_ribbon_fill_identity()</code></summary>


[Full reference](reference/scale_ribbon_fill_stepsn.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |

</details>

<span id="api-scale_ribbon_fill_manual"></span>
<details><summary><code>scale_ribbon_fill_manual()</code></summary>


[Full reference](reference/scale_ribbon_fill_stepsn.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Manual values, or positions for gradient colours where supported by the corresponding ggplot2 scale. |

</details>

<span id="api-scale_ribbon_fill_stepsn"></span>
<details><summary><code>scale_ribbon_fill_stepsn()</code></summary>


[Full reference](reference/scale_ribbon_fill_stepsn.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |
| `colours` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Gradient colours. |
| `values` | `NULL` | Vector/scale rule; interpreted by scale type | Manual values, or positions for gradient colours where supported by the corresponding ggplot2 scale. |
| `colors` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Gradient colours. |
| `guide` | `ggplot2::guide_coloursteps(available_aes = "ribbon_fill")` | Guide object or name | A guide object; defaults accept the ribbon_fill aesthetic. |

</details>

<span id="api-scale_ribbon_linetype_manual"></span>
<details><summary><code>scale_ribbon_linetype_manual()</code></summary>


[Full reference](reference/scale_ribbon_alpha_continuous.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to the corresponding ggplot2 scale. |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Manual aesthetic values. |

</details>

<span id="api-scale_seq_color_manual"></span>
<details><summary><code>scale_seq_color_manual()</code></summary>


[Full reference](reference/scale_seq_colour_manual.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to ggplot2::scale_colour_manual(). |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | A set of aesthetic values to map data values to. |

</details>

<span id="api-scale_seq_colour_manual"></span>
<details><summary><code>scale_seq_colour_manual()</code></summary>


[Full reference](reference/scale_seq_colour_manual.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to ggplot2::scale_colour_manual(). |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | A set of aesthetic values to map data values to. |

</details>

<span id="api-scale_seq_position_continuous"></span>
<details><summary><code>scale_seq_position_continuous()</code></summary>


[Full reference](reference/scale_seq_position_continuous.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `name` | `ggplot2::waiver()` | R object; constraints below | Scale name; position guides are disabled by default. |
| `breaks` | `ggplot2::waiver()` | Vector/scale rule; interpreted by scale type | Standard continuous-scale controls. |
| `minor_breaks` | `ggplot2::waiver()` | Vector/scale rule; interpreted by scale type | Standard continuous-scale controls. |
| `labels` | `ggplot2::waiver()` | R object; constraints below | Standard continuous-scale controls. |
| `limits` | `NULL` | Vector/scale rule; interpreted by scale type | Standard continuous-scale controls. |
| `expand` | `ggplot2::waiver()` | R object; constraints below | Standard continuous-scale controls. |
| `oob` | `scales::censor` | R object; constraints below | Standard continuous-scale controls. |
| `transform` | `"identity"` | Character/enum; choices below | Standard continuous-scale controls. |
| `n.breaks` | `6` | Numeric/integer; range in description | Standard continuous-scale controls. |
| `coordinates` | `c("relative", "source")` | Vector; choices versus multiple selection below | In linear layouts, use relative locus coordinates or original source coordinates retained by focus_ggchord_data(). This affects axis breaks and labels, not the underlying genomic geometry. Breaks returned in transformed scale space are converted back to genomic positions in every coordinate mode. Explicit breaks = NULL and minor_breaks = NULL suppress the corresponding ticks. |
| `sequences` | `NULL` | Nullable R object; type/inheritance below | Optional character vector selecting sequence axes by their displayed accver. NULL keeps all axes; character() hides every axis. Focused or segmented views use window IDs, not parent source accessions. Unknown or undisplayed identifiers produce an error. This controls axis lines, ticks and labels together without dropping sequence geometry. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional arguments passed to ggplot2::continuous_scale(). |

</details>

<span id="api-scale_seq_ring_manual"></span>
<details><summary><code>scale_seq_ring_manual()</code></summary>


[Full reference](reference/scale_seq_ring_manual.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Arguments passed to ggplot2::discrete_scale(). |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | A named or unnamed numeric vector of positive radii. |
| `name` | `"Ring"` | Character/enum; choices below | Scale title. The guide is hidden by default because the scale controls geometry rather than a visible legend aesthetic. |
| `limits` | `NULL` | Vector/scale rule; interpreted by scale type | Optional ring order. |
| `guide` | `"none"` | Character/enum; choices below | Guide specification, default "none". |

</details>

<span id="api-scale_track_colour_gradientn"></span>
<details><summary><code>scale_track_colour_gradientn()</code></summary>


[Full reference](reference/scale_track_fill_gradientn.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `track` | `Required/alternative entry; see description` | R object; constraints below | Track ID. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Standard ggplot2 scale arguments. |
| `colours` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Gradient colours. |
| `guide` | `NULL` | Guide object or name | Optional standard guide for the fill/colour scale. |

</details>

<span id="api-scale_track_colour_manual"></span>
<details><summary><code>scale_track_colour_manual()</code></summary>


[Full reference](reference/scale_track_fill_gradientn.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `track` | `Required/alternative entry; see description` | R object; constraints below | Track ID. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Standard ggplot2 scale arguments. |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Manual category colours. |

</details>

<span id="api-scale_track_fill_gradientn"></span>
<details><summary><code>scale_track_fill_gradientn()</code></summary>


[Full reference](reference/scale_track_fill_gradientn.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `track` | `Required/alternative entry; see description` | R object; constraints below | Track ID. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Standard ggplot2 scale arguments. |
| `colours` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Gradient colours. |
| `guide` | `NULL` | Guide object or name | Optional standard guide for the fill/colour scale. |

</details>

<span id="api-scale_track_fill_manual"></span>
<details><summary><code>scale_track_fill_manual()</code></summary>


[Full reference](reference/scale_track_fill_gradientn.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `track` | `Required/alternative entry; see description` | R object; constraints below | Track ID. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Standard ggplot2 scale arguments. |
| `values` | `Required/alternative entry; see description` | Vector/scale rule; interpreted by scale type | Manual category colours. |

</details>

<span id="api-scale_track_value_continuous"></span>
<details><summary><code>scale_track_value_continuous()</code></summary>


[Full reference](reference/scale_track_value_continuous.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `track` | `Required/alternative entry; see description` | R object; constraints below | Track ID from position_track(). |
| `name` | `track` | R object; constraints below | Axis title, default the track ID. |
| `limits` | `NULL` | Vector/scale rule; interpreted by scale type | Standard continuous-scale controls. |
| `breaks` | `ggplot2::waiver()` | Vector/scale rule; interpreted by scale type | Standard continuous-scale controls. |
| `labels` | `ggplot2::waiver()` | R object; constraints below | Standard continuous-scale controls. |
| `transform` | `"identity"` | Character/enum; choices below | Standard continuous-scale controls. |
| `oob` | `scales::censor` | R object; constraints below | Standard continuous-scale controls. |
| `unit` | `NULL` | Nullable R object; type/inheritance below | Optional data-unit label retained in track exports. |
| `guide` | `c("axis", "none")` | Vector; choices versus multiple selection below | "axis" or "none" for the numerical track axis, or a native ggplot2::guide_axis() with local theme (axis.line/text/title), angle and check.overlap. Other guide_axis options fail explicitly because this sequence-local numerical reference is not a Cartesian panel axis. Local guide styling overrides global ggchord.track.axis.* theme elements. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional ggplot2 continuous-scale arguments. |

</details>


## Legends and guides {#guides}

<span id="api-guide_ggchord_colorbar"></span>
<details><summary><code>guide_ggchord_colorbar()</code></summary>


[Full reference](reference/guide_ggchord_legend.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `title` | `ggplot2::waiver()` | R object; constraints below | Guide title. |
| `theme` | `NULL` | Nullable R object; type/inheritance below | Optional guide-specific theme. |
| `nbin` | `NULL` | Nullable R object; type/inheritance below | Number of colourbar bins. |
| `display` | `"raster"` | Character/enum; choices below | Colourbar display method. |
| `alpha` | `NA` | R object; constraints below | Colourbar alpha override. |
| `draw.ulim` | `TRUE` | Logical | Whether to draw upper/lower limit ticks. |
| `draw.llim` | `TRUE` | Logical | Whether to draw upper/lower limit ticks. |
| `angle` | `NULL` | Nullable R object; type/inheritance below | Tick-label angle. |
| `position` | `NULL` | Position object or documented string | Guide position and direction. |
| `direction` | `NULL` | Nullable R object; type/inheritance below | Guide position and direction. |
| `reverse` | `FALSE` | Logical | Whether to reverse the key order. |
| `order` | `0` | Numeric/integer; range in description | Guide order. |
| `available_aes` | `c("colour", "color", "fill", "ribbon_fill", "gene_fill", "feature_fill",      "region_fill")` | Vector; choices versus multiple selection below | Aesthetics accepted by the guide. |
| `size_scale` | `1` | Numeric/integer; range in description | Positive multiplier for keys, spacing and typography. ggchord's automatically generated guides derive this from the active export device; direct calls default to 1. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional arguments passed to the corresponding ggplot2 guide. |

</details>

<span id="api-guide_ggchord_colourbar"></span>
<details><summary><code>guide_ggchord_colourbar()</code></summary>


[Full reference](reference/guide_ggchord_legend.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `title` | `ggplot2::waiver()` | R object; constraints below | Guide title. |
| `theme` | `NULL` | Nullable R object; type/inheritance below | Optional guide-specific theme. |
| `nbin` | `NULL` | Nullable R object; type/inheritance below | Number of colourbar bins. |
| `display` | `"raster"` | Character/enum; choices below | Colourbar display method. |
| `alpha` | `NA` | R object; constraints below | Colourbar alpha override. |
| `draw.ulim` | `TRUE` | Logical | Whether to draw upper/lower limit ticks. |
| `draw.llim` | `TRUE` | Logical | Whether to draw upper/lower limit ticks. |
| `angle` | `NULL` | Nullable R object; type/inheritance below | Tick-label angle. |
| `position` | `NULL` | Position object or documented string | Guide position and direction. |
| `direction` | `NULL` | Nullable R object; type/inheritance below | Guide position and direction. |
| `reverse` | `FALSE` | Logical | Whether to reverse the key order. |
| `order` | `0` | Numeric/integer; range in description | Guide order. |
| `available_aes` | `c("colour", "color", "fill", "ribbon_fill", "gene_fill", "feature_fill",      "region_fill")` | Vector; choices versus multiple selection below | Aesthetics accepted by the guide. |
| `size_scale` | `1` | Numeric/integer; range in description | Positive multiplier for keys, spacing and typography. ggchord's automatically generated guides derive this from the active export device; direct calls default to 1. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional arguments passed to the corresponding ggplot2 guide. |

</details>

<span id="api-guide_ggchord_legend"></span>
<details><summary><code>guide_ggchord_legend()</code></summary>


[Full reference](reference/guide_ggchord_legend.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `title` | `ggplot2::waiver()` | R object; constraints below | Guide title. |
| `theme` | `NULL` | Nullable R object; type/inheritance below | Optional guide-specific theme. |
| `position` | `NULL` | Position object or documented string | Guide position and direction. |
| `direction` | `NULL` | Nullable R object; type/inheritance below | Guide position and direction. |
| `override.aes` | `list()` | Named list | A list of legend-key aesthetic overrides. |
| `nrow` | `NULL` | Nullable R object; type/inheritance below | Legend key layout. |
| `ncol` | `NULL` | Nullable R object; type/inheritance below | Legend key layout. |
| `reverse` | `FALSE` | Logical | Whether to reverse the key order. |
| `order` | `0` | Numeric/integer; range in description | Guide order. |
| `size_scale` | `1` | Numeric/integer; range in description | Positive multiplier for keys, spacing and typography. ggchord's automatically generated guides derive this from the active export device; direct calls default to 1. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional arguments passed to the corresponding ggplot2 guide. |

</details>


## Themes and physical appearance {#themes}

<span id="api-theme_ggchord"></span>
<details><summary><code>theme_ggchord()</code></summary>


[Full reference](reference/theme_ggchord.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `axis` | `NULL` | Nullable R object; type/inheritance below | NULL to show automatic sequence axes, or ggplot2::element_blank() to hide the complete axis decoration. |
| `axis.line` | `NULL` | Nullable R object; type/inheritance below | Axis line elements. |
| `axis.ticks` | `NULL` | Nullable R object; type/inheritance below | Axis line elements. |
| `axis.minor.ticks` | `NULL` | Nullable R object; type/inheritance below | Axis line elements. |
| `axis.text` | `NULL` | Nullable R object; type/inheritance below | Axis text element. |
| `axis.gap` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.ticks.length` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.minor.ticks.length` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.text.offset` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.text.orientation` | `NULL` | Nullable R object; type/inheritance below | One of "horizontal", "parallel", "perpendicular", or a finite numeric angle. In circular layouts, parallel labels follow the circle glyph by glyph. Horizontal and numeric orientations retain straight text. |
| `axis.text.check.overlap` | `NULL` | Nullable R object; type/inheritance below | Whether overlapping axis labels are hidden. Defaults to TRUE. |
| `seq.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `gene.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `gene.label.segment` | `NULL` | Nullable R object; type/inheritance below | Line elements for gene, feature and restriction-site leader lines. |
| `feature.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `feature.label.segment` | `NULL` | Nullable R object; type/inheritance below | Line elements for gene, feature and restriction-site leader lines. |
| `restriction.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `restriction.label.segment` | `NULL` | Nullable R object; type/inheritance below | Line elements for gene, feature and restriction-site leader lines. |
| `seq.center.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `seq.center.name` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `seq.center.length` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `legend.seq` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.seq.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.seq.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.seq.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.seq.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.seq.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.seq.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.seq.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.seq.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.seq.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.seq.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.seq.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.seq.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.seq.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.link` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.position` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.direction` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.background` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key.size` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key.width` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key.height` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.text` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.text.position` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.title` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.title.position` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.margin` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.spacing` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.ribbon` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.ribbon.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.ribbon.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.ribbon.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.ribbon.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.ribbon.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.ribbon.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.ribbon.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.ribbon.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.ribbon.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.ribbon.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.ribbon.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.ribbon.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.ribbon.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.ribbon.ticks` | `NULL` | Nullable R object; type/inheritance below | Ribbon colourbar lines. |
| `legend.ribbon.ticks.length` | `NULL` | Nullable R object; type/inheritance below | Ribbon colourbar tick length. |
| `legend.ribbon.axis.line` | `NULL` | Nullable R object; type/inheritance below | Ribbon colourbar lines. |
| `legend.gene` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.gene.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.gene.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.gene.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.gene.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.gene.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.gene.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.gene.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.gene.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.gene.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.gene.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.gene.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.gene.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.gene.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.feature` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.feature.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.feature.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.feature.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.feature.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.feature.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.feature.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.feature.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.feature.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.feature.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.feature.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.feature.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.feature.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.feature.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.region` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.region.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.region.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.region.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.region.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.region.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.region.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.region.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.region.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.region.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.region.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.region.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.region.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.region.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Standard ggplot2 theme arguments, including common text, plot and legend settings. Common legend settings also supply role inheritance. Removed underscore-style arguments produce a migration error. |

</details>

<span id="api-theme_ggchord_dark"></span>
<details><summary><code>theme_ggchord_dark()</code></summary>


[Full reference](reference/theme_ggchord.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `axis` | `NULL` | Nullable R object; type/inheritance below | NULL to show automatic sequence axes, or ggplot2::element_blank() to hide the complete axis decoration. |
| `axis.line` | `NULL` | Nullable R object; type/inheritance below | Axis line elements. |
| `axis.ticks` | `NULL` | Nullable R object; type/inheritance below | Axis line elements. |
| `axis.minor.ticks` | `NULL` | Nullable R object; type/inheritance below | Axis line elements. |
| `axis.text` | `NULL` | Nullable R object; type/inheritance below | Axis text element. |
| `axis.gap` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.ticks.length` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.minor.ticks.length` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.text.offset` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.text.orientation` | `NULL` | Nullable R object; type/inheritance below | One of "horizontal", "parallel", "perpendicular", or a finite numeric angle. In circular layouts, parallel labels follow the circle glyph by glyph. Horizontal and numeric orientations retain straight text. |
| `axis.text.check.overlap` | `NULL` | Nullable R object; type/inheritance below | Whether overlapping axis labels are hidden. Defaults to TRUE. |
| `seq.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `gene.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `gene.label.segment` | `NULL` | Nullable R object; type/inheritance below | Line elements for gene, feature and restriction-site leader lines. |
| `feature.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `feature.label.segment` | `NULL` | Nullable R object; type/inheritance below | Line elements for gene, feature and restriction-site leader lines. |
| `restriction.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `restriction.label.segment` | `NULL` | Nullable R object; type/inheritance below | Line elements for gene, feature and restriction-site leader lines. |
| `seq.center.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `seq.center.name` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `seq.center.length` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `legend.seq` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.seq.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.seq.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.seq.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.seq.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.seq.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.seq.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.seq.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.seq.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.seq.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.seq.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.seq.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.seq.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.seq.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.link` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.position` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.direction` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.background` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key.size` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key.width` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key.height` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.text` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.text.position` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.title` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.title.position` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.margin` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.spacing` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.ribbon` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.ribbon.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.ribbon.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.ribbon.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.ribbon.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.ribbon.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.ribbon.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.ribbon.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.ribbon.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.ribbon.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.ribbon.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.ribbon.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.ribbon.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.ribbon.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.ribbon.ticks` | `NULL` | Nullable R object; type/inheritance below | Ribbon colourbar lines. |
| `legend.ribbon.ticks.length` | `NULL` | Nullable R object; type/inheritance below | Ribbon colourbar tick length. |
| `legend.ribbon.axis.line` | `NULL` | Nullable R object; type/inheritance below | Ribbon colourbar lines. |
| `legend.gene` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.gene.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.gene.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.gene.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.gene.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.gene.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.gene.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.gene.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.gene.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.gene.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.gene.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.gene.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.gene.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.gene.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.feature` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.feature.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.feature.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.feature.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.feature.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.feature.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.feature.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.feature.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.feature.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.feature.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.feature.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.feature.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.feature.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.feature.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.region` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.region.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.region.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.region.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.region.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.region.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.region.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.region.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.region.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.region.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.region.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.region.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.region.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.region.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Standard ggplot2 theme arguments, including common text, plot and legend settings. Common legend settings also supply role inheritance. Removed underscore-style arguments produce a migration error. |

</details>

<span id="api-theme_ggchord_minimal"></span>
<details><summary><code>theme_ggchord_minimal()</code></summary>


[Full reference](reference/theme_ggchord.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `axis` | `NULL` | Nullable R object; type/inheritance below | NULL to show automatic sequence axes, or ggplot2::element_blank() to hide the complete axis decoration. |
| `axis.line` | `NULL` | Nullable R object; type/inheritance below | Axis line elements. |
| `axis.ticks` | `NULL` | Nullable R object; type/inheritance below | Axis line elements. |
| `axis.minor.ticks` | `NULL` | Nullable R object; type/inheritance below | Axis line elements. |
| `axis.text` | `NULL` | Nullable R object; type/inheritance below | Axis text element. |
| `axis.gap` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.ticks.length` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.minor.ticks.length` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.text.offset` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.text.orientation` | `NULL` | Nullable R object; type/inheritance below | One of "horizontal", "parallel", "perpendicular", or a finite numeric angle. In circular layouts, parallel labels follow the circle glyph by glyph. Horizontal and numeric orientations retain straight text. |
| `axis.text.check.overlap` | `NULL` | Nullable R object; type/inheritance below | Whether overlapping axis labels are hidden. Defaults to TRUE. |
| `seq.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `gene.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `gene.label.segment` | `NULL` | Nullable R object; type/inheritance below | Line elements for gene, feature and restriction-site leader lines. |
| `feature.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `feature.label.segment` | `NULL` | Nullable R object; type/inheritance below | Line elements for gene, feature and restriction-site leader lines. |
| `restriction.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `restriction.label.segment` | `NULL` | Nullable R object; type/inheritance below | Line elements for gene, feature and restriction-site leader lines. |
| `seq.center.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `seq.center.name` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `seq.center.length` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `legend.seq` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.seq.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.seq.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.seq.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.seq.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.seq.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.seq.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.seq.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.seq.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.seq.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.seq.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.seq.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.seq.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.seq.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.link` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.position` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.direction` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.background` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key.size` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key.width` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key.height` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.text` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.text.position` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.title` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.title.position` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.margin` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.spacing` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.ribbon` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.ribbon.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.ribbon.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.ribbon.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.ribbon.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.ribbon.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.ribbon.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.ribbon.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.ribbon.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.ribbon.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.ribbon.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.ribbon.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.ribbon.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.ribbon.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.ribbon.ticks` | `NULL` | Nullable R object; type/inheritance below | Ribbon colourbar lines. |
| `legend.ribbon.ticks.length` | `NULL` | Nullable R object; type/inheritance below | Ribbon colourbar tick length. |
| `legend.ribbon.axis.line` | `NULL` | Nullable R object; type/inheritance below | Ribbon colourbar lines. |
| `legend.gene` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.gene.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.gene.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.gene.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.gene.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.gene.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.gene.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.gene.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.gene.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.gene.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.gene.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.gene.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.gene.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.gene.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.feature` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.feature.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.feature.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.feature.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.feature.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.feature.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.feature.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.feature.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.feature.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.feature.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.feature.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.feature.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.feature.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.feature.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.region` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.region.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.region.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.region.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.region.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.region.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.region.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.region.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.region.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.region.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.region.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.region.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.region.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.region.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Standard ggplot2 theme arguments, including common text, plot and legend settings. Common legend settings also supply role inheritance. Removed underscore-style arguments produce a migration error. |

</details>

<span id="api-theme_ggchord_plasmid"></span>
<details><summary><code>theme_ggchord_plasmid()</code></summary>


[Full reference](reference/theme_ggchord.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `axis` | `NULL` | Nullable R object; type/inheritance below | NULL to show automatic sequence axes, or ggplot2::element_blank() to hide the complete axis decoration. |
| `axis.line` | `NULL` | Nullable R object; type/inheritance below | Axis line elements. |
| `axis.ticks` | `NULL` | Nullable R object; type/inheritance below | Axis line elements. |
| `axis.minor.ticks` | `NULL` | Nullable R object; type/inheritance below | Axis line elements. |
| `axis.text` | `NULL` | Nullable R object; type/inheritance below | Axis text element. |
| `axis.gap` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.ticks.length` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.minor.ticks.length` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.text.offset` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.text.orientation` | `NULL` | Nullable R object; type/inheritance below | One of "horizontal", "parallel", "perpendicular", or a finite numeric angle. In circular layouts, parallel labels follow the circle glyph by glyph. Horizontal and numeric orientations retain straight text. |
| `axis.text.check.overlap` | `NULL` | Nullable R object; type/inheritance below | Whether overlapping axis labels are hidden. Defaults to TRUE. |
| `seq.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `gene.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `gene.label.segment` | `NULL` | Nullable R object; type/inheritance below | Line elements for gene, feature and restriction-site leader lines. |
| `feature.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `feature.label.segment` | `NULL` | Nullable R object; type/inheritance below | Line elements for gene, feature and restriction-site leader lines. |
| `restriction.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `restriction.label.segment` | `NULL` | Nullable R object; type/inheritance below | Line elements for gene, feature and restriction-site leader lines. |
| `seq.center.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `seq.center.name` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `seq.center.length` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `legend.seq` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.seq.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.seq.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.seq.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.seq.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.seq.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.seq.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.seq.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.seq.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.seq.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.seq.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.seq.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.seq.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.seq.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.link` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.position` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.direction` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.background` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key.size` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key.width` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key.height` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.text` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.text.position` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.title` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.title.position` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.margin` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.spacing` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.ribbon` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.ribbon.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.ribbon.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.ribbon.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.ribbon.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.ribbon.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.ribbon.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.ribbon.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.ribbon.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.ribbon.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.ribbon.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.ribbon.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.ribbon.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.ribbon.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.ribbon.ticks` | `NULL` | Nullable R object; type/inheritance below | Ribbon colourbar lines. |
| `legend.ribbon.ticks.length` | `NULL` | Nullable R object; type/inheritance below | Ribbon colourbar tick length. |
| `legend.ribbon.axis.line` | `NULL` | Nullable R object; type/inheritance below | Ribbon colourbar lines. |
| `legend.gene` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.gene.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.gene.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.gene.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.gene.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.gene.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.gene.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.gene.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.gene.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.gene.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.gene.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.gene.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.gene.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.gene.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.feature` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.feature.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.feature.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.feature.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.feature.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.feature.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.feature.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.feature.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.feature.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.feature.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.feature.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.feature.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.feature.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.feature.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.region` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.region.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.region.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.region.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.region.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.region.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.region.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.region.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.region.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.region.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.region.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.region.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.region.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.region.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Standard ggplot2 theme arguments, including common text, plot and legend settings. Common legend settings also supply role inheritance. Removed underscore-style arguments produce a migration error. |

</details>

<span id="api-theme_ggchord_publication"></span>
<details><summary><code>theme_ggchord_publication()</code></summary>


[Full reference](reference/theme_ggchord.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `axis` | `NULL` | Nullable R object; type/inheritance below | NULL to show automatic sequence axes, or ggplot2::element_blank() to hide the complete axis decoration. |
| `axis.line` | `NULL` | Nullable R object; type/inheritance below | Axis line elements. |
| `axis.ticks` | `NULL` | Nullable R object; type/inheritance below | Axis line elements. |
| `axis.minor.ticks` | `NULL` | Nullable R object; type/inheritance below | Axis line elements. |
| `axis.text` | `NULL` | Nullable R object; type/inheritance below | Axis text element. |
| `axis.gap` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.ticks.length` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.minor.ticks.length` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.text.offset` | `NULL` | Nullable R object; type/inheritance below | Scalar grid::unit() objects controlling physical axis spacing. A negative axis.gap places the circular axis inside the backbone; a positive value places it outside. Zero attaches ticks to the inner backbone edge and is the circular default, including the plasmid preset. Circular major/minor tick lengths default to 1.5/0.7 mm, or 2.2/0.8 mm with the plasmid preset. Their default line widths are 0.45/0.3 mm, or 0.6/0.35 mm with that preset. Explicit length and line settings take priority. |
| `axis.text.orientation` | `NULL` | Nullable R object; type/inheritance below | One of "horizontal", "parallel", "perpendicular", or a finite numeric angle. In circular layouts, parallel labels follow the circle glyph by glyph. Horizontal and numeric orientations retain straight text. |
| `axis.text.check.overlap` | `NULL` | Nullable R object; type/inheritance below | Whether overlapping axis labels are hidden. Defaults to TRUE. |
| `seq.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `gene.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `gene.label.segment` | `NULL` | Nullable R object; type/inheritance below | Line elements for gene, feature and restriction-site leader lines. |
| `feature.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `feature.label.segment` | `NULL` | Nullable R object; type/inheritance below | Line elements for gene, feature and restriction-site leader lines. |
| `restriction.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `restriction.label.segment` | `NULL` | Nullable R object; type/inheritance below | Line elements for gene, feature and restriction-site leader lines. |
| `seq.center.label` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `seq.center.name` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `seq.center.length` | `NULL` | Nullable R object; type/inheritance below | Text elements for sequence, gene, feature, restriction-site and circular centre labels. |
| `legend.seq` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.seq.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.seq.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.seq.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.seq.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.seq.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.seq.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.seq.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.seq.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.seq.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.seq.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.seq.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.seq.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.seq.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.link` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.position` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.direction` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.background` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key.size` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key.width` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.key.height` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.text` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.text.position` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.title` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.title.position` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.margin` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.link.spacing` | `NULL` | Nullable R object; type/inheritance below | Point-link guide appearance; follows the same role-specific rules. |
| `legend.ribbon` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.ribbon.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.ribbon.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.ribbon.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.ribbon.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.ribbon.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.ribbon.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.ribbon.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.ribbon.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.ribbon.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.ribbon.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.ribbon.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.ribbon.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.ribbon.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.ribbon.ticks` | `NULL` | Nullable R object; type/inheritance below | Ribbon colourbar lines. |
| `legend.ribbon.ticks.length` | `NULL` | Nullable R object; type/inheritance below | Ribbon colourbar tick length. |
| `legend.ribbon.axis.line` | `NULL` | Nullable R object; type/inheritance below | Ribbon colourbar lines. |
| `legend.gene` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.gene.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.gene.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.gene.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.gene.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.gene.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.gene.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.gene.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.gene.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.gene.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.gene.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.gene.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.gene.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.gene.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.feature` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.feature.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.feature.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.feature.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.feature.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.feature.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.feature.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.feature.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.feature.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.feature.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.feature.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.feature.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.feature.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.feature.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `legend.region` | `NULL` | Nullable R object; type/inheritance below | NULL to show and inherit common settings, or ggplot2::element_blank() to hide every guide for that role. |
| `legend.region.position` | `NULL` | Nullable R object; type/inheritance below | Per-role guide positions. |
| `legend.region.direction` | `NULL` | Nullable R object; type/inheritance below | Per-role directions. |
| `legend.region.background` | `NULL` | Nullable R object; type/inheritance below | Per-role backgrounds. |
| `legend.region.key` | `NULL` | Nullable R object; type/inheritance below | Per-role key elements. |
| `legend.region.key.size` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.region.key.width` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.region.key.height` | `NULL` | Nullable R object; type/inheritance below | Per-role key dimensions. Automatic Identity rulers measure their boundary labels independently of generic legend swatch sizes. Set a Ribbon-specific dimension or an explicit guide for a manual ruler length. |
| `legend.region.text` | `NULL` | Nullable R object; type/inheritance below | Per-role text elements. |
| `legend.region.text.position` | `NULL` | Nullable R object; type/inheritance below | Per-role text positions. |
| `legend.region.title` | `NULL` | Nullable R object; type/inheritance below | Per-role title elements. |
| `legend.region.title.position` | `NULL` | Nullable R object; type/inheritance below | Per-role title positions. |
| `legend.region.margin` | `NULL` | Nullable R object; type/inheritance below | Per-role margins. |
| `legend.region.spacing` | `NULL` | Nullable R object; type/inheritance below | Per-role spacing. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Standard ggplot2 theme arguments, including common text, plot and legend settings. Common legend settings also supply role inheritance. Removed underscore-style arguments produce a migration error. |

</details>


## Input, conversion, diagnostics and output {#workflows}

<span id="api-adapt_ggchord"></span>
<details><summary><code>adapt_ggchord()</code></summary>


[Full reference](reference/adapt_ggchord.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `mode` | `c("adaptive", "fixed", "scale")` | Vector; choices versus multiple selection below | One of adaptive, fixed, or scale. |
| `reference_size` | `c(8, 6)` | Vector; choices versus multiple selection below | Positive width and height in inches. Scale modes use the fourth root of target/reference area, with a bounded text multiplier. |
| `text_scale_range` | `c(0.85, 1.35)` | Vector; choices versus multiple selection below | Lower and upper positive text multipliers containing 1. Lines and default spacing use the square root of the text multiplier. |
| `min_text_size` | `2.5` | Numeric/integer; range in description | Lower bound in mm for automatically reduced theme text. Defaults already smaller than this value are preserved rather than enlarged. Explicit text sizes and existing label-fit policies keep their own contracts. |

</details>

<span id="api-as_ggchord_alignment"></span>
<details><summary><code>as_ggchord_alignment()</code></summary>


[Full reference](reference/as_ggchord_alignment.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `sequences` | `Required/alternative entry; see description` | R object; constraints below | Named character vector of aligned strings, or a data frame with unique accver and character sequence. All strings must have the same positive width and each record must contain a nongap symbol. |
| `alphabet` | `Required/alternative entry; see description` | R object; constraints below | Explicit "DNA", "RNA" or "AA". IUPAC nucleotide ambiguity codes and protein ambiguity/rare-residue codes are accepted; * is accepted only for AA. The alphabet is never guessed. |
| `source_data` | `NULL` | Nullable R object; type/inheritance below | Optional data frame with one row per alignment accver, start (first nongap source coordinate) and strand ("+" or "-"). Optional source_accver declares a different source ID; optional character sequence supplies the full forward source for verification. Coordinates must remain positive exact integers. Input metadata is retained in the ggchord_conversion attribute. |

</details>

<span id="api-as_ggchord_feature"></span>
<details><summary><code>as_ggchord_feature()</code></summary>


[Full reference](reference/as_ggchord_feature.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `x` | `Required/alternative entry; see description` | R object; constraints below | A data frame or GRanges object. |
| `mapping` | `NULL` | aes mapping | Optional named character vector mapping output roles to source columns, as in as_ggchord_seq(). |

</details>

<span id="api-as_ggchord_link"></span>
<details><summary><code>as_ggchord_link()</code></summary>


[Full reference](reference/as_ggchord_link.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `x` | `Required/alternative entry; see description` | R object; constraints below | A data frame or GInteractions object. |
| `mapping` | `NULL` | aes mapping | Optional named character vector mapping canonical link roles to source columns, as in as_ggchord_seq(). |

</details>

<span id="api-as_ggchord_seq"></span>
<details><summary><code>as_ggchord_seq()</code></summary>


[Full reference](reference/as_ggchord_seq.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `x` | `Required/alternative entry; see description` | R object; constraints below | A data frame, GRanges, Seqinfo, GInteractions, or a Biostrings XStringSet/XString object. |
| `mapping` | `NULL` | aes mapping | Optional named character vector output_role = "input_column" for data-frame inputs. For example c(accver = "chrom", length = "size"). |
| `lengths` | `NULL` | Nullable R object; type/inheritance below | Optional named numeric vector of known sequence lengths. Fills unknown lengths and must agree with any known source lengths. |
| `ids` | `NULL` | Character path/identifier/vector; see constraints | Optional character vector of explicit identifiers for unnamed XStringSet/XString inputs only. Existing source names must agree. |

</details>

<span id="api-build_ggchord_cds"></span>
<details><summary><code>build_ggchord_cds()</code></summary>


[Full reference](reference/build_ggchord_cds.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `cds_data` | `Required/alternative entry; see description` | R object; constraints below | Data frame with feature_id, accver, ordered one-based inclusive start/end, strand (+/-) and GFF3 phase (0–2, numeric or literal text as returned by read_gff3()). . remains unknown; unknown initial phase requires explicit codon_start for nested CDS. If type is present, select CDS rows explicitly before building a map. Alternatively supply one row per CDS with a segments list-column and explicit initial phase or codon_start (1–3). Nested genomic pieces retain explicit segment_order when supplied and are never replaced by the parent's outer envelope. |
| `genetic_code` | `NULL` | Nullable R object; type/inheritance below | Optional positive integer code ID, scalar or named by CDS feature_id. Must agree with a supplied genetic_code column. Missing IDs remain unknown; they do not prevent coordinate-only projection. |
| `partial` | `c("error", "trim")` | Vector; choices versus multiple selection below | "error" (default) or explicit "trim". |

</details>

<span id="api-bundle_ggchord_ribbons"></span>
<details><summary><code>bundle_ggchord_ribbons()</code></summary>


[Full reference](reference/bundle_ggchord_ribbons.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `ribbon_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | Alignment data in ggchord ribbon format. |
| `seq_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | Sequence data containing accver and length. |
| `bins` | `80L` | Numeric/integer; range in description | Positive integer number of normalized midpoint bins per sequence, default 80. |
| `min_bundle` | `2L` | Numeric/integer; range in description | Minimum number of rows required for aggregation, default 2. Smaller groups remain as individual ribbons. |
| `weight` | `c("length", "count", "pident")` | Vector; choices versus multiple selection below | Bundle weight: "length", "count", or "pident". Identity weight is aligned length multiplied by identity proportion. |
| `group_by` | `NULL` | Nullable R object; type/inheritance below | Optional character columns that must agree before ribbons can be bundled. |

</details>

<span id="api-clean_ggchord_data"></span>
<details><summary><code>clean_ggchord_data()</code></summary>


[Full reference](reference/clean_ggchord_data.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `seq_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | data.frame/tibble, required. Must contain accver and length; used as the coordinate reference for clipping. Reference IDs must always be non-missing, non-empty and unique, regardless of unknown_id. |
| `ribbon_data` | `NULL` | Table/supported object; data may be a function | data.frame/tibble or NULL, optional. Alignment results. Other input types are rejected rather than treated as absent. |
| `gene_data` | `NULL` | Table/supported object; data may be a function | data.frame/tibble or NULL, optional. Gene annotation data. Other input types are rejected rather than treated as absent. |
| `unknown_id` | `c("drop", "error", "keep")` | Vector; choices versus multiple selection below | Character, default "drop". Policy for rows whose sequence ID in ribbon_data or gene_data is missing, empty or unknown: "drop" removes them, "error" stops with a message, "keep" leaves them unchanged. A kept row without a reference sequence is reported, but its coordinates and annotations are not cleaned and it is not guaranteed to be drawable. IDs are never guessed or replaced. Other structural checks, including finite numeric coordinates and valid gene strands, still apply. |
| `out_of_range` | `c("clip", "drop", "error", "keep")` | Vector; choices versus multiple selection below | Character, default "clip". Policy for coordinates outside [1, sequence length]: "clip" clamps them, "drop" removes the row, "error" stops, "keep" leaves them unchanged. Coordinates are one-based inclusive, so start == end is a valid single-base interval and is retained. With "clip", intervals wholly outside their sequence are removed and reported; partially overlapping intervals retain their intersection, even a single base. |
| `reversed_interval` | `c("sort", "drop", "error", "keep")` | Vector; choices versus multiple selection below | Character, default "sort". Policy for intervals with start > end: "sort" swaps the endpoints so the feature draws stably (the original direction is recorded in the report), "drop" removes the row, "error" stops, "keep" leaves them unchanged. |
| `invalid_pident` | `c("clip", "drop", "error", "keep")` | Vector; choices versus multiple selection below | Character, default "clip". Policy for pident outside [0, 100]: "clip" clamps them, "drop" removes the row, "error" stops, "keep" leaves them unchanged. |
| `empty_annotation` | `c("keep", "drop", "replace")` | Vector; choices versus multiple selection below | Character, default "keep". Policy for missing/empty anno in gene data: "replace" fills them with replacement_annotation, "drop" removes the row, "keep" leaves them unchanged. |
| `replacement_annotation` | `"unannotated"` | Character/enum; choices below | Character, default "unannotated". Replacement annotation used when empty_annotation = "replace". |
| `feature_data` | `NULL` | Table/supported object; data may be a function | Unified annotation input; mutually exclusive with gene_data. |

</details>

<span id="api-deduplicate_ggchord_ribbons"></span>
<details><summary><code>deduplicate_ggchord_ribbons()</code></summary>


[Full reference](reference/deduplicate_ggchord_ribbons.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `ribbon_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | data.frame in ribbon_data format. |
| `tolerance` | `0` | Numeric/integer; range in description | Numeric, default 0. Maximum absolute difference (in bp) allowed on each of qstart/qend/sstart/send for two rows to count as coordinate duplicates (used with by = "coordinates"). |
| `by` | `c("exact", "coordinates", "overlap")` | Vector; choices versus multiple selection below | Character, default "exact". "exact" removes rows with identical pair and coordinates; "coordinates" additionally treats rows within tolerance bp as duplicates; "overlap" treats blocks whose reciprocal overlap is at least min_reciprocal_overlap on both sequences as duplicates. Each interval must meet the threshold, not only the shorter one. Opposite relative endpoint directions are kept as distinct relationships by all methods. Single-base endpoints use supplied strand/direction evidence; otherwise their direction is unknown. |
| `keep` | `c("best_pident", "longest", "first")` | Vector; choices versus multiple selection below | Character, default "best_pident". Which representative to keep: "best_pident" (highest pident), "longest" (longest length) or "first" (first occurrence in the input). |
| `min_reciprocal_overlap` | `0.9` | Numeric/integer; range in description | Numeric (0-1), default 0.9. Reciprocal overlap threshold used with by = "overlap". |

</details>

<span id="api-diagnose_ggchord"></span>
<details><summary><code>diagnose_ggchord()</code></summary>


[Full reference](reference/diagnose_ggchord.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `plot` | `Required/alternative entry; see description` | R object; constraints below | A ggchord plot. |
| `width` | `7` | Numeric/integer; range in description | Positive target canvas dimensions. |
| `height` | `5` | Numeric/integer; range in description | Positive target canvas dimensions. |
| `units` | `c("in", "cm", "mm", "px")` | Vector; choices versus multiple selection below | Dimensions in "in", "cm", "mm", or "px". |
| `device` | `c("png", "pdf", "svg")` | Vector; choices versus multiple selection below | Target device: "png", "pdf", or "svg". It must be available in the current R environment. Only this device is inspected. |
| `dpi` | `144` | Numeric/integer; range in description | Raster resolution; also defines pixel dimensions. |
| `min_size` | `2` | Numeric/integer; range in description | Advisory minimum text size and line width, both in millimetres. These do not modify the supplied style. |
| `min_linewidth` | `0.2` | Numeric/integer; range in description | Advisory minimum text size and line width, both in millimetres. These do not modify the supplied style. |
| `colour_vision` | `FALSE` | Logical | Request optional simulated colour-vision hints. |

</details>

<span id="api-export_ggchord_fusion"></span>
<details><summary><code>export_ggchord_fusion()</code></summary>


[Full reference](reference/export_ggchord_fusion.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `plot` | `Required/alternative entry; see description` | R object; constraints below | A Fusion object. |
| `width` | `10` | Numeric/integer; range in description | Positive canvas dimensions in inches. |
| `height` | `8` | Numeric/integer; range in description | Positive canvas dimensions in inches. |

</details>

<span id="api-export_ggchord_layout"></span>
<details><summary><code>export_ggchord_layout()</code></summary>


[Full reference](reference/export_ggchord_layout.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `plot` | `Required/alternative entry; see description` | R object; constraints below | A ggchord plot. |
| `include` | `c("seq", "ribbon", "gene", "feature", "labels", "axis", "track")` | Vector; choices versus multiple selection below | Character vector selecting any of "seq", "ribbon", "gene", "feature", "labels", "alignment", "link", "restriction", "track", and "axis". Point links and restriction sites require explicit selection. Dotplots default to alignment and axis; axis_layout contains reversible axis offsets, directions and source origins. Alignment exports retain all four endpoints, layer identity and source row; omissions appear in connection_report. |
| `original_data` | `FALSE` | Logical | Logical. Include the resolved per-layer input tables under original_data, default FALSE. |

</details>

<span id="api-filter_ggchord_ribbons"></span>
<details><summary><code>filter_ggchord_ribbons()</code></summary>


[Full reference](reference/filter_ggchord_ribbons.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `ribbon_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | data.frame with at least qaccver and saccver columns (normally the full ribbon_data format). |
| `seq_ids` | `NULL` | Nullable R object; type/inheritance below | Optional character vector. Keep only rows where both the query and the subject are in seq_ids. |
| `min_pident` | `NULL` | Nullable R object; type/inheritance below | Optional numeric. Lower/upper bounds on pident. |
| `max_pident` | `NULL` | Nullable R object; type/inheritance below | Optional numeric. Lower/upper bounds on pident. |
| `min_length` | `NULL` | Nullable R object; type/inheritance below | Optional numeric. Lower bound on length. |
| `max_evalue` | `NULL` | Nullable R object; type/inheritance below | Optional numeric. Upper bound on evalue (requires an evalue column). |
| `min_bitscore` | `NULL` | Nullable R object; type/inheritance below | Optional numeric. Lower bound on bitscore (requires a bitscore column). |
| `min_query_coverage` | `NULL` | Nullable R object; type/inheritance below | Optional numeric (0-100). Lower bound on the query coverage. Uses the qcovs column when present; otherwise computed from qlen and the query interval (requires qlen). |
| `min_subject_coverage` | `NULL` | Nullable R object; type/inheritance below | Optional numeric (0-100). Lower bound on the subject coverage, computed from slen and the subject interval (requires an slen column). |
| `keep_pairs` | `NULL` | Nullable R object; type/inheritance below | Optional data.frame/list/matrix describing an undirected set of sequence pairs. A data.frame or matrix with the first two columns used as query/subject IDs, or a list of length-2 character vectors. A row is kept when its query/subject pair matches any pair in either direction. |
| `drop_self_links` | `TRUE` | Logical | Logical, default TRUE. Remove rows where qaccver == saccver. |
| `sort_by` | `NULL` | Nullable R object; type/inheritance below | Optional character vector of column names. Prefix a name with "-" or "desc:" for descending order (e.g. c("pident", "-evalue")). |

</details>

<span id="api-filter_restriction_sites"></span>
<details><summary><code>filter_restriction_sites()</code></summary>


[Full reference](reference/filter_restriction_sites.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `sites` | `Required/alternative entry; see description` | R object; constraints below | Result from find_restriction_sites(). |
| `set` | `c("all", "unique", "unique_dual", "six_plus", "unique_6plus",      "commercial", "reference")` | Vector; choices versus multiple selection below | Display preset. "reference" reuses the enzyme-set profile stored for an exactly recognized reference sequence. It errors when no such profile is available. A saved "None" profile returns no sites. |
| `enzymes` | `NULL` | Nullable R object; type/inheritance below | Additional filters. |
| `min_site_length` | `NULL` | Nullable R object; type/inheritance below | Additional filters. |
| `cuts` | `NULL` | Nullable R object; type/inheritance below | Additional filters. |
| `window` | `NULL` | Nullable R object; type/inheritance below | Additional filters. |
| `commercial_only` | `FALSE` | Logical | Additional filters. |
| `parent_set` | `c("all", "commercial_all", "commercial_nonredundant")` | Vector; choices versus multiple selection below | Enzyme catalogue subset. "commercial_nonredundant" keeps one stable commercial representative per equivalent recognition group and site, preferring the database's preferred enzyme when present. |

</details>

<span id="api-find_common_features"></span>
<details><summary><code>find_common_features()</code></summary>


[Full reference](reference/find_common_features.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `sequence` | `Required/alternative entry; see description` | Sequence object/string/complete-sequence table; see description | DNA text, a named character vector, or a data frame with accver and sequence. |
| `database` | `NULL` | Nullable R object; type/inheritance below | Optional normalized database list or custom data frame. |
| `types` | `NULL` | Nullable R object; type/inheritance below | Optional feature-type and feature-name/ID filters. |
| `features` | `NULL` | Nullable R object; type/inheritance below | Optional feature-type and feature-name/ID filters. |
| `circular` | `TRUE` | Logical | Search across the sequence origin. |
| `mode` | `c("auto", "dna", "protein")` | Vector; choices versus multiple selection below | Matching mode. "auto" first uses curated annotations for an exactly known sequence, then follows each record's detection mode. Records marked exactProteinMatch use exact protein followed by near-exact DNA; other records use exact DNA. "protein" explicitly enables exact and approximate protein matching. |
| `min_protein_identity` | `0.9` | Numeric/integer; range in description | Approximate protein thresholds, each one finite number in [0, 1]. |
| `min_protein_coverage` | `0.8` | Numeric/integer; range in description | Approximate protein thresholds, each one finite number in [0, 1]. |
| `resolve` | `c("best", "all")` | Vector; choices versus multiple selection below | Return deterministic best annotations or every candidate. |

</details>

<span id="api-find_ggchord_loci"></span>
<details><summary><code>find_ggchord_loci()</code></summary>


[Full reference](reference/find_ggchord_loci.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `seq_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | Sequence table containing unique accver and positive integer length values. Optional genome and genome_id columns are retained in the returned windows. |
| `feature_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | Feature table containing feature_id, accver, start, end, and optionally strand. Missing strand means unknown. |
| `targets` | `Required/alternative entry; see description` | R object; constraints below | Non-empty character vector of distinct feature_id values. Selection is exact; no annotation text or accession versions are guessed. |
| `upstream` | `0` | Numeric/integer; range in description | Single non-negative integer numbers of bases added before and after each target in its biological direction. |
| `downstream` | `0` | Numeric/integer; range in description | Single non-negative integer numbers of bases added before and after each target in its biological direction. |

</details>

<span id="api-find_primer_bindings"></span>
<details><summary><code>find_primer_bindings()</code></summary>


[Full reference](reference/find_primer_bindings.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `sequence` | `Required/alternative entry; see description` | Sequence object/string/complete-sequence table; see description | DNA text, a named character vector, or a data frame with accver and sequence. |
| `primers` | `NULL` | Nullable R object; type/inheritance below | Optional custom primers with sequence and optional name and primer_id (one ID per distinct sequence), or a named character vector. |
| `set` | `c("universal", "all", "reference")` | Vector; choices versus multiple selection below | Bundled catalogue subset. reference returns authoritative primer annotations from a recognized bundled reference plasmid. |
| `match` | `c("exact", "three_prime")` | Vector; choices versus multiple selection below | Full exact matching or longest exact 3-prime annealing. |
| `min_annealed_bases` | `15L` | Numeric/integer; range in description | One finite integer of at least 8, the minimum suffix length for three_prime matching. Fractional values are rejected. |
| `circular` | `TRUE` | Logical | Allow matches to cross the sequence origin. |
| `unique_only` | `FALSE` | Logical | Keep only primer sequences with one binding site per target. |

</details>

<span id="api-find_restriction_sites"></span>
<details><summary><code>find_restriction_sites()</code></summary>


[Full reference](reference/find_restriction_sites.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `sequence` | `Required/alternative entry; see description` | Sequence object/string/complete-sequence table; see description | DNA text, a named character vector, or a data frame with accver and sequence. |
| `enzymes` | `NULL` | Nullable R object; type/inheritance below | Optional enzyme names. A data frame is a compatibility spelling for patterns. |
| `patterns` | `NULL` | Nullable R object; type/inheritance below | Optional custom definitions. |
| `circular` | `TRUE` | Logical | Allow recognition and cleavage coordinates to wrap origin. |
| `database` | `NULL` | Nullable R object; type/inheritance below | Optional parsed pattern database. |
| `methylation` | `c("auto", "none", "dam_dcm")` | Vector; choices versus multiple selection below | Methylation context. "auto" uses the audited context for recognized reference plasmids and otherwise makes no blocked-site claim; "none" disables blocking; "dam_dcm" evaluates source-backed Dam/Dcm sensitivities against the supplied sequence. |

</details>

<span id="api-focus_ggchord_data"></span>
<details><summary><code>focus_ggchord_data()</code></summary>


[Full reference](reference/focus_ggchord_data.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `seq_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | Sequence table containing accver and length. |
| `ribbon_data` | `NULL` | Table/supported object; data may be a function | Optional ggchord alignment table. |
| `gene_data` | `NULL` | Table/supported object; data may be a function | Optional gene or feature table. |
| `loci` | `Required/alternative entry; see description` | Explicit window/anchor table or documented equivalent | Data frame containing accver, start, and end, and optionally locus_id. |
| `expand` | `0` | Numeric/integer; range in description | Non-negative number of bases added on both sides, scalar or one value per locus. |
| `boundary` | `c("trim", "drop")` | Vector; choices versus multiple selection below | Boundary policy. "trim" clips partially overlapping genes and alignments; "drop" keeps only rows fully contained in selected loci. |
| `feature_data` | `NULL` | Table/supported object; data may be a function | Unified annotation input; mutually exclusive with gene_data. |

</details>

<span id="api-focus_ggchord_track"></span>
<details><summary><code>focus_ggchord_track()</code></summary>


[Full reference](reference/focus_ggchord_track.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | Quantitative input table. Metadata and mapped values are retained. |
| `seq_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | Sequence reference table. |
| `loci` | `Required/alternative entry; see description` | Explicit window/anchor table or documented equivalent | Explicit accver, start, end and optional locus_id, as in focus_ggchord_data(). |
| `type` | `c("point", "interval")` | Vector; choices versus multiple selection below | "point" requires pos; "interval" requires start/end. |
| `boundary` | `c("trim", "drop")` | Vector; choices versus multiple selection below | "trim" (default) or "drop" for partly retained intervals. |

</details>

<span id="api-fusion_panel"></span>
<details><summary><code>fusion_panel()</code></summary>


[Full reference](reference/fusion_panel.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `plot` | `Required/alternative entry; see description` | R object; constraints below | An independent ggchord genomic plot, another Fusion, or a Cartesian detail from register_ggchord_alignment(). Dotplot, facets and unregistered external ggplots are not yet registered here. Character panels require explicit cross-panel fill; backbone-colour gradient is unavailable. |
| `id` | `Required/alternative entry; see description` | R object; constraints below | Unique nonempty panel name; "center" is reserved. |
| `side` | `c("right", "top", "left", "bottom")` | Vector; choices versus multiple selection below | One of right, top, left, bottom, topright, topleft, bottomright or bottomleft; one panel per named side. Placement never rotates the child Coord. |
| `size` | `NULL` | Nullable R object; type/inheritance below | Optional fractional width and height. Out-of-frame allocations fail; overlap fails unless explicitly allowed by the composition. |
| `box` | `NULL` | Nullable R object; type/inheritance below | Optional explicit rectangle (xmin, xmax, ymin, ymax), optionally named in any order. Numeric values are fractions of the area excluding the composition legend; absolute grid units retain their physical size. Use instead of side/size/angle/distance. Out-of-frame boxes fail at build. |
| `angle` | `NULL` | Nullable R object; type/inheritance below | Optional finite placement angle in degrees counterclockwise from the right, instead of side. Does not rotate the plot or its labels. |
| `distance` | `NULL` | Nullable R object; type/inheritance below | Positive radial distance from the centre allocation, as a fraction of half the shorter physical allocation dimension; default .8. Angles describe physical canvas directions even on rectangular devices. The allocation must fit without implicit clamping. |

</details>

<span id="api-fusion_ribbon"></span>
<details><summary><code>fusion_ribbon()</code></summary>


[Full reference](reference/fusion_ribbon.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | A relationship data frame, preserving original endpoints. |
| `mapping` | `NULL` | aes mapping | Standard ggplot2 aesthetic mapping: fill, colour, alpha, linewidth and linetype. With pident and no explicit fill, Identity is mapped by default with 0/50/80/90/95/100 breaks. Composition uses ordinary scale_fill_* etc.; panel role scales remain local to each panel. |
| `qpanel` | `NULL` | Nullable R object; type/inheritance below | Optional fixed panel IDs. Do not also supply their columns. |
| `spanel` | `NULL` | Nullable R object; type/inheritance below | Optional fixed panel IDs. Do not also supply their columns. |
| `gradient` | `FALSE` | Logical | Opt in to query-to-subject sequence-colour interpolation. Cannot be combined with an explicit fill or fill mapping. |
| `tip` | `c("none", "subject")` | Vector; choices versus multiple selection below | "none" or a decorative "subject" tip. Source fronts remain exported in full even when the tip is drawn. |
| `allow_repeated` | `FALSE` | Logical | Explicitly allow repeated source_hsp_id display. |
| `curvature` | `0.3` | Numeric/integer; range in description | Control distance relative to endpoint separation, in [0,1]. |
| `fill` | `NULL` | Nullable R object; type/inheritance below | Fixed styles; explicit values take priority over mappings. The default fill without pident is steelblue. |
| `colour` | `NA` | R object; constraints below | Fixed styles; explicit values take priority over mappings. The default fill without pident is steelblue. |
| `alpha` | `0.35` | Numeric/integer; range in description | Fixed styles; explicit values take priority over mappings. The default fill without pident is steelblue. |
| `linewidth` | `0.25` | Numeric/integer; range in description | Fixed styles; explicit values take priority over mappings. The default fill without pident is steelblue. |
| `linetype` | `1` | Numeric/integer; range in description | Fixed styles; explicit values take priority over mappings. The default fill without pident is steelblue. |
| `show.legend` | `NA` | R object; constraints below | Standard ggplot2 legend selection: NA follows mappings; TRUE/FALSE or a named logical vector controls guide participation. |

</details>

<span id="api-get_chord_layout"></span>
<details><summary><code>get_chord_layout()</code></summary>


[Full reference](reference/get_chord_layout.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `plot` | `Required/alternative entry; see description` | R object; constraints below | A ggchord plot. |
| `build` | `TRUE` | Logical | Logical. Build the layout when it is not cached, default TRUE. |

</details>

<span id="api-ggchord"></span>
<details><summary><code>ggchord()</code></summary>


[Full reference](reference/ggchord.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `seq_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | data.frame/tibble with accver and length, or a GRanges, Seqinfo or GInteractions object with known full sequence lengths, or a named nucleotide or protein Biostrings XStringSet. Unnamed sequence objects must first use as_ggchord_seq() with explicit ids. Optional S4 inputs use as_ggchord_seq(). The sequence_unit column is "bp" (default when omitted) or "aa". All rows use the same unit. Protein lengths come from complete proteins. Native aa plots accept links and generic features; gene, primer and restriction layers require bp. Arbitrary character widths are not biological sequence lengths. |
| `ribbon_data` | `NULL` | Table/supported object; data may be a function | Optional alignment table with qaccver, saccver, qstart, qend, sstart and send. length and pident are optional unless requested by a mapping, filter or statistic. A GInteractions object uses as_ggchord_link() without inventing identity or alignment length. |
| `gene_data` | `NULL` | Table/supported object; data may be a function | Compatible legacy annotation input. The third positional argument retains its meaning. Use feature_data in new code; a non-NULL legacy input emits one message at construction, never during rendering. Untyped legacy tables retain their directional gene semantics. No removal version is scheduled. Do not supply both annotation inputs. |
| `debug` | `FALSE` | Logical | Logical. Whether to output debug information, default FALSE |
| `validate` | `c("warn", "error", "none")` | Vector; choices versus multiple selection below | Character, default "warn". How to run the structured input-data validation (see validate_ggchord_data): "warn" emits a single summary warning when the data has problems and caches the full report on the plot object (p$ggchord$validation); "error" stops on severe problems; "none" skips the diagnostic validation (the cheap structural checks that prevent crashes are still performed). |
| `feature_data` | `NULL` | Table/supported object; data may be a function | Optional unified annotation table or GRanges. Requires accver, start, end; absent strand becomes unknown ("."). Generic feature layers use the full table. Default gene layers select gene/CDS from feature_class, type, or role (in that priority order), and require known strands. tRNA/rRNA remain generic features. Explicit layer data overrides the constructor input. Use disjoint subsets when combining gene and feature layers; each layer keeps its own data independently. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Reserved for clear migration errors from removed constructor arguments. Use labs(), coord_chord() and theme() for presentation. |

</details>

<span id="api-ggchord_fusion"></span>
<details><summary><code>ggchord_fusion()</code></summary>


[Full reference](reference/ggchord_fusion.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `center` | `Required/alternative entry; see description` | R object; constraints below | A ggchord genomic plot (Circular, Chord, Collinear, Sankey or Radial), or another Fusion. Nested biological endpoints use full panel paths. |
| `center_size` | `NULL` | Nullable R object; type/inheritance below | Optional width and height as fractions of the composition area, excluding its cross-panel legend. NULL allocates available space. |
| `layout` | `c("fixed", "adaptive")` | Vector; choices versus multiple selection below | "fixed" uses fixed panel allocations; "adaptive" also adjusts the central allocation to the physical aspect ratio. Both rebuild labels and guides at actual panel sizes. Neither rescales explicit fonts. |
| `allow_overlap` | `FALSE` | Logical | Explicitly allow container overlap. Default FALSE fails on overlap. When TRUE, later panels paint over earlier panels. Connections use disjoint native grid clips to avoid accumulating translucent fills. This is deliberate placement, not global collision avoidance. |
| `center_box` | `NULL` | Nullable R object; type/inheritance below | Optional centre allocation: xmin, xmax, ymin, ymax as four fractions of the composition area excluding its legend, or four absolute grid units. Use instead of center_size. See fusion_panel(). |

</details>

<span id="api-link_ggchord_clusters"></span>
<details><summary><code>link_ggchord_clusters()</code></summary>


[Full reference](reference/link_ggchord_clusters.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `feature_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | Feature data frame with feature_id, accver, start, end, and optional strand, accepted by as_ggchord_feature(). |
| `membership` | `Required/alternative entry; see description` | R object; constraints below | Data frame with feature_id and cluster_id. A feature can belong to one cluster; extra columns are retained in the report. |
| `pairs` | `Required/alternative entry; see description` | R object; constraints below | Data frame with explicit qaccver and saccver comparisons. |
| `paralogs` | `c("drop", "all")` | Vector; choices versus multiple selection below | How to handle multiple members per cluster and sequence: "drop" (default) or explicit "all". |
| `segment_policy` | `c("error", "span")` | Vector; choices versus multiple selection below | "error" (default) rejects repeated feature IDs; "span" explicitly uses a validated discontinuous feature's outer bounds. |
| `max_links` | `10000` | Numeric/integer; range in description | Positive integer upper limit for generated rows, default 10000. Exceeding it stops before Cartesian products are allocated. |

</details>

<span id="api-magnify_ggchord"></span>
<details><summary><code>magnify_ggchord()</code></summary>


[Full reference](reference/magnify_ggchord.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `plot` | `Required/alternative entry; see description` | R object; constraints below | A registered ggchord genomic source plot or nested Fusion. |
| `region` | `Required/alternative entry; see description` | R object; constraints below | Nonempty data frame with accver, start, end in original one-based inclusive source coordinates. Bounds must be ordered positive integers. Optional instance and detail_instance select exact display instances in the source and detail panels for every row. Ambiguity fails; origin-wrapping intervals must be split before using this interface. With a list of details, a detail_panel column selects the named target for each row. Nested leaves use full paths, for example detail/center; optional source_panel similarly selects a source leaf under center. Ambiguous container selection fails. |
| `detail_plot` | `Required/alternative entry; see description` | R object; constraints below | An explicit registered ggchord genomic plot, nested Fusion, a character plot from register_ggchord_alignment(), or a uniquely named list of such plots. Every panel must have at least one source region; the name center is reserved. The selected source interval must be fully visible in both panels, in the same bp/aa unit. Character plots require explicit source/column registration; ordinary unregistered plots are rejected. |
| `side` | `c("right", "top", "left", "bottom")` | Vector; choices versus multiple selection below | A named side from fusion_panel(), including diagonal sides. For a list of plots, supply a character vector named by panel ID, with distinct sides. Multiple panels require explicit routing and placement. |
| `center_size` | `NULL` | Nullable R object; type/inheritance below | Optional panel width/height fractions, following ggchord_fusion() and fusion_panel(). Defaults reserve a smaller source panel and a wider detail container for a single detail plot. With multiple details, use a named list of sizes or NULL for Fusion defaults; center_size also follows Fusion defaults. Out-of-frame allocations fail; overlap requires explicit allow_overlap. |
| `detail_size` | `NULL` | Nullable R object; type/inheritance below | Optional panel width/height fractions, following ggchord_fusion() and fusion_panel(). Defaults reserve a smaller source panel and a wider detail container for a single detail plot. With multiple details, use a named list of sizes or NULL for Fusion defaults; center_size also follows Fusion defaults. Out-of-frame allocations fail; overlap requires explicit allow_overlap. |
| `connector` | `c("leader", "none")` | Vector; choices versus multiple selection below | "leader" draws two decorative boundary leaders per row; "none" retains the highlights only. These are not biological ribbons or alignments. |
| `colour` | `"#B45309"` | Character/enum; choices below | Fixed highlight/leader colour. |
| `alpha` | `0.2` | Numeric/integer; range in description | Highlight opacity in [0,1]. Leaders use the opaque colour. |
| `highlight_width` | `1.5` | Numeric/integer; range in description | Positive highlight band width in physical mm. |
| `linewidth` | `0.35` | Numeric/integer; range in description | Nonnegative leader width using the usual ggplot2 convention. |
| `linetype` | `2` | Numeric/integer; range in description | Leader line type, default dashed. |
| `detail_box` | `NULL` | Nullable R object; type/inheritance below | Explicit detail allocation instead of side/detail_size, using the box convention in fusion_panel(). For multiple plots, supply a named list covering every detail. Physical boxes are rechecked at target size. |
| `allow_overlap` | `FALSE` | Logical | Deliberately permit overlapping containers, following ggchord_fusion(). Later panels cover earlier ones; default FALSE. |
| `center_box` | `NULL` | Nullable R object; type/inheritance below | Explicit source allocation instead of center_size. |

</details>

<span id="api-map_alignment_position"></span>
<details><summary><code>map_alignment_position()</code></summary>


[Full reference](reference/map_alignment_position.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | Long alignment frame from as_ggchord_alignment(). |
| `accver` | `Required/alternative entry; see description` | R object; constraints below | One alignment-row ID, even if its source_accver differs. |
| `position` | `Required/alternative entry; see description` | Position object or documented string | Numeric positions, positive integral or NA. |
| `inverse` | `FALSE` | Logical | FALSE maps alignment columns to source positions; TRUE maps source positions to alignment columns. |

</details>

<span id="api-map_ggchord_position"></span>
<details><summary><code>map_ggchord_position()</code></summary>


[Full reference](reference/map_ggchord_position.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `plot` | `Required/alternative entry; see description` | R object; constraints below | A ggchord plot using coord_collinear(), coord_radial_synteny() or coord_dotplot(). |
| `accver` | `Required/alternative entry; see description` | R object; constraints below | Sequence identifiers, scalar or one per position. |
| `position` | `Required/alternative entry; see description` | Position object or documented string | Numeric genomic positions, or layout x when inverse = TRUE. Radial inverse mapping instead requires a matrix/data frame with x and y columns on the specified backbone. Off-backbone positions are rejected. |
| `inverse` | `FALSE` | Logical | Convert layout x back to native sequence positions. |
| `source` | `FALSE` | Logical | Use original source coordinates after focusing. |
| `axis` | `NULL` | Nullable R object; type/inheritance below | Dotplot axis to map, either "x" (query) or "y" (subject). Required for Dotplots; omitted for linear and radial layouts. Dotplot gaps return NA on inverse mapping, with the chosen axis recorded in the result. |

</details>

<span id="api-merge_ggchord_ribbons"></span>
<details><summary><code>merge_ggchord_ribbons()</code></summary>


[Full reference](reference/merge_ggchord_ribbons.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `ribbon_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | data.frame in ribbon_data format. |
| `max_gap` | `0` | Numeric/integer; range in description | Numeric, default 0. Maximum gap (in bp) allowed between two blocks on both sequences for them to be merged. |
| `min_pident_difference` | `0` | Numeric/integer; range in description | Numeric, default 0. When > 0, two blocks are only merged when their pident difference is <= this value. |
| `require_same_orientation` | `TRUE` | Logical | Logical, default TRUE. Only merge blocks with the same relative query/subject direction and compatible progression of their bounds. Reversing both endpoints preserves relative direction. FALSE explicitly permits a coarser mixed-direction summary. |
| `group_by` | `c("qaccver", "saccver")` | Vector; choices versus multiple selection below | Character vector, default c("qaccver", "saccver"). Additional grouping columns; query/subject IDs always remain grouping keys. |
| `extra_columns` | `c("consistent", "first")` | Vector; choices versus multiple selection below | How to handle columns other than the required ribbon fields when multiple rows are merged. "consistent" (default) keeps a value only when all source rows agree and otherwise stores a typed missing value; "first" retains the first source row's value. |

</details>

<span id="api-optimize_ggchord_layout"></span>
<details><summary><code>optimize_ggchord_layout()</code></summary>


[Full reference](reference/optimize_ggchord_layout.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `seq_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | Sequence data containing accver and length. |
| `ribbon_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | Alignment data in ggchord ribbon format. |
| `seq_order` | `NULL` | Geometry number/vector/documented formats | Optional initial complete sequence order. The default uses the row order of seq_data. |
| `seq_orientation` | `NULL` | Geometry number/vector/documented formats | Optional initial orientation (1 or -1), scalar, complete vector, or named vector. |
| `optimize` | `c("order", "orientation")` | Vector; choices versus multiple selection below | Any of "order" and "orientation". |
| `weight` | `c("length", "count", "pident")` | Vector; choices versus multiple selection below | Objective weight: "length", "count", or "pident". |
| `max_iter` | `50L` | Numeric/integer; range in description | Positive integer maximum number of improving search rounds. |

</details>

<span id="api-prepare_ggchord_data"></span>
<details><summary><code>prepare_ggchord_data()</code></summary>


[Full reference](reference/prepare_ggchord_data.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `seq_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | Sequence input accepted by as_ggchord_seq(). |
| `ribbon_data` | `NULL` | Table/supported object; data may be a function | Optional paired intervals accepted by as_ggchord_link(). |
| `gene_data` | `NULL` | Table/supported object; data may be a function | Optional directional gene input accepted by as_ggchord_feature(). Nondirectional intervals belong in features. |
| `features` | `list()` | Named list | Named lists of feature or paired-interval inputs. |
| `links` | `list()` | Named list | Named lists of feature or paired-interval inputs. |
| `mapping` | `list()` | aes mapping | Named list with optional seq, ribbon, gene column mappings, and features/links lists of mappings named like the respective collections. Each column mapping is a named character vector. |
| `steps` | `list()` | Named list | Named list of explicit helper argument lists. Supported names: clean, filter, deduplicate, merge, bundle, optimize. Data arguments are supplied by preparation and cannot be overridden. |
| `feature_data` | `NULL` | Table/supported object; data may be a function | Unified annotation table, mutually exclusive with gene_data. |

</details>

<span id="api-primer_catalog"></span>
<details><summary><code>primer_catalog()</code></summary>


[Full reference](reference/primer_catalog.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `set` | `c("all", "universal")` | Vector; choices versus multiple selection below | Return all records or only the commonly used universal subset. |
| `aliases` | `FALSE` | Logical | Return source alias records instead of the sequence catalogue. |

</details>

<span id="api-project_ggchord_features"></span>
<details><summary><code>project_ggchord_features()</code></summary>


[Full reference](reference/project_ggchord_features.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `feature_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | Domain/feature data frame with feature_id, parent_id (or a single-valued parent_ids list-column), start and end. |
| `cds` | `Required/alternative entry; see description` | R object; constraints below | A map from build_ggchord_cds(). |
| `from` | `c("protein", "genome")` | Vector; choices versus multiple selection below | "protein" (default) or "genome". |
| `codon_boundary` | `c("error", "expand")` | Vector; choices versus multiple selection below | Inverse projection policy: "error" or "expand". |

</details>

<span id="api-project_ggchord_links"></span>
<details><summary><code>project_ggchord_links()</code></summary>


[Full reference](reference/project_ggchord_links.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `link_data` | `Required/alternative entry; see description` | R object; constraints below | Data frame with qfeature_id, sfeature_id identifying CDS models and ordered inclusive protein qstart, qend, sstart, send. Query and subject residue spans must be equal. If present, gapopen must be zero and source length must equal this residue span. |
| `cds` | `Required/alternative entry; see description` | R object; constraints below | A map from build_ggchord_cds(). |
| `ungapped` | `FALSE` | Logical | Must explicitly be TRUE, confirming an ungapped residue correspondence rather than inferring it from equal endpoint spans. |

</details>

<span id="api-read_aligned_fasta"></span>
<details><summary><code>read_aligned_fasta()</code></summary>


[Full reference](reference/read_aligned_fasta.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `file` | `Required/alternative entry; see description` | Character path/identifier/vector; see constraints | A single local aligned FASTA path; no multi-alignment concatenation. |
| `alphabet` | `Required/alternative entry; see description` | R object; constraints below | Explicit "DNA", "RNA" or "AA". IUPAC nucleotide ambiguity codes and protein ambiguity/rare-residue codes are accepted; * is accepted only for AA. The alphabet is never guessed. |
| `source_data` | `NULL` | Nullable R object; type/inheritance below | Optional data frame with one row per alignment accver, start (first nongap source coordinate) and strand ("+" or "-"). Optional source_accver declares a different source ID; optional character sequence supplies the full forward source for verification. Coordinates must remain positive exact integers. Input metadata is retained in the ggchord_conversion attribute. |
| `source_file` | `FALSE` | Logical | Logical. Also retain .source_file on each character row. |

</details>

<span id="api-read_bed"></span>
<details><summary><code>read_bed()</code></summary>


[Full reference](reference/read_bed.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `file` | `NULL` | Character path/identifier/vector; see constraints | Optional path vector or glob patterns; read before files. |
| `files` | `NULL` | Character path/identifier/vector; see constraints | Optional paths or wildcard patterns for multiple BED files. |
| `source_file` | `FALSE` | Logical | Logical. Add a .source_file column, default FALSE. |

</details>

<span id="api-read_bedpe"></span>
<details><summary><code>read_bedpe()</code></summary>


[Full reference](reference/read_bedpe.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `file` | `NULL` | Character path/identifier/vector; see constraints | Optional local path vector or glob patterns; read before files. |
| `files` | `NULL` | Character path/identifier/vector; see constraints | Optional additional local paths or glob patterns. |
| `source_file` | `FALSE` | Logical | Logical. Add a .source_file column, default FALSE. |
| `unmapped` | `c("error", "keep")` | Vector; choices versus multiple selection below | "error" (default) or "keep" for incomplete anchors. |

</details>

<span id="api-read_blast"></span>
<details><summary><code>read_blast()</code></summary>


[Full reference](reference/read_blast.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `file` | `NULL` | Character path/identifier/vector; see constraints | Optional path vector or glob patterns; read before files. |
| `files` | `NULL` | Character path/identifier/vector; see constraints | Optional character vector of BLAST tabular output files (literal paths and/or wildcard patterns). All matched files are read and combined into one data.frame. |
| `format` | `c("auto", "outfmt6", "outfmt7", "custom")` | Vector; choices versus multiple selection below | Character. "auto" (default) detects the column layout from the file; "outfmt6" requires the standard 12-column layout, "outfmt7" parses its # Fields: declaration, and "custom" requires col_names. |
| `col_names` | `NULL` | Nullable R object; type/inheritance below | Optional character vector naming the columns in the file, used with format = "custom" or to override auto-detection. |
| `comment` | `"#"` | Character/enum; choices below | Character comment character, default "#" (BLAST outfmt 7 header lines start with #). |
| `source_file` | `FALSE` | Logical | Logical. Add a .source_file column when reading one or more files, default FALSE. |
| `...` | `Required/alternative entry; see description` | Forwarded arguments | Additional arguments passed to utils::read.table() (e.g. na.strings). |

</details>

<span id="api-read_fasta_lengths"></span>
<details><summary><code>read_fasta_lengths()</code></summary>


[Full reference](reference/read_fasta_lengths.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `file` | `NULL` | Character path/identifier/vector; see constraints | Optional path vector or glob patterns; read before files. |
| `files` | `NULL` | Character path/identifier/vector; see constraints | Optional character vector of FASTA files (literal paths and/or wildcard patterns). All matched files are read and combined. |
| `header_delim` | `NULL` | Nullable R object; type/inheritance below | Optional character. When given, each header is split at every occurrence of this delimiter and only the first piece is kept. |
| `source_file` | `FALSE` | Logical | Logical. Add a .source_file column when reading one or more files, default FALSE. |
| `sequence_unit` | `NULL` | Nullable R object; type/inheritance below | Optional "bp" or "aa" declaration. When supplied, adds that column to the result; the alphabet is not inferred. Omitted units retain the existing output columns and ggchord interprets them as bp. |
| `keep_sequence` | `FALSE` | Logical | Logical, default FALSE. Also retain the actual FASTA characters in sequence, removing whitespace only and preserving case. No alphabet, reverse complement or alignment columns are inferred. Gap characters, if present, are counted in length; gapped input needs an explicit alignment-to-source map before genomic use. |

</details>

<span id="api-read_genbank"></span>
<details><summary><code>read_genbank()</code></summary>


[Full reference](reference/read_genbank.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `file` | `NULL` | Character path/identifier/vector; see constraints | Optional path vector or glob patterns; read before files. |
| `files` | `NULL` | Character path/identifier/vector; see constraints | Optional paths or wildcard patterns for several files. |
| `feature_types` | `NULL` | Nullable R object; type/inheritance below | Optional feature types to select. NULL retains all. |
| `lengths` | `NULL` | Nullable R object; type/inheritance below | Optional named numeric vector of known full sequence lengths, following as_ggchord_seq(). |
| `python` | `"python3"` | Character path/identifier/vector; see constraints | Interpreter name or path, default "python3". It must already have the optional Biopython package installed. No environment is modified. |
| `unsupported` | `c("error", "report")` | Vector; choices versus multiple selection below | "error" or "report" for locations without a safe local interval representation. |
| `source_file` | `FALSE` | Logical | Add .source_file columns to output tables, default FALSE. File provenance is retained in the report regardless. |

</details>

<span id="api-read_gff3"></span>
<details><summary><code>read_gff3()</code></summary>


[Full reference](reference/read_gff3.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `file` | `NULL` | Character path/identifier/vector; see constraints | Optional path vector or glob patterns; read before files. |
| `files` | `NULL` | Character path/identifier/vector; see constraints | Optional character vector of GFF3 files (literal paths and/or wildcard patterns). All matched files are read and combined. |
| `feature_types` | `"CDS"` | Character/enum; choices below | Character vector of feature types to keep, default "CDS". |
| `anno_from` | `c("product", "Name", "gene", "ID")` | Vector; choices versus multiple selection below | Character vector of GFF3 attribute keys, tried in order to fill the anno column. |
| `unstranded` | `c("plus", "drop", "keep")` | Vector; choices versus multiple selection below | "plus" (the default) assigns unknown strands to +; "drop" removes them. Use "keep" with geom_feature() to retain them as nondirectional . intervals. source_strand always records the input. |
| `source_file` | `FALSE` | Logical | Logical. Add a .source_file column when reading one or more files, default FALSE. |

</details>

<span id="api-read_paf"></span>
<details><summary><code>read_paf()</code></summary>


[Full reference](reference/read_paf.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `file` | `NULL` | Character path/identifier/vector; see constraints | Optional path vector or glob patterns; read before files. |
| `files` | `NULL` | Character path/identifier/vector; see constraints | Optional paths or wildcard patterns for multiple PAF files. |
| `source_file` | `FALSE` | Logical | Logical. Add a .source_file column, default FALSE. |

</details>

<span id="api-register_ggchord_alignment"></span>
<details><summary><code>register_ggchord_alignment()</code></summary>


[Full reference](reference/register_ggchord_alignment.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `plot` | `Required/alternative entry; see description` | R object; constraints below | A single-panel Cartesian ggplot with geom_residue layers mapping x = alignment_column, y = accver, label = symbol. Identity/reverse x scales and a discrete y scale are supported. Additional independent annotations are allowed. Facets, nonlinear scales, flipped coordinates and position adjustments moving registered cells are rejected. |
| `data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | The displayed rows of as_ggchord_alignment(), read_aligned_fasta() or translate_ggchord_cds()$nucleotides/residues. One cell per accver/alignment_column; source coordinates must be unique within each row. Keep original columns when filtering a window. Every selected source base must be present; split spliced/origin-wrapping regions explicitly. Gap columns remain visible but have no source coordinate. |

</details>

<span id="api-segment_ggchord_data"></span>
<details><summary><code>segment_ggchord_data()</code></summary>


[Full reference](reference/segment_ggchord_data.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `seq_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | Full sequence table with unique accver and integer lengths. |
| `ribbon_data` | `NULL` | Table/supported object; data may be a function | Optional paired interval table. |
| `gene_data` | `NULL` | Table/supported object; data may be a function | Optional gene or generic feature table, including segments. |
| `windows` | `Required/alternative entry; see description` | R object; constraints below | Table with accver, start, end, unique locus_id, optional row_id (defaults to source accver) and orientation (1 or -1, defaults to 1). Bounds must be ordered, integral and within the sequence. |
| `point_data` | `NULL` | Nullable R object; type/inheritance below | Optional point links with qaccver, saccver, qpos, spos. |
| `cross_gap` | `c("split", "drop")` | Vector; choices versus multiple selection below | How to handle interval links clipped by selected windows. |
| `max_links` | `10000` | Numeric/integer; range in description | Positive integer upper bound on candidate ribbon/point copies. Increase explicitly for larger, intentionally overlapping views. |
| `feature_data` | `NULL` | Table/supported object; data may be a function | Unified annotation input; mutually exclusive with gene_data. |

</details>

<span id="api-translate_ggchord_cds"></span>
<details><summary><code>translate_ggchord_cds()</code></summary>


[Full reference](reference/translate_ggchord_cds.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `cds` | `Required/alternative entry; see description` | R object; constraints below | A map from build_ggchord_cds(). |
| `sequences` | `Required/alternative entry; see description` | R object; constraints below | Named character vector of complete forward source strings, or a data frame with unique accver and sequence (e.g. read_fasta_lengths() with keep_sequence = TRUE). DNAStringSet/RNAStringSet are also accepted. Every CDS source must exist and contain its declared intervals. |
| `alphabet` | `c("DNA", "RNA")` | Vector; choices versus multiple selection below | Explicit "DNA" (default) or "RNA"; ungapped IUPAC only. |
| `initiation` | `c("literal", "start")` | Vector; choices versus multiple selection below | "literal" (default) translates every codon literally. "start" uses the genetic code's alternative initiation codons for the first codon only; it requires a complete biological 5-prime end. |
| `ambiguous` | `c("error", "solve")` | Vector; choices versus multiple selection below | "error" (default) rejects ambiguous codons; "solve" uses Biostrings to resolve them, emitting X when not uniquely resolvable. |

</details>

<span id="api-validate_ggchord_data"></span>
<details><summary><code>validate_ggchord_data()</code></summary>


[Full reference](reference/validate_ggchord_data.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `seq_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | data.frame/tibble, required. Basic sequence information (columns accver, length). |
| `ribbon_data` | `NULL` | Table/supported object; data may be a function | data.frame/tibble, optional. Alignment results (columns qaccver, saccver, qstart, qend, sstart, send). |
| `gene_data` | `NULL` | Table/supported object; data may be a function | data.frame/tibble, optional. Gene annotation data (columns accver, start, end, strand; anno is optional). |
| `strict` | `FALSE` | Logical | Logical. When TRUE, collect issues and then stop if any severe problem is found. When FALSE (default), return the full diagnostic report without stopping. |
| `check_coordinates` | `TRUE` | Logical | Logical, default TRUE. Whether to check that ribbon/gene coordinates stay inside [1, sequence length]. |
| `check_duplicates` | `TRUE` | Logical | Logical, default TRUE. Whether to look for fully duplicated, near-duplicated and highly overlapping records. Reciprocal overlap requires the intersection to cover at least 90 percent of each interval, including descending ribbon endpoints. Reversed feature bounds may cross a circular origin and are excluded from linear overlap diagnostics; exact duplicates are still reported. |
| `check_self_links` | `TRUE` | Logical | Logical, default TRUE. Whether to flag alignment rows where qaccver == saccver. |
| `feature_data` | `NULL` | Table/supported object; data may be a function | Unified annotation table, mutually exclusive with gene_data. Requires accver, start, end. Missing strand is unknown; +, -, ., and +/- are accepted. Reports use table "feature" and data_summary$n_feature; the input is never modified. |

</details>

<span id="api-validate_ggchord_hierarchy"></span>
<details><summary><code>validate_ggchord_hierarchy()</code></summary>


[Full reference](reference/validate_ggchord_hierarchy.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `seq_data` | `Required/alternative entry; see description` | Table/supported object; data may be a function | A data frame containing accver and genome foreign keys. |
| `genomes` | `Required/alternative entry; see description` | R object; constraints below | A data frame containing a unique genome primary key. |
| `strict` | `FALSE` | Logical | Logical. Stop if relational errors are present; warnings about genomes without contigs do not cause failure. |

</details>

<span id="api-view_ggchord"></span>
<details><summary><code>view_ggchord()</code></summary>


[Full reference](reference/view_ggchord.html)


| Parameter | Signature default | Type/value rule | Meaning, range and constraints |
| --- | --- | --- | --- |
| `plot` | `ggplot2::last_plot()` | R object; constraints below | A ggchord or ggplot object, default last_plot(). |
| `width` | `12.39` | Numeric/integer; range in description | Positive output width. With no explicit dimensions, coord_chord() uses 12.39 inches while coord_circular() derives its width from the annotation envelope. |
| `height` | `9.71` | Numeric/integer; range in description | Positive output height. With no explicit dimensions, coord_chord() uses 9.71 inches while coord_circular() derives both dimensions. Use NULL to fit the sequence, labels and legends while preserving equal coordinate units. Explicit width/height values are always respected. |
| `units` | `c("in", "cm", "mm", "px")` | Vector; choices versus multiple selection below | Output units: "in", "cm", "mm", or "px". |
| `device` | `c("png", "svg")` | Vector; choices versus multiple selection below | Preview device, "png" or "svg". SVG output requires the optional svglite package or a working base Cairo SVG device. |
| `dpi` | `144` | Numeric/integer; range in description | Positive raster resolution, default 144. |
| `bg` | `NULL` | Nullable R object; type/inheritance below | Optional background colour passed to ggsave(). |
| `viewer` | `c("auto", "ide", "browser", "none")` | Vector; choices versus multiple selection below | How to open the preview: "auto" prefers getOption("viewer") and falls back to a browser in interactive sessions; "ide", "browser", and "none" select an explicit destination. "none" is useful in scripts and tests. |

</details>


## Upstream reexports and objects {#upstream}

These retain upstream ggplot2/grid semantics rather than defining new ggchord layers. After_stat and after_scale differ by stage even if identical compares their R objects equal.

<span id="api-aes"></span>
<details><summary><code>aes</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |
| `x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `y` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `...` | `Required/alternative entry` | Forwarded arguments; Follow upstream help for types, units and accepted values. |

</details>

<span id="api-after_scale"></span>
<details><summary><code>after_scale</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |
| `x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |

</details>

<span id="api-after_stat"></span>
<details><summary><code>after_stat</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |
| `x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |

</details>

<span id="api-annotate"></span>
<details><summary><code>annotate</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |
| `geom` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `x` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `y` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `xmin` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `xmax` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `ymin` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `ymax` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `xend` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `yend` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `...` | `Required/alternative entry` | Forwarded arguments; Follow upstream help for types, units and accepted values. |
| `na.rm` | `FALSE` | Logical; Follow upstream help for types, units and accepted values. |

</details>

<span id="api-arrow"></span>
<details><summary><code>arrow</code> · grid</summary>

Upstream entry reference: [Rd](reference/grid-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |
| `angle` | `30` | Numeric/integer; range in description; Follow upstream help for types, units and accepted values. |
| `length` | `unit(0.25, "inches")` | grid physical-unit object; Follow upstream help for types, units and accepted values. |
| `ends` | `"last"` | Character/enum; choices below; Follow upstream help for types, units and accepted values. |
| `type` | `"open"` | Character/enum; choices below; Follow upstream help for types, units and accepted values. |

</details>

<span id="api-element_blank"></span>
<details><summary><code>element_blank</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |

</details>

<span id="api-element_line"></span>
<details><summary><code>element_line</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |
| `colour` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `linewidth` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `linetype` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `lineend` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `color` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `linejoin` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `arrow` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `arrow.fill` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `inherit.blank` | `FALSE` | Logical; Follow upstream help for types, units and accepted values. |
| `size` | `deprecated()` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `...` | `Required/alternative entry` | Forwarded arguments; Follow upstream help for types, units and accepted values. |

</details>

<span id="api-element_rect"></span>
<details><summary><code>element_rect</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |
| `fill` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `colour` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `linewidth` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `linetype` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `color` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `linejoin` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `inherit.blank` | `FALSE` | Logical; Follow upstream help for types, units and accepted values. |
| `size` | `deprecated()` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `...` | `Required/alternative entry` | Forwarded arguments; Follow upstream help for types, units and accepted values. |

</details>

<span id="api-element_text"></span>
<details><summary><code>element_text</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |
| `family` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `face` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `colour` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `size` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `hjust` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `vjust` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `angle` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `lineheight` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `color` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `margin` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `debug` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `inherit.blank` | `FALSE` | Logical; Follow upstream help for types, units and accepted values. |
| `...` | `Required/alternative entry` | Forwarded arguments; Follow upstream help for types, units and accepted values. |

</details>

<span id="api-expansion"></span>
<details><summary><code>expansion</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |
| `mult` | `0` | Numeric/integer; range in description; Follow upstream help for types, units and accepted values. |
| `add` | `0` | Numeric/integer; range in description; Follow upstream help for types, units and accepted values. |

</details>

<span id="api-ggsave"></span>
<details><summary><code>ggsave</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |
| `filename` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `plot` | `get_last_plot()` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `device` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `path` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `scale` | `1` | Numeric/integer; range in description; Follow upstream help for types, units and accepted values. |
| `width` | `NA` | Geometry number/vector/documented formats; Follow upstream help for types, units and accepted values. |
| `height` | `NA` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `units` | `c("in", "cm", "mm", "px")` | Vector; choices versus multiple selection below; Follow upstream help for types, units and accepted values. |
| `dpi` | `300` | Numeric/integer; range in description; Follow upstream help for types, units and accepted values. |
| `limitsize` | `TRUE` | Logical; Follow upstream help for types, units and accepted values. |
| `bg` | `NULL` | Nullable R object; type/inheritance below; Follow upstream help for types, units and accepted values. |
| `create.dir` | `FALSE` | Logical; Follow upstream help for types, units and accepted values. |
| `...` | `Required/alternative entry` | Forwarded arguments; Follow upstream help for types, units and accepted values. |

</details>

<span id="api-ggtitle"></span>
<details><summary><code>ggtitle</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |
| `label` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `subtitle` | `waiver()` | R object; constraints below; Follow upstream help for types, units and accepted values. |

</details>

<span id="api-guides"></span>
<details><summary><code>guides</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |
| `...` | `Required/alternative entry` | Forwarded arguments; Follow upstream help for types, units and accepted values. |

</details>

<span id="api-labs"></span>
<details><summary><code>labs</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |
| `...` | `Required/alternative entry` | Forwarded arguments; Follow upstream help for types, units and accepted values. |
| `title` | `waiver()` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `subtitle` | `waiver()` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `caption` | `waiver()` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `tag` | `waiver()` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `dictionary` | `waiver()` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `alt` | `waiver()` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `alt_insight` | `waiver()` | R object; constraints below; Follow upstream help for types, units and accepted values. |

</details>

<span id="api-last_plot"></span>
<details><summary><code>last_plot</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |

</details>

<span id="api-margin"></span>
<details><summary><code>margin</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |
| `t` | `0` | Numeric/integer; range in description; Follow upstream help for types, units and accepted values. |
| `r` | `0` | Numeric/integer; range in description; Follow upstream help for types, units and accepted values. |
| `b` | `0` | Numeric/integer; range in description; Follow upstream help for types, units and accepted values. |
| `l` | `0` | Numeric/integer; range in description; Follow upstream help for types, units and accepted values. |
| `unit` | `"pt"` | Character/enum; choices below; Follow upstream help for types, units and accepted values. |
| `...` | `Required/alternative entry` | Forwarded arguments; Follow upstream help for types, units and accepted values. |

</details>

<span id="api-position_identity"></span>
<details><summary><code>position_identity</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |

</details>

<span id="api-PositionIdentity"></span>
<details><summary><code>PositionIdentity</code> · object</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |

</details>

<span id="api-rel"></span>
<details><summary><code>rel</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |
| `x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |

</details>

<span id="api-theme"></span>
<details><summary><code>theme</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |
| `...` | `Required/alternative entry` | Forwarded arguments; Follow upstream help for types, units and accepted values. |
| `line` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `rect` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `text` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `title` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `point` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `polygon` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `geom` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `spacing` | `Required/alternative entry` | Geometry number/vector/documented formats; Follow upstream help for types, units and accepted values. |
| `margins` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `aspect.ratio` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.title` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.title.x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.title.x.top` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.title.x.bottom` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.title.y` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.title.y.left` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.title.y.right` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.text` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.text.x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.text.x.top` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.text.x.bottom` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.text.y` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.text.y.left` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.text.y.right` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.text.theta` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.text.r` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks.x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks.x.top` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks.x.bottom` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks.y` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks.y.left` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks.y.right` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks.theta` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks.r` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.minor.ticks.x.top` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.minor.ticks.x.bottom` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.minor.ticks.y.left` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.minor.ticks.y.right` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.minor.ticks.theta` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.minor.ticks.r` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks.length` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks.length.x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks.length.x.top` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks.length.x.bottom` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks.length.y` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks.length.y.left` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks.length.y.right` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks.length.theta` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.ticks.length.r` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.minor.ticks.length` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.minor.ticks.length.x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.minor.ticks.length.x.top` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.minor.ticks.length.x.bottom` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.minor.ticks.length.y` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.minor.ticks.length.y.left` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.minor.ticks.length.y.right` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.minor.ticks.length.theta` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.minor.ticks.length.r` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.line` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.line.x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.line.x.top` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.line.x.bottom` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.line.y` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.line.y.left` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.line.y.right` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.line.theta` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `axis.line.r` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.background` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.margin` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.spacing` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.spacing.x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.spacing.y` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.key` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.key.size` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.key.height` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.key.width` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.key.spacing` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.key.spacing.x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.key.spacing.y` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.key.justification` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.frame` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.ticks` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.ticks.length` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.axis.line` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.text` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.text.position` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.title` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.title.position` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.position` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.position.inside` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.direction` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.byrow` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.justification` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.justification.top` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.justification.bottom` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.justification.left` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.justification.right` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.justification.inside` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.location` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.box` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.box.just` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.box.margin` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.box.background` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `legend.box.spacing` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `panel.background` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `panel.border` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `panel.spacing` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `panel.spacing.x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `panel.spacing.y` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `panel.grid` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `panel.grid.major` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `panel.grid.minor` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `panel.grid.major.x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `panel.grid.major.y` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `panel.grid.minor.x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `panel.grid.minor.y` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `panel.ontop` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `panel.widths` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `panel.heights` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `plot.background` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `plot.title` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `plot.title.position` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `plot.subtitle` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `plot.caption` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `plot.caption.position` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `plot.tag` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `plot.tag.position` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `plot.tag.location` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `plot.margin` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `strip.background` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `strip.background.x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `strip.background.y` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `strip.clip` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `strip.placement` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `strip.text` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `strip.text.x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `strip.text.x.bottom` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `strip.text.x.top` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `strip.text.y` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `strip.text.y.left` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `strip.text.y.right` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `strip.switch.pad.grid` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `strip.switch.pad.wrap` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `complete` | `FALSE` | Logical; Follow upstream help for types, units and accepted values. |
| `validate` | `TRUE` | Logical; Follow upstream help for types, units and accepted values. |

</details>

<span id="api-unit"></span>
<details><summary><code>unit</code> · grid</summary>

Upstream entry reference: [Rd](reference/grid-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |
| `x` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `units` | `Required/alternative entry` | R object; constraints below; Follow upstream help for types, units and accepted values. |
| `data` | `NULL` | Table/supported object; data may be a function; Follow upstream help for types, units and accepted values. |

</details>

<span id="api-waiver"></span>
<details><summary><code>waiver</code> · ggplot2</summary>

Upstream entry reference: [Rd](reference/ggplot2-helpers.html)

| Parameter | Signature default | Type/contract |
| --- | --- | --- |

</details>
