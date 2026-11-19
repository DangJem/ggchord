<a id="api"></a>

## API migration {#api}

Since v0.50.0, sequence layout, feature lanes, tracks and guides use the revised API. Update existing scripts one layer at a time. Use `layout_sequences()` for sequence order and orientation, `position_feature_stack()` for feature lanes and `track_band()` to share a band across quantitative layers.

[Full migration table and replacement calls](../../MIGRATION.md) · [Browse functions by purpose](../../man/)

<a id="pages"></a>

## Previous documentation addresses {#pages}

Bookmarks for plotting entries, cases and figure anchors remain usable. Browse the [plotting guide](modes.md) by figure or [data and operations](formats.md) by task.

## Removed layers and positions {#removed-layers}

- Replace `geom_gene()` and `geom_feature_plasmid()` with `geom_feature(preset = "gene")` and `geom_feature(preset = "plasmid")` respectively. Pass the same layer data. [Feature presets](../../man/geom_feature.Rd#arg-preset)
- Replace `geom_feature_label_repel()` with `geom_feature_label(placement = "repel")`. Set the placement policy with `label_layout`, control long labels with `overflow`, and select the feature layer with `anchor_layer` when needed. [Feature labels](../../man/geom_feature_label.Rd)
- Replace `position_plasmid(offset)` with `position_offset(offset)`. For lanes, use `position_feature_stack()` to set the starting position and direction. [Offsets](../../man/position_offset.Rd) · [Feature lanes](../../man/position_feature_stack.Rd)
- Declare the physical band with `track_band("depth", width = w)` and pass `track = "depth"` on each quantitative layer. The former `position_track()` placed the band. The layer’s `position` now adjusts its data by stacking, dodging or nudging. [Track bands](../../man/track_band.Rd) · [Quantitative layers](../../man/geom_track.Rd)
