<a id="shared-band"></a>

## Overlay bars and a line on one scale {#shared-band}

The installed bedGraph contains GC percentages in successive pUC19 windows. Bars span the observed intervals and the line joins their midpoints. Declare the physical band once with `track_band("gc")`. Both geometry layers select it with `track = "gc"`. A height of 60 has the same meaning in both layers.

{{example:gc-shared-band}}

<a id="domain"></a>

## Separate physical width, value and colour {#domain}

Use the band's `width` for space on the page. Use `scale_track_value_continuous()` with explicit `limits` for the numerical axis. Changing the width leaves each GC percentage intact. A tile with `aes(fill = value)` represents a value by colour across an interval. For a second measure with different units, declare another band and give it its own value scale.

To compute coverage from intervals and compare groups, see [window statistics and facets](tracks-statistics.md).
