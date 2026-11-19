Use a chord diagram to see which sequence intervals align across several references. The backbones follow arcs. Each end of a ribbon occupies the interval recorded on that sequence.

<a id="draw"></a>

## One alignment on four backbones {#draw}

The first row of `ribbon_data_example` connects two of the four references. Its query spans 26,298–51,139 bp on MT108731.1. Its subject spans 7,121–31,959 bp on MT118296.1. The two other arcs provide context.

{{example:first-chord}}

<a id="adapt"></a>

## Add your own alignments {#adapt}

Supply `accver` and complete `length` for each backbone. Each ribbon needs `qaccver`, `saccver` and the four endpoints `qstart`, `qend`, `sstart` and `send`. Keep descending endpoints when the alignment reports them. To compare identities across many ribbons, retain the `pident` column in the link table.

[Four-genome chord map](case-chord.md) · [Ribbon layer](../../man/geom_link_ribbon.Rd) · [Chord coordinates](../../man/coord_chord.Rd)
