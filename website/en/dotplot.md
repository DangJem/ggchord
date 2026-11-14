## HSP relationships on two sequence axes {#task}

This matrix uses all 32 original HSP records across four complete sequence inputs. Both axes preserve accession identity and local positions. Identity controls colour; observed HSP length controls stroke width.

{{manualfigure:flagdotplot}}

## What the module supports {#path}

- Multiple query and subject sequences, explicit axis order and native-unit gaps.
- Forward and reverse HSP endpoints, independent display orientation on either axis.
- Identity and other supplied metrics mapped through link scales.
- Exact source-to-display conversion, visible windows and traceable input rows.

[Data, axes and mapping](dotplot-start.html) · [Alignment imports](alignments.html) · [Protein coordinates](protein-alignments.html)

## Required inputs {#inputs}

Supply complete lengths and query/subject IDs, plus qstart/qend/sstart/send for each HSP. Axis concatenation does not fuse biological sequence identities. bp and aa need separate coordinate domains or an explicit projection.

## Boundaries {#support}

Each mark is an HSP endpoint segment, not a base-by-base reconstruction. Mismatches and gaps need observed aligned characters or an exact alignment trace. Backbone gene/feature/track layers do not apply here. Dotplot is not a registered Fusion endpoint.
