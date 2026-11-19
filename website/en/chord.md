A chord diagram wraps the sequence coordinates around a circle. Follow an arc to read positions along one reference. Follow a ribbon across the circle to find the corresponding interval on another reference. Each example below runs in a fresh R session with ggchord and dplyr installed.

<a id="backbones"></a>

## Draw the two references {#backbones}

The first two phage genomes in the bundled sequence table are MT108731.1 (64,323 bp) and MT118296.1 (32,090 bp). Select those accessions so the figure contains exactly the sequences you will compare. The longer genome takes up about twice as much of the circle. At this stage the empty centre tells you nothing about how well the genomes match.

{{example:first-backbones}}

<a id="draw"></a>

## Join one pair of aligned intervals {#draw}

The first alignment for this pair is a local match, also called a high-scoring segment pair (HSP). It joins positions 26,298–51,139 bp on the query MT108731.1 to 7,121–31,959 bp on the subject MT118296.1. The two ends of the ribbon span those intervals on their respective arcs. The rest of each genome stays visible.

{{example:first-chord}}

<a id="features"></a>

## Locate the coding sequences {#features}

Keep the same references and the same ribbon. Add the selected coding sequences from the bundled annotation table to see which features lie inside or outside each ribbon end. Arrow direction comes from the original annotation in `source_strand`. The first arrow on MT108731.1 spans 13,064–14,554 bp on the plus strand, before the aligned interval begins.

{{example:first-features}}

<a id="adapt"></a>

## Draw your own alignments {#adapt}

Replace the reference table with one row per sequence and its full length. Give the ribbon layer the original endpoints on both references, retaining their order for reverse alignments. The [data guide](data.md#alignments) describes all six required columns and the optional identity field. To read the same pair left to right, continue to the [collinear guide](collinear.md). For a complete figure with multiple genomes, see the [four-genome chord map](case-chord.md).

[Sequence layer](../../man/geom_seq.Rd) · [Ribbon layer](../../man/geom_link_ribbon.Rd) · [Chord coordinates](../../man/coord_chord.Rd)
