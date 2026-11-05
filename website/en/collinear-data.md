The [basic input tables](chord-data.html) also apply here. This page focuses on display orientation and relationships between rows.

## Read the three tables as rows {#tables}

One sequence row gives one backbone, and each link row describes two intervals. Annotations retain the source `strand` when a Position separates them above and below the backbone. Display gaps in a local window do not measure missing DNA.

## Read an observed alignment {#example}

{{code:frontlinktable}}

{{output:frontlinktable}}

Both intervals in the first row increase. The query covers 26,298–51,139 bp and the subject covers 7,121–31,959 bp. Locate them on the first two rows of the [collinear example](collinear-build.html#links). Preserve decreasing subject endpoints when present, and configure display flips separately.

## Obtain the inputs {#import}

Use FASTA for full lengths, GFF3 for gene intervals, and BLAST or PAF for the two aligned intervals. IDs must match exactly. Gene intervals are one-based and inclusive, while decreasing link endpoints may encode reverse alignments. With your own files, follow the [source-file workflow](real-data.html#read) through to the [linear plot](real-data.html#linear). Use the [format guide](formats.html) for other inputs and [CDS projection](projection.html) for protein coordinates.

## Current conversion boundaries {#contract}

You can set backbone display direction manually or let the coordinate system choose it from alignment evidence. Choose one approach for a plot. Display flips retain the source strand. Adjacency follows the displayed row order. Enable non-adjacent links explicitly when you want to include them.

When reusing layers across layouts, keep circular radius, curvature, ring and radial sequence-label controls with their applicable coordinates.

Continue to [cross-row links](collinear-routing.html).
