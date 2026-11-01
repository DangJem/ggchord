<div class="home-intro">
<img class="package-logo" src="../assets/logo.webp" alt="ggchord" width="150" height="150">
<p class="eyebrow">GENOMIC VISUALIZATION · R / GGPLOT2</p>
<p class="home-lead">Make sequence relationships visible.</p>
<p>Compose <strong>chord comparisons, circular maps and collinear plots</strong> with ggplot2 layers. Turn existing sequences, annotations and alignments into interpretable, reproducible genomic figures.</p>
<div class="home-actions"><a class="primary-action" href="installation.html">Draw your first plot →</a><a href="articles.html">Find a workflow</a><a href="reference/index.html">Look up a function</a></div>
</div>

## Start with the question {#choose}
| What do you want to understand? | Choose | Start here |
| --- | --- | --- |
| Which aligned regions connect several sequences? | A chord plot for the overall relationship pattern | [Real sequence comparisons](real-data.html) |
| Which annotations lie on one circular sequence? | A circular map retaining its full length and direction | [The pUC19 workflow](plasmids.html) |
| How do gene positions and directions compare locally? | A focused collinear plot bringing corresponding intervals together | [Compare real CDS](homology.html) |

A circular layout does not suit every question. The [plot selection guide](modes.html) explains inputs, interpretation and limits.

## Chord · See the relationship pattern {#chord}
{{figure:realchord}}

Existing alignments among four phage sequences. Categorical colours identify sequences; a continuous scale encodes reported `pident`. Outer arrows are explicitly selected CDS. Read endpoints, then identity; crossings alone do not establish rearrangements.
[Import the original FASTA, GFF3 and BLAST records and reproduce this figure →](real-data.html#read)

## Circular · Read one sequence's annotations {#circular}
{{figure:plasmid}}

CDS, replication origin and promoters on pUC19 have distinct category colours. Primers and restriction sites have their own lines and labels. Three enzymes are deliberately displayed; all matches remain in the source result.
[From bundled DNA to a complete plasmid map →](plasmids.html)

## Collinear · Compare local positions {#collinear}
{{figure:homologs}}

AmpR CDS with identical coding DNA in pUC19 and pETDuet-1, shown with 200 bp flanks. Axes retain source positions; matching arrows identify the sequence-supported CDS relationship. The band connects the declared intervals; it is not a newly computed alignment.
[Check CDS, establish membership and focus the neighbourhood →](homology.html)

## Continue from where you are {#learn}
<div class="route-grid">
<a class="route-card" href="installation.html"><strong>New to ggchord</strong><span>Install → minimal inputs → first plot → interpret → save</span></a>
<a class="route-card" href="formats.html"><strong>Bring your own data</strong><span>Choose a reader; check IDs, lengths and coordinate units</span></a>
<a class="route-card" href="styling.html"><strong>Improve an existing plot</strong><span>Colours, guides, labels, local views and export dimensions</span></a>
<a class="route-card" href="reference/index.html"><strong>Find an exact parameter</strong><span>Complete signatures, arguments, return values and examples</span></a>
</div>

## Begin with bundled data {#quick-start}
Draw the four bundled sequence backbones to check IDs and complete lengths:

```r
library(ggchord)
data(seq_data_example)
sequences <- seq_data_example[c("accver", "length")]
ggchord(seq_data = sequences) +
  geom_seq() +
  geom_seq_label() +
  coord_chord()
```

[Get started](installation.html#input) provides independent preparation code: just two sequence fields, four gene fields and six link fields enter the plot. You can draw before learning every layer or handling historical metadata.

Updating older code? Read the [migration guide](migration.html). These pages describe local 0.16.5 development; the formal release remains 0.15.0. See the [version and installation notes](installation.html#source).
