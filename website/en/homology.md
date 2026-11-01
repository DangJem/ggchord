<div class="task-intro"><p>Compare AmpR neighbourhoods in two plasmids. Check real coding DNA, establish correspondence, then focus to each CDS plus 200 bp on either side.</p></div>

## Prepare this page {#prepare}
Only bundled data is needed; no download or new alignment. Load the package and optional display settings first.

{{prepare:docstyle}}

## 1. Extract actual CDS and coding DNA {#membership}
Select AmpR annotations from bundled pUC19 and pETDuet-1, retaining IDs, coordinates, strand and text. Reverse-complement negative-strand DNA so comparisons use the same biological orientation rather than annotation names alone.

{{code:homologinput}}

`cluster_features` should contain two actual CDS; `cds_dna` retains their oriented coding sequences. This strict example checks exact DNA equality. General homology analysis can use external alignments and clustering; homologs need not be identical.

## 2. Declare members and sequence pairs {#links}
Build membership from exact coding-DNA equality. Different DNA does not enter one group merely because both annotations say AmpR.

{{code:homologmembership}}

Inspect `membership` and `homologs$report$comparisons` and `$members`. The linking utility consumes existing groups; it does not compute a new alignment. Unknown members and ambiguous paralogs are reported; only explicit `paralogs = "all"` requests all combinations within the documented bounds.

## 3. Focus on comparable neighbourhoods {#focus}
Full plasmid lengths would shrink the individual CDS. Construct windows from each real CDS source interval and add 200 bp on either side, retaining the original sequence and annotation tables.

{{code:homologs}}
{{figure:homologs}}

Gold arrows identify the checked identical coding sequences; the purple band connects their CDS ranges and internal text names AmpR. Axes use each original sequence's positions, so one horizontal location can have different source labels. Break markers indicate that only local windows remain.

`local_homologs$report` records retention and clipping. Genes outside the windows are omitted, not absent. The band represents declared interval correspondence, not a base alignment; envelopes of segmented annotations likewise do not establish exon-level correspondence.

## 4. Project supported protein intervals to both CDS {#projection}
{{code:linkprojection}}
{{results:projection}}

The first 20 residues of identical CDS provide sequence evidence for `ungapped = TRUE`. Equal spans alone do not justify that claim for ordinary BLASTp HSPs. The utility retains protein statistics rather than presenting them as nucleotide identity.

## Check your result {#check-result}
Both internal labels should read AmpR and arrows retain the actual negative strand. Coding DNA supports correspondence; source axes and focus reports explain the local display. Save at [explicit 10 × 5 in dimensions](output.html).

Next: [read real BLASTp query intervals](alignments.html#protein), understand [CDS projection and segments](projection.html), or create [independent local windows](views.html).

## Save this result {#save-result}
Save the plot just created to the working directory; use matching dimensions for preview.

```r
ggsave("ampR-neighbourhood.png", p_homologs, width = 10, height = 5, units = "in", dpi = 144)
```
