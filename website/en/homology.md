<div class="task-intro"><p>Compare AmpR neighbourhoods in two plasmids. Check real coding DNA, establish correspondence, then focus to each CDS plus 200 bp on either side.</p></div>

## Extract actual CDS and coding DNA {#membership}
Select AmpR annotations from bundled pUC19 and pETDuet-1, retaining IDs, coordinates, strand and text. Reverse-complement negative-strand DNA so comparisons use the same biological orientation rather than annotation names alone.

{{code:homologinput}}

{{output:homologinput}}

Read the negative-strand CDS in coding orientation, then compare the two DNA strings.

{{code:homologdna}}
{{output:homologdna}}

`cluster_features` should contain two actual CDS; `cds_dna` retains their oriented coding sequences. This strict example checks exact DNA equality. General homology analysis can use external alignments and clustering; homologs need not be identical.

## Declare members and sequence pairs {#links}
Build membership from exact coding-DNA equality. Different DNA does not enter one group merely because both annotations say AmpR.

{{code:homologmembership}}

{{output:homologmembership}}

Inspect `membership` and `homologs$report$comparisons` and `$members`. The linking utility consumes existing groups; it does not compute a new alignment. Unknown members and ambiguous paralogs are reported; only explicit `paralogs = "all"` requests all combinations within the documented bounds.

## Focus on comparable neighbourhoods {#focus}
Full plasmid lengths would shrink the individual CDS. Construct windows from each real CDS source interval and add 200 bp on either side, retaining the original sequence and annotation tables.

{{code:homologs}}
{{figure:homologs}}

Blue arrows show the reverse-strand AmpR CDS and the grey band connects the verified intervals. Axes use each original sequence's positions, so one horizontal location can have different source labels. Break markers indicate that only local windows remain.

`local_homologs$report` records retention and clipping. Genes outside the windows are omitted, not absent. The band represents declared interval correspondence, not a base alignment; envelopes of segmented annotations likewise do not establish exon-level correspondence.

## Project supported protein intervals to both CDS {#projection}
{{code:linkprojection}}

{{console:linkprojection}}
{{results:projection}}

The first 20 residues of identical CDS provide sequence evidence for `ungapped = TRUE`. Equal spans alone do not justify that claim for ordinary BLASTp HSPs. The utility retains protein statistics rather than presenting them as nucleotide identity.

## Check your result {#check-result}
Both AmpR arrows show the reverse strand. Coding DNA supports correspondence; source axes and focus reports explain the local display. Save at [explicit 10 × 5 in dimensions](output.html).

Next: [read real BLASTp query intervals](alignments.html#protein), understand [CDS projection and segments](projection.html), or create [independent local windows](views.html).

## Save this result {#save-result}
Save the plot just created to the working directory; use matching dimensions for preview.

{{code:courseblocke}}
