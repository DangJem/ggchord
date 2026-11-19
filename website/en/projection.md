<a id="cds"></a>

## Locate a residue span on the genome {#cds}

To locate residues 2–5 of pUC19 lacZα, select its annotated CDS and keep the source strand. This complete CDS starts in phase 0 and uses bacterial genetic code 11. The selected residue span is an illustration. It projects to bases 455–466 on the reverse strand.

{{example:project-residues}}

<a id="segments"></a>

## Handle joined CDS segments and aligned pairs {#segments}

For a split CDS, build the model from all its original segments in their annotated reading frame. A projected span can include multiple genomic pieces. Inspect the returned segments to distinguish the coding bases from bases between exons.

`project_ggchord_links()` projects both sides of an established correspondence. A gapped protein HSP gives outer boundaries but may not establish a residue-to-residue match for every column. Use the available alignment evidence when choosing which coordinates to project. If the intended figure is a protein comparison, provide complete protein lengths and retain native aa coordinates.
