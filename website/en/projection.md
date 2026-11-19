<a id="prepare"></a>

## Choose the input {#prepare}

Start with a CDS table and a known amino-acid interval. A complete CDS supplies source coordinates, strand and reading frame.

<a id="cds"></a>

## Locate an amino-acid interval on the nucleotide axis {#cds}

The bundled pUC19 annotation supplies the complete lacZα CDS. Residues 2–5 map to 455–466 bp on the original negative strand.

```r
record <- flagship_circular_example$references[["pUC19"]]
coding <- record$features |>
  dplyr::filter(anno == "lacZα", type == "CDS") |>
  dplyr::transmute(accver, start, end, strand, feature_id = "lacZ_alpha", phase = 0L)
cds <- build_ggchord_cds(coding, genetic_code = 11)
span <- data.frame(feature_id = "display", parent_id = "lacZ_alpha", start = 2, end = 5)
mapped <- project_ggchord_features(span, cds)
mapped[c("accver", "start", "end", "strand")]
```

<a id="boundaries"></a>

## Segment boundaries and phase {#boundaries}

Several CDS segments can share one parent ID. Projection orders them in transcription direction. The returned segments column preserves intron gaps that the outer interval spans.

<a id="links"></a>

## Project a known pair of intervals {#links}

Each schematic axis below contains ten codons. The example specifies a one-to-one correspondence between three residues and shows endpoint order on both strands. Real data require their own evidence for that correspondence.

```r
models <- data.frame(
  feature_id = c("A", "B"), accver = c("row1", "row2"),
  start = c(101, 301), end = c(130, 330), strand = c("+", "-"), phase = 0L
)
pair <- data.frame(qfeature_id = "A", sfeature_id = "B",
  qstart = 2, qend = 4, sstart = 5, send = 7)
project_ggchord_links(pair, build_ggchord_cds(models), ungapped = TRUE)
```

<a id="check-result"></a>

## Check the projection {#check-result}

Check parent CDS IDs, native units and the full segment list. For a gapped BLASTp hit, determine whether the evidence supports a single query span or a residue correspondence at both ends.
