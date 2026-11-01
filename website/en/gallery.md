<p class="lesson-intro">Choose an example by the question it answers. Thumbnails help you navigate; click for full size. Every figure links to preparation, code and interpretation.</p>

## Find a figure by task {#tasks}
| Overall relationships | One annotated sequence | Local positions | Quantities and text |
| --- | --- | --- | --- |
| [Chord](#realchord) | [Plasmid](#plasmid) | [CDS neighbourhood](#homologs) | [GC tracks](#multitrack) |
| [Point links](#pointlinks) | [Feature stacking](#featurestack) | [Focus](#focus) | [Label layouts](#labelauto) |

## One source for both languages {#source}
All figures use bundled source-backed data or attributed real input records, defined in the shared
[documentation example source](../assets/documentation-examples.R).
No raw FASTA/GFF3/BLAST record was changed to improve appearance. The same image
bytes are used by both languages, README and vignette; captions provide localized
interpretation. The source and asset hashes are recorded in the [asset manifest](../assets/plots.csv).

Display filters are stated in the relevant articles. GC is computed from actual
DNA; exact-CDS membership is checked from sequence content. Deliberate corruption
is confined to diagnostic examples and is never presented as biological data.

## Figures and dimensions {#figures}
{{gallery}}

## Reproduce {#reproduce}
From the source repository, with its suggested development dependencies available:

```sh
Rscript tools/documentation-plots.R
python3 tools/sync-tutorials.py --check
python3 tools/build-website.py
python3 tools/check-documentation.py
```

The PNG device is used only to update these documentation figures. This is not
an output-device matrix or release validation. Source code and dimensions live
in `tools/documentation-plots.R`; changing the figure parameters requires
regenerating the manifest. See [output](output.html) for your own export.

## Asset sources {#assets}
[Artwork, template and dependency licenses](../assets/ASSET-SOURCES.md).
