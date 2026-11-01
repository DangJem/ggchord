## One source for both languages {#source}
All six figures use **synthetic teaching inputs**, defined in the shared
[documentation example source](../assets/documentation-examples.R).
No raw FASTA/GFF3/BLAST record was changed to improve appearance. The same image
bytes are used by both languages, README and vignette; captions provide localized
interpretation. The source and asset hashes are recorded in the [asset manifest](../assets/plots.csv).

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
