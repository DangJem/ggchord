# Maintaining the ggchord documentation

[English](README.md) · [现代汉语](README-Hans.md)

The website documents the public API in this checkout. Check the package version and documentation channel when comparing a local build with the deployed site. [Package metadata](../DESCRIPTION) · [Documentation channel](version.json)

## Update an article or a navigation label

Edit article prose in `website/en/`. Edit English titles, navigation labels, captions and other structured text in `website/content-en/`. A localizable field looks like `{"en": "Get started"}`. The content policy lists the English source, the saved Han translation in `website/translations/zh/` and the generated output for each item. The renderer reads Han Markdown from `website/zh/` and uses the paired JSON for structured text. [Content policy](content-policy.json)

After changing the English, prepare translation drafts from the repository root.

```sh
python3 tools/sync-website-languages.py --prepare
```

Review each `needs_review` entry against the current English. The draft records the English text and its SHA-256 digest. Update the Han text, set its status to `reviewed`, then generate and verify the outputs.

```sh
python3 tools/sync-website-languages.py
python3 tools/sync-website-languages.py --check
```

The check reports missing translations, stale reviews and direct edits to generated outputs. Both editions share code, directives, figure identifiers, navigation identifiers, link destinations and explicit section anchors. Match heading levels and table columns when translating.

## Import complete reviewed drafts

For an existing article, use the English source path from the policy and the path to the reviewed Han Markdown. To add several pages, place the reviewed Han drafts in one directory with the same filenames as their English sources. The directory command validates every supplied pair before registering new pages and writing their outputs.

```sh
python3 tools/sync-website-languages.py --reviewed-page en/installation.md /tmp/installation-reviewed.md
python3 tools/sync-website-languages.py --reviewed-directory /tmp/reviewed-han
```

For structured content, edit the English JSON first. The reviewed draft must have the same structure, English strings and shared metadata, with both `en` and `zh` in each localized field. Register new structured files in the content policy before importing them. Run the following command to record a reviewed homepage draft.

```sh
python3 tools/sync-website-languages.py --reviewed-content home.json /tmp/home-reviewed.json
```

Article links point to actual repository files, such as a sibling `.md` page, a help topic in `../../man/` or an image in `../../man/figures/`. The renderer converts these links to website URLs. When retiring an article, remove its policy entry, source, translation and generated file together. Update incoming links and record destinations for its old URL and section anchors.

## Update the API reference

The API reference documents the public inventory and categories for the package source in this checkout. Edit English page prose in `website/reference-manual/en.json`, then review the corresponding Han text in `website/reference-manual/zh.json`. Record the SHA-256 of the reviewed English source in the Han file. The shared `structure.json` holds the invocation, parameter and role keys, examples and cross-links. Both language editions expose the same code and anchors. Native facet pages link to the upstream ggplot2 documentation.

Give each multi-function family one article with shared usage and argument sections. Use its prose to explain when to choose each function. In `structure.json`, `featured_examples` selects the examples shown first. Other member scripts remain directly linkable in the same example section. Include only the inputs and calls needed to demonstrate the result. Declare any shared input blocks in the example’s `example_fixtures` in setup order. The displayed code, downloads and rendered results use that same script.

Run these commands to check and build the API reference independently.

```sh
python3 tools/check-documentation.py --reference-only
python3 tools/build-website.py --reference-only
```

The local entry is `.site/en/reference/index.html`, with the matching Han entry under `.site/zh/reference/`. The API reference uses the same pkgdown template, navigation and appearance controls as the handbook. The builder runs the downloadable R scripts, adds their plots and console output to both language editions, and refreshes `.site/assets/ggchord-source.zip` from the package sources in the same checkout. The example runner uses downlit, ragg, png and the optional reader dependencies rtracklayer, VariantAnnotation, Biopython and MUMmer’s `show-coords`. Set `GGCHORD_PYTHON` or `GGCHORD_SHOW_COORDS` if those tools are installed outside the default search path.

## Update function help, package READMEs or offline tutorials

- **Native R help.** Edit the English roxygen source and regenerate its Rd. Review the corresponding API reference topic when the public contract changes.
- **Package READMEs.** Both editions share layout and code in `website/readme-template.md`. Edit `website/content-en/readme-text.json`, review the Han translation and regenerate the paired JSON. Run `python3 tools/sync-tutorials.py` to update the two READMEs, then run it with `--check` to verify them.
- **Offline tutorials.** Edit `vignettes/ggchord_vignette.Rmd` and update the Han edition in `vignettes/ggchord_guide_hans.Rmd`. After review, record the SHA-256 of the English file in the Han header. The README sync command also checks the tutorial digest, R chunks, section anchors and local asset references.
- **Offline figures.** Run `python3 tools/prepare-vignette-figures.py` to prepare the images embedded in the standalone tutorial HTML distributed by CRAN. This maintenance step uses Pillow. During release checks, `tools/check-vignette-resources.py` verifies the local resources referenced by the tutorial HTML in the source archive.

## Build and inspect the local site

Design the website for desktop computers. On article and reference pages, the chapter list and the outline sit together on the left. Keep tables and code within the article and allow wide values to scroll inside their own block. Use source, build and static HTML checks for website validation. Check the generated navigation, stylesheets, fragment links and example inputs from the build output.

The core build requires Python 3.9 or newer, R, Pandoc, the package dependencies and the dependencies of the examples. The site uses `pkgdown`, Bootstrap 5, `pkgload`, `rmarkdown` and `jsonlite`. The full build renders both languages, replays the handbook and API examples, and checks the output before updating the local preview. Run the following commands from the repository root.

```sh
python3 tools/sync-website-languages.py --check
python3 tools/check-documentation.py --source-only
python3 tools/build-website.py --render-jobs 2
python3 tools/check-documentation.py
python3 -m http.server 8765 --directory .site --bind 127.0.0.1
```

Open `http://127.0.0.1:8765/` to inspect the files in `.site/`. Both editions use the same images for the four groups of homepage previews. The case scripts under `website/examples/` supply the linked full figures. The homepage previews use the same cases and are regenerated with `Rscript tools/render-site-flagships.R`. For changes to a figure in the website, package README or offline tutorial, keep a record of the input data, filtering rules, plotting code and image source. [Artwork, third-party assets and licenses](ASSET-SOURCES.md)

## CI and deployment

Pushes and pull requests that affect documentation inputs run source checks. A manual GitHub Actions run can produce a preview artifact. Deployment requires explicit authorization. The `--publication` flag controls the build label. The historical site is archived in `docs/`. GitHub Releases records package releases and GitHub Actions records website deployments.
