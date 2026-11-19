Start with the [reference and layer tables](data.md). A pair of sequence lengths and one observed interval pair are enough to compare two genomes. Use the bundled records to [draw your first chord diagram](chord.md#backbones). The guide begins with two labelled backbones, adds a single alignment and then locates the coding sequences on those same arcs.

<a id="reading-order"></a>

## Read the same alignment in several views {#reading-order}

1. Begin with the [chord diagram](chord.md#backbones). Read the length of each reference and find the two ends of the ribbon.
2. Move the same records onto [parallel backbones](collinear.md). Then compare a [radial arrangement](radial-synteny.md) and a [two-column Sankey view](sankey.md). The same source positions can be easier to read in a different arrangement.
3. Follow the question into a [single circular map](circular.md), a [two-axis dot plot](dotplot.md), or an [overview with a detail panel](fusion.md). These use different inputs or panel arrangements.

<a id="modules"></a>

## Choose a view for your question {#modules}

Six coordinate systems arrange data inside a panel. Fusion combines independently drawn panels and connects their corresponding intervals.

| Question | Guide |
| --- | --- |
| Which intervals match among sequences arranged around a circle? | [Chord diagram](chord.md) |
| Where does an alignment start and end along two parallel sequences? | [Collinear plot](collinear.md) |
| Which intervals match when the sequences radiate out from a centre? | [Radial synteny](radial-synteny.md) |
| Which query interval on the left matches which subject interval on the right? | [Sankey view](sankey.md) |
| Where are features, primers and cut sites on one circular sequence? | [Circular map](circular.md) |
| What do the slopes of local alignments show about their relative direction? | [Alignment dot plot](dotplot.md) |
| Which region of an overview is enlarged in another panel? | [Fusion composition](fusion.md) |

<a id="adjust"></a>

## Adjust an existing plot {#adjust}

Use [sequence order and orientation](multi-sequence-layout.md) to change how you read the backbones, [features and labels](annotations.md) to identify selected annotations, and [quantitative tracks](tracks.md) to add measured values. For a geometric setting that varies by sequence or strand, consult [flexible parameters](parameters.md).
