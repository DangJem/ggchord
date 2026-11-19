<a id="fields"></a>

## Match records to the reference {#fields}

| Table | Fields to start with | Meaning |
| --- | --- | --- |
| Reference | `accver, length` | One row per source ID, with its full sequence length. Declare `sequence_unit = "aa"` for protein coordinates. |
| Features | `accver, start, end` | One interval on a reference. Preserve `strand` where known. A stable `feature_id` lets subset labels find the corresponding feature. |
| Ribbons | `qaccver, saccver, qstart, qend, sstart, send` | Two intervals on named references. Endpoint order carries alignment direction. |
| Signal | `accver, start, end` plus the value column | Map the value to height or fill for the chosen geometry. |

Positions are one-based and interval ends are inclusive. The interval 1–100 spans 100 bases on a bp reference. The same numeric interval on an aa reference spans 100 residues. Preserve the identifiers in the source files, including accession versions.

<a id="prepared"></a>

## Pass tables to the right component {#prepared}

Start with `ggchord(reference)`. The second argument is `mapping`. Give other tables to the layers that draw them, for example `geom_feature(data = features)` or `geom_link_ribbon(data = hits)`. Alternatively, prepare tables once. The example reads back the feature input. The expression `ggchord(prepared) + geom_seq() + geom_feature()` draws it on the reference.

{{example:prepared-tables}}

The teaching annotation table retains its original strand in `source_strand`. Select that field explicitly. When an interval is measured on a protein and the target is genomic, use its [CDS coordinate model](projection.md).
