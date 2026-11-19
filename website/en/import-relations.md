<a id="alignments"></a>

## Import HSP endpoints {#alignments}

Read an installed BLAST table for MT108731.1 against MT118296.1. An HSP is a local aligned segment. Its four endpoints identify the two spans and `pident` is percent identity. Reference lengths come from a sequence table when the alignment file omits them.

{{example:installed-blast}}

For a PAF file, `read_paf()` also returns query and subject lengths. The reader converts file boundaries to one-based inclusive endpoints and preserves reverse orientation. Pass the resulting intervals directly to a [ribbon layer](multi-sequence-links.md).

<a id="members"></a>

## Read a list of group members {#members}

A group file records membership. Here the first named group contains four plant features. Keep the identifiers and line number if you later join them to a GFF3 annotation. Pairwise connections require explicit pairs or a declared choice of which members to compare.

{{example:installed-members}}
