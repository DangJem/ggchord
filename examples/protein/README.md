# Complete proteins and observed BLASTp comparisons

Four unchanged complete proteins from the public EMALE/Mavirus FASTA archive:
BVI_029_0020 (302 aa), BVI_058_0011 (302 aa), E4-10_023_0015 (306 aa),
and ADZ16414.1 (310 aa). The three EMALE records share the homologous group
containing BVI_023A_0009; the Mavirus protein is the corresponding homolog.
Selection includes close and more distant comparisons, without modifying residues.

`all-vs-all-blastp.tsv` is actual NCBI BLASTp 2.17.0+ output, using the same
complete FASTA as query and subject, `-outfmt 6 -evalue 1e-5 -max_hsps 1`.
It retains all 16 ordered pairs, including self hits. `sources.json` records
source and output SHA256 values and the exact command. The two bundled R tables
retain the four complete proteins and one direction for each of six unordered pairs.

Rebuild from the source archive with `data-raw/prepare-bundled-proteins.py`, then
run `Rscript data-raw/generate_protein_data.R`. BLASTp is a development-only tool;
loading and plotting the bundled tables does not require it.

Source archive: https://raw.githubusercontent.com/thackl/gggenomes/master/data-raw/emales.tgz
Repository: https://github.com/thackl/gggenomes (GPL-3). Records and source headers
are retained; this selection does not copy gggenomes implementation code.
