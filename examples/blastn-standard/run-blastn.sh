#!/usr/bin/env bash
# Optional fresh computation. Bundled .tsv files are historical row-preserving
# column conversions; this script does not claim byte-identical new BLAST runs.
set -euo pipefail
seqs=(MT108731.1 MT118296.1 OQ646790.1 OR222515.1)
for ((i=0; i<${#seqs[@]}-1; i++)); do
  for ((j=i+1; j<${#seqs[@]}; j++)); do
    blastn -query "examples/fasta/${seqs[$i]}.fna" \
      -subject "examples/fasta/${seqs[$j]}.fna" -outfmt 6 \
      -out "examples/blastn-standard/${seqs[$i]}__${seqs[$j]}.tsv"
  done
done
