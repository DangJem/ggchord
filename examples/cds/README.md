# Source-backed CDS / 有来源的 CDS

`NC_000932.gb` is the unchanged Arabidopsis thaliana chloroplast GenBank
NC_000932.1 record from the Biopython test archive, retrieved 2026-11-12:
<https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_000932.gb>.
It is a historical record (15-APR-2009), not a claim of the latest reference.
Upstream repository: <https://github.com/biopython/biopython>.
Biopython's permission notice is retained at
[../alignment/BIOPYTHON-LICENSE.rst](../alignment/BIOPYTHON-LICENSE.rst).
NCBI sequence reuse policy: <https://www.ncbi.nlm.nih.gov/home/about/policies/>.
No upstream implementation or plotting code is copied.

`python3 data-raw/prepare-cds-inputs.py` extracts the unchanged full forward
sequence to FASTA, and preserves the annotated psbK (+), rps16 (spliced -)
and petD (spliced +) CDS locations, biological segment order, codon_start,
translation table and original protein translations in expected.json.
Hashes identify the original record and generated FASTA. No alignment,
translation prediction or annotation analysis is performed by this script.
The three original source translations are compared after removing the single
terminal stop retained by ggchord; internal stop characters are never removed.

本目录保留真实源记录、注释链向、分段和预期翻译。使用其明确记录的表11，
不根据图形猜测phase或基因功能；分段间phase仅按已声明codon_start与累积长度推导。

Translation delegates to the public Biostrings API (Artistic-2.0):
<https://bioconductor.org/packages/Biostrings>. No translation code is copied.
