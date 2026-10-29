# Run from the source repository root after loading ggchord.
# Reuse original FASTA/GFF3/BLAST records; no source records are rewritten.
library(ggchord)
seqs <- read_fasta_lengths(files = "examples/fasta/*.fna", source_file = TRUE)
genes <- read_gff3(files = "examples/gff3/*.gff3", unstranded = "keep", source_file = TRUE)
alignments <- read_blast(files = "examples/blastn/*.o7", source_file = TRUE)

# Genome membership is supplied explicitly by the analyst. These assemblies
# represent four distinct original sequence records, rather than merged contigs.
seqs$genome <- seqs$accver
genomes <- data.frame(genome = seqs$genome)
hierarchy <- validate_ggchord_hierarchy(seqs, genomes, strict = TRUE)
prepared <- prepare_ggchord_data(hierarchy$seq_data, alignments,
  features = list(genes = genes))
prepared$report$validation$after$summary
prepared$report$pair_counts

# A user-supplied cluster table can separately connect homologous features.
# Cluster assignments cannot be inferred from gene names or this GFF3 table.
# membership <- read.delim("my_membership.tsv") # feature_id, cluster_id
# pairs <- data.frame(qaccver = "explicit query", saccver = "explicit subject")
# homologs <- link_ggchord_clusters(genes, membership, pairs)
# geom_link_ribbon(data = homologs$data) uses interval homology anchors;
# it does not represent a newly computed nucleotide alignment.

ggchord(prepared$seq_data, prepared$ribbon_data) + geom_seq() +
  geom_link_ribbon() + geom_feature(data = prepared$features$genes) +
  coord_collinear()
