library(ggchord)
hit <- ribbon_data_example[1, ]
ggchord(seq_data_example[1:2, ]) +
  geom_alignment_segment(data = hit, linewidth = 1) +
  coord_dotplot(x_seqs = hit$qaccver, y_seqs = hit$saccver)
