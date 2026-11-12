# Optional dependency: dplyr. One real HSP is deliberately shown 12 times.
library(ggchord)
library(dplyr)
observed <- ribbon_data_example |> slice_head(n = 1) |>
  mutate(source_hsp_id = 'explicit-one-HSP-twelve-views')
source_seqs <- seq_data_example |> select(accver, length)
target_seqs <- source_seqs |> filter(accver == observed$saccver)
center <- ggchord(source_seqs) + geom_seq(show.legend = FALSE) +
  geom_seq_label(size = 2.7) +
  coord_chord() + theme(ggchord.axis = element_blank())
auto_fusion <- ggchord_fusion(center, layout = 'adaptive')
for (i in seq_len(12)) {
  id <- paste0('view-', i)
  view <- ggchord(target_seqs) + geom_seq(show.legend = FALSE) +
    coord_collinear(seq_axis = if (i %% 2) 'horizontal' else 'vertical') +
    theme(ggchord.axis = element_blank()) + ggplot2::labs(title = paste(id, observed$saccver, sep = '\n'))
  auto_fusion <- auto_fusion +
    fusion_panel(view, id, side = 'auto') +
    fusion_ribbon(observed, qpanel = 'center', spanel = id, allow_repeated = TRUE,
      fill = '#548CA0', alpha = .2)
}
auto_fusion
