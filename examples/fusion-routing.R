# Optional dependency: dplyr. These are real bundled HSP records.
library(ggchord)
library(dplyr)
first_hit <- ribbon_data_example |> slice_head(n = 1)
channel_hits <- ribbon_data_example |>
  filter(qaccver == first_hit$qaccver, saccver == first_hit$saccver) |>
  mutate(route_group = 'explicit-display-channel')
# Shared corridors describe drawing; they do not aggregate biological evidence.
panel <- ggchord(seq_data_example) +
  geom_seq(show.legend = FALSE) +
  geom_seq_label(size = 3, seq_label_rotation = 90) +
  coord_collinear(seq_axis = 'vertical') +
  theme(ggchord.axis = element_blank())
composition <- ggchord_fusion(panel, center_box = c(.02, .28, .12, .88)) +
  fusion_panel(panel, 'target', box = c(.72, .98, .12, .88)) +
  fusion_panel(panel, 'obstacle', box = c(.4, .6, .22, .78))
routed <- composition +
  fusion_ribbon(channel_hits, qpanel = 'center', spanel = 'target',
    route = 'bundle', route_width = .8, fill = '#548CA0', alpha = .3)
routed
