library(ggchord)
refs <- head(seq_data_example, 2)
pair <- subset(ribbon_data_example, qaccver == refs$accver[1] & saccver == refs$accver[2])
hit <- pair[which.min(pair$length), ]
overview <- ggchord(refs[1, ]) + geom_seq() + coord_circular()
detail <- ggchord(refs) +
  geom_seq() + geom_seq_label() +
  coord_collinear(seq_axis = "vertical")
ggchord_fusion(overview) +
  fusion_panel(detail, id = "detail", side = "right") +
  fusion_ribbon(hit, qpanel = "center", spanel = "detail",
    fill = "#269B82", colour = "#237F70", linewidth = 0.2, alpha = 0.55)
