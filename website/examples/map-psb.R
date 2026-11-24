library(ggchord)
library(ggplot2)
library(dplyr)
library(scales)

draw_reference_map <- function(reference, overview = FALSE) {
  sequence <- as_ggchord_seq(reference$seq_data, topology = "circular")
  features <- reference$features
  # Keep the full segmented feature table; the label layer follows its lanes.
  lanes <- position_feature_stack(
    spacing = 1.3, direction = "negative", order = "priority"
  )
  p <- ggchord(sequence) +
    geom_seq(geometry = "double", thickness = 1.4, linewidth = 0.52, colour = "#292826") +
    geom_feature(
      aes(fill = feature_color, shape = feature_shape),
      data = features, id = "annotations", thickness = 5.3, position = lanes,
      colour = "#393635", linewidth = 0.27, segment_boundaries = TRUE
    ) +
    scale_fill_identity() +
    scale_shape_identity() +
    geom_seq_label_center(
      aes(label = paste0(accver, "\n", comma(length), " bp")),
      size = 3.55, fontface = "bold"
    ) +
    scale_seq_continuous(
      name = NULL, labels = label_comma(),
      guide = guide_seq_axis(side = "negative")
    ) +
    coord_circular(fit = "labels", clip = "off") +
    theme_ggchord_plasmid()
  label_data <- if (overview) filter(features, type %in% c("CDS", "rep_origin")) else features
  p <- p +
    geom_feature_label(
      aes(label = label, fill = feature_color),
      data = label_data, anchor_layer = "annotations",
      size = if (overview) 3.6 else 3.0,
      segment_colour = "#77736E", segment_linewidth = 0.20,
      label_layout = layout_label_repel(
        arrangement = "local", side = if (overview) "positive" else "both",
        offset = if (overview) 3 else 0, fit = "none",
        outer_size = if (overview) NULL else 3.3
      )
    )
  if (overview) return(p)

  if (nrow(reference$sites)) {
    # Label each documented boundary from its preceding residue on the circle.
    site_labels <- reference$sites |>
      mutate(cut_boundary = coalesce(cut_top, cut_bottom)) |>
      filter(!is.na(cut_boundary)) |>
      group_by(accver, cut_boundary) |>
      summarise(
        enzymes = paste(unique(enzyme), collapse = " – "),
        source_members = list(site_id),
        .groups = "drop"
      ) |>
      mutate(
        pos = if_else(cut_boundary == 0, sequence$length[[1L]], cut_boundary),
        enzyme_label = enzymes,
        coordinate_label = paste0("(", cut_boundary, ")"),
        label = paste(enzymes, coordinate_label)
      )
    p <- p +
      geom_site_label(
        aes(label = label), data = site_labels, anchor_layer = NULL, size = 3, fontface = "bold", colour = "#383535",
        segment_colour = "#8A8177", segment_linewidth = 0.19,
        label_layout = layout_label_repel(
          arrangement = "local", orientation = "horizontal", side = "positive",
          offset = 2.6, fit = "none"
        )
      )
  }
  if (nrow(reference$primers)) {
    primers <- reference$primers |>
      mutate(label = paste0(label, " (", start, "–", end, ")"))
    p <- p +
      geom_primer(
        data = primers, id = "primers",
        shape = shape_feature("arrow", head_length = 0.8), fill = "#7852A2",
        colour = "#7852A2", position = position_offset(0), thickness = 0.76
      ) +
      geom_primer_label(
        aes(label = label), anchor_layer = "primers", size = 3, fontface = "bold", colour = "#7852A2",
        segment_colour = "#8D65BB", segment_linewidth = 0.26,
        label_layout = layout_label_repel(
          arrangement = "local", orientation = "horizontal", side = "positive",
          offset = 2.7, fit = "none"
        )
      )
  }
  p
}

reference <- flagship_circular_example$references[["pSB1C3"]]
draw_reference_map(reference)
