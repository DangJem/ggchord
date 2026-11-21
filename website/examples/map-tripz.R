library(ggchord)
library(ggplot2)
library(dplyr)
library(scales)

draw_reference_map <- function(reference, overview = FALSE) {
  sequence <- as_ggchord_seq(reference$seq_data, topology = "circular")
  features <- reference$features |>
    mutate(
      source_strand = strand,
      feature_shape = if_else(strand == "+/-", "block", feature_shape),
      strand = if_else(strand == "+/-", ".", strand)
    )
  # Keep the full segmented feature table; the label layer follows its lanes.
  lanes <- position_feature_stack(
    spacing = 1.8, direction = "negative", base_position = position_offset(-3.4)
  )
  p <- ggchord(sequence) +
    geom_seq(geometry = "double", thickness = 1.25, linewidth = 0.52, colour = "#292826") +
    geom_feature(
      aes(fill = feature_color, shape = feature_shape),
      data = features, id = "annotations", thickness = 3.25, position = lanes,
      colour = "#393635", linewidth = 0.22, segment_boundaries = TRUE
    ) +
    scale_fill_identity() +
    scale_shape_identity() +
    geom_seq_center_label(
      aes(label = paste0(accver, "\n", comma(length), " bp")),
      size = 3.55, fontface = "bold"
    ) +
    scale_seq_position(
      breaks = if (overview) breaks_pretty(3) else waiver(), labels = label_comma(),
      guide = guide_seq_axis(side = "negative")
    ) +
    coord_circular(start_angle = 90, fit = "labels", clip = "off") +
    theme_ggchord_plasmid()
  label_data <- if (overview) filter(features, type %in% c("CDS", "rep_origin")) else features
  p <- p +
    geom_feature_label_repel(
      aes(label = label, fill = feature_color),
      data = label_data, anchor_layer = "annotations",
      size = if (overview) 2.8 else 2.55,
      segment_colour = "#77736E", segment_linewidth = 0.20,
      label_layout = label_layout_repel(
        arrangement = "local", side = "both",
        offset = 0, fit = "shrink", min_size = 2.15
      )
    )
  if (overview) return(p)

  if (nrow(reference$sites)) {
    # Use the documented cut boundary, not the start of the recognition motif.
    site_labels <- reference$sites |>
      mutate(pos = coalesce(cut_top, cut_bottom)) |>
      group_by(accver, pos) |>
      summarise(enzymes = paste(unique(enzyme), collapse = " – "), .groups = "drop") |>
      mutate(label = paste0(enzymes, " (", pos, ")"))
    p <- p +
      geom_restriction_site(
        data = reference$sites, size = 1.45, colour = "#514E4C"
      ) +
      geom_site_label_repel(
        aes(label = label), data = site_labels, size = 2.2, colour = "#383535",
        segment_colour = "#8A8177", segment_linewidth = 0.19,
        label_layout = label_layout_repel(
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
        data = primers, id = "primers", shape = "primer_arrow", fill = "#7852A2",
        colour = "#7852A2", position = position_offset(-0.65), thickness = 0.8
      ) +
      geom_primer_label_repel(
        aes(label = label), anchor_layer = "primers", size = 2.35, colour = "#7852A2",
        segment_colour = "#8D65BB", segment_linewidth = 0.26,
        label_layout = label_layout_repel(
          arrangement = "local", orientation = "horizontal", side = "positive",
          offset = 2.7, fit = "none"
        )
      )
  }
  p
}

reference <- flagship_circular_example$references[["pTRIPZ"]]
draw_reference_map(reference)
