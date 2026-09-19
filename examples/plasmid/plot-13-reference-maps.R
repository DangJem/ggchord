# Reproduce the 13 circular plasmid examples used by the v0.13.0 visual
# acceptance script. Run from the package root:
#
#   Rscript examples/plasmid/plot-13-reference-maps.R [output-directory]
#
# The returned PNGs use view_ggchord() so each map receives a content-derived
# physical size. In particular, the dense pETDuet-1 perimeter uses a wide
# canvas instead of compressing its text into a square device.

devtools::load_all(quiet = TRUE)

reference_map_names <- c(
  "pBR322",
  "pUC19",
  "pBluescript II SK(+)",
  "pSB1C3",
  "pET-28a(+)",
  "pETDuet-1",
  "pcDNA3.1(+)",
  "pTRE-Tight-BI",
  "pSpCas9(BB)-2A-GFP (PX458)",
  "pDONR221",
  "pCAMBIA1300",
  "pEarleyGate 201",
  "pTRIPZ"
)

read_reference_fasta <- function(name) {
  path <- file.path("examples", "plasmid", paste0(name, ".fna"))
  lines <- readLines(path, warn = FALSE)
  dna <- paste0(lines[!grepl("^>", lines)], collapse = "")
  data.frame(
    accver = make.names(name),
    label = name,
    length = nchar(dna),
    sequence = dna,
    stringsAsFactors = FALSE
  )
}

reference_restriction_sites <- function(sequence, name) {
  sites <- find_restriction_sites(sequence)
  if (identical(name, "pETDuet-1")) {
    # The checked-in pETDuet-1 reference PNG uses a dense unique-6+ display.
    return(filter_restriction_sites(
      sites,
      set = "unique_6plus",
      parent_set = "commercial_nonredundant"
    ))
  }
  filter_restriction_sites(
    sites,
    set = "reference",
    parent_set = "commercial_nonredundant"
  )
}

build_reference_plot <- function(name) {
  sequence <- read_reference_fasta(name)
  features <- find_common_features(sequence)
  sites <- reference_restriction_sites(sequence, name)
  primers <- find_primer_bindings(sequence, set = "reference")

  feature_tracks <- position_feature_stack(
    spacing = .085,
    base_position = position_plasmid()
  )
  primer_tracks <- position_feature_stack(
    spacing = .085,
    base_position = position_plasmid()
  )

  primer_layers <- NULL
  if (nrow(primers)) {
    primers$feature_label <- primers$name
    primer_layers <- list(
      geom_primer(data = primers, position = primer_tracks),
      geom_feature_label_repel(
        data = primers,
        position = primer_tracks,
        external = TRUE,
        max_overlaps = 0
      )
    )
  }

  ggchord(sequence, validate = "none") +
    geom_seq(seq_style = "double") +
    geom_feature_plasmid(data = features, position = feature_tracks) +
    geom_feature_label_repel(
      data = features,
      position = feature_tracks,
      external = TRUE,
      max_overlaps = 0
    ) +
    geom_restriction_site(data = sites) +
    primer_layers +
    geom_seq_center_label() +
    coord_circular(rotation = 90) +
    theme_ggchord_plasmid()
}

render_reference_plot <- function(name, output_dir) {
  plot <- build_reference_plot(name)
  preview_file <- view_ggchord(plot, viewer = "none")
  output_file <- file.path(output_dir, paste0(name, ".png"))
  copied <- file.copy(preview_file, output_file, overwrite = TRUE)
  if (!copied) stop("Could not copy preview to: ", output_file)
  output_file
}

arguments <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(arguments) && nzchar(arguments[1L])) {
  arguments[1L]
} else {
  file.path(tempdir(), "ggchord-13-reference-maps")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# The 13 explicit examples. Use build_reference_plot("pBR322") when an
# editable ggplot object is wanted instead of a rendered PNG.
output_files <- vapply(
  reference_map_names,
  render_reference_plot,
  character(1L),
  output_dir = output_dir,
  USE.NAMES = TRUE
)

writeLines(c(
  "Rendered 13 reference maps with view_ggchord():",
  unname(output_files)
))
