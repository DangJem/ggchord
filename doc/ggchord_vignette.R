## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.width = 7,
  fig.height = 7
)

## ----install, eval=FALSE------------------------------------------------------
#  # Install from local
#  devtools::install_local("path/to/ggchord_0.2.0.tar.gz")
#  
#  # Load the package
#  library(ggchord)
#  library(dplyr)  # For data processing

## ----load, echo=FALSE---------------------------------------------------------
library(ggchord)
library(dplyr)

## ----load-data----------------------------------------------------------------
# Load built-in example data
data(seq_data_example)    # Sequence data
data(ribbon_data_example) # BLAST alignment data
data(gene_data_example)   # Gene annotation data (short genes filtered out)

# Inspect data structure
head(seq_data_example)
head(ribbon_data_example)
head(gene_data_example)

## ----part1-1------------------------------------------------------------------
# Basic chord diagram: sequences only
p1 <- ggchord(
  seq_data = seq_data_example
)
print(p1)

## ----part1-2------------------------------------------------------------------
p2 <- ggchord(
  seq_data = seq_data_example,
  seq_order = c("MT118296.1", "OR222515.1", "MT108731.1", "OQ646790.1"),  # Sequence order
  seq_orientation = c(1, -1, 1, -1),  # Direction (1 = forward, -1 = reverse)
  seq_curvature = c(0, 2, -2, 6),     # Curvature (0 = straight line, 1 = standard arc)
  seq_colors = c("steelblue", "orange", "pink", "yellow")  # Sequence colors
)
print(p2)

## ----part2-1------------------------------------------------------------------
# Show sequence alignment relationships
p3 <- ggchord(
  seq_data = seq_data_example,
  ribbon_data = ribbon_data_example
)
print(p3)

## ----part2-2------------------------------------------------------------------
# Color by query sequence
p4 <- ggchord(
  seq_data = seq_data_example,
  ribbon_data = ribbon_data_example,
  ribbon_color_scheme = "query"  # Ribbon color matches query sequence
)
print(p4)

# Single color
p5 <- ggchord(
  seq_data = seq_data_example,
  ribbon_data = ribbon_data_example,
  ribbon_color_scheme = "single",  # Uniform color
  ribbon_colors = "orange"         # Custom color
)
print(p5)

## ----part2-4------------------------------------------------------------------
p6 <- ggchord(
  seq_data = seq_data_example,
  ribbon_data = ribbon_data_example,
  seq_orientation = c(1, -1, 1, -1),  # Sequence orientation
  seq_curvature = c(0, 2, -2, 6),     # Curvature
  seq_gap = c(0.1, 0.05, 0.09, 0.05), # Inter-sequence gaps
  seq_radius = c(1, 5, 1, 1)          # Sequence radii
)
print(p6)

## ----part3-1------------------------------------------------------------------
# Show gene annotations
p7 <- ggchord(
  seq_data = seq_data_example,
  gene_data = gene_data_example
)
print(p7)

## ----part3-2------------------------------------------------------------------
# Color by gene annotation category
p8 <- ggchord(
  seq_data = seq_data_example,
  gene_data = gene_data_example,
  gene_color_scheme = "manual"  # Color by gene annotation (anno)
)
print(p8)

# Display gene labels and adjust styles
p9 <- ggchord(
  seq_data = seq_data_example,
  gene_data = gene_data_example,
  gene_label_show = TRUE,        # Show labels
  gene_label_rotation = 45,      # Label rotation angle
  gene_label_radial_offset = 0.1, # Radial offset of labels
  panel_margin = list(l = 0.2)    # Left margin
)
print(p9)

## ----part4-2------------------------------------------------------------------
p10 <- ggchord(
  seq_data = seq_data_example,
  ribbon_data = ribbon_data_example,
  gene_data = gene_data_example,
  title = "Multiple Sequence Alignment Chord Diagram (with Gene Annotations)",  # Title
  
  # Sequence parameters
  seq_gap = 0.03,
  seq_radius = c(3, 2, 2, 1),
  seq_orientation = c(-1, -1, -1, 1),
  seq_curvature = c(0, 1, -1, 1.5),
  
  # Gene parameters
  gene_offset = list(             # Gene arrow offset (strand-specific)
    c("+" = 0.2, "-" = -0.2),
    c("+" = 0.2, "-" = -0.2),
    c("+" = 0.2, "-" = 0),
    c("+" = 0.2, "-" = 0.1)
  ),
  gene_label_rotation = list(     # Gene label rotation (strand-specific)
    c("+" = 45, "-" = -45),
    c("+" = 0.2, "-" = -0.2),
    c("+" = 0.2, "-" = 0),
    c("+" = 0.2, "-" = 0.1)
  ),
  gene_width = 0.08,              # Gene arrow width
  gene_label_show = TRUE,         # Show gene labels
  
  # Ribbon parameters
  ribbon_gap = 0.1,
  ribbon_color_scheme = "pident", # Color by identity
  
  # Axis parameters
  axis_label_orientation = c(0, 45, 80, 130), # Tick label angles
  axis_tick_major_length = 0.03,
  axis_label_size = 2,
  
  # Overall style
  rotation = 45,  # Plot rotation angle
  show_axis = TRUE,
  panel_margin = list(t = 0.1)
)
print(p10)

