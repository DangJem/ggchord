# geom-seq.R - sequence arc layer
# Fetches pre-computed sequence arc data from the package environment and renders it with geom_path
# Sequence layout parameters are specified in this layer and stored for use at print time

# ---------------------------------------------------------------------------
# Internal geom used by the optional sequence-group label layer.  It exposes
# text colour under an internal aesthetic name ("zcolour") so that the group
# labels can have their own identity scale without colliding with the Seq ID
# colour scale used by geom_seq().
# ---------------------------------------------------------------------------
rename_text_colour_geom <- function(geom = ggplot2::GeomText,
                                    old_aes = "colour",
                                    new_aes = "group_colour") {
  new_geom <- ggplot2::ggproto(
    paste0("GeomChord", sub("^Geom", "", class(geom)[1])), geom
  )

  aes_names <- names(new_geom$default_aes)
  aes_names[aes_names == old_aes] <- new_aes
  names(new_geom$default_aes) <- aes_names

  old_handle_na <- geom$handle_na
  new_geom$handle_na <- function(self, data, params) {
    colnames(data)[colnames(data) == new_aes] <- old_aes
    old_handle_na(data, params)
  }
  old_draw_key <- geom$draw_key
  new_geom$draw_key <- function(data, params, size) {
    colnames(data)[colnames(data) == new_aes] <- old_aes
    old_draw_key(data, params, size)
  }
  new_geom
}

seq_group_label_geom <- rename_text_colour_geom()
seq_geom <- rename_geom_aes(GeomPath, renames = c(colour = "seq_colour"))

#' Add a sequence arc layer
#'
#' Draws arcs (or straight lines, depending on the curvature setting) representing sequences in the chord diagram.
#' Sequence layout parameters (order, orientation, radius, curvature, colors, etc.) are specified here.
#' Sequences can be visually grouped with an extra inter-group gap and optional group labels.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param seq_order Optional character vector. Specifies the drawing order of sequences
#' @param seq_labels Optional character vector or named vector. Sequence labels
#' @param seq_orientation Optional numeric (1 or -1). Sequence orientation, default 1
#' @param seq_gap Optional numeric. Gap proportion between sequences, default 0.03
#' @param seq_radius Optional numeric (> 0). Sequence arc radius, default 1.0
#' @param seq_curvature Optional numeric. Arc curvature (0=straight, 1=standard arc, >1=more curved), default 1.0
#' @param seq_colors Optional color vector or named vector. Sequence colors
#' @param seq_group Optional group specification. NULL disables grouping unless
#'   \code{seq_data} contains a \code{seq_group} column; a single string naming
#'   such a column uses that column; otherwise a single value, named vector,
#'   unnamed vector matching the sequences, or a list (see
#'   \code{process_sequence_param()}) is accepted.
#' @param seq_group_gap Numeric, default 0.08. Extra gap proportion inserted
#'   between consecutive groups (in addition to \code{seq_gap}).
#' @param seq_group_labels Logical or character, default TRUE. When TRUE the
#'   group names are drawn at the angular midpoint of each group; a named
#'   character vector can be used to override the group label text.
#' @param seq_group_label_radius Numeric, default 1.35. Radial position of the
#'   group labels as a multiplier of the group's outermost sequence radius
#'   (same convention as \code{seq_label_radius}: \code{1} on the arc,
#'   \code{> 1} outside).
#' @param seq_group_colors Optional named color vector by group name (or an
#'   unnamed vector recycled positionally). Colours the group labels.
#' @param linewidth Arc line width, default 0.9
#' @param show_legend Whether to show the legend for this layer, default TRUE
#' @param legend_position Position of this layer's legend (the Seq ID legend):
#'   one of "left", "right", "top", "bottom" or "inside", default "right". Pass
#'   NULL to let the legend follow \code{theme(legend.position = ...)} together
#'   with the other legends. Can also be set with \code{theme(legend.position.seq = ...)}.
#' @param ... Additional arguments passed to \code{geom_path()}
#'
#' @return A list of ggplot2 layers
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' p <- ggchord(seq_data_example) + geom_seq()
#' p
geom_seq <- function(mapping = NULL, data = NULL,
                     seq_order = NULL,
                     seq_labels = NULL,
                     seq_orientation = NULL,
                     seq_gap = NULL,
                     seq_radius = NULL,
                     seq_curvature = NULL,
                     seq_colors = NULL,
                     seq_group = NULL,
                     seq_group_gap = 0.08,
                     seq_group_labels = TRUE,
                     seq_group_label_radius = 1.35,
                     seq_group_colors = NULL,
                     linewidth = 0.9,
                     show_legend = TRUE,
                     legend_position = "right",
                     ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  if (!missing(seq_labels) && !is.null(seq_labels)) {
    ggchord_deprecate_once(
      "geom_seq(seq_labels)",
      "scale_seq_colour_manual(labels = ...) and geom_seq_label(labels = ...)"
    )
  }
  if (!missing(legend_position)) {
    ggchord_deprecate_once(
      "geom_seq(legend_position)",
      "guides(seq_colour = guide_ggchord_legend(position = ...))"
    )
  }

  if (!is.numeric(seq_group_gap) || length(seq_group_gap) != 1 ||
      !is.finite(seq_group_gap) || seq_group_gap < 0) {
    ggchord_stop("seq_group_gap must be a finite non-negative number")
  }
  if (!is.numeric(seq_group_label_radius) || length(seq_group_label_radius) != 1 ||
      !is.finite(seq_group_label_radius)) {
    ggchord_stop("seq_group_label_radius must be a finite number")
  }

  # The layout is computed at build time (ggplot_build.ggchord). The
  # parameters are attached to the layer itself so that the plot object is
  # fully self-contained.
  lyr <- ggplot2::layer(
    data        = data.frame(x = numeric(0), y = numeric(0),
                             seq_id = character(0)),
    mapping     = aes(x = x, y = y, group = seq_id, seq_colour = seq_id),
    stat        = "identity",
    geom        = seq_geom,
    position    = "identity",
    show.legend = if (identical(show_legend, TRUE)) {
                    c(seq_colour = TRUE, fill = FALSE)
                  } else show_legend,
    inherit.aes = FALSE,
    check.param = FALSE,
    key_glyph   = key_glyph_seq,
    params      = list(
      linewidth = linewidth,
      arrow = grid::arrow(type = "closed", length = grid::unit(2.4, "mm")),
      ...
    )
  )
  lyr$ggchord_type <- "seq"
  lyr$ggchord_params <- list(
    type                  = "seq",
    seq_order             = seq_order,
    seq_labels            = seq_labels,
    seq_orientation       = seq_orientation,
    seq_gap               = seq_gap,
    seq_radius            = seq_radius,
    seq_curvature         = seq_curvature,
    seq_colors            = seq_colors,
    seq_group             = seq_group,
    seq_group_gap         = seq_group_gap,
    seq_group_labels      = seq_group_labels,
    seq_group_label_radius = seq_group_label_radius,
    seq_group_colors      = seq_group_colors,
    legend_position       = legend_position
  )
  lyr <- ggchord_capture_layer_input(
    lyr, data, mapping, c("seq_id", "length", "seq_group")
  )
  lyr <- ggchord_add_legacy_scale(
    lyr, !missing(seq_colors) && !is.null(seq_colors), "seq_colors",
    "seq_colour", "scale_seq_colour_manual(values = ...)"
  )
  lyr <- ggchord_add_legacy_scale(
    lyr, !missing(seq_group_colors) && !is.null(seq_group_colors),
    "seq_group_colors", "group_colour",
    "scale_group_colour_manual(values = ...)"
  )

  list(lyr)
}

#' Add sequence-group labels
#'
#' Draws labels for contiguous sequence groups defined by `seq_group` in
#' `geom_seq()`. This dedicated layer lets group-label typography be
#' controlled independently through `ggchord.group.label`.
#'
#' @param mapping,data Optional layer mapping and sequence data.
#' @param labels Logical or character labels. `TRUE` uses group names;
#'   named values replace selected group names.
#' @param radius Radial position relative to the outermost sequence radius.
#' @param size Optional text size in millimetres. The theme element is used
#'   when `NULL`.
#' @param show_legend Whether to show a group-colour legend.
#' @param ... Additional fixed arguments passed to `geom_text()`.
#'
#' @return A list containing one ggplot2 layer.
#' @export
geom_seq_group_label <- function(mapping = NULL, data = NULL,
                                 labels = TRUE, radius = 1.35, size = NULL,
                                 show_legend = FALSE, ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  if (!is.numeric(radius) || length(radius) != 1 || !is.finite(radius)) {
    ggchord_stop("geom_seq_group_label(): radius must be a finite number")
  }
  if (!is.null(size) && (!is.numeric(size) || length(size) != 1 ||
      !is.finite(size) || size <= 0)) {
    ggchord_stop("geom_seq_group_label(): size must be NULL or a positive number")
  }

  lyr <- ggplot2::layer(
    data = data.frame(
      text_x = numeric(0), text_y = numeric(0), label = character(0),
      text_angle = numeric(0), hjust = numeric(0), vjust = numeric(0),
      size = numeric(0), group_id = character(0)
    ),
    mapping = aes(
      x = text_x, y = text_y, label = label, angle = text_angle,
      hjust = hjust, vjust = vjust, size = I(size),
      group_colour = group_id
    ),
    stat = "identity", geom = seq_group_label_geom, position = "identity",
    show.legend = show_legend, inherit.aes = FALSE, check.param = FALSE,
    params = list(...)
  )
  lyr$ggchord_type <- "seq_group_label"
  lyr$ggchord_theme_element <- "ggchord.group.label"
  lyr$ggchord_params <- list(
    type = "seq_group_label", labels = labels, radius = radius, size = size
  )
  lyr <- ggchord_capture_layer_input(
    lyr, data, mapping, c("seq_id", "length", "seq_group")
  )
  list(lyr)
}

#' Build a concrete sequence-group label layer from computed layout data
#'
#' Group labels are appended at build time (not when the user calls
#' \code{geom_seq()}), which keeps \code{geom_seq()} backward compatible: it
#' always returns a single sequence-arc layer.
#' @keywords internal
ggchord_group_label_layer <- function(group_labels) {
  if (is.null(group_labels) || nrow(group_labels) == 0) return(NULL)
  lyr <- ggplot2::layer(
    data = group_labels,
    mapping = aes(x = text_x, y = text_y, label = label,
                  angle = text_angle, hjust = hjust, vjust = vjust,
                  size = I(size), group_colour = group_id),
    stat = "identity",
    geom = seq_group_label_geom,
    position = "identity",
    show.legend = FALSE,
    inherit.aes = FALSE,
    check.param = FALSE,
    params = list()
  )
  lyr$ggchord_theme_element <- "ggchord.group.label"
  lyr
}
