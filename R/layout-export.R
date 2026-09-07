# layout-export.R - stable public layout export

#' Export computed ggchord geometry
#'
#' Builds a ggchord plot and returns selected computed geometry in a stable,
#' layer-aware structure. Every exported row records its source layer and,
#' where applicable, its source input row. Original mapped columns already
#' attached to the geometry are preserved.
#'
#' @param plot A ggchord plot.
#' @param include Character vector selecting any of \code{"seq"},
#'   \code{"ribbon"}, \code{"gene"}, \code{"feature"}, \code{"labels"},
#'   \code{"link"}, \code{"restriction"}, and \code{"axis"}. Point links and
#'   restriction sites require explicit selection.
#' @param original_data Logical. Include the resolved per-layer input tables
#'   under \code{original_data}, default \code{FALSE}.
#'
#' @return An object of class \code{ggchord_layout_export}. Selected component
#'   names contain data frames. The \code{metadata} member records rotation,
#'   aspect ratio, limits, coordinate units and transformation state.
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' data(ribbon_data_example)
#' p <- ggchord(seq_data_example, ribbon_data_example) +
#'   geom_seq() + geom_link_ribbon()
#' exported <- export_ggchord_layout(p, include = c("seq", "ribbon"))
#' names(exported)
export_ggchord_layout <- function(
    plot,
    include = c("seq", "ribbon", "gene", "feature", "labels", "axis"),
    original_data = FALSE) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  if (!inherits(plot, "ggchord") || is.null(plot$ggchord)) {
    ggchord_stop("export_ggchord_layout(): plot must be a ggchord object")
  }
  allowed <- c("seq", "link", "ribbon", "gene", "feature", "labels",
               "restriction", "axis")
  if (!is.character(include) || anyNA(include) ||
      any(!include %in% allowed)) {
    ggchord_stop(
      "export_ggchord_layout(): include must contain only: ",
      paste(allowed, collapse = ", ")
    )
  }
  include <- unique(include)
  if (!is.logical(original_data) || length(original_data) != 1L ||
      is.na(original_data)) {
    ggchord_stop("export_ggchord_layout(): original_data must be TRUE or FALSE")
  }

  layout <- get_chord_layout(plot, build = TRUE)
  exported <- stats::setNames(vector("list", length(include)), include)
  for (nm in include) exported[[nm]] <- data.frame()
  kept_inputs <- list()

  classify_component <- function(component, layer) {
    switch(component,
      seq = "seq",
      ribbon = "ribbon",
      link = "link",
      gene_poly = if (isTRUE(layer$ggchord_params$is_feature)) {
        "feature"
      } else {
        "gene"
      },
      seq_region = "feature",
      restriction_site = "restriction",
      gene_text = "labels",
      gene_text_repel = "labels",
      gene_label_segment = "labels",
      gene_label_repel = "labels",
      seq_label = "labels",
      seq_center_label = "labels",
      axis_line = "axis",
      axis_seg = "axis",
      axis_text = "axis",
      axis = "axis",
      NA_character_
    )
  }

  seen_layer_ids <- character(0)
  for (layer in plot$layers) {
    layer_id <- layer$ggchord_layer_id %||% ""
    if (!nzchar(layer_id)) next
    if (layer_id %in% seen_layer_ids) next
    seen_layer_ids <- c(seen_layer_ids, layer_id)
    matching <- Filter(function(x) {
      identical(x$ggchord_layer_id %||% "", layer_id)
    }, plot$layers)
    if (any(vapply(matching, function(x) {
      isTRUE(x$ggchord_params$is_feature)
    }, logical(1)))) {
      layer$ggchord_params$is_feature <- TRUE
    }
    registry <- layout$layer_geometry[[layer_id]]
    if (is.null(registry)) next

    for (component in names(registry)) {
      category <- classify_component(component, layer)
      if (is.na(category) || !category %in% include) next
      geometry <- registry[[component]]
      if (!is.data.frame(geometry) || nrow(geometry) == 0L) next
      geometry <- as.data.frame(geometry, stringsAsFactors = FALSE)
      input <- layout$layer_inputs[[layer_id]][[component]]
      if (!"source_row" %in% names(geometry)) {
        geometry$source_row <- if (!is.null(input) &&
            all(c("accver") %in% names(geometry)) &&
            "accver" %in% names(input)) {
          as.integer(match(
            as.character(geometry$accver), as.character(input$accver)
          ))
        } else {
          rep(NA_integer_, nrow(geometry))
        }
      }
      geometry$layer_id <- rep(layer_id, nrow(geometry))
      geometry$component <- rep(component, nrow(geometry))
      exported[[category]] <- if (nrow(exported[[category]]) == 0L &&
          ncol(exported[[category]]) == 0L) {
        geometry
      } else {
        ggchord_rbind_fill(list(exported[[category]], geometry))
      }

      if (isTRUE(original_data)) {
        if (!is.null(input)) {
          if (is.null(kept_inputs[[layer_id]])) kept_inputs[[layer_id]] <- list()
          kept_inputs[[layer_id]][[component]] <- input
        }
      }
    }
  }

  # The sequence axis is now an automatic decoration rather than a public
  # layer, so export it directly from the computed layout under a stable
  # synthetic layer id.
  if ("axis" %in% include && nrow(exported$axis) == 0L) {
    automatic_axis <- ggchord_axis_geometry(layout)
    if (nrow(automatic_axis) > 0L) {
      automatic_axis$source_row <- NA_integer_
      automatic_axis$layer_id <- ".automatic-axis"
      automatic_axis$component <- automatic_axis$.component
      exported$axis <- automatic_axis
    }
  }

  coord <- plot$coordinates
  fitted <- ggchord_adaptive_limits(layout)
  metadata <- list(
    coordinate_space = "ggchord_layout",
    units = "data",
    rotation = layout$rotation %||%
      if (isTRUE(coord$ggchord_coord)) coord$rotation else 0,
    rotation_applied = TRUE,
    ratio = coord$ratio %||% 1,
    coord_transform_applied = FALSE,
    fit = coord$fit %||% "labels",
    xlim = coord$user_xlim %||% fitted$xlim,
    ylim = coord$user_ylim %||% fitted$ylim,
    clip = coord$clip %||% "off"
  )
  if (isTRUE(coord$ggchord_circular)) {
    metadata$coordinate <- "circular"
    metadata$gap <- coord$circular_gap
    metadata$direction <- coord$circular_direction
  } else {
    metadata$coordinate <- "chord"
  }

  exported$metadata <- metadata
  if (isTRUE(original_data)) exported$original_data <- kept_inputs
  class(exported) <- c("ggchord_layout_export", "list")
  exported
}
