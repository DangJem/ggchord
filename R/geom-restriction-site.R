GeomRestrictionSite <- ggplot2::ggproto(
  "GeomRestrictionSite", ggplot2::Geom,
  required_aes = c("x", "y"),
  default_aes = ggplot2::aes(
    xend = NA_real_, yend = NA_real_, label = NA_character_,
    enzyme_label = NA_character_, coordinate_label = NA_character_,
    label_order = NA_character_, enzyme_fontface = "plain",
    .component = NA_character_, group = NA_integer_, colour = "#202020",
    alpha = 1, linewidth = .28, linetype = 1, size = 2.9, angle = 0,
    hjust = .5, vjust = .5, family = "", fontface = 1, lineheight = 1.2
  ),
  draw_key = ggplot2::draw_key_path,
  draw_panel = function(data, panel_params, coord, na.rm = FALSE,
                        segment_params = list(), text_params = list(),
                        composite_labels = TRUE) {
    paths <- data[data$.component == "path", , drop = FALSE]
    labels <- data[data$.component == "label", , drop = FALSE]
    for (nm in names(segment_params)) {
      if (nm %in% names(paths)) paths[[nm]] <- segment_params[[nm]]
    }
    for (nm in names(text_params)) {
      if (nm %in% names(labels)) labels[[nm]] <- text_params[[nm]]
    }
    grobs <- list()
    if (nrow(paths)) {
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPath$draw_panel(
        paths, panel_params, coord, lineend = "round", linejoin = "round",
        na.rm = na.rm
      )
    }
    if (nrow(labels)) {
      composite <- isTRUE(composite_labels) &&
        all(c("enzyme_label", "coordinate_label", "label_order") %in%
          names(labels))
      if (composite) {
        placed <- coord$transform(labels, panel_params)
        fragments <- vector("list", nrow(placed))
        for (i in seq_len(nrow(placed))) {
          enzyme_face <- if ("enzyme_fontface" %in% names(placed) &&
              !is.na(placed$enzyme_fontface[i])) {
            placed$enzyme_fontface[i]
          } else "plain"
          parts <- ggchord_site_label_parts(placed$enzyme_label[i], placed$coordinate_label[i],
            placed$size[i], placed$family[i], enzyme_face, placed$lineheight[i],
            scales::alpha(placed$colour[i], placed$alpha[i]))
          enzyme <- parts$enzyme; coordinate <- parts$coordinate; gap <- parts$gap
          enzyme$y <- coordinate$y <- grid::unit(placed$y[i], 'native')
          enzyme$vjust <- coordinate$vjust <- placed$vjust[i]
          has_enzyme <- !is.na(placed$enzyme_label[i]) &&
            nzchar(placed$enzyme_label[i])
          has_coordinate <- !is.na(placed$coordinate_label[i]) &&
            nzchar(placed$coordinate_label[i])
          centre <- grid::unit(placed$x[i], "native")
          if (has_enzyme && has_coordinate) {
            total <- grid::grobWidth(enzyme) + gap +
              grid::grobWidth(coordinate)
            if (placed$label_order[i] == "enzyme_position") {
              enzyme$x <- centre - total / 2
              enzyme$hjust <- 0
              coordinate$x <- enzyme$x + grid::grobWidth(enzyme) + gap
              coordinate$hjust <- 0
            } else {
              enzyme$x <- centre + total / 2
              enzyme$hjust <- 1
              coordinate$x <- enzyme$x - grid::grobWidth(enzyme) - gap
              coordinate$hjust <- 1
            }
            fragments[[i]] <- grid::grobTree(enzyme, coordinate)
          } else if (has_enzyme) {
            enzyme$x <- centre
            enzyme$hjust <- .5
            fragments[[i]] <- enzyme
          } else if (has_coordinate) {
            coordinate$x <- centre
            coordinate$hjust <- .5
            fragments[[i]] <- coordinate
          } else {
            fragments[[i]] <- grid::nullGrob()
          }
        }
        grobs[[length(grobs) + 1L]] <- do.call(grid::grobTree, fragments)
      } else {
        grobs[[length(grobs) + 1L]] <- ggplot2::GeomText$draw_panel(
          labels, panel_params, coord, parse = FALSE,
          check_overlap = FALSE, na.rm = na.rm
        )
      }
    }
    do.call(grid::grobTree, grobs)
  }
)

#' Draw restriction sites, labels and deterministic leaders
#' @param mapping,data Standard layer inputs. Positions refer to the displayed
#'   nucleotide backbone. An explicit `sequence_unit` must be `"bp"`.
#' @param label Draw labels.
#' @param label_style Label enzyme and position, enzyme only, or position only.
#' @param label_side Draw labels outside, inside, or automatically.
#' @param leader Leader routing style. The default `"radial"` measures actual
#'   text boxes, separates site anchors from final radial label positions,
#'   and routes either a direct connector or an independent radial stub plus
#'   fan segment.
#'   All enzyme-facing text edges follow one close circular contour in genomic
#'   order. Dense lateral clusters use evenly spaced contour fans, while sparse
#'   clusters retain their natural angular positions.
#'   `"straight"` connects the tick tip directly to the final label edge,
#'   without detours around the sequence or other leaders. `"elbow"` uses
#'   one bend, starting with a local-normal departure where possible.
#'   Label positions and grouping are shared by all routing styles.
#' @param min_label_gap Minimum genomic fraction between label slots. `NULL`
#'   derives a compact device-aware default from the rendered labels.
#' @param tick_length,label_offset Local-normal distances.
#' @param colour,linewidth,label_size,fontface,family Fixed appearance. With
#'   the default `fontface = NULL`, enzyme names are bold only when that enzyme
#'   cuts once in the supplied sequence; coordinates remain regular. An
#'   explicit `fontface` disables this automatic distinction.
#' @param label_order Automatic mirrored label order, enzyme first, or
#'   position first.
#' @param id Optional unique layer ID used in layout and source exports.
#' @param position,show.legend,inherit.aes Standard layer arguments.
#' @param ... Additional geom parameters.
#' @return A ggplot2 layer.
#' @seealso \code{\link{ggchord_data_contract}}
#' @examples
#' sites <- data.frame(accver = "A", position = 10, enzyme = "E")
#' ggchord(data.frame(accver = "A", length = 100)) +
#'   geom_seq() +
#'   geom_restriction_site(data = sites, id = "sites") +
#'   coord_circular()
#' @export
geom_restriction_site<-function(mapping=NULL,data=NULL,label=TRUE,label_style=c("enzyme_position","enzyme","position"),label_order=c("auto","enzyme_position","position_enzyme"),label_side=c("outside","inside","auto"),leader=c("radial","elbow","straight"),min_label_gap=NULL,tick_length=.04,label_offset=.18,colour="#202020",linewidth=.28,label_size=3.2,fontface=NULL,family=NULL,position="identity",show.legend=FALSE,inherit.aes=FALSE,...,id=NULL){
  ggchord_check_argument_names(sys.call(), sys.function())
  colour_supplied<-!missing(colour);linewidth_supplied<-!missing(linewidth);label_size_supplied<-!missing(label_size);fontface_supplied<-!missing(fontface)&&!is.null(fontface);family_supplied<-!missing(family)&&!is.null(family)
  dots <- list(...)
  colour_supplied <- colour_supplied || "color" %in% names(dots)
  alias <- ggchord_colour_alias(colour, dots, "geom_restriction_site()", sys.call())
  colour <- alias$colour
  if(is.null(data))ggchord_stop("geom_restriction_site(): supply restriction-site data")
  if (identical(leader, "trunk")) ggchord_stop(
    "geom_restriction_site(): leader = 'trunk' was removed in v0.16.1; use leader = 'radial'")
  label_style<-match.arg(label_style);label_order<-match.arg(label_order);label_side<-match.arg(label_side);leader<-match.arg(leader)
  vals<-c(tick_length,label_offset);if(!is.numeric(vals)||any(!is.finite(vals))||any(vals<0)||(!is.null(min_label_gap)&&(!is.numeric(min_label_gap)||length(min_label_gap)!=1L||!is.finite(min_label_gap)||min_label_gap<0)))ggchord_stop("geom_restriction_site(): invalid gaps or offsets")
  segment_params<-list();text_params<-list()
  if(colour_supplied){segment_params$colour<-colour;text_params$colour<-colour}
  if(linewidth_supplied)segment_params$linewidth<-linewidth
  if(label_size_supplied)text_params$size<-label_size
  if(fontface_supplied)text_params$fontface<-fontface
  if(family_supplied)text_params$family<-family
  ggchord_check_geom_params(alias$dots, GeomRestrictionSite, "geom_restriction_site()", position)
  lyr<-ggplot2::layer(data=data.frame(x=numeric(),y=numeric()),mapping=ggplot2::aes(x=x,y=y,group=group,label=label,.component=I(.component),enzyme_label=I(enzyme_label),coordinate_label=I(coordinate_label),label_order=I(label_order),enzyme_fontface=I(enzyme_fontface)),stat="identity",geom=GeomRestrictionSite,position=position,show.legend=show.legend,inherit.aes=inherit.aes,check.aes=FALSE,check.param=FALSE,params=c(list(colour=colour,linewidth=linewidth,size=label_size,segment_params=segment_params,text_params=text_params,composite_labels=!fontface_supplied),alias$dots))
  lyr$ggchord_type<-"restriction_site";lyr$ggchord_theme_components<-c(segment_params="ggchord.restriction.label.segment",text_params="ggchord.restriction.label")
  lyr$ggchord_params<-list(id=id,type="restriction_site",label=label,label_style=label_style,label_order=label_order,label_side=label_side,leader=leader,min_label_gap=min_label_gap,tick_length=tick_length,label_offset=label_offset,label_size=label_size)
  ggchord_capture_layer_input(lyr,data,mapping,c("accver","position","enzyme"))
}
