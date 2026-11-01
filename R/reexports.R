#' Selected ggplot2 helpers for ggchord plots
#'
#' ggchord re-exports a deliberately small set of ggplot2 helpers that are
#' directly useful when constructing, annotating, styling, or saving a chord
#' plot. Cartesian geoms, coordinates, facets, generic scales, and complete
#' theme presets remain in ggplot2 so IDE completion stays focused.
#'
#' Mapping helpers are \code{aes()}, \code{after_stat()} and \code{after_scale()};
#' \code{ggplot2::stage()} combines stages of one aesthetic. Use \code{labs()} and
#' \code{ggtitle()} for titles, \code{guides()} for guide selection, and \code{theme()} with
#' the \code{element_*()} helpers for appearance. \code{annotate()} adds fixed annotations.
#' \code{margin()} and \code{rel()} specify theme sizes; \code{waiver()} and \code{expansion()}
#' specify automatic or expanded scale settings. \code{position_identity()} and
#' \code{PositionIdentity} are the upstream identity placement helpers.
#' \code{ggsave()} saves a supplied plot; \code{last_plot()} retrieves ggplot2's last plot.
#' Each symbol is the original ggplot2 function or object, and can also be
#' called with its \code{ggplot2::} namespace. \code{unit()} and \code{arrow()} come from grid.
#'
#' @param x,y Aesthetic expressions for mappings; x is also the staged expression or relative size, depending on the helper.
#' @param ... Additional named mappings, guide specifications, theme elements, fixed annotation aesthetics or device arguments, as appropriate for the called helper.
#' @param geom,xmin,xmax,ymin,ymax,xend,yend,na.rm For annotate(), the geometry name, positional bounds and endpoints, and whether missing values are removed silently. In theme(), geom instead sets geometry defaults.
#' @param title,subtitle,caption,tag,dictionary,alt,alt_insight,label Plot labels, title, subtitle, caption and tag; dictionary replaces variable labels and alt/alt_insight provide accessible descriptions. In theme(), title supplies inherited title styling.
#' @param line,rect,text,point,polygon,spacing,margins,aspect.ratio,axis.title,axis.title.x,axis.title.x.top,axis.title.x.bottom,axis.title.y,axis.title.y.left,axis.title.y.right,axis.text,axis.text.x,axis.text.x.top,axis.text.x.bottom,axis.text.y,axis.text.y.left,axis.text.y.right,axis.text.theta,axis.text.r,axis.ticks,axis.ticks.x,axis.ticks.x.top,axis.ticks.x.bottom,axis.ticks.y,axis.ticks.y.left,axis.ticks.y.right,axis.ticks.theta,axis.ticks.r,axis.minor.ticks.x.top,axis.minor.ticks.x.bottom,axis.minor.ticks.y.left,axis.minor.ticks.y.right,axis.minor.ticks.theta,axis.minor.ticks.r,axis.ticks.length,axis.ticks.length.x,axis.ticks.length.x.top,axis.ticks.length.x.bottom,axis.ticks.length.y,axis.ticks.length.y.left,axis.ticks.length.y.right,axis.ticks.length.theta,axis.ticks.length.r,axis.minor.ticks.length,axis.minor.ticks.length.x,axis.minor.ticks.length.x.top,axis.minor.ticks.length.x.bottom,axis.minor.ticks.length.y,axis.minor.ticks.length.y.left,axis.minor.ticks.length.y.right,axis.minor.ticks.length.theta,axis.minor.ticks.length.r,axis.line,axis.line.x,axis.line.x.top,axis.line.x.bottom,axis.line.y,axis.line.y.left,axis.line.y.right,axis.line.theta,axis.line.r,legend.background,legend.margin,legend.spacing,legend.spacing.x,legend.spacing.y,legend.key,legend.key.size,legend.key.height,legend.key.width,legend.key.spacing,legend.key.spacing.x,legend.key.spacing.y,legend.key.justification,legend.frame,legend.ticks,legend.ticks.length,legend.axis.line,legend.text,legend.text.position,legend.title,legend.title.position,legend.position,legend.position.inside,legend.direction,legend.byrow,legend.justification,legend.justification.top,legend.justification.bottom,legend.justification.left,legend.justification.right,legend.justification.inside,legend.location,legend.box,legend.box.just,legend.box.margin,legend.box.background,legend.box.spacing,panel.background,panel.border,panel.spacing,panel.spacing.x,panel.spacing.y,panel.grid,panel.grid.major,panel.grid.minor,panel.grid.major.x,panel.grid.major.y,panel.grid.minor.x,panel.grid.minor.y,panel.ontop,panel.widths,panel.heights,plot.background,plot.title,plot.title.position,plot.subtitle,plot.caption,plot.caption.position,plot.tag,plot.tag.position,plot.tag.location,plot.margin,strip.background,strip.background.x,strip.background.y,strip.clip,strip.placement,strip.text,strip.text.x,strip.text.x.bottom,strip.text.x.top,strip.text.y,strip.text.y.left,strip.text.y.right,strip.switch.pad.grid,strip.switch.pad.wrap,complete,validate Theme elements and layout settings. Line, rectangle and text elements control appearance; spacing and margins use physical units, positions and justifications control placement. More specific axis, legend, panel, plot and strip elements inherit from their parents. complete marks a complete theme; validate checks element types. Consult ggplot2::theme for each element and its accepted values.
#' @param colour,linewidth,linetype,lineend,color,linejoin,arrow,arrow.fill,inherit.blank,size,fill,family,face,hjust,vjust,angle,lineheight,margin,debug Element styling: colours, fill, line width and type, line ends and joins, arrows, font family/face/size, alignment, angle, line spacing and text margins. inherit.blank follows blank parent elements; debug shows text bounds. size for line/rectangle elements is deprecated upstream; use linewidth.
#' @param t,r,b,l,unit Top, right, bottom and left margins and their physical unit.
#' @param filename,plot,device,path,scale,width,height,units,dpi,limitsize,bg,create.dir Save settings: destination filename/path, plot, device, scale multiplier, width/height and units, raster resolution, size limit, background and permission to create the directory.
#' @param mult,add Multiplicative and additive scale expansion, each a single value or lower/upper pair.
#' @name ggplot2-helpers
NULL

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 aes
#' @export
ggplot2::aes

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 after_stat
#' @export
ggplot2::after_stat

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 after_scale
#' @export
ggplot2::after_scale

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 annotate
#' @export
ggplot2::annotate

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 labs
#' @export
ggplot2::labs

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 ggtitle
#' @export
ggplot2::ggtitle

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 guides
#' @export
ggplot2::guides

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 theme
#' @export
ggplot2::theme

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 element_blank
#' @export
ggplot2::element_blank

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 element_line
#' @export
ggplot2::element_line

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 element_rect
#' @export
ggplot2::element_rect

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 element_text
#' @export
ggplot2::element_text

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 margin
#' @export
ggplot2::margin

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 rel
#' @export
ggplot2::rel

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 ggsave
#' @export
ggplot2::ggsave

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 last_plot
#' @export
ggplot2::last_plot

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 waiver
#' @export
ggplot2::waiver

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 expansion
#' @export
ggplot2::expansion

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 position_identity
#' @export
ggplot2::position_identity

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 PositionIdentity
#' @export
ggplot2::PositionIdentity

#' Physical units and direction arrows
#'
#' These are the original grid functions, available after loading ggchord.
#' Use physical units for theme spacing and arrows for point links.
#' @param x Numeric unit values.
#' @param units Unit names, such as mm, inches or npc.
#' @param data Optional supplementary data for special grid units.
#' @param angle Arrow-head angle in degrees.
#' @param length Arrow-head length as a unit object.
#' @param ends Which ends receive arrows: first, last or both.
#' @param type Arrow-head style: open or closed.
#' @return A grid unit or arrow object.
#' @name grid-helpers
#' @examples
#' unit(2, "mm")
#' arrow(length = unit(2, "mm"), type = "closed")
NULL

#' @rdname grid-helpers
#' @importFrom grid unit
#' @export
grid::unit

#' @rdname grid-helpers
#' @importFrom grid arrow
#' @export
grid::arrow
