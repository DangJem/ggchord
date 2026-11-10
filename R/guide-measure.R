# Per-role text/key requests supplement the actually measured guide-box gtable.
ggchord_measure_guides <- function(plot,boxes) {
  if (!nrow(boxes)) return(data.frame())
  positions <- sub("^guide-box-","",boxes$guide)
  active <- character()
  for(layer in plot$layers) {
    if(identical(layer$show.legend,FALSE))next
    mapping <- setdiff(names(layer$mapping),names(layer$show.legend)[layer$show.legend %in% FALSE])
    if(isTRUE(layer$inherit.aes))mapping<-union(mapping,names(plot$mapping))
    active<-union(active,setdiff(mapping,names(layer$aes_params)))
  }
  result <- list()
  for(scale in plot$scales$scales) {
    aes <- intersect(scale$aesthetics,active)[1]
    if(is.na(aes))next
    guide <- plot$guides$guides[[aes]] %||% scale$guide
    if(identical(guide,'none') || identical(guide,FALSE))next
    role <- if(grepl('^track_',aes))'track' else switch(aes,seq_colour='seq',ribbon_fill='ribbon',
      ribbon_colour='ribbon',ribbon_alpha='ribbon',gene_fill='gene',feature_fill='feature',region_fill='region',NULL)
    if(is.null(role))next
    labels <- scale$get_labels()
    labels <- labels %||% character()
    params <- if(inherits(guide,'Guide'))guide$params else list()
    own_theme <- params$theme
    text <- own_theme$legend.text %||% plot$theme$legend.text
    title_style <- own_theme$legend.title %||% plot$theme$legend.title
    size <- function(style,element,fallback) {
      value <- if(inherits(style,'element_text'))style@size else NULL
      if(is.numeric(value) && !inherits(value,'rel') && length(value)==1 && is.finite(value))value else
        ggchord_theme_point_size(plot,element,fallback)
    }
    family <- if(inherits(text,'element_text'))text@family %||% '' else ''
    face <- if(inherits(text,'element_text'))text@face %||% 'plain' else 'plain'
    font_size <- size(text,'legend.text',8)
    title <- params$title
    if(is.null(title)||inherits(title,'waiver'))title<-scale$name
    if(is.null(title)||inherits(title,'waiver'))title<-plot$labels[[aes]] %||% aes
    measure <- function(value,points,face='plain') {
      g<-grid::textGrob(as.character(value),gp=grid::gpar(fontsize=points,fontfamily=family,fontface=face))
      c(width=grid::convertWidth(grid::grobWidth(g),'mm',TRUE),height=grid::convertHeight(grid::grobHeight(g),'mm',TRUE))
    }
    metrics<-if(!length(labels) || inherits(text,"element_blank")) matrix(c(0,0),nrow=1) else do.call(rbind,lapply(labels,measure,points=font_size,face=face))
    title_metrics<-measure(title,size(title_style,'legend.title',9),
      if(inherits(title_style,'element_text'))title_style@face %||% 'plain' else 'plain')
    requested <- function(field) {
      value<-own_theme[[field]] %||% plot$theme[[field]]
      if(grid::is.unit(value))grid::convertUnit(value,'mm',valueOnly=TRUE) else NA_real_
    }
    spec<-ggchord_role_guide_spec(plot,role,identical(aes,'ribbon_fill'))
    where <- params$position %||% spec$position %||% NA_character_
    if (is.numeric(where)) where <- "inside"
    if (!is.na(where) && !where %in% positions) next
    result[[length(result)+1L]]<-data.frame(role=role,aesthetic=aes,title=as.character(title),n_items=length(scale$get_breaks()),
      max_label_width_mm=max(metrics[,1]),max_label_height_mm=max(metrics[,2]),
      title_width_mm=title_metrics[1],title_height_mm=title_metrics[2],font_size_pt=font_size,font_family=family,
      requested_key_width_mm=requested('legend.key.width'),requested_key_height_mm=requested('legend.key.height'),
      position=where,direction=params$direction %||% spec$direction %||% NA_character_,
      ncol=params$ncol %||% NA_integer_,nrow=params$nrow %||% NA_integer_,
      shared_height_budget_mm=spec$height_budget %||% NA_real_,
      measurement=if(is.expression(labels))'estimated_expression' else 'device_text_metrics',
      font_fallback='not_inspected',labels=I(list(as.character(labels))),stringsAsFactors=FALSE)
  }
  ggchord_rbind_fill(result) %||% data.frame()
}
