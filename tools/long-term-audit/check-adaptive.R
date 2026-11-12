# Targeted development checks; synthetic category labels are explicitly labelled.
# No testthat, package checks, vignette rebuild or release device matrix.
pkgload::load_all(quiet=TRUE)
library(dplyr)
output <- 'tools/long-term-audit'
policy_plot <- function(p, width, height, file=NULL) {
  path <- file %||% tempfile(fileext='.png')
  owned <- is.null(file)
  previous <- grDevices::dev.cur()
  ggchord:::ggchord_diagnostic_device(path,'png',width,height,120)
  on.exit({grDevices::dev.off(); if(previous %in% grDevices::dev.list()) grDevices::dev.set(previous)
    if(owned) unlink(path)},add=TRUE)
  table <- ggplot2::ggplotGrob(p)
  if (!owned) {grid::grid.newpage(); grid::grid.draw(table)}
  list(context=attr(table,'ggchord_render_context'),layout=attr(table,'ggchord_render_layout'))
}
`%||%` <- function(x,y) if(is.null(x)) y else x
# These are synthetic contract inputs, not biological observations.
seqs <- data.frame(accver=paste0('Sequence-',1:10),length=1000)
features <- data.frame(accver=seqs$accver,start=100,end=300,strand='+',type='gene',
  category=paste('Feature category',1:10))
hits <- data.frame(qaccver=seqs$accver[1:9],saccver=seqs$accver[2:10],
  qstart=100,qend=200,sstart=300,send=400,pident=seq(80,100,length.out=9))
base <- ggchord(seqs,hits,feature_data=features)+geom_seq()+
  geom_feature(aes(feature_fill=category))+geom_link_ribbon()
p <- base+adapt_ggchord()
joint <- policy_plot(p,6,4,file=file.path(output,'v0.24.0-adaptive-guides.png'))
guides <- joint$context$guide_layout
stopifnot(guides$status=='feasible',!any(guides$bounds$overflow),
  identical(sort(guides$guides$n_items),c(6L,10L,10L)),
  guides$guides$type[guides$guides$title=='Identity (%)']=='continuous',
  identical(p$ggchord$data$seq_data,seqs),identical(p$ggchord$data$ribbon_data,hits),
  is.null(base$ggchord$layout_policy))
# A same-target redraw must choose the same public layout; it must not rescale
# an already scaled build or mutate the user plot when another target is drawn.
again <- policy_plot(p,6,4)
stopifnot(identical(guides$guides,again$context$guide_layout$guides))
real <- ggchord(seq_data_example,ribbon_data_example)+geom_seq()+geom_link_ribbon()
fixed <- policy_plot(real+adapt_ggchord('fixed'),6,4)
scaled <- policy_plot(real+adapt_ggchord('scale'),10,8,
  file=file.path(output,'v0.24.0-adaptive-real.png'))
stopifnot(fixed$context$adaptive_layout$text_factor==1,
  scaled$context$adaptive_layout$text_factor>1,
  scaled$context$adaptive_layout$text_factor<=1.35,
  identical(export_ggchord_layout(real,include='ribbon',original_data=TRUE)$original_data,
    export_ggchord_layout(real+adapt_ggchord(),include='ribbon',original_data=TRUE)$original_data))
locked <- real+geom_seq_label(size=4)+
  theme_ggchord(seq.label=ggplot2::element_text(size=18))+
  list(ggplot2::theme(legend.text=ggplot2::element_text(size=12),
      title=ggplot2::element_text(size=13)),
    ggplot2::theme(legend.text=ggplot2::element_text(colour='red')))+
  ggplot2::guides(seq_colour=ggplot2::guide_legend(ncol=2,position='bottom',direction='horizontal'))+
  adapt_ggchord()
d <- diagnose_ggchord(locked,width=6,height=4,device='png')
g <- filter(d$guide_layout$guides,aesthetic=='seq_colour')
s <- filter(d$adaptive_layout$styles,element %in% c('ggchord.seq.label','legend.text','legend.title'))
labels <- filter(d$text, label %in% seq_data_example$accver)
stopifnot(g$final_position=='bottom',g$columns==2,g$direction=='horizontal',
  nrow(s)==3L,nrow(labels)==nrow(seq_data_example),
  all(s$user_locked),all(s$effective==s$requested),
  all(labels$size_mm==4))
# Shared rendering entry points: a single representative per backbone/axis type.
views <- list(collinear=coord_collinear(),sankey=coord_sankey(),radial=coord_radial_synteny())
view_results <- lapply(views,function(coord) {
  x <- policy_plot(real+coord+adapt_ggchord(),8,6)
  x$context$guide_layout$status
})
dot <- ggchord(seq_data_example,ribbon_data_example)+geom_alignment_segment()+coord_dotplot()+adapt_ggchord()
view_results$dotplot <- policy_plot(dot,8,6)$context$guide_layout$status
circle <- ggchord(plasmid_example_pUC19)+geom_seq()+coord_circular()+adapt_ggchord()
view_results$circular <- policy_plot(circle,6,4)$context$guide_layout$status
stopifnot(all(unlist(view_results) %in% c('feasible','no_guides')))
# Fusion rebuilds the policy inside each allocated panel, in its own mm context.
region <- data.frame(accver='pUC19',start=101,end=1900)
focused <- focus_ggchord_data(plasmid_example_pUC19,loci=region)
detail <- ggchord(focused$seq_data)+geom_seq()+coord_collinear()+adapt_ggchord()
fusion <- export_ggchord_fusion(magnify_ggchord(circle,region,detail),12,8)
stopifnot(nrow(fusion$leaders)==2L,
  all(vapply(fusion$diagnostics,function(x) length(x$adaptive_layout$canvas_mm)==2L,logical(1))),
  all(fusion$focus_regions$start==101),all(fusion$focus_regions$end==1900))
# A locked guide wider than the target must be reported as infeasible.
impossible <- ggchord(seq_data_example)+geom_seq()+
  ggplot2::guides(seq_colour=ggplot2::guide_legend(nrow=1,position='bottom',
    theme=ggplot2::theme(legend.key.width=grid::unit(100,'mm'))))+adapt_ggchord()
small <- suppressWarnings(policy_plot(impossible,4,4))
stopifnot(small$context$guide_layout$status=='insufficient_space',
  any(small$context$guide_layout$bounds$overflow))
jsonlite::write_json(list(joint=guides,scale=scaled$context$adaptive_layout,
  explicit_constraints=g,coordinate_views=view_results,fusion_panel_sizes=lapply(fusion$diagnostics,
    function(x)x$adaptive_layout$canvas_mm),infeasible_reported=TRUE,
  source_data_preserved=TRUE,same_target_stable=TRUE,no_release_matrix=TRUE),
  file.path(output,'v0.24.0-adaptive-results.json'),auto_unbox=TRUE,pretty=TRUE)
cat('PASS: joint guides, explicit policies, source preservation and panel contexts\n')
