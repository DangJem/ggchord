# Direct F2/M3 checks; not the publication-device or package test matrix.
pkgload::load_all(quiet=TRUE)
library(dplyr)
out <- 'tools/long-term-audit'
fail <- function(expr,pattern) {
  e <- tryCatch({force(expr);NULL},error=identity)
  stopifnot(inherits(e,'error'),grepl(pattern,conditionMessage(e)))
}
data(seq_data_example);data(ribbon_data_example)
seqs <- seq_data_example |> select(accver,length)
hit <- ribbon_data_example |>
  filter(pmin(qstart,qend,sstart,send)>1) |> slice_head(n=1)
stopifnot(nrow(hit)==1)
# This is one real observed HSP displayed repeatedly, never eight observations.
hit$source_hsp_id <- 'explicit-real-HSP-1'
base <- ggchord(seqs)+geom_seq(show.legend=FALSE)
whole <- base+coord_chord()
targets <- list(right=base+coord_collinear(seq_axis='vertical'),
  top=base+coord_sankey(),left=base+coord_radial_synteny(start_angle=35),
  bottom=base+coord_collinear(),topright=base+coord_sankey(),
  topleft=base+coord_radial_synteny(),bottomright=base+coord_collinear(),
  bottomleft=base+coord_collinear(seq_axis='vertical'))
fusion <- ggchord_fusion(whole)
source_records <- list()
for(id in names(targets)) {
  fusion <- fusion+fusion_panel(targets[[id]],id,side=id)
  row <- hit
  if(isTRUE(targets[[id]]$coordinates$ggchord_sankey)) {
    instances <- get_chord_layout(targets[[id]])$display_instances
    row$s_instance <- instances$instance_id[instances$accver==row$saccver & instances$side=='right']
  }
  fusion <- fusion+fusion_ribbon(row,qpanel='center',spanel=id,fill='#548CA0',alpha=.35,
    allow_repeated=TRUE,show.legend=FALSE)
  source_records[[id]] <- row
}
advanced <- export_ggchord_fusion(fusion,16,13)
stopifnot(nrow(advanced$panels)==9,nrow(advanced$connections)==8,
  nrow(advanced$source_observations)==1,advanced$source_observations$displayed_connections==8)
for(i in seq_along(source_records)) stopifnot(identical(advanced$relationship_data[[i]],source_records[[i]]))
# Independent public source-position mapping and final measured native scales.
for(id in names(targets)) {
  front <- filter(advanced$fronts,panel_id==id,endpoint=='s')
  panel <- filter(advanced$panels,panel_id==id)
  context <- advanced$diagnostics[[id]]$render_context
  position <- map_ggchord_position(targets[[id]],unique(front$instance_id),front$position_boundary,source=TRUE)
  x <- panel$panel_xmin+(position$x-context$x_range[1])/diff(context$x_range)*(panel$panel_xmax-panel$panel_xmin)
  y <- panel$panel_ymin+(position$y-context$y_range[1])/diff(context$y_range)*(panel$panel_ymax-panel$panel_ymin)
  stopifnot(all(abs(front$x-x)<.02),all(abs(front$y-y)<.02))
}
# Missing the exact mirrored instance must not select the first side.
ambiguous <- ggchord_fusion(whole)+fusion_panel(targets$top,'mirror')+
  fusion_ribbon(hit,qpanel='center',spanel='mirror',fill='grey')
fail(export_ggchord_fusion(ambiguous),'ambiguous')
# Child orientation stays its own Coord choice; angle only allocates a container.
angle <- ggchord_fusion(whole)+fusion_panel(base+coord_collinear(),'angle',angle=31,distance=.8,size=c(.16,.22))
angle_export <- export_ggchord_fusion(angle,12,8)
centre <- function(d) c(mean(c(d$container_xmin,d$container_xmax)),mean(c(d$container_ymin,d$container_ymax)))
v <- centre(angle_export$panels[2,])-centre(angle_export$panels[1,])
stopifnot(angle_export$panels$requested_angle[2]==31,angle_export$panels$coordinate[2]=='collinear',
  abs(atan2(v[2],v[1])*180/pi-31)<.01)
physical <- ggchord_fusion(whole,center_box=unit(c(5,85,10,100),'mm'))+
  fusion_panel(base+coord_collinear(),'physical',box=unit(c(100,210,10,100),'mm'))
pe <- export_ggchord_fusion(physical,12,6)
stopifnot(all(abs(c(pe$panels$container_xmin[2],pe$panels$container_xmax[2])-c(100,210))<.01))
fail(ggchord_fusion(whole)+fusion_panel(base,'bad',box=c(.8,1.2,.1,.3)),'ordered nonnegative')
fail(ggchord_fusion(whole)+fusion_panel(base,'bad',box=c(.1,.2,.3,.4),side='right'),'box alone')
# Root/child containers, explicit original metadata and exact nested routes.
genome <- read_fasta_lengths('examples/format-inputs/emales-BVI_023A.fna',keep_sequence=TRUE)
coding <- read_gff3('examples/format-inputs/emales-BVI_023A.gff3') |> filter(feature_id=='BVI_023A_0009')
translation <- translate_ggchord_cds(build_ggchord_cds(coding,genetic_code=11),genome,initiation='start')
letters <- translation$nucleotides |> filter(alignment_column<=60)
region <- letters |> summarise(accver=first(source_accver),start=min(source_position),end=max(source_position))
local_region <- coding |> select(accver,start,end)
focused <- focus_ggchord_data(genome,loci=local_region)
local <- ggchord(focused$seq_data)+geom_seq(show.legend=FALSE)+coord_collinear()
character <- ggplot2::ggplot(letters,aes(alignment_column,accver,label=symbol))+
  geom_residue(aes(fill=symbol),overflow='error',size=2.7,show.legend=FALSE)+
  ggplot2::scale_fill_manual(values=c(A='#B8D8BA',C='#A6CBE3',G='#F2D28C',T='#EAA8A8'))+
  ggplot2::theme_minimal()+ggplot2::labs(x='Observed coding base',y=NULL,title='60 real negative-strand CDS bases')
character <- register_ggchord_alignment(character,letters)
inner <- magnify_ggchord(local,region,character,
  center_box=c(.02,.2,.2,.8),detail_box=c(.27,.99,.15,.85))
macro <- ggchord(genome |> select(accver,length),feature_data=coding)+
  geom_seq(show.legend=FALSE)+geom_gene()+coord_collinear()
outer_region <- local_region |> mutate(detail_panel='detail/center',purpose='Complete source CDS')
nested <- magnify_ggchord(macro,outer_region,inner,
  center_box=c(.02,.22,.2,.8),detail_box=c(.27,.99,.1,.9))
nested_export <- export_ggchord_fusion(nested,16,6)
stopifnot(nrow(nested_export$leaders)==4,
  setequal(nested_export$focus_regions$composition_id,c('root','detail/root')),
  identical(nested_export$focus_region_data_by_composition$root,outer_region),
  identical(nested_export$focus_region_data_by_composition[['detail/root']],region))
# Both scopes' source boundaries and exact source origins survive nesting.
for(scope in c('root','detail/root')) {
  original <- nested_export$focus_region_data_by_composition[[scope]]
  resolved <- filter(nested_export$focus_regions,composition_id==scope)
  stopifnot(isTRUE(all.equal(resolved$start,original$start)),isTRUE(all.equal(resolved$end,original$end)))
}
# Repeated declared evidence across nesting requires all layers to allow it.
one <- hit; one$source_hsp_id <- 'same-across-nesting'
child <- ggchord_fusion(whole)+fusion_panel(base+coord_collinear(),'child')+
  fusion_ribbon(one,qpanel='center',spanel='child',fill='grey')
parent <- ggchord_fusion(whole)+fusion_panel(child,'nested')+
  fusion_ribbon(one,qpanel='center',spanel='nested/child',fill='grey',allow_repeated=TRUE)
fail(export_ggchord_fusion(parent,14,9),'across nested')
# Explicit overlaps use disjoint clips; alpha is not accumulated by repainting.
cover <- ggchord_fusion(whole,allow_overlap=TRUE,center_box=c(.05,.8,.1,.9))+
  fusion_panel(base+coord_collinear(),'inset',box=c(.45,.95,.35,.65))
covered <- export_ggchord_fusion(cover,12,8)
stopifnot(nrow(covered$overlaps)==1,covered$overlaps$over=='inset')
fail(export_ggchord_fusion(ggchord_fusion(whole,center_box=c(.05,.8,.1,.9))+
  fusion_panel(base+coord_collinear(),'inset',box=c(.45,.95,.35,.65)),12,8),'overlap')
# Three required M3 sizes use real source windows; character panels retain their own overflow policy.
r <- hit |> transmute(accver=qaccver,start=pmin(qstart,qend),end=pmax(qstart,qend))
w <- focus_ggchord_data(seqs,loci=r)
detail <- ggchord(w$seq_data)+geom_seq(show.legend=FALSE)+coord_collinear()
m <- magnify_ggchord(whole,r,detail)
for(size in list(c(4,4),c(8,8),c(12,8))) {
  e <- export_ggchord_fusion(m,size[1],size[2])
  stopifnot(nrow(e$leaders)==2,e$focus_regions$start==r$start,e$focus_regions$end==r$end)
}
ggplot2::ggsave(file.path(out,'v0.27.0-fusion-eight.png'),fusion,width=16,height=13,dpi=110)
ggplot2::ggsave(file.path(out,'v0.27.0-m3-nested.png'),nested,width=16,height=6,dpi=120)
jsonlite::write_json(list(status='passed directed F2/M3 checks; global routing and release matrices remain separate',
  eight_panels=advanced$panels |> select(-requested_box),
  source_observations=advanced$source_observations,
  nested_regions=nested_export$focus_regions,nested_scope_names=names(nested_export$focus_region_data_by_composition),
  original_nested_region=region,overlaps=covered$overlaps,
  source='bundled real HSP and unchanged EMALE genome/CDS; no new inferred alignment'),
  file.path(out,'v0.27.0-fusion-advanced-results.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
cat('PASS exact derived-view anchors, eight directions/free placement, nested source identity and explicit overlap\n')
