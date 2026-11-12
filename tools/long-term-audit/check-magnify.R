# Direct v0.22.0 source/physical checks, not a package/release test suite.
pkgload::load_all(quiet=TRUE)
library(dplyr)
out <- 'tools/long-term-audit'
started <- proc.time()[['elapsed']]
expect_error <- function(expr,pattern) {
  e <- tryCatch({force(expr);NULL},error=identity)
  stopifnot(inherits(e,'error'),grepl(pattern,conditionMessage(e)))
}
plasmid <- plasmid_example_pETDuet_1
features <- find_common_features(plasmid)
region <- data.frame(accver=plasmid$accver,start=101,end=1900)
local <- focus_ggchord_data(plasmid,feature_data=filter(features,start<=end),loci=region)
primers <- find_primer_bindings(plasmid,set='universal') |> filter(start<=end)
primers <- focus_ggchord_data(plasmid,feature_data=primers,loci=region)$feature_data
sites <- find_restriction_sites(plasmid) |> filter_restriction_sites(enzymes=c('EcoRI','BamHI','HindIII'))
sites <- focus_ggchord_data(plasmid,feature_data=sites,loci=region)$feature_data
site_names <- sites |> summarise(accver=first(accver),start=min(start),end=max(end),strand='.',anno=paste(unique(enzyme),collapse=' / '))
labels <- local$feature_data |> filter(type %in% c('CDS','rep_origin'),anno!='ATG')
whole <- ggchord(plasmid,feature_data=features)+geom_seq(show.legend=FALSE)+
  geom_feature(position=position_feature_stack(base_position=position_plasmid()),show.legend=FALSE)+
  geom_seq_center_label()+coord_circular(rotation=-30)
detail <- ggchord(local$seq_data,feature_data=local$feature_data)+
  geom_seq(seq_orientation=-1,show.legend=FALSE)+
  geom_feature(position=position_feature_stack(base_position=position_strand(.1)),show.legend=FALSE)+
  geom_feature_label(data=labels,label_orientation='horizontal',label_overlap='allow',colour='#202020',size=3)+
  geom_feature(data=primers,position=position_strand(.22),width=.035,colour='#7C3A94',show.legend=FALSE)+
  geom_feature_label(aes(feature_label=name),data=primers,position=position_strand(.22),label_orientation='horizontal',label_overlap='allow',colour='#202020',size=3)+
  geom_feature(data=sites,shape='block',position=position_strand(.34),width=.035,colour='#B45309',show.legend=FALSE)+
  geom_feature_label(aes(feature_label=anno),data=site_names,position=position_strand(.34),label_orientation='horizontal',label_overlap='allow',label_radial_offset=.3,colour='#202020',size=3)+
  coord_collinear(seq_axis='vertical')+
  scale_seq_position_continuous(coordinates='source',breaks=c(101,500,1000,1500,1900))
view <- magnify_ggchord(whole,region,detail)
before <- local$seq_data
check <- function(e,bounds=c(100,1900)) {
  stopifnot(nrow(e$focus_regions)==1L,nrow(e$leaders)==2L,
    all(e$leaders$semantic=='decorative_focus_leader'),
    all(is.finite(unlist(e$leaders[c('x','y','xend','yend')]))),
    e$focus_regions$stat_scope=='explicit_detail_plot',
    is.finite(e$focus_regions$length_ratio),e$focus_regions$length_ratio>0)
  for(id in c('center','detail')) {
    f <- e$focus_fronts[e$focus_fronts$panel_id==id,]
    stopifnot(isTRUE(all.equal(f$position_boundary[c(1,nrow(f))],bounds,check.attributes=FALSE)))
  }
}
e <- export_ggchord_fusion(view,12,8); check(e)
stopifnot(nrow(e$connections)==0L,length(e$relationship_data)==0L,
  e$focus_regions$detail_origin==101,nchar(before$sequence)==1800,
  identical(e$focus_region_data,region))
e2 <- export_ggchord_fusion(view,10,7);check(e2)
stopifnot(identical(e$focus_regions[c('accver','start','end')],e2$focus_regions[c('accver','start','end')]),
  !isTRUE(all.equal(e$leaders$x,e2$leaders$x)),identical(before,local$seq_data))
reflected <- export_ggchord_fusion(magnify_ggchord(whole+ggplot2::scale_x_reverse(),region,detail),12,8)
check(reflected)
a <- e$focus_fronts[e$focus_fronts$panel_id=='center',]
b <- reflected$focus_fronts[reflected$focus_fronts$panel_id=='center',]
stopifnot(sum(diff(a$x)*diff(b$x))<0)
# Original accession and origin survive a renamed focus window.
renamed <- focus_ggchord_data(plasmid,loci=transform(region,locus_id='local-window'))
p_renamed <- ggchord(renamed$seq_data)+geom_seq(show.legend=FALSE)+coord_collinear()
er <- export_ggchord_fusion(magnify_ggchord(whole,region,p_renamed,side='bottom'),10,8)
check(er);stopifnot(er$focus_regions$detail_instance=='local-window')
# Four actual sequences, an explicit repeated display instance and a real HSP
# remain distinct from the decorative source-region mapping.
seqs <- seq_data_example[c('accver','length')]
r <- data.frame(accver=seqs$accver[1],start=101,end=1900)
f <- focus_ggchord_data(seqs,loci=r)
p_detail <- ggchord(f$seq_data)+geom_seq(show.legend=FALSE)+coord_collinear()
p_chord <- ggchord(seqs)+geom_seq(show.legend=FALSE)+coord_chord()
chord <- magnify_ggchord(p_chord,r,p_detail)
ec <- export_ggchord_fusion(chord,10,8);check(ec)
p_other <- ggchord(seqs)+geom_seq(show.legend=FALSE)+coord_collinear()
combined <- chord+fusion_panel(p_other,'other','top',size=c(.1,.07))+
  fusion_ribbon(ribbon_data_example[1,,drop=FALSE],qpanel='center',spanel='other',fill='grey70')
# Use a compact center to keep the third explicit panel separate.
combined$panels$center$size <- c(.4,.7)
eboth <- export_ggchord_fusion(combined,12,8);check(eboth)
stopifnot(nrow(eboth$connections)==1L,length(eboth$relationship_data)==1L,
  identical(eboth$relationship_data[[1]],ribbon_data_example[1,,drop=FALSE]))
repeat_source <- ggchord(seqs)+geom_seq(seq_order=setNames(c(seqs$accver[1],seqs$accver[2],seqs$accver[1]),c('A1','B','A2')),show.legend=FALSE)+coord_chord()
expect_error(export_ggchord_fusion(magnify_ggchord(repeat_source,r,p_detail)),'ambiguous')
selected <- magnify_ggchord(repeat_source,transform(r,instance='A2'),p_detail)
stopifnot(export_ggchord_fusion(selected)$focus_regions$source_instance=='A2')
# A single real source position retains closed-cell width.
single <- data.frame(accver=plasmid$accver,start=300,end=300)
f1 <- focus_ggchord_data(plasmid,loci=single)
p1 <- ggchord(f1$seq_data)+geom_seq(show.legend=FALSE)+coord_collinear()
es <- export_ggchord_fusion(magnify_ggchord(whole,single,p1),10,8);check(es,c(299,300))
proteins <- read_fasta_lengths('examples/protein/four-proteins.faa',sequence_unit='aa')
pr <- data.frame(accver=proteins$accver[1],start=21,end=60)
pf <- focus_ggchord_data(proteins,loci=pr)
pc <- ggchord(proteins)+geom_seq(show.legend=FALSE)+coord_chord()
pd <- ggchord(pf$seq_data)+geom_seq(show.legend=FALSE)+coord_collinear()
ep <- export_ggchord_fusion(magnify_ggchord(pc,pr,pd),10,8);check(ep,c(20,60))
stopifnot(ep$focus_regions$sequence_unit=='aa')
no_leader <- export_ggchord_fusion(magnify_ggchord(whole,region,detail,connector='none'),12,8)
stopifnot(nrow(no_leader$leaders)==0L,nrow(no_leader$focus_fronts)>0)
expect_error(magnify_ggchord(whole,transform(region,start=1900,end=101),detail),'ordered positive')
expect_error(export_ggchord_fusion(magnify_ggchord(whole,transform(region,accver='unknown'),detail)),'no fully visible')
expect_error(magnify_ggchord(whole,region,ggplot2::ggplot()),'ggchord plots')
expect_error(export_ggchord_fusion(magnify_ggchord(whole,region,p1)),'no fully visible')
ggplot2::ggsave(file.path(out,'v0.22.0-magnify.png'),view,width=12,height=8,dpi=120,bg='white')
utils::write.csv(e$focus_regions,file.path(out,'v0.22.0-focus-regions.csv'),row.names=FALSE)
utils::write.csv(e$focus_fronts,file.path(out,'v0.22.0-focus-fronts.csv'),row.names=FALSE)
utils::write.csv(e$leaders,file.path(out,'v0.22.0-leaders.csv'),row.names=FALSE)
jsonlite::write_json(list(version='0.22.0.9000',region=region,feature_records=nrow(local$feature_data),
  primer_records=nrow(primers),recognition_records=nrow(sites),source_bounds=c(100,1900),
  exact_source_mapping=TRUE,independent_detail_and_statistics=TRUE,biological_links=nrow(e$connections),
  decorative_leaders=nrow(e$leaders),sampled_length_ratio=e$focus_regions$length_ratio,
  measured_sizes_inches=list(c(12,8),c(10,7)),native_scale_reflection=TRUE,
  renamed_focus_and_instance_resolution=TRUE,real_hsp_and_focus_kept_separate=TRUE,
  native_aa_and_single_position=TRUE,source_data_unchanged=TRUE,
  elapsed_seconds=unname(proc.time()[['elapsed']]-started),
  limits=c('explicit detail plot and statistics','Circular/Chord source and Collinear detail only',
    'one ordered source region','no global collision guarantee','no release device matrix')),
  file.path(out,'v0.22.0-magnify-results.json'),pretty=TRUE,auto_unbox=TRUE)
message('Magnify: real source windows, exact fronts, scale reflection, independent HSPs, units and failures passed')
