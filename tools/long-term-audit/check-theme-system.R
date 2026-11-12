# Targeted theme/Track checks, not a package/release test matrix.
pkgload::load_all(quiet=TRUE)
library(dplyr)
out <- 'tools/long-term-audit'
fail <- function(expr,pattern) {
  e <- tryCatch({force(expr);NULL},error=identity)
  stopifnot(inherits(e,'error'),grepl(pattern,conditionMessage(e)))
}
data(plasmid_example_pUC19)
sequence <- plasmid_example_pUC19
# Two genuine statistics of the unchanged sequence, not fabricated coverage.
starts <- seq.int(1L,sequence$length,100L)
bins <- tibble(accver=sequence$accver,start=starts,end=pmin(starts+99L,sequence$length)) |>
  rowwise() |>
  mutate(fragment=substr(sequence$sequence,start,end),
    gc=mean(strsplit(fragment,'',fixed=TRUE)[[1]] %in% c('G','C')),
    adenine=sum(strsplit(fragment,'',fixed=TRUE)[[1]]=='A'),pos=(start+end)/2) |>
  ungroup()
gc_band <- position_track('GC',width=unit(5,'mm'))
a_band <- position_track('adenine',order=2,width=unit(5,'mm'))
p <- ggchord(sequence) + geom_seq() +
  geom_track(aes(x=pos,y=gc,colour=gc),data=bins,geom='line',track=gc_band) +
  geom_track(aes(x=pos,y=gc),data=bins,geom='point',track=gc_band,colour='orange',size=.8) +
  geom_track(aes(x=pos,y=adenine,fill=adenine),data=bins,geom='col',track=a_band) +
  scale_track_value_continuous('GC',limits=c(0,1),unit='fraction',
    guide=ggplot2::guide_axis(theme=theme(axis.text=element_text(colour='#E5A5C0',size=12),
      axis.title=element_text(colour='#E5A5C0',face='bold')))) +
  scale_track_value_continuous('adenine',limits=c(0,100),unit='bases per bin') +
  scale_track_colour_gradientn('GC',colours=c('#AECEE8','#59B59C'),name='GC fraction') +
  scale_track_fill_gradientn('adenine',colours=c('#CCAE65','#F9E4AB'),name='A bases per bin') +
  coord_circular() + theme_ggchord_dark() +
  theme(ggchord.axis.text.orientation='horizontal',
    ggchord.track.axis.text=element_text(size=10),
    ggchord.track.axis.title=element_text(face='italic'),legend.position='bottom')
original <- p$ggchord$data
report <- export_ggchord_layout(p)
axes <- filter(report$axis,component=='track_axis')
stopifnot(setequal(unique(axes$track_id),c('GC','adenine')),
  all(filter(axes,track_id=='GC',.component=='track_axis_text')$.axis_colour=='#E5A5C0'),
  all(filter(axes,track_id=='adenine',.component=='track_axis_text')$.axis_colour=='#F1F3F5'),
  all(filter(axes,track_id=='GC',.component=='track_axis_text')$.axis_size >
    filter(axes,track_id=='adenine',.component=='track_axis_text')$.axis_size[1]),
  nrow(filter(axes,.component=='track_axis_title'))==2)
# A standard theme override and the old convenience spelling reach one setting.
a <- p+theme_ggchord(axis.gap=unit(3,'mm'),axis.text.orientation='horizontal')
b <- p+theme_ggchord()+theme(ggchord.axis.gap=unit(3,'mm'),ggchord.axis.text.orientation='horizontal')
stopifnot(identical(ggchord:::ggchord_plot_settings(a)$axis,ggchord:::ggchord_plot_settings(b)$axis))
c <- a+theme(ggchord.axis.gap=unit(1,'mm'))
stopifnot(as.numeric(ggchord:::ggchord_plot_settings(c)$axis$gap)==1,
  as.numeric(ggchord:::ggchord_plot_settings(a)$axis$gap)==3)
# Later global styles preserve the first Track's local Guide override.
q <- p+theme(ggchord.track.axis.text=element_text(colour='cyan',size=18))
qr <- export_ggchord_layout(q)
qaxes <- filter(qr$axis,component=='track_axis',.component=='track_axis_text')
stopifnot(all(filter(qaxes,track_id=='GC')$.axis_colour=='#E5A5C0'),
  all(filter(qaxes,track_id=='adenine')$.axis_colour=='cyan'))
hidden <- p+theme(ggchord.track.axis.title=element_blank())
h <- export_ggchord_layout(hidden)$axis
# GC has an explicit local title; the other Track inherits the blank.
stopifnot(all(filter(h,component=='track_axis',track_id=='adenine',.component=='track_axis_title')$.axis_hidden),
  !any(filter(h,component=='track_axis',track_id=='GC',.component=='track_axis_title')$.axis_hidden))
fail(scale_track_value_continuous('GC',guide=ggplot2::guide_axis(n.dodge=2)),'unsupported option')
fail(export_ggchord_layout(p+theme(ggchord.axis.gap=unit(1,'npc'))),'absolute physical')
# Real target-device draw: styles and measured text bounds share the final font.
ggplot2::ggsave(file.path(out,'v0.27.0-theme-tracks.png'),p,width=9,height=7,dpi=120)
stopifnot(identical(original,p$ggchord$data))
# All presets reuse the same vocabulary; source tables and statistics stay identical.
presets <- list(default=theme_ggchord,minimal=theme_ggchord_minimal,
  publication=theme_ggchord_publication,dark=theme_ggchord_dark,plasmid=theme_ggchord_plasmid)
for (name in names(presets)) {
  x <- export_ggchord_layout(p+presets[[name]]())
  stopifnot(identical(x$track_statistics,report$track_statistics),
    identical(x$track_source_data,report$track_source_data))
}
# Representative shared decoration in existing multi-sequence views.
data(seq_data_example)
seqs <- seq_data_example[1:2,c('accver','length')]
d <- data.frame(accver=rep(seqs$accver,each=2),pos=rep(c(1,100),2),value=c(1,2,2,3))
view <- ggchord(seqs)+geom_seq()+geom_track(aes(x=pos,y=value),data=d,geom='line',
  track=position_track('observed'))+scale_track_value_continuous('observed',name='Declared control')+
  theme_ggchord_dark()+theme(ggchord.track.axis.text=element_text(colour='cyan'))
views <- list(chord=coord_chord(),collinear=coord_collinear(),radial=coord_radial_synteny(),sankey=coord_sankey())
for (name in names(views)) {
  # Values in this last block are explicitly synthetic style/coordinate controls.
  r <- export_ggchord_layout(view+views[[name]])
  stopifnot(nrow(filter(r$axis,component=='track_axis'))>0)
}
jsonlite::write_json(list(status='passed targeted theme and Track checks; not release acceptance',
  presets=names(presets),views=names(views),track_ids=unique(axes$track_id),
  source='unchanged pUC19 sequence; GC fraction and adenine count per 100bp bin',
  axes=axes[,intersect(c('track_id','label','.component','.axis_colour','.axis_size','.axis_hidden'),names(axes))],
  semantics='theme changes decoration; independent named scales preserve domains and source statistics'),
  file.path(out,'v0.27.0-theme-results.json'),pretty=TRUE,auto_unbox=TRUE)
cat('PASS theme inheritance, scoped Track guides, source preservation and representative views\n')
