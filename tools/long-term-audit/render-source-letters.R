# Four independent source excerpts, not a multiple sequence alignment.
source("tools/long-term-audit/semantic-prototypes.R")
d <- utils::read.csv("tools/long-term-audit/v0.20.0-source-letters.csv",stringsAsFactors=FALSE)
d$row <- paste0(d$accver," (",d$sequence_unit,")")
p <- ggplot2::ggplot(d,ggplot2::aes(source_position,row)) +
  ggplot2::geom_tile(fill="#E5EDF2",colour="white",height=.7) +
  ggplot2::geom_text(ggplot2::aes(label=symbol),family="mono",size=3) +
  ggplot2::facet_wrap(~alphabet,ncol=1,scales="free") +
  ggplot2::scale_x_continuous(breaks=scales::breaks_pretty(n=8),expand=ggplot2::expansion(add=.5)) +
  ggplot2::labs(x="Original source position (1-based)",y=NULL,
    title="Actual source characters: independent excerpts",
    subtitle="Two 60 bp DNA excerpts and two 40 aa protein excerpts; columns do not imply homology") +
  ggplot2::theme_minimal(base_size=11) +
  ggplot2::theme(panel.grid=ggplot2::element_blank(),strip.text=ggplot2::element_text(face="bold"))
ggplot2::ggsave("tools/long-term-audit/v0.20.0-source-letters.png",p,width=12,height=4,dpi=120,bg="white")
