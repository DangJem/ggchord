# Targeted v0.12 visual acceptance. Outputs never replace website figures.
pkgload::load_all(quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else tempfile("ggchord-v012-")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
s <- data.frame(accver = c("A", "B", "C"), length = 1000)
line <- data.frame(qaccver = "A", saccver = c("B", "C"), qpos = 200, spos = 500)
ribbon <- transform(line, qstart = 100, qend = 400, sstart = 600, send = 400)
feature <- data.frame(accver = "A", start = c(100,300,500,700),
  end = c(200,400,600,800), strand = c("-","-","+","+"),
  type = c("CDS","tRNA","repeat","promoter"))
shape <- c(CDS="arrow",tRNA="block","repeat"="chevron",promoter="lollipop")
base <- ggchord(s) + geom_seq()
plots <- list(
  line_query = base + geom_link_line(data=line,link_branch="query",arrow=arrow(length=unit(2,"mm"))),
  line_subject = base + geom_link_line(data=transform(line,qaccver=saccver,saccver=qaccver,qpos=spos,spos=qpos),link_branch="subject",arrow=arrow(ends="both",length=unit(2,"mm"))),
  ribbon_query = ggchord(s,ribbon)+geom_seq()+geom_link_ribbon(fill="#3B90B6",link_branch="query"),
  feature_keys = base+geom_feature(aes(feature_shape=type),data=feature)+scale_feature_shape_manual(values=shape),
  smooth = ggchord(s,ribbon)+geom_seq()+geom_feature(aes(feature_shape=type),data=feature)+scale_feature_shape_manual(values=shape)+geom_link_ribbon(fill="#3B90B6",link_avoid="smooth"),
  uniform = ggchord(s,ribbon)+geom_seq()+geom_feature(aes(feature_shape=type),data=feature)+scale_feature_shape_manual(values=shape)+geom_link_ribbon(fill="#3B90B6",link_avoid="uniform")
)
for (name in names(plots)) {
  for (ext in c("png","pdf","svg")) {
    ggplot2::ggsave(file.path(out,paste0(name,".",ext)),plots[[name]],width=8,height=6,dpi=120)
  }
}
writeLines(capture.output(sessionInfo()),file.path(out,"session.txt"))
message("Acceptance outputs: ",normalizePath(out))
