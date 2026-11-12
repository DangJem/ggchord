# Reuse original canonical preparations without replacing data or geometry.
load_linear_lesson <- function() {
  lines <- readLines('tools/documentation-examples.R',warn=FALSE)
  markers <- grep('^# ---- [a-z]+ ----$',lines)
  chunks <- setNames(lapply(seq_len(length(markers)-1L),function(i)
    lines[seq.int(markers[i]+1L,markers[i+1L]-1L)]),
    sub('^# ---- ([a-z]+) ----$','\\1',lines[markers[-length(markers)]]))
  environment <- new.env(parent=globalenv())
  for(key in c('setup','sequencepalette','featurepalette','docstyle','showcasedata',
      'learnseq','learnlinks','learngenes','learnlabels','learnlinear'))
    eval(parse(text=chunks[[key]]),environment)
  environment
}
