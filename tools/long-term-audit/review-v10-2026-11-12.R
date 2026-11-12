# Focused audit probes, with explicitly synthetic inputs; not a package test suite.
# Records observations rather than accepting the current defects as a contract.
pkgload::load_all(quiet = TRUE)
seqs <- data.frame(accver = c('A', 'B'), length = 1000)
args <- commandArgs(trailingOnly = TRUE)
output <- args[!startsWith(args, '--')]
if (!length(output)) output <- 'tools/long-term-audit/review-v10-2026-11-12-results.json'
results <- list(baseline = trimws(system2('git', c('rev-parse', '--short', 'HEAD'), stdout = TRUE)), R = R.version.string,
  ggplot2 = as.character(utils::packageVersion('ggplot2')))

observations <- data.frame(accver = 'A', pos = c(100, 100), category = c('one', 'two'))
p <- ggchord(seqs) + geom_seq() + geom_track(
  aes(x = pos, fill = category, group = 1), data = observations,
  geom = 'bar', track = position_track('counts')) + coord_collinear()
track <- get_chord_layout(p)$track_statistics[[1]]
reference <- ggplot2::ggplot(observations, aes(x = pos, fill = category, group = 1)) +
  ggplot2::geom_bar()
native <- suppressWarnings(ggplot2::ggplot_build(reference)$data[[1]])
results$explicit_group <- list(track_count = track$count, native_count = native$count,
  defect_observed = !identical(track$count, native$count))

points <- data.frame(accver = 'A', pos = c(100, 150, 200))
p <- ggchord(seqs) + geom_seq() + geom_track(aes(x = pos), data = points,
  geom = 'col', stat = 'bin', breaks = c(100, 200),
  track = position_track('bins')) + coord_collinear()
binned <- get_chord_layout(p)$track_statistics[[1]]
results$bin_source_membership <- list(count = binned$count,
  source_rows = unname(binned$source_rows), expected_source_rows = list(1:3),
  defect_observed = !identical(unname(binned$source_rows[[1]]), 1:3))

hits <- data.frame(qaccver = 'A', saccver = 'B', qstart = 100, qend = 200,
  sstart = 300, send = 400, pident = 90, source_hsp_id = 'observed-1')
p <- ggchord(seqs, hits) + geom_seq() + geom_link_ribbon() + geom_link_ribbon()
acceptance <- lapply(list(chord = coord_chord(), collinear = coord_collinear(),
  radial = coord_radial_synteny()), function(coord) {
  tryCatch({ get_chord_layout(p + coord); list(accepted = TRUE) },
    error = function(e) list(accepted = FALSE, message = conditionMessage(e)))
})
results$declared_duplicate_hsp <- list(views = acceptance,
  defect_observed = any(vapply(acceptance, function(x)
    isTRUE(x$accepted) || !grepl('allow_repeated', x$message), logical(1))))

source <- ggchord(seqs) + geom_seq(seq_order = c(A1 = 'A', B = 'B', A2 = 'A')) + coord_chord()
detail <- ggchord(seqs) + geom_seq() + coord_collinear()
region <- data.frame(accver = 'A', start = 100, end = 200, instance_note = 'A2')
selected <- tryCatch(export_ggchord_fusion(magnify_ggchord(source, region, detail))$focus_regions,
  error = conditionMessage)
plain <- region[c('accver', 'start', 'end')]
without_note <- tryCatch({export_ggchord_fusion(magnify_ggchord(source, plain, detail)); 'accepted'},
  error = conditionMessage)
results$metadata_selector <- list(with_metadata = selected,
  without_metadata = without_note,
  defect_observed = !is.character(selected) || !grepl('ambiguous', selected) ||
    !grepl('ambiguous', without_note))

if ('--check' %in% args) {
  seqs <- data.frame(accver=LETTERS[1:3],length=1000)
  h <- data.frame(qaccver='A',saccver=c('B','C'),qstart=100,qend=200,sstart=300,send=400,pident=90,source_hsp_id=c('ab','ac'))
  base <- ggchord(seqs,h)+geom_seq()
  # One layer keeps AB, another keeps AC: omitted duplicates must not cause errors.
  p <- base+geom_link_ribbon(pairs=data.frame(qaccver='A',saccver='B'))+
   geom_link_ribbon(pairs=data.frame(qaccver='A',saccver='C'))+coord_collinear()
  out <- get_chord_layout(p)
  stopifnot(sum(out$connection_report$retained)==2L)
  # Explicit names activate instance routing, while default adjacency filters AC.
  p <- ggchord(seqs,h)+geom_seq(seq_order=c(one='A',two='B',three='C'))+
   geom_link_ribbon()+geom_link_ribbon(data=h[2,],adjacent_only=FALSE)+coord_collinear()
  out <- get_chord_layout(p)
  stopifnot(sum(out$instance_routing$retained)==2L,!any(out$instance_routing$repeated_display))
  p <- ggchord(seqs,h[1,])+geom_seq()+geom_link_ribbon(allow_repeated=TRUE)+
   geom_link_ribbon(allow_repeated=TRUE)+coord_collinear()
  stopifnot(nrow(get_chord_layout(p)$ribbon_polys)>0)
  # Bundled derived rows still refer to declared identities of contributing HSPs.
  h2 <- rbind(h[1,],transform(h[1,],qstart=110,qend=210,source_hsp_id='ab2'))
  bundle <- ggchord(seqs,h2,validate='none')+geom_seq()+stat_ribbon_bundle(bins=1,min_bundle=2,weight='count')+geom_link_ribbon()+coord_collinear()
  err <- tryCatch(get_chord_layout(bundle),error=conditionMessage)
  stopifnot(is.character(err),grepl('allow_repeated',err))
  # Same Stat memberships survive padding/drop and a data Position; callbacks run once/group.
  d <- data.frame(accver='A',pos=c(100,150+1e-8,200),weight=c(2,3,4))
  for(closed in c('right','left')) {
    calls <- 0L
    breaks <- function(x) { calls <<- calls+1L; c(100,150,200) }
    p <- ggchord(seqs)+geom_seq()+geom_track(aes(x=pos,weight=weight),data=d,
      geom='col',stat='bin',breaks=breaks,closed=closed,pad=TRUE,drop='all',
      position='nudge',track=position_track('bins'))+coord_collinear()
    b <- get_chord_layout(p)$track_statistics[[1]]
    stopifnot(calls==1L, setequal(unlist(b$source_rows),1:3),
      identical(b$count, vapply(b$source_rows,function(i) sum(d$weight[i]),numeric(1))))
  }
  invisible(ggplot2::ggplotGrob(p))
  results$related_checks <- list(display_selection = TRUE, explicit_consent = TRUE,
    bundle_source_members = TRUE, native_bin_boundaries = TRUE, bin_grob = TRUE)
}

results$runtime_md5 <- as.list(tools::md5sum(c('R/track-grammar.R', 'R/display-instances.R',
  'R/ggchord-build.R', 'R/coord-radial-synteny.R', 'R/magnify.R')))
results$code_state <- if (length(system2('git', c('diff', '--name-only', '--', 'R'), stdout = TRUE)))
  'working tree changes on baseline' else 'committed runtime'
jsonlite::write_json(results, output[1],
  pretty = TRUE, auto_unbox = TRUE)
defects <- sum(vapply(results[4:7], function(x) isTRUE(x$defect_observed), logical(1)))
cat('Audit probes completed; observed defects:', defects, '\n')
if ('--check' %in% args) stopifnot(defects == 0L)
