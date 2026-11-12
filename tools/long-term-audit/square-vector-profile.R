source('tools/long-term-audit/route-vector-profile-candidate.R')
helper <- readLines('tools/long-term-audit/check-label-quality.R', warn = FALSE)
end <- grep('^# Explicit geometry controls', helper)[1L]-1L
eval(parse(text = helper[seq_len(end)]))
lines <- readLines('tools/documentation-examples.R', warn = FALSE)
marks <- grep('^# ---- [a-z]+ ----$', lines)
chunks <- setNames(lapply(seq_len(length(marks)-1L), function(i)
  lines[seq.int(marks[i]+1L, marks[i+1L]-1L)]),
  sub('^# ---- ([a-z]+) ----$', '\\1', lines[marks[-length(marks)]]))
for (name in c('docstyle', 'showcaseplasmiddata', 'showcaseplasmid')) eval(parse(text = chunks[[name]]))
namespace <- asNamespace('ggchord')
unlockBinding('ggchord_route_shortest', namespace)
assign('ggchord_route_shortest', candidate, envir = namespace)
lockBinding('ggchord_route_shortest', namespace)
prefix <- 'v0.27.0-flagship-vector-prototype-'
Rprof('/tmp/ggchord-square-vector-profile.out', interval = .01)
result <- tryCatch({
  setTimeLimit(elapsed = 120, transient = FALSE)
  measure(p_show_plasmid, 'square', 8, 8)
}, error = function(e) list(interrupted = conditionMessage(e)), finally = {
  setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  Rprof(NULL)
})
profile <- summaryRprof('/tmp/ggchord-square-vector-profile.out')
print(head(profile$by.total, 18))
print(head(profile$by.self, 18))
jsonlite::write_json(list(scope = 'Bounded diagnostic sample of the in-process candidate; not a quality pass or a performance threshold',
  production_override = TRUE, candidate_sha256 = digest::digest(deparse(body(candidate)), algo = 'sha256'),
  outcome = result, total_samples_seconds = profile$sampling.time,
  by_total = head(profile$by.total, 25), by_self = head(profile$by.self, 25)),
  'tools/long-term-audit/v0.27.0-square-profile-results.json', pretty = TRUE, auto_unbox = TRUE)
cat('Square bounded profile recorded\n')
