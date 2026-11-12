# Pure process/format controls: no biological datasets or analysis engines.
pkgload::load_all(quiet = TRUE)
python <- Sys.getenv('PARSER_PYTHON', 'python3')
probe <- r"---(
import sys,json
with open(sys.argv[1],encoding='utf-8') as h: config=json.load(h)
print('format-control process message',file=sys.stderr)
with open(sys.argv[2],'w',encoding='utf-8') as h:
    json.dump({'prefix':sys.prefix,'base':sys.base_prefix,'configuration':config},h)
)---"
paths <- c('literal $(untouched).txt', 'format control with spaces.txt')
checked <- ggchord_parse_external_source(probe, list(files = as.list(paths)), python, 'process check')
stopifnot(identical(unlist(checked$data$configuration$files), paths),
  'format-control process message' %in% checked$messages)
# Real decoder APIs on explicitly synthetic format-only records. No claimed
# molecule, observation, homology, source characters, support or new analysis.
temporary <- tempfile('ggchord-parser-control-'); dir.create(temporary)
on.exit_cleanup <- function() unlink(temporary, recursive = TRUE)
first <- file.path(temporary, 'format control one.gb')
second <- file.path(temporary, 'format control two.gb')
binary <- file.path(temporary, 'format control.bcf')
make_controls <- r"---(
import json,sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pysam
with open(sys.argv[1],encoding='utf-8') as h: config=json.load(h)
for number,path in enumerate(config['genbank'],1):
    record=SeqRecord(Seq(None,length=8*number),id='FORMAT_ONLY_'+str(number),description='Synthetic format control, no biological claim')
    record.annotations['molecule_type']='DNA'
    SeqIO.write(record,path,'genbank')
header=pysam.VariantHeader()
header.contigs.add('FORMAT_ONLY',length=100)
with pysam.VariantFile(config['bcf'],'wb',header=header) as output:
    output.write(output.new_record(contig='FORMAT_ONLY',start=4,stop=5,alleles=('A','C')))
with open(sys.argv[2],'w') as h: json.dump({'made':True},h)
)---"
ggchord_parse_external_source(make_controls, list(genbank = as.list(c(first, second)), bcf = binary),
  python, 'format fixture preparation')
tryCatch({
  genbank <- read_genbank(c(first, second), python = python, source_file = TRUE)
  stopifnot(identical(genbank$seq_data$accver, c('FORMAT_ONLY_1', 'FORMAT_ONLY_2')),
    identical(as.numeric(genbank$seq_data$length), c(8, 16)))
  single <- read_genbank(first, python = python)
  stopifnot(single$seq_data$length == 8)
  variant <- read_vcf(binary, python = python, source_file = TRUE)
  stopifnot(nrow(variant$features$variant) == 1,
    variant$features$variant$accver == 'FORMAT_ONLY', variant$features$variant$start == 5,
    variant$report$records[[1L]]$binary_bridge$package == 'pysam')
}, finally = on.exit_cleanup())
jsonlite::write_json(list(status = 'passed', json_paths_unchanged = TRUE,
  stderr_message_retained = TRUE, selected_environment_prefix = checked$data$prefix,
  virtual_environment_preserved = checked$data$prefix != checked$data$base,
  genbank_single_file_control = TRUE, genbank_two_file_control = TRUE, bcf_control = TRUE,
  fixtures = 'explicitly synthetic format controls, no biological claims'),
  'tools/long-term-audit/v0.27.0-source-parser-results.json', pretty = TRUE, auto_unbox = TRUE)
cat('Optional parser process boundary: passed\n')
