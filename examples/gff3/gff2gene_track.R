library(tidyverse)

gff3FilesPath <- list.files(path = ".", pattern = "*.gff3")

gff3Table <- map_df(gff3FilesPath,~read_tsv(.x,show_col_types = F,comment = "#",col_names = F) %>% set_names(c("seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes")))

geneTrackTable <- gff3Table %>% filter(type=="CDS") %>% mutate(anno=str_extract(attributes,"(?<=product=)[^;]+(?=;)")) %>% select(seq_id,start,end,strand,anno)

write_tsv(geneTrackTable,"gene_track.tsv")
