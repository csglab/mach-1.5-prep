setwd('/scratch/asabe/projects/mach/mach-2')

library(data.table)
library(stringr)
library(magrittr)
library(purrr, warn.conflicts = F)

set.seed(42)
num_threads <- 64

setDTthreads(num_threads)

### Reading inputs for all organisms
processed_sequence_files = list.files('data/refseq/sequences', '.seqs.csv.gz$', full.names = T)

seq_dt = data.table(
  file = processed_sequence_files,
  meta = processed_sequence_files %>% str_replace('.seqs.csv.gz$', '.metadata.csv.gz'),
  prefix = processed_sequence_files %>% basename() %>% str_remove('.seqs.csv.gz$'),
  accession = processed_sequence_files %>% basename() %>% str_split('_', simplify = T) %>% extract(, 2) %>% str_c('GCF_', .))

### Filtering the accessions based on the new phylogenic codes
refseq_phylogeny_file = 'data/refseq/phylogeny/refseq.metadata.passed_vertebrates.euarchontoglires.commontree.parsed.csv'
phylogeny = fread(refseq_phylogeny_file)

seq_dt = merge(seq_dt,
               phylogeny[, .(organism_name, taxid, accessions, phylogenic_code)],
               by.x = 'accession',
               by.y = 'accessions')

seq_dt = seq_dt[, .(file, meta, phylogenic_code)]
colnames(seq_dt) = c('file_path', 'meta_path', 'updated_code')

.update_phylogenic_code = function(file_path, meta_path, updated_code) {
  
  cat(file_path, ' — ', updated_code, '\n')
  seqs = fread(file_path)
  seqs[, text := text %>% str_remove('^.*I') %>% str_c(updated_code, 'I', .)]
  seqs[, seq_len := nchar(text)]
  
  metadata = fread(meta_path)
  metadata[, phylogenic_code := as.character(phylogenic_code)]
  metadata[, phylogenic_code := updated_code]
  
  .rename = function(x) str_replace(x, '.csv.gz$', '.processed.csv.gz')
  fwrite(seqs, .rename(file_path))
  fwrite(metadata, .rename(meta_path))
  
}

pmap(seq_dt, .update_phylogenic_code)

