setwd('/scratch/asabe/projects/mach/mach-2')

### Loading libs
library(data.table)
library(stringr)
library(magrittr)

splitted_data_dir = 'data/refseq/splitted'
csvs = list.files(splitted_data_dir, pattern = '*.csv.gz', full.names = TRUE)
csvs_info = csvs %>% basename() %>% str_remove('.csv.gz') %>% str_split_fixed('.transcript_seqs_', 2)
organisms = csvs_info %>% extract(, 1)
splits = csvs_info %>% extract(, 2)

get_stats <- function(csv, organism, split) {
  print(str_glue('Processing {organism} {split}'))
  dat = fread(csv)
  num_genes = dat[, gene_id %>% unique() %>% length()]
  num_trs = dat[, transcript_id %>% unique() %>% length()]
  num_tokens = dat[, sum(seq_len)] %>% sum()
  num_exonic_tokens = dat[, str_count(text, '[ACGU]')] %>% sum()
  num_intronic_tokens = dat[, str_count(text, '[acgu]')] %>% sum()
  num_unk_tokens = dat[, str_count(text, '[Nn]')] %>% sum()
  num_contains_unk = dat[, sum(contains_unk)]
  
  data.table(organism = organism,
             split = split,
             num_genes = num_genes,
             num_trs = num_trs,
             num_tokens = num_tokens,
             num_exonic_tokens = num_exonic_tokens,
             num_intronic_tokens = num_intronic_tokens,
             num_unk_tokens = num_unk_tokens,
             num_contains_unk = num_contains_unk)
}

stats_list <- mapply(get_stats, 
                    csvs, 
                    organisms, 
                    splits,
                    SIMPLIFY = FALSE)

stats_dt <- rbindlist(stats_list)
fwrite(stats_dt, 'data/refseq/splitted/stats.csv')
