setwd('/scratch/asabe/projects/mach/mach-2')

library(data.table)
library(stringr)
library(magrittr)
library(purrr, warn.conflicts = F)

set.seed(42)
num_threads <- 64

setDTthreads(num_threads)

num_samples = 128
full_dataset_file = 'data/refseq/sequences/GCF_000001405.40_GRCh38.p14.seqs.processed.csv.gz'
output_prefix = full_dataset_file %>% str_replace('.csv.gz$', '.sampled') %>% str_replace('/sequences/', '/overfit/')

dat = fread(full_dataset_file)
long_sampled_dat = dat[order(-seq_len)][seq(num_samples)]
short_sampled_dat = dat[order(seq_len)][seq(num_samples * 8)]

fwrite(long_sampled_dat, str_glue("{output_prefix}.long.csv.gz"))
fwrite(short_sampled_dat, str_glue("{output_prefix}.short.csv.gz"))

