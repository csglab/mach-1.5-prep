setwd('/scratch/asabe/projects/mach/mach-2')

### Loading libs
library(data.table)
library(stringr)
library(magrittr)
shhh = suppressPackageStartupMessages
library(IRanges) %>% shhh()
library(GenomicRanges) %>% shhh()
library(scales)

set.seed(42)

num_threads <- 64
max_seq_len <- 65536 ## 64K

setDTthreads(num_threads)

### Reading inputs for all organisms
processed_gtf_files = list.files('data/refseq/passed-transcripts', '.csv', full.names = T)
prefixes = processed_gtf_files %>% str_remove('.passed_transcripts.csv')
genome_files = processed_gtf_files %>% str_replace('/passed-transcripts/', '/genome/') %>% str_replace('.passed_transcripts.csv', '.fna.gz')
accessions = processed_gtf_files %>% basename() %>% str_split('_', simplify = T) %>% extract(, 2) %>% str_c('GCF_', .)

refseq_metadata_file = 'data/refseq/metadata/assembly_summary_refseq.txt'
metadata = fread(refseq_metadata_file)
setnames(metadata, '#assembly_accession', 'assembly_accession')
metadata = metadata[assembly_accession %in% accessions]

### Handling the species token based on phylogeny
refseq_phylogeny_file = 'data/refseq/phylogeny/refseq.metadata.passed_vertebrates.commontree.parsed.csv'
phylogeny = fread(refseq_phylogeny_file)

### From NCBI Taxonomy:
phylogeny = phylogeny[accessions, on = .(accessions), .(phylogenic_code, organism, taxid)]
tax_ids = phylogeny$taxid
phylogenic_codes = phylogeny$phylogenic_code
species = phylogeny$organism
phylogenic_code_len = phylogenic_codes %>% nchar() %>% max()

### Processing each organism
for(i in seq_along(processed_gtf_files)) {
  i = 26
  processed_gtf_file = processed_gtf_files[i]
  prefix = prefixes[i]
  genome_file = genome_files[i]
  accession = accessions[i]
  species_ = species[i]
  phylogenic_code = phylogenic_codes[i]
  tax_id = tax_ids[i]

  output_prefix = prefix %>% str_replace('/passed-transcripts/', '/sequences/')
  output_prefix %>% dirname() %>% dir.create(showWarnings = FALSE)
  output_file = str_glue("{output_prefix}.{tax_id}_{phylogenic_code}.transcript_seqs.csv.gz")

  if (file.exists(output_file)) {
    cat("Output file already exists:", output_file, "\n")
    cat("Skipping...\n")
    next
  }
  
  cat(i, processed_gtf_file, prefix, genome_file, accession, species_, phylogenic_code, sep = " | ", "\n")
  
  gtf = fread(processed_gtf_file)
  
  gtf[, exon_number := str_extract(attribute, 'exon_number ""(\\d+)""') %>% str_extract("\\d+") %>% as.integer()]
  gtf[, num_exons := max(exon_number, na.rm = T), transcript_id]
  gtf[, contain_cds := feature %>% str_c(collapse = ',') %>% str_detect('CDS'), transcript_id]
  gtf = gtf[contain_cds == TRUE]
  
  exons = gtf[feature == 'exon', .(seqname, start, end, strand, gene_id, transcript_id, exon_number)]
  cds = gtf[feature == 'CDS', .(seqname, start, end, strand, gene_id, transcript_id, exon_number)]
  
  exons_gr = exons[, GRanges(
    seqnames = str_c(seqname, gene_id, transcript_id, sep = '@'),
    ranges = IRanges(
      start = start,
      end = end),
    strand = strand)]
  
  cds_gr = cds[, GRanges(
    seqnames = str_c(seqname, gene_id, transcript_id, sep = '@'),
    ranges = IRanges(
      start = start,
      end = end),
    strand = strand)]
  
  utr_gr = GenomicRanges::subtract(exons_gr, cds_gr)
  
  .gr2dt = function(gr, feat)
    gr %>% 
    { if(class(gr) %like% 'List') unlist(.) else .} %>% 
    as.data.table() %>% 
    tidyr::separate(seqnames, c('seqname', 'gene_id', 'transcript_id'), '@', convert = T) %>% 
    dplyr::mutate(feature = feat) %>% 
    as.data.table()
  
  utr = .gr2dt(utr_gr, 'UTR')
  exons = .gr2dt(exons_gr, 'exon')
  cds = .gr2dt(cds_gr, 'CDS')
  
  # trs = rbindlist(list(exons, utr, cds))
  cds_utrs = rbind(utr, cds)
  
  cds_utrs = merge(cds_utrs, gtf[, .(seqname, start, end, strand, gene_id, transcript_id, feature, exon_number)],
                   by = c('seqname', 'start', 'end', 'strand', 'gene_id', 'transcript_id', 'feature'),
                   all.x = T)
  
  ### Find 5' or 3' of UTR regions
  cds_utrs[order(transcript_id, start, end), region_index := fifelse(strand == '+', 1:.N, .N:1), transcript_id]
  
  cds_utrs[feature == 'CDS', median_cds_index := median(region_index), transcript_id]
  cds_utrs[, median_cds_index := median_cds_index %>% na.omit() %>% unique(), transcript_id]
  cds_utrs[feature == 'UTR', feature := fifelse(region_index < median_cds_index, '5UTR', '3UTR')]
  
  cds_utrs[, num_exons := max(exon_number, na.rm = TRUE), transcript_id]
  cds_utrs[feature == '5UTR', exon_number := 0]
  cds_utrs[feature == '3UTR', exon_number := num_exons + 1]
  
  cds_utrs = cds_utrs[, .(seqname, strand, gene_id, transcript_id, width, start, end, feature, exon_number)]
  
  ### Appending exons
  exons = merge(exons, gtf[, .(seqname, start, end, strand, gene_id, transcript_id, feature, exon_number)],
                by = c('seqname', 'start', 'end', 'strand', 'gene_id', 'transcript_id', 'feature'),
                all.x = T)
  exons = exons[, .(seqname, strand, gene_id, transcript_id, width, start, end, feature, exon_number)]
  
  trs = rbind(exons, cds_utrs)
  
  ### Removing transcripts with 1 exon only
  # gtf = gtf[num_exons > 1]
  
  gtf[, tr_start := min(start), transcript_id]
  gtf[, tr_end := max(end), transcript_id]
  
  ### Adding flanking regions to the gene coords
  gtf[, tokenized_tr_len := (tr_end - tr_start + 1) + 1 + phylogenic_code_len + 1 + 1 + 1] ## CLS, [SPECIES], S, P, SEP 
  
  gtf <- gtf[tokenized_tr_len <= max_seq_len]
  gtf[, tokenized_tr_len := NULL]
  passed_trs = gtf[, transcript_id %>% unique()]
  
  trs = trs[transcript_id %in% passed_trs]
  
  ### Reversing exon_number on the negative strand
  ### Not needed for RefSeq annotation
  ## gtf[strand == '-', exon_number := (num_exons - exon_number + 1)]
  
  # setnames(gtf, 'feature', 'type')
  
  introns <- trs[feature == 'exon']
  introns <- introns[, start_ := start]
  
  introns[order(end), start := c(end[-.N] + 1, NA), transcript_id]
  introns[order(end), end := c(start_[-1] - 1, NA), transcript_id]
  introns[, start_ := NULL]
  introns[, feature := 'intron']
  introns <- introns[!is.na(start)]
  
  introns[strand == '-', exon_number := exon_number - 1]
  trs <- rbind(trs, introns)
  
  trs <- trs[strand != '*']
  
  invalid_trs <- trs[(start <= 0) | (end <= 0) | (start > end), transcript_id]
  print(str_glue('{length(invalid_trs)} invalid transcripts (due to negative coordinates)'))
  trs <- trs[!(transcript_id %in% invalid_trs)]
  
  tr_info <- gtf[feature == 'transcript', 
                 .(seqname, start, end, gene_id, transcript_id, num_exons)]
  
  colnames(tr_info) <- c('chr', 'tr_start', 'tr_end', 'gene_id', 'transcript_id', 'num_exons')
  
  trs = trs[feature %in% c('5UTR', 'CDS', 'intron', '3UTR')]
  bed = copy(trs)
  bed[, name := str_c(transcript_id, feature, exon_number, sep = '@')]
  bed[, score := '.']
  bed = bed[, .(seqname, start, end, name, score, strand)]
  bed[, start := start - 1]
  
  fwrite(bed, sep='\t', row.names = F, col.names = F, quote = F, file = str_glue("{prefix}.bed"))
  
  # bedtools_getfasta_script <- str_glue("
  #     gunzip --keep --force --stdout {genome_file} > {genome_file}.fasta && \\
  #     bedtools getfasta -s -bedOut -fi {genome_file}.fasta -bed {prefix}.bed > {prefix}.with_seqs.bed && \\
  #     rm --force {genome_file}.fasta*
  #   ")

  bedtools_getfasta_script <- str_glue("
  gunzip --keep --force --stdout {genome_file} > {genome_file}.fasta && \\
  /scratch/asabe/envs/preprocess_seqs/bin/bedtools getfasta -s -bedOut -fi {genome_file}.fasta -bed {prefix}.bed > {prefix}.with_seqs.bed  && \\
  rm --force {genome_file}.fasta*")
  
  # Execute the script
  system(bedtools_getfasta_script)
  
  seqs = fread(str_glue("{prefix}.with_seqs.bed"), col.names = c('chr', 'start', 'end', 'name', 'score', 'strand', 'seq'))
  seqs %<>% tidyr::separate(name, c('transcript_id', 'feature', 'exon_number'), convert = T, sep = '@') %>% as.data.table()
  seqs = seqs[, .(chr, strand, start, end, feature, transcript_id, exon_number, seq)]
  setnames(seqs, 'chr', 'seqname')
  seqs[, start := start + 1] ## undo bed[, start := start - 1]
  
  trs = merge(trs, seqs, by = c('seqname', 'strand', 'start', 'end', 'feature', 'transcript_id', 'exon_number'))
  
  invalid_trs <- trs[str_count(seq, 'N') > 0, transcript_id %>% unique()]
  print(str_glue('{length(invalid_trs)} invalid_rna_seqs (due to N)'))
  # gtf <- gtf[!(transcript_id %in% invalid_trs)]
  
  trs[feature == 'CDS', seq := str_to_upper(seq)]
  trs[feature == 'intron', seq := str_to_lower(seq)]
  
  utr_transformer <- c('A' = 'W', 'C' = 'X', 'G' = 'Y', 'T' = 'Z', 'N' = 'U')
  transform_to_utr <- function(seq) {
    utr_transformer[
      seq %>% str_split('') %>% unlist()
    ] %>% str_c(collapse = '')
  }
  str_to_utr <- Vectorize(transform_to_utr, 'seq')
  
  trs[feature %like% 'UTR', seq := seq %>% str_to_upper() %>% str_to_utr()]
  
  trs[order(transcript_id, start, end), region_index := fifelse(strand == '+', 1:.N, .N:1), transcript_id]
  
  trs[feature == 'CDS', is_first_cds := region_index == min(region_index), transcript_id]
  trs[feature == 'CDS', is_last_cds := region_index == max(region_index), transcript_id]
  
  trs[is_first_cds == TRUE, seq := str_c('S', seq)]
  trs[is_last_cds == TRUE, seq := str_c(seq, 'E')]
  trs[, c('is_first_cds', 'is_last_cds') := NULL]
  
  trs[, feature := factor(feature, levels = c('5UTR', 'CDS', 'intron', '3UTR'), ordered = TRUE)]
  
  tr_seqs <- trs[
    order(transcript_id, region_index),
    .(seq = str_c(seq, collapse = '')),
    by = .(transcript_id, strand)]
  
  tr_seqs[, seq := str_c(phylogenic_code, 'I', seq, 'P')]
  
  # tr_seqs = tr_seqs[transcript_id == 'XM_006201594.3']
  
  ### Handling stop_codon. in RefSeq GTF's stop_codon is not part of the last CDS
  cds_transformer <- c('W' = 'A', 'X' = 'C', 'Y' = 'G', 'Z' = 'T', 'U' = 'N',
                       'a' = 'a', 'c' = 'c', 'g' = 'g', 't' = 't', 'n' = 'n')
  
  transform_to_cds <- function(seq) {
    cds_transformer[
      seq %>% str_split('') %>% unlist()
    ] %>% str_c(collapse = '')
  }
  str_to_cds <- Vectorize(transform_to_cds, 'seq')
  
  tr_seqs[, c('before_E', 'after_E') := tstrsplit(seq, 'E')]
  tr_seqs[, stop_codon_end := after_E %>% gregexpr('W|X|Y|Z', .) %>% sapply(extract, 3)]
  tr_seqs[, before_E := after_E %>% str_sub(end = stop_codon_end) %>% str_to_cds() %>% str_c(before_E, .)]
  tr_seqs[, after_E := after_E %>% str_sub(start = stop_codon_end + 1)]
  tr_seqs[, seq := str_c(before_E, 'E', after_E)]
  tr_seqs[, seq_len := nchar(seq)]
  tr_seqs[, contains_unk := str_detect(seq, 'N|n|U')]
  
  tr_seqs = tr_seqs[, .(transcript_id, seq, seq_len, contains_unk)]
  tr_seqs <- merge(tr_seqs, tr_info, by = 'transcript_id', all.x = TRUE)
  
  tr_metadata = tr_seqs[, .(transcript_id, gene_id, chr, tr_start, tr_end, num_exons, seq_len)]
  tr_metadata[, species := species_]
  tr_metadata[, tax_id := tax_id]
  tr_metadata[, phylogenic_code := phylogenic_code]
  tr_metadata[, accession := accession]
  tr_metadata[, prefix := prefix]
  tr_metadata[, processed_gtf_file := processed_gtf_file]
  tr_metadata[, genome_file := genome_file]
  
  tr_seqs = tr_seqs[, .(transcript_id, gene_id, chr, seq, seq_len, contains_unk)]
  setorderv(tr_seqs, c('transcript_id', 'gene_id', 'chr'))
  setnames(tr_seqs, 'seq', 'text')
  
  ### Decided not to tokenized UTRs as different tokens
  
  tr_seqs[, c('before_S', 'after_S') := tstrsplit(text, 'S')]
  tr_seqs[, c('before_S_before_I', 'before_S_after_I') := tstrsplit(before_S, 'I')]
  tr_seqs[, c('after_S_before_E', 'after_S_after_E') := tstrsplit(after_S, 'E')]
  tr_seqs[, after_S_after_E := str_remove(after_S_after_E, 'P')]
  tr_seqs[, text2 := str_c(
    before_S_before_I,
    'I',
    str_to_cds(before_S_after_I),
    after_S_before_E,
    str_to_cds(after_S_after_E),
    'P')]
    
  fwrite(tr_seqs, str_glue("{output_prefix}.{tax_id}_{phylogenic_code}.transcript_seqs.csv.gz"))
  fwrite(tr_metadata, str_glue("{output_prefix}.{tax_id}_{phylogenic_code}.transcript_metadata.csv.gz"))
  
  num_genes = tr_metadata[, gene_id %>% unique() %>% length()] %>% comma()
  num_trs = tr_metadata[, transcript_id %>% unique() %>% length()] %>% comma()
  num_tokens = tr_metadata[, seq_len %>% sum()] %>% comma()

  cat(i, species_, num_genes, num_trs, num_tokens, sep = " | ", "\n")
  print('======================')
  
}