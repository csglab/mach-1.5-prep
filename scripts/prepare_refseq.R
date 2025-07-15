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

### Resolving comma and slash issue in the phylogenic code
# refseq_phylogeny_with_comma_and_slash_file = 'data/refseq/phylogeny/refseq.metadata.passed_vertebrates.commontree.parsed.with_comma_and_slash.csv'
# phylogeny_with_comma_and_slash = fread(refseq_phylogeny_with_comma_and_slash_file)
# setnames(phylogeny_with_comma_and_slash, 'phylogenic_code', 'phylogenic_code_with_comma_and_slash')
# phylogeny = merge(phylogeny,
#                   phylogeny_with_comma_and_slash[, .(taxid, phylogenic_code_with_comma_and_slash)],
#                   by = 'taxid')
# phylogeny = phylogeny[phylogenic_code != phylogenic_code_with_comma_and_slash]
# accessions = accessions[accessions %in% phylogeny$accessions]
# processed_gtf_files = processed_gtf_files[processed_gtf_files %like% str_c(accessions, collapse = '|')]
# 
# print(processed_gtf_files)

### From NCBI Taxonomy:
phylogeny = phylogeny[accessions, on = .(accessions), .(phylogenic_code, organism, taxid)]
tax_ids = phylogeny$taxid
phylogenic_codes = phylogeny$phylogenic_code
species = phylogeny$organism
phylogenic_code_len = phylogenic_codes %>% nchar() %>% max()

### Processing each organism
for(i in seq_along(processed_gtf_files)) {
  
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
  
  print('===========================')
  cat(i, processed_gtf_file, prefix, genome_file, accession, species_, phylogenic_code, sep = " | ", "\n")
  
  gtf = fread(processed_gtf_file)
  
  gtf[, exon_number := str_extract(attribute, 'exon_number ""(\\d+)""') %>% str_extract("\\d+") %>% as.integer()]
  gtf[, num_exons := max(exon_number, na.rm = T), transcript_id]
  gtf[, contain_cds := feature %>% str_c(collapse = ',') %>% str_detect('CDS'), transcript_id]
  gtf = gtf[contain_cds == TRUE]
  
  exons = gtf[feature == 'exon', .(seqname, start, end, strand, gene_id, transcript_id, exon_number)]
  
  exons_gr = exons[, GRanges(
    seqnames = str_c(seqname, gene_id, transcript_id, sep = '@'),
    ranges = IRanges(
      start = start,
      end = end),
    strand = strand)]
  
  .gr2dt = function(gr, feat)
    gr %>%
    as.data.table() %>% 
    tidyr::separate(seqnames, c('seqname', 'gene_id', 'transcript_id'), '@', convert = T) %>% 
    dplyr::mutate(feature = feat) %>% 
    as.data.table()
  
  exons = .gr2dt(exons_gr, 'exon')
  
  ### Appending exons
  exons = merge(exons, gtf[, .(seqname, start, end, strand, gene_id, transcript_id, feature, exon_number)],
                by = c('seqname', 'start', 'end', 'strand', 'gene_id', 'transcript_id', 'feature'),
                all.x = T)
  exons = exons[, .(seqname, strand, gene_id, transcript_id, width, start, end, feature, exon_number)]
  
  ### Removing transcripts with 1 exon only
  # gtf = gtf[num_exons > 1]
  
  gtf[, tr_start := min(start), transcript_id]
  gtf[, tr_end := max(end), transcript_id]
  
  ### Adding flanking regions to the gene coords
  gtf[, tokenized_tr_len := (tr_end - tr_start + 1) + 1 + phylogenic_code_len + 1 + 1 + 1] ## CLS, [SPECIES], S, P, SEP 
  
  gtf <- gtf[tokenized_tr_len <= max_seq_len]
  gtf[, tokenized_tr_len := NULL]
  passed_trs = gtf[, transcript_id %>% unique()]
  
  exons = exons[transcript_id %in% passed_trs]
  
  ### Reversing exon_number on the negative strand
  ### Not needed for RefSeq annotation
  ## gtf[strand == '-', exon_number := (num_exons - exon_number + 1)]
  # setnames(gtf, 'feature', 'type')
  
  introns <- copy(exons)
  introns <- introns[, start_ := start]
  
  introns[order(end), start := c(end[-.N] + 1, NA), transcript_id]
  introns[order(end), end := c(start_[-1] - 1, NA), transcript_id]
  introns[, start_ := NULL]
  introns[, feature := 'intron']
  introns <- introns[!is.na(start)]
  
  introns[strand == '-', exon_number := exon_number - 1]
  trs <- rbind(exons, introns)
  
  trs <- trs[strand != '*']
  
  invalid_trs <- trs[(start <= 0) | (end <= 0) | (start > end), transcript_id]
  print(str_glue('{length(invalid_trs)} invalid transcripts (due to negative coordinates)'))
  trs <- trs[!(transcript_id %in% invalid_trs)]
  
  tr_info <- gtf[feature == 'transcript', 
                 .(seqname, start, end, gene_id, transcript_id, num_exons)]
  
  colnames(tr_info) <- c('chr', 'tr_start', 'tr_end', 'gene_id', 'transcript_id', 'num_exons')
  
  trs = trs[feature %in% c('exon', 'intron')]
  bed = copy(trs)
  bed[, name := str_c(transcript_id, feature, exon_number, sep = '@')]
  bed[, score := '.']
  bed = bed[, .(seqname, start, end, name, score, strand)]
  bed[, start := start - 1]
  
  fwrite(bed, sep='\t', row.names = F, col.names = F, quote = F, file = str_glue("{prefix}.bed"))
  
  bedtools_getfasta_script <- str_glue("
      gunzip --keep --force --stdout {genome_file} > {genome_file}.fasta && \\
      bedtools getfasta -s -bedOut -fi {genome_file}.fasta -bed {prefix}.bed > {prefix}.with_seqs.bed && \\
      rm --force {genome_file}.fasta*
    ")
  
  # bedtools_getfasta_script <- str_glue("
  # gunzip --keep --force --stdout {genome_file} > {genome_file}.fasta && \\
  # /scratch/asabe/envs/preprocess_seqs/bin/bedtools getfasta -s -bedOut -fi {genome_file}.fasta -bed {prefix}.bed > {prefix}.with_seqs.bed  && \\
  # rm --force {genome_file}.fasta*")
  
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
  
  trs[feature == 'exon', seq := str_to_upper(seq)]
  trs[feature == 'intron', seq := str_to_lower(seq)]
  
  trs[order(transcript_id, start, end), region_index := fifelse(strand == '+', 1:.N, .N:1), transcript_id]
  
  trs[, feature := factor(feature, levels = c('exon', 'intron'), ordered = TRUE)]
  
  tr_seqs <- trs[
    order(transcript_id, region_index),
    .(seq = str_c(seq, collapse = '')),
    by = .(transcript_id, strand)]
  
  tr_seqs[, seq := str_replace_all(seq, 'T', 'U')]
  tr_seqs[, seq := str_replace_all(seq, 't', 'u')]
  tr_seqs[, seq := str_c(phylogenic_code, 'I', seq, 'P')]
  
  tr_seqs[, contains_unk := str_detect(seq, '[RYSWKMBDHVNryswhmbdhvn.-]')] ## from infer_phylogenic_code.R (IUPAC_codes)
  tr_seqs[, seq_len := nchar(seq)]
  
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
  
  if(str_detect(phylogenic_code, '/'))
    phylogenic_code = str_replace_all(phylogenic_code, '/', '-slash-')
  
  fwrite(tr_seqs, str_glue("{output_prefix}.{tax_id}.transcript_seqs.csv.gz"))
  fwrite(tr_metadata, str_glue("{output_prefix}.{tax_id}.transcript_metadata.csv.gz"))
  
  num_genes = tr_metadata[, gene_id %>% unique() %>% length()] %>% comma()
  num_trs = tr_metadata[, transcript_id %>% unique() %>% length()] %>% comma()
  num_tokens = tr_metadata[, seq_len %>% sum()] %>% comma()
  
  cat(i, species_, num_genes, num_trs, num_tokens, sep = " | ", "\n")
  
}