setwd('/scratch/asabe/projects/mach/mach-2')

library(data.table)
library(stringr)
library(magrittr)

set.seed(42)
num_threads = 64
setDTthreads(num_threads)

### Reading inputs for all organisms
processed_gtf_files = list.files('data/refseq/passed-transcripts', '.csv', full.names = T)
prefixes = processed_gtf_files %>% str_remove('.passed_transcripts.csv')
genome_files = processed_gtf_files %>% str_replace('/passed-transcripts/', '/genome/') %>% str_replace('.passed_transcripts.csv', '.fna.gz')
accessions = processed_gtf_files %>% basename() %>% str_split('_', simplify = T) %>% extract(, 2) %>% str_c('GCF_', .)

refseq_metadata_file = 'data/refseq/metadata/assembly_summary_refseq.txt'
metadata = fread(refseq_metadata_file) %>% suppressWarnings()
setnames(metadata, '#assembly_accession', 'assembly_accession')

# vertebrates = metadata[group %like% '^vertebrate', .(taxid, organism_name, assembly_accession)]
# vertebrates[order(taxid), .(organism_name)] %>% unique() %>% fwrite('data/refseq/metadata/refseq.metadata.vertebrates.organism_name.txt', sep = '\n', quote = F, col.names = F, row.names = F)
# vertebrates[order(taxid), .(taxid)] %>% unique() %>% fwrite('data/refseq/metadata/refseq.metadata.vertebrates.tax_id.txt', sep = '\n', quote = F, col.names = F, row.names = F)

metadata = metadata[assembly_accession %in% accessions]
# metadata[order(taxid), .(taxid)] %>% unique() %>% fwrite('data/refseq/metadata/refseq.metadata.passed_vertebrates.tax_id.txt', sep = '\n', quote = F, col.names = F, row.names = F)

commontree_file = 'data/refseq/phylogeny/refseq.metadata.passed_vertebrates.euarchontoglires.commontree.txt'

lines = readLines(commontree_file)

tree = data.table(
  line = lines,
  rank_char = lines %>% 
    str_remove('[A-Z].*') %>% 
    str_replace_all("\\+", "\\\\"),
  content = lines %>% 
    str_extract_all('[A-Za-z].*') %>% 
    unlist()
)

tree[, rank := nchar(rank_char)]
tree[, rank := rank / 2]
tree[, rank := rank + 1]

tree = tree[, .(line, rank, content)]
tree[, trend := rank %>% subtract(c(rank[2:.N], 0)) %>% sign()]
tree[, is_organism := fifelse(trend >= 0, T, F)]
tree[, index := .I]

taxa_levels = c('class', 'order', 'family', 'genus', 'species', 'subspecies')
rank_taxa = data.table(
  rank = taxa_levels %>% length() %>% seq(),
  taxa = taxa_levels)

org_index = tree[is_organism == TRUE, index]

organisms_ = lapply(org_index, function(org_i) {
  org_tree = tree[.(1:org_i), on = .(index)]
  org_rank = org_tree[.(org_i), on = .(index), rank]
  org_phylo = org_tree[rank %>% rev() %>% duplicated() %>% not() %>% rev()][
    rank <= org_rank][, .(content, rank)]
  org_phylo = merge(org_phylo, rank_taxa, by = 'rank', all.x = T, sort = F)
  org_phylo = org_phylo[, .(taxa, content)]
  org_phylo = org_phylo[.(taxa_levels), on = .(taxa)]
  org_phylo = org_phylo %>% dcast(. ~ taxa, value.var = 'content')
  org_phylo = org_phylo[, ..taxa_levels]
  if(is.na(org_phylo$species))
    org_phylo[, c('family', 'genus', 'species') := org_phylo[, .(order, family, genus)]]
   
  org_phylo
})

organisms = rbindlist(organisms_)
organisms[, organism_name := fifelse(is.na(subspecies), species, subspecies)]
organisms[, subspecies := NULL]
taxa_levels = taxa_levels[taxa_levels != 'subspecies']

tokenizer_chars = c(0:9, LETTERS, letters)

IUPAC_codes = c(
    'A',    # Adenine
    'C',    # Cytosine
    'G',    # Guanine
    'T',    # Thymine (DNA)
    'U',    # Uracil (RNA)
    'R',    # puRine (A or G)
    'Y',    # pYrimidine (C or U)
    'S',    # Strong bond (G or C) - 3 hydrogen bonds
    'W',    # Weak bond (A or U) - 2 hydrogen bonds
    'K',    # Keto (G or U)
    'M',    # aMino (A or C)
    'B',    # not A (C or G or U)
    'D',    # not C (A or G or U)
    'H',    # not G (A or C or U)
    'V',    # not U (A or C or G)
    'N',    # aNy base (A or C or G or U)
    '-',    # gap in alignment
    '.'     # alternative gap symbol
)

reserved_chars = c('I' ,'P', 'S', 'E', 'W', 'X', 'Y', 'Z')
reserved_chars = c(reserved_chars, str_to_lower(reserved_chars))
IUPAC_codes = c(IUPAC_codes, str_to_lower(IUPAC_codes))
reserved_chars = c(reserved_chars, IUPAC_codes)
tokenizer_chars = tokenizer_chars[!(tokenizer_chars %in% reserved_chars)]
num_chars = length(tokenizer_chars)

group2code = data.table(
  GRP = seq(0, num_chars - 1),
  code = tokenizer_chars)

setorder(group2code, GRP)

for(taxa_i in seq_along(taxa_levels)) {
  
  parent_taxas = taxa_levels[1:taxa_i-1]
  taxa = taxa_levels[taxa_i]
  
  if(length(parent_taxas) == 0) {
    organisms[, code := '1']
    setnames(organisms, 'code', str_c(taxa, 'code', sep = '_'))
    next
  }
  
  taxa_dt = organisms[, c(parent_taxas, taxa), with = F] %>% unique() %>% na.omit()
  taxa_dt[, GRP := 1:.N, by = parent_taxas]
  
  if(max(taxa_dt$GRP) > max(group2code$GRP)) {
    print(str_glue('taxa_i : {taxa_i}'))
    print(str_glue('taxa: {taxa}'))
    print(str_glue('max(taxa_dt$GRP): {max(taxa_dt$GRP)}'))
    print(str_glue('max(group2code$GRP): {max(group2code$GRP)}'))
    
    stop('Group2Code more characters is needed')
  }
  
  taxa2code = taxa_dt[group2code, on = .(GRP), nomatch = NULL][, c(taxa, 'code'), with = F]
  setnames(taxa2code, 'code', str_c(taxa, 'code', sep = '_'))
  organisms = merge(organisms, taxa2code, by = taxa, sort = F, all.x = T)
  
}

organisms[is.na(organisms)] <- 0

code_cols = taxa_levels %>% str_c(., 'code', sep = '_')
organisms[, phylogenic_code := .SD %>%
            as.matrix() %>%
            apply(., 1, str_c, collapse = '', simplify = F) %>%
            unlist(),
          .SDcols = code_cols]

organisms = unique(organisms)
organisms = merge(organisms, metadata[, .(organism_name, taxid, accessions)], by = 'organism_name', sort = F)
setcolorder(organisms, c('organism_name', 'taxid', 'accessions', taxa_levels))

fwrite(organisms, str_replace(commontree_file, '.txt$', '.parsed.csv'))
