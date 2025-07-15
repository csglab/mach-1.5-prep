import os
import urllib.request
import subprocess
import pandas as pd
from tqdm import tqdm

# import warnings
pd.options.mode.chained_assignment = None
# warnings.filterwarnings('ignore', category=FutureWarning, message='.*Downcasting object dtype arrays.*')

##############################
#### Download RefSeq Data ####
##############################

### Create directories
os.makedirs('data/refseq/passed-transcripts', exist_ok=True)

refseq_assembly_summary_file = 'data/refseq/metadata/assembly_summary_refseq.txt'
asm = pd.read_csv(refseq_assembly_summary_file, sep='\t', skiprows=1, header=0, low_memory=False)

# refseq_assembly_summary_file = 'data/refseq/metadata/assembly_summary_refseq_primates.txt'
# asm = pd.read_csv(refseq_assembly_summary_file, header=0, low_memory=False)

asm = asm[asm.group.isin(['vertebrate_mammalian', 'vertebrate_other'])]
# asm = asm[asm.group.isin(['vertebrate_mammalian'])]
asm = asm[asm.refseq_category.isin(['reference genome', 'representative genome'])]

ftp_path = asm.ftp_path
accession = [filename[9] for filename in ftp_path.str.split('/')]

gtf_urls = [ftp + '/' + acc + '_genomic.gtf.gz' for ftp, acc in zip(ftp_path, accession)]
fna_urls = [ftp + '/' + acc + '_genomic.fna.gz' for ftp, acc in zip(ftp_path, accession)]

gtf_files = ['data/refseq/annotation/' + acc + '_genomic.gtf.gz' for acc in accession]

### Extract splice sites regions from each species
pbar = tqdm(total = len(gtf_files), smoothing=0)
for gtf_file in gtf_files:

    print('=' * 30)
    accession_id = os.path.basename(gtf_file).replace('.gtf.gz', '')
    pbar.set_description(accession_id)

    output_csv = 'data/refseq/passed-transcripts/' + accession_id + '.passed_transcripts.csv'
    if os.path.exists(output_csv):
        print(f'{output_csv} already exists.')
        pbar.update(1)
        continue

    ### Read GTF file
    try:
        gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)
    except FileNotFoundError:
        print(f"Error: Could not find the GTF file at path: {gtf_file}")
        pbar.update(1)
        continue

    gtf.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    ### Extract transcript annotations from GTF file
    trs_annot = gtf[gtf.feature == 'transcript']
    trs_annot = pd.concat([
        trs_annot.start,
        trs_annot.end,
        trs_annot.attribute.str.extract('transcript_id "(?P<transcript_id>[^"]+)"', expand=True),
        trs_annot.attribute.str.extract('transcript_biotype "(?P<transcript_biotype>[^"]+)"', expand=True),
        trs_annot.attribute.str.extract('gene_id "(?P<gene_id>[^"]+)"', expand=True),
        trs_annot.attribute.str.extract('model_evidence "(?P<model_evidence>[^"]+)"', expand=True)], axis=1)

    ### Filter on protein-coding transcripts only
    trs_annot = trs_annot[trs_annot.transcript_biotype == 'mRNA']

    ## Some species, there is no 'transcript' annotation in their GTF file,
    ## We decided to ignore these species for now
    num_trs = trs_annot.shape[0]
    if num_trs == 0:
        print('No transcripts found in GTF file for ' + gtf_file)
        pbar.update(1)
        continue

    curated_trs = trs_annot.loc[trs_annot.transcript_id.str.contains('NM_'), ['transcript_id', 'gene_id']]
    print(f'Curated transcripts: {curated_trs.shape[0]}')
    
    predicted_trs = trs_annot.loc[trs_annot.transcript_id.str.contains('XM_'), ['transcript_id', 'model_evidence', 'gene_id']]
    print(f'Predicted transcripts: {predicted_trs.shape[0]}')

    ### Extrat Gnomon model evidence for predicted transcripts
    predicted_trs = pd.concat([
        predicted_trs,
        predicted_trs['model_evidence'].str.extract('(?P<perc_covered_features>[0-9]+)% coverage of the annotated genomic feature by RNAseq alignments', expand=True).infer_objects(copy=False).fillna(0).astype(int),
        predicted_trs['model_evidence'].str.extract(r'including (?P<num_samples_support>[0-9]+) sample(?:s)? with support for all annotated introns', expand=True).infer_objects(copy=False).fillna(0).astype(int)],
        axis=1)    
    ## Passed predicted transcripts are those with:
    ##  100% coverage of the annotated genomic feature by RNAseq alignments, and
    ##  at least 1 sample with support for all annotated introns
    passed_predicted_trs = predicted_trs.loc[(predicted_trs.perc_covered_features == 100) & (predicted_trs.num_samples_support >= 1), ['transcript_id', 'gene_id']]
    print(f'Passed predicted transcripts: {passed_predicted_trs.shape[0]}')

    all_passed_trs = pd.concat([curated_trs, passed_predicted_trs], ignore_index=True)
    print(f'All passed transcripts: {all_passed_trs.shape[0]}')

    if all_passed_trs.shape[0] == 0:
        print('No curated or passed predicted transcripts for ' + gtf_file)
        pbar.update(1)
        continue

    gtf['transcript_id'] = gtf.attribute.str.extract('transcript_id "(?P<transcript_id>[^"]+)"', expand=True)
    gtf['gene_id'] = gtf.attribute.str.extract('gene_id "(?P<gene_id>[^"]+)"', expand=True)

    gtf = gtf.loc[(gtf.transcript_id.isin(all_passed_trs.transcript_id)) & (gtf.gene_id.isin(all_passed_trs.gene_id))]

    gtf.to_csv(output_csv, index=False)
    pbar.update(1)
