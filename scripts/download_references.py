import os
import urllib.request
import subprocess
import pandas as pd
from tqdm import tqdm

pd.options.mode.chained_assignment = None

##############################
#### Download RefSeq Data ####
##############################

### Create directories
os.makedirs('data/refseq', exist_ok=True)
os.makedirs('data/refseq/annotation', exist_ok=True)
os.makedirs('data/refseq/genome', exist_ok=True)
os.makedirs('data/refseq/passed-transcripts', exist_ok=True)

### Download RefSeq assembly summary to get ftp paths to genomes and annotations
refseq_assembly_summary_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt'
refseq_assembly_summary_file = 'data/refseq/metadata/assembly_summary_refseq.txt'
urllib.request.urlretrieve(refseq_assembly_summary_url, refseq_assembly_summary_file)
# refseq_assembly_summary_file = 'data/refseq/assembly_summary_refseq_primates.txt'

### Filter species to vertebrates
asm = pd.read_csv(refseq_assembly_summary_file, sep='\t', skiprows=1, header=0, low_memory=False)
# asm = pd.read_csv(refseq_assembly_summary_file, header=0, low_memory=False)

asm = asm[asm.group.isin(['vertebrate_mammalian', 'vertebrate_other'])]
# asm = asm[asm.group.isin(['vertebrate_mammalian'])]


### Filter refseq genome category
## From README_assembly_summary.txt:
## Column  5: "refseq_category"
##    RefSeq Category: whether the assembly is a reference or representative genome
##    in the NCBI Reference Sequence (RefSeq) project classification. 
##    Values:
##            reference genome      - a manually selected high quality genome 
##                                    assembly that NCBI and the community have 
##                                    identified as being important as a standard 
##                                    against which other data are compared
##            representative genome - a genome computationally or manually selected
##                                    as a representative from among the best 
##                                    genomes available for a species or clade that
##                                    does not have a designated reference genome
##            na                    - no RefSeq category assigned to this assembly
##    Prokaryotes may have more than one reference or representative genome per 
##    species.
asm = asm[asm.refseq_category.isin(['reference genome', 'representative genome'])]

ftp_path = asm.ftp_path
accession = [filename[9] for filename in ftp_path.str.split('/')]

gtf_urls = [ftp + '/' + acc + '_genomic.gtf.gz' for ftp, acc in zip(ftp_path, accession)]
fna_urls = [ftp + '/' + acc + '_genomic.fna.gz' for ftp, acc in zip(ftp_path, accession)]

gtf_files = ['data/refseq/annotation/' + acc + '_genomic.gtf.gz' for acc in accession]
fna_files = ['data/refseq/genome/' + acc + '_genomic.fna.gz' for acc in accession]

def download_file(url, output_file):
    try:
        urllib.request.urlretrieve(url, output_file)
        return True
    except urllib.error.HTTPError as e:
        print(f"HTTP Error {e.code}: {url}")
        return False
    except (urllib.error.URLError, Exception) as e:
        print(f"Error downloading {url}: {e}")
        return False

for gtf_file, gtf_url in zip(gtf_files, gtf_urls):
    if not os.path.exists(gtf_file):
        print(gtf_url)
        download_file(gtf_url, gtf_file)

for fna_file, fna_url in zip(fna_files, fna_urls):
    if not os.path.exists(fna_file):
        print(fna_url)
        download_file(fna_url, fna_file)
