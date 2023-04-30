import os
import subprocess

# Download genome
#url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz'
#filename = 'GRCh38.fa.gz'
#if not os.path.isfile(filename):
#    print('Downloading genome...')
#    subprocess.run(['wget', url])
#    subprocess.run(['gunzip', filename + '.gz'])

#os.system('mkdir -p ./bwa_mem2')
#os.system('mkdir -p ./bwa_mem2_ert')
#os.system('mkdir -p ./bwa_sent')
#os.system('cp GRCh38.fa.gz {./bwa_mem2,./bwa_mem2_ert,./bwa_sent/}')


import sys
fa=sys.argv[1]

os.chdir('bwa_mem2')
# Build BWA-MEM2 index
print('Building BWA-MEM2 index...')
subprocess.run(['../../../resources/bwa-mem2/bwa-mem2', 'index', '-t', '10', '-p', f'{fa}', '-a', 'bwtsw', f'{fa}'])

os.chdir('../bwa_mem2_ert')
print('Building BWA-MEM2 ert index...')
subprocess.run(['../../../resources/ert/bwa-mem2', 'index', '-t', '10', '-p', f'{fa}', '-a', 'ert', f'{fa}'])

os.chdir('../bwa_mem_sent')
# Build Sentieon BWA mem index
print('Building Sentieon BWA mem index...')
subprocess.run(['sentieon', 'bwa', 'index',  f'{fa}'])
