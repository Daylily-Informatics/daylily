import urllib.request

url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz'
filename = 'GRCh38.fa.gz'

print('Downloading genome...')
urllib.request.urlretrieve(url, filename)
print('Genome downloaded and saved to', filename)
