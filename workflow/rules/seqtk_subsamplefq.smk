
#  https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/seqtk/subsample/pe.html

if 'subsample' in config and 'subsample_n' in config and 'subsample_seed' in config:

    rule seqtk_subsample_pe:
        input:
            f1="{sample}.1.fastq.gz",
            f2="{sample}.2.fastq.gz"
        output:
            f1="{sample}.1.subsampled.fastq.gz",
            f2="{sample}.2.subsampled.fastq.gz"
        params:
            n=3,
            seed=12345
        log:
            "logs/seqtk_subsample/{sample}.log"
        threads: 4 if 'seqtk' not in conifig else config['seqtk']['threads']
        wrapper:
            "v1.1.0/bio/seqtk/subsample/pe"
