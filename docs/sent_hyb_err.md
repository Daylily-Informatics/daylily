


# Sentieon Hybrid Pipeline Notes

I am encountering an error towards the end(?) of the pipeline.
```
ERROR:sentieon_cli.executor:Error running command, 'sentieon driver --input /fsx/scratch/tmpszj20ekz/hybrid_stage1.bam --input /fsx/scratch/tmpszj20ekz/stage1_hap.bam --reference /fsx/data/genomic_data/organism_references/H_sapiens/hg38/broad_hg38/Homo_sapiens_assembly38.fasta --thread_count 120 --interval_padding 0 --algo HybridStage2 --model /fsx/data/cached_envs/sentieon-genomics-202503/bundles/HybridIlluminaONT1.0.bundle/HybridStage2.model --mode 0 --all_bed /fsx/scratch/tmpszj20ekz/hybrid_stage2.bed --hap_bed /fsx/scratch/tmpszj20ekz/stage1_hap.bed --unmap_bam /fsx/scratch/tmpszj20ekz/hybrid_stage2_unmap.bam --alt_bam /fsx/scratch/tmpszj20ekz/hybrid_stage2_alt.bam --min_kmer_size 80 --prune_factor 5 --target_extension 30'
```

## Spot Instance Details

Start Datetime: 2025-04-29 05:29:28  
Region: us-west-2 
AZ: us-west-2d
Type: c6i.32xlarge
Spot price: 1.732800 USD/hour

## Sentieon CLI Conda Env

### sentHyb.yaml

```yaml
channels:
  - conda-forge
  - bioconda
  - anaconda
  - r
  - defaults
dependencies:
  - mbuffer
  - curl
  - samtools
  - jemalloc
  - htslib
  - pigz
  - isa-l==2.31.0
  - bedtools
  - bcftools>=1.10
  - multiqc==1.25.2
  - mosdepth==0.3.2
  - pip
  - pip:
      - https://github.com/Sentieon/sentieon-cli/releases/download/v1.2.1/sentieon_cli-1.2.1.tar.gz      
```

### Create Conda Env

```bash
conda env create -n sentHyb -f sentHyb.yaml
```


## Pipeline Inputs

### ILMN 30x HG001 fastqs
```bash
ls -lth /fsx/data/genomic_data/organism_reads/H_sapiens/giab/google-brain/novaseq/pcr_free/30x/30x/HG001.novaseq.pcr-free.30x.R*
-rwxr-xr-x 1 root root 25G Feb 19 06:51 /fsx/data/genomic_data/organism_reads/H_sapiens/giab/google-brain/novaseq/pcr_free/30x/30x/HG001.novaseq.pcr-free.30x.R1.fastq.gz
-rwxr-xr-x 1 root root 26G Feb 19 06:51 /fsx/data/genomic_data/organism_reads/H_sapiens/giab/google-brain/novaseq/pcr_free/30x/30x/HG001.novaseq.pcr-free.30x.R2.fastq.gz
```

### ONT HG001 cram

```bash
ls -lth  /fsx/data/cram_data/ont_crams/HG001-sup-PAW81754.haplotagged.cram*
-rwxr-xr-x 1 root root 434K Apr 28 07:24 /fsx/data/cram_data/ont_crams/HG001-sup-PAW81754.haplotagged.cram.crai
-rwxr-xr-x 1 root root  72G Apr 27 23:37 /fsx/data/cram_data/ont_crams/HG001-sup-PAW81754.haplotagged.cram
```

## Pipeline Invocation

### hyb_test.sh

```bash

export SENTIEON_LICENSE=/fsx/data/cached_envs/Life_Sciences_Manufacturing_Corporation_eval.lic

export PATH=$PATH:/fsx/data/cached_envs/sentieon-genomics-202503/bin/

jemalloc_path=$(find "$CONDA_PREFIX" -name "libjemalloc*" | grep -E '\.so|\.dylib' | head -n 1);
echo "JEMALLOC=$jemalloc_path"

LD_PRELOAD="$jemalloc_path";

ulimit -n 65536

echo "starting hybrid test"

LD_PRELOAD=$LD_PRELOAD sentieon-cli -v dnascope-hybrid \
             -t 120 \
             -r /fsx/data/genomic_data/organism_references/H_sapiens/hg38/broad_hg38/Homo_sapiens_assembly38.fasta \
             --sr_r1_fastq /fsx/data/genomic_data/organism_reads/H_sapiens/giab/google-brain/novaseq/pcr_free/30x/30x/HG001.novaseq.pcr-free.30x.R1.fastq.gz \
             --sr_r2_fastq /fsx/data/genomic_data/organism_reads/H_sapiens/giab/google-brain/novaseq/pcr_free/30x/30x/HG001.novaseq.pcr-free.30x.R2.fastq.gz \
             --sr_readgroups "@RG\tID:hg001-1\tSM:hg001\tLB:hg001-LB-1\tPL:ILLUMINA" \
             --lr_aln /fsx/data/cram_data/ont_crams/HG001-sup-PAW81754.haplotagged.cram \
             -m /fsx/data/cached_envs/sentieon-genomics-202503/bundles/HybridIlluminaONT1.0.bundle \
             -b /fsx/data/genomic_data/organism_annotations/H_sapiens/hg38/controls/giab/snv/v4.2.1/HG001/wgsHC/HG001.bed \
             --longread_tech ONT \
             HG001_hybrid.vcf.gz

echo "hybrid test complete"

```
Calling the pipeline
```

conda activate sentHyb
source hyb.sh > hyb_test.log 2>&1


  