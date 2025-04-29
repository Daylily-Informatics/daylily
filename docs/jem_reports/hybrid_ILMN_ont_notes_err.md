


# Sentieon Hybrid ILMN+ONT Pipeline Notes

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

### HG001 GIAB HC Regions.bed
```bash
ls -lth /fsx/data/genomic_data/organism_annotations/H_sapiens/hg38/controls/giab/snv/v4.2.1/HG001/wgsHC/HG001.bed
-rwxr-xr-x 1 root root 15M Feb 19 06:37 /fsx/data/genomic_data/organism_annotations/H_sapiens/hg38/controls/giab/snv/v4.2.1/HG001/wgsHC/HG001.bed
```

### Sentieon Bundle
```bash
ls -lth /fsx/data/cached_envs/sentieon-genomics-202503/bundles/HybridIlluminaONT1.0.bundle
-rwxr-xr-x 1 root root 158M Apr 16 08:24 /fsx/data/cached_envs/sentieon-genomics-202503/bundles/HybridIlluminaONT1.0.bundle
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

### Calling the pipeline

```bash
conda activate sentHyb
source hyb.sh > hyb_test.log 2>&1
```

## Pipeline Output

### Log
- [./docs/jem_reports/hyb_test.log](hyb_test.log)
  

#### Log Tail

```bash

tail -n 30 hyb_test.log

  features: 1f8bfbff fffa320b
  extended: f1bf07ab 00407f5a
  amd bits: 2c100800 00000121
     brand: Intel(R) Xeon(R) Platinum 8375C CPU @ 2.90GHz
threads: 120 max 128
algo: util-sort
license: sentieon:util=1
output file size: 3430975612
output reads: 30048671
bam_mem_sort: 51 calls 24.019 user 0.853 sys 24.973 real
bam_write: 51 calls 18.342 user 2.891 sys 64.560 real
execute: 1 calls 9.854 user 3.944 sys 452.427 real
merge_files: 1 calls 9.840 user 3.933 sys 34.611 real
parse_chunk: 802 calls 30.888 user 12.719 sys 45.921 real
read_chunk: 2299 calls 14.907 user 14.134 sys 389.124 real
sort_block: 1 calls 11.242 user 13.759 sys 357.040 real
write_chunk: 1495 calls 9.589 user 2.572 sys 12.293 real
overall: 5441556480 mem 760.045 user 88.229 sys 452.444 real
INFO:sentieon_cli.executor:Running: sentieon driver --input /fsx/scratch/tmpszj20ekz/hybrid_stage1.bam --input /fsx/scratch/tmpszj20ekz/stage1_hap.bam --reference /fsx/data/genomic_data/organism_references/H_sapiens/hg38/broad_hg38/Homo_sapiens_assembly38.fasta --thread_count 120 --interval_padding 0 --algo HybridStage2 --model /fsx/data/cached_envs/sentieon-genomics-202503/bundles/HybridIlluminaONT1.0.bundle/HybridStage2.model --mode 0 --all_bed /fsx/scratch/tmpszj20ekz/hybrid_stage2.bed --hap_bed /fsx/scratch/tmpszj20ekz/stage1_hap.bed --unmap_bam /fsx/scratch/tmpszj20ekz/hybrid_stage2_unmap.bam --alt_bam /fsx/scratch/tmpszj20ekz/hybrid_stage2_alt.bam --min_kmer_size 80 --prune_factor 5 --target_extension 30
INFO:sentieon_cli.executor:Running: rm /fsx/scratch/tmpszj20ekz/stage1_nocoor.fa /fsx/scratch/tmpszj20ekz/stage1_nocoor.bed /fsx/scratch/tmpszj20ekz/stage1_ins.fa /fsx/scratch/tmpszj20ekz/stage1_ins.bed /fsx/scratch/tmpszj20ekz/stage1_hap.vcf
rm: cannot remove '/fsx/scratch/tmpszj20ekz/stage1_nocoor.fa': No such file or directory
rm: cannot remove '/fsx/scratch/tmpszj20ekz/stage1_nocoor.bed': No such file or directory
cmdline: /fsx/data/cached_envs/sentieon-genomics-202503/libexec/driver --input /fsx/scratch/tmpszj20ekz/hybrid_stage1.bam --input /fsx/scratch/tmpszj20ekz/stage1_hap.bam --reference /fsx/data/genomic_data/organism_references/H_sapiens/hg38/broad_hg38/Homo_sapiens_assembly38.fasta --thread_count 120 --interval_padding 0 --algo HybridStage2 --model /fsx/data/cached_envs/sentieon-genomics-202503/bundles/HybridIlluminaONT1.0.bundle/HybridStage2.model --mode 0 --all_bed /fsx/scratch/tmpszj20ekz/hybrid_stage2.bed --hap_bed /fsx/scratch/tmpszj20ekz/stage1_hap.bed --unmap_bam /fsx/scratch/tmpszj20ekz/hybrid_stage2_unmap.bam --alt_bam /fsx/scratch/tmpszj20ekz/hybrid_stage2_alt.bam --min_kmer_size 80 --prune_factor 5 --target_extension 30
rm: cannot remove '/fsx/scratch/tmpszj20ekz/stage1_hap.vcf': No such file or directory
This software is licensed to johnm@lsmc.com by Sentieon Inc.
Error: Failed to open /fsx/scratch/tmpszj20ekz/stage1_hap.bam: No such file or directory

ERROR:sentieon_cli.executor:Error running command, 'sentieon driver --input /fsx/scratch/tmpszj20ekz/hybrid_stage1.bam --input /fsx/scratch/tmpszj20ekz/stage1_hap.bam --reference /fsx/data/genomic_data/organism_references/H_sapiens/hg38/broad_hg38/Homo_sapiens_assembly38.fasta --thread_count 120 --interval_padding 0 --algo HybridStage2 --model /fsx/data/cached_envs/sentieon-genomics-202503/bundles/HybridIlluminaONT1.0.bundle/HybridStage2.model --mode 0 --all_bed /fsx/scratch/tmpszj20ekz/hybrid_stage2.bed --hap_bed /fsx/scratch/tmpszj20ekz/stage1_hap.bed --unmap_bam /fsx/scratch/tmpszj20ekz/hybrid_stage2_unmap.bam --alt_bam /fsx/scratch/tmpszj20ekz/hybrid_stage2_alt.bam --min_kmer_size 80 --prune_factor 5 --target_extension 30'
CommandError: Execution failed

```


### Hyb Pipeline Run Dir Contents

```bash
ls -lt
total 14247188
-rw-rw-r-- 1 ubuntu ubuntu     2075073 Apr 29 07:10 hyb_test.log
-rw-rw-r-- 1 ubuntu ubuntu       22693 Apr 29 06:33 HG001_hybrid.cnv.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu         157 Apr 29 06:33 HG001_hybrid.cnv.vcf.gz.tbi
drwxrwxr-x 3 ubuntu ubuntu       33280 Apr 29 06:33 HG001_hybrid_metrics
-rw-rw-r-- 1 ubuntu ubuntu 14484932828 Apr 29 06:29 HG001_hybrid_deduped.cram
-rw-rw-r-- 1 ubuntu ubuntu     2745264 Apr 29 06:29 HG001_hybrid_deduped.cram.bai
-rw-rw-r-- 1 ubuntu ubuntu     1328307 Apr 29 06:29 HG001_hybrid_deduped.cram.crai
-rw-rw-r-- 1 ubuntu ubuntu     1366922 Apr 29 06:26 HG001_hybrid.sv.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu       56294 Apr 29 06:26 HG001_hybrid.sv.vcf.gz.tbi
-rw-rw-r-- 1 ubuntu ubuntu     1301811 Apr 29 05:59 HG001_hybrid_mosdepth_0.mosdepth.global.dist.txt
-rw-rw-r-- 1 ubuntu ubuntu      470175 Apr 29 05:59 HG001_hybrid_mosdepth_0.mosdepth.region.dist.txt
-rw-rw-r-- 1 ubuntu ubuntu       13727 Apr 29 05:59 HG001_hybrid_mosdepth_0.mosdepth.summary.txt
-rw-rw-r-- 1 ubuntu ubuntu    34783653 Apr 29 05:59 HG001_hybrid_mosdepth_0.regions.bed.gz
-rw-rw-r-- 1 ubuntu ubuntu        9229 Apr 29 05:59 HG001_hybrid_mosdepth_0.regions.bed.gz.csi
-rw-rw-r-- 1 ubuntu ubuntu    49252564 Apr 29 05:59 HG001_hybrid_mosdepth_0.thresholds.bed.gz
-rw-rw-r-- 1 ubuntu ubuntu        9459 Apr 29 05:59 HG001_hybrid_mosdepth_0.thresholds.bed.gz.csi
-rw-rw-r-- 1 ubuntu ubuntu        1318 Apr 29 05:54 hyb.sh
-rw-rw-r-- 1 ubuntu ubuntu        1221 Apr 29 05:52 hyb.sh~
drwxrwxr-x 2 ubuntu ubuntu       33280 Apr 29 05:47 tmp

```


# Other Notes

- The disk had not filled up, there is a ton of open space.

```bash
df -h .
10.0.1.237@tcp:/w55azb4v  6.6T  2.2T  4.5T  33% /fsx
```

