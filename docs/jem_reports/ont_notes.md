
# ONT Aligned CRAMs, Variant Called w/Sentieon DNAscope (notes)

## Results In Brief
- The sentieon ONT-solo SNV calling pipeline is stable and we can run a very large number of CRAMs in a reproducible and scalable AWS environment.
  - a large number being on the order of 1000's of crams an hour (more honestly).
- Avg runtime for SNV calling: `28.9m`.
- Avg EC2 spot cost per cram: `$1.57`.
- Avg Coverage (and other cov stats): _rolling in this AM_
  - [ONT LSMC MultiQC Reports -coming shortly-](./)
- Avg SNP Fscore (GIAB high confidence region bed): `0.9985`.  **!!! wow !!!**
- Avg SNP Fscore (Ultima high confidence region bed): `0.964`.

### Analysis Results Data
I'll be processing these in more detail, but for now, here are the raw data files.

#### Concordance

> I expected the SNP Fscore performance in the ultima bed regions to be higher. This needs to be investigated.

- [docs/jem_reports/ont_solo_concordance.tsv](ont_solo_concordance.tsv)
 

| CmpFootprint   | SNPClass   |   Mean_Fscore |   Min_Fscore |   Max_Fscore |
|:---------------|:-----------|--------------:|-------------:|-------------:|
| ultima         | DEL_50     |      0.881195 |     0.852995 |     0.913246 |
| ultima         | INS_50     |      0.895766 |     0.879083 |     0.920608 |
| ultima         | Indel_50   |      0.864182 |     0.797539 |     0.922895 |
| ultima         | SNPts      |      0.964536 |     0.95021  |     0.99902  |
| ultima         | SNPtv      |      0.951728 |     0.933225 |     0.998589 |
| wgsHC          | DEL_50     |      0.898143 |     0.87696  |     0.94355  |
| wgsHC          | INS_50     |      0.89924  |     0.879083 |     0.940202 |
| wgsHC          | Indel_50   |      0.877399 |     0.841214 |     0.951348 |
| wgsHC          | SNPts      |      0.99883  |     0.998553 |     0.999274 |
| wgsHC          | SNPtv      |      0.998351 |     0.997778 |     0.998882 |

- I include the ultima bed segmented data for comparison to the [ultima solo results](ultima_notes.md).

> `python bin/calc_concordance_stats.py docs/jem_reports/ont_solo_concordance.tsv`

#### Benchmarking (Runtime, memory, IO, costs)
- [docs/jem_reports/ultima_benchmarks.tsv](ultima_benchmarks.tsv)

| rule                      |   avg_task_cost |   min_task_cost |   max_task_cost |   avg_minutes |   min_minutes |   max_minutes |   avg_cpu_efficiency |   min_cpu_efficiency |   max_cpu_efficiency |
|:--------------------------|----------------:|----------------:|----------------:|--------------:|--------------:|--------------:|---------------------:|---------------------:|---------------------:|
| sent.sentdont.1-24        |          1.5680 |          1.1503 |          2.3880 |       28.9300 |       25.2662 |       33.1148 |             113.9016 |             105.6283 |             126.4068 |
| sent.sentdont.concat.fofn |          0.0001 |          0.0000 |          0.0012 |        0.2191 |        0.0062 |        2.7626 |               0.0000 |               0.0000 |               0.0000 |
| sent.sentdont.concordance |          0.3626 |          0.3443 |          0.3815 |       16.0118 |       15.6672 |       16.7594 |               0.3591 |               0.0056 |               0.5569 |
| sent.sentdont.merge       |          0.0065 |          0.0063 |          0.0068 |        0.4303 |        0.4133 |        0.4513 |               1.8715 |               1.7309 |               1.9498 |

> `python bin/calc_benchmark_stats.py docs/jem_reports/ont_benchmarks.tsv`

#### VCF Files

- Data will migrate from FSx to S3 today, I'll drop the links in here.

## Job Details

### Spot Instances Jobs May End Up On
4 instance types are defined in the partition running these jobs. All are high IO, Mem 350GB-1TB, 196vcpu. Generally spot prices are $1.5 to $3.5/hr (dedicated ~$13/hr+).
 - m7i.metal-48xl 
 - m7i.48xlarge 
 - r7i.48xlarge 
 - r7i.metal-48xl 

> _note!_: the 2 base 192 vcpu instance types: `c7i.48xlarge` &  `c7i.metal-48xl`, if used, will wedge IO (so are not available to high IO demanding jobs).

#### A Representative Spot Instance
One of the spots which processed the ONT CRAMs:
- m7i.48xlarge
- 196vcpu
- $2.73/hr
- Ran for ~29min (cost of $1.365)

### ONT Input crams
I ran DNAscope on the following subset of ONT crams, I wanted at least one each for every GIAB sample (I'll be re-running the entire set once the initial data review of whats here checks out).

```bash
ls -1 /fsx/data/cram_data/ont_crams/*cram
/fsx/data/cram_data/ont_crams/HG001-sup-PAW79146.haplotagged.cram
/fsx/data/cram_data/ont_crams/HG001-sup-PAW81754.haplotagged.cram
/fsx/data/cram_data/ont_crams/HG002-sup-PAW70337.haplotagged.cram
/fsx/data/cram_data/ont_crams/HG002-sup-PAW71238.haplotagged.cram
/fsx/data/cram_data/ont_crams/HG003-sup-PAY87794.haplotagged.cram
/fsx/data/cram_data/ont_crams/HG003-sup-PAY87954.haplotagged.cram
/fsx/data/cram_data/ont_crams/HG004-sup-PAY87778.haplotagged.cram
/fsx/data/cram_data/ont_crams/HG004-sup-PAY88428.haplotagged.cram
/fsx/data/cram_data/ont_crams/HG005-sup-PAW87816.haplotagged.cram
/fsx/data/cram_data/ont_crams/HG005-sup-PAW88001.haplotagged.cram
/fsx/data/cram_data/ont_crams/HG006-sup-PAY77227.haplotagged.cram
/fsx/data/cram_data/ont_crams/HG006-sup-PBA16846.haplotagged.cram
/fsx/data/cram_data/ont_crams/HG007-sup-PAY12990.haplotagged.cram
# /fsx/data/cram_data/ont_crams/HG007-sup-PBA20413.haplotagged.cram # This sample was not processed in this batch
```

### CRAM Coverage + Other QC Stats 

> **running**

### DNAscope Commands

```bash
timestamp=$(date +%Y%m%d%H%M%S);
export TMPDIR=/fsx/scratch/sentdug_tmp_$timestamp;
mkdir -p $TMPDIR;

trap "rm -rf "$TMPDIR" || echo '$TMPDIR rm fails' >>  results/day/hg38/ONT1_ANA1-HG004_DBC0_0/align/sent/snv/sentdont/log/vcfs/ONT1_ANA1-HG004_DBC0_0.sent.sentdont.1-24.snv.log  2>&1" EXIT;

ulimit -n 65536;

# Find the jemalloc library in the active conda environment
jemalloc_path=$(find "$CONDA_PREFIX" -name "libjemalloc*" | grep -E '\.so|\.dylib' | head -n 1);

# Check if jemalloc was found and set LD_PRELOAD accordingly
if [[ -n "$jemalloc_path" ]]; then
    LD_PRELOAD="$jemalloc_path";
    echo "LD_PRELOAD set to: $LD_PRELOAD" 
else
    echo "libjemalloc not found in the active conda environment $CONDA_PREFIX.";
    exit 3;
fi

LD_PRELOAD=$LD_PRELOAD /fsx/data/cached_envs/sentieon-genomics-202503/bin/sentieon driver \
        -t 182             \
        -r /fsx/data/genomic_data/organism_references/H_sapiens/hg38/broad_hg38/Homo_sapiens_assembly38.fasta \
        -i results/day/hg38/ONT1_ANA1-HG004_DBC0_0/align/sent/ONT1_ANA1-HG004_DBC0_0.cram \            
        --interval chr1:,chr2:,chr3:,chr4:,chr5:,chr6:,chr7:,chr8:,chr9:,chr10:,chr11:,chr12:,chr13:,chr14:,chr15:,chr16:,chr17:,chr18:,chr19:,chr20:,chr21:,chr22:,chrX:,chrY:   \          
        --algo DNAscope \
        --model /fsx/data/cached_envs/sentieon-genomics-202503/bundles/DNAscopeONT2.1/diploid_model \            
        --emit_mode variant  \
                   results/day/hg38/ONT1_ANA1-HG004_DBC0_0/align/sent/snv/sentdont/vcfs/1-24/ONT1_ANA1-HG004_DBC0_0.sent.sentdont.1-24.snv.vcf.tmp \
                   >> results/day/hg38/ONT1_ANA1-HG004_DBC0_0/align/sent/snv/sentdont/log/vcfs/ONT1_ANA1-HG004_DBC0_0.sent.sentdont.1-24.snv.log 2>&1;

LD_PRELOAD=$LD_PRELOAD /fsx/data/cached_envs/sentieon-genomics-202503/bin/sentieon driver \
        -t 192             \
        -r /fsx/data/genomic_data/organism_references/H_sapiens/hg38/broad_hg38/Homo_sapiens_assembly38.fasta \
        --algo DNAModelApply  \
        --model /fsx/data/cached_envs/sentieon-genomics-202503/bundles/DNAscopeONT2.1/diploid_model           \
        -v results/day/hg38/ONT1_ANA1-HG004_DBC0_0/align/sent/snv/sentdont/vcfs/1-24/ONT1_ANA1-HG004_DBC0_0.sent.sentdont.1-24.snv.vcf.tmp \
        results/day/hg38/ONT1_ANA1-HG004_DBC0_0/align/sent/snv/sentdont/vcfs/1-24/ONT1_ANA1-HG004_DBC0_0.sent.sentdont.1-24.snv.vcf \
        >> results/day/hg38/ONT1_ANA1-HG004_DBC0_0/align/sent/snv/sentdont/log/vcfs/ONT1_ANA1-HG004_DBC0_0.sent.sentdont.1-24.snv.log 2>&1;

```

