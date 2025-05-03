
# Ultima Aligned CRAMs, Variant Called w/Sentieon DNAscope (notes)

## Results In Brief
- The sentieon Ultima-solo SNV calling pipeline is stable and we can run a very large number of CRAMs in a reproducible and scalable AWS environment.
- We can process 1000's of crams an hour (more honestly).
- Avg runtime for SNV calling: `38.4m`.
- Avg EC2 spot cost per cram: `$1.98`.
- Avg Coverage (and other cov stats): _rolling in this AM_
  - [Ultima LSMC MultiQC Reports -coming shortly-](./)
- Avg SNP Fscore (GIAB high confidence region bed): `0.994`.
- Avg SNP Fscore (Ultima high confidence region bed): `0.97`.

### Analysis Results Data
I'll be processing these in more detail, but for now, here are the raw data files.

#### Concordance

> I expected the SNP Fscore performance in the ultima bed regions to be higher. This needs to be investigated.

- [docs/jem_reports/ultima_solo_concordance.tsv](ultima_solo_concordance.tsv)
 
| CmpFootprint   | SNPClass   |   Mean_Fscore |   Min_Fscore |   Max_Fscore |
|:---------------|:-----------|--------------:|-------------:|-------------:|
| ultima         | DEL_50     |      0.886204 |     0.855814 |     0.91461  |
| ultima         | INS_50     |      0.908542 |     0.885621 |     0.923475 |
| ultima         | Indel_50   |      0.877434 |     0.801254 |     0.94875  |
| ultima         | SNPts      |      0.971477 |     0.957412 |     0.994657 |
| ultima         | SNPtv      |      0.965741 |     0.94875  |     0.992529 |
| wgsHC          | DEL_50     |      0.904299 |     0.881632 |     0.943329 |
| wgsHC          | INS_50     |      0.930903 |     0.912714 |     0.959161 |
| wgsHC          | Indel_50   |      0.847819 |     0.801254 |     0.933683 |
| wgsHC          | SNPts      |      0.994608 |     0.993956 |     0.995441 |
| wgsHC          | SNPtv      |      0.992171 |     0.990928 |     0.99303  |

> `python bin/calc_concordance_stats.py docs/jem_reports/ont_solo_concordance.tsv`

#### Benchmarking (Runtime, memory, IO, costs)
- [docs/jem_reports/ultima_benchmarks.tsv](ultima_benchmarks.tsv)

| rule                     |   avg_task_cost |   min_task_cost |   max_task_cost |   avg_minutes |   min_minutes |   max_minutes |   avg_cpu_efficiency |   min_cpu_efficiency |   max_cpu_efficiency |
|:-------------------------|----------------:|----------------:|----------------:|--------------:|--------------:|--------------:|---------------------:|---------------------:|---------------------:|
| sent.sentdug.1-24        |          1.9834 |          1.7239 |          2.4139 |       38.4436 |       34.9796 |       42.7305 |              21.4038 |              19.6914 |              23.1942 |
| sent.sentdug.concat.fofn |          0.0000 |          0.0000 |          0.0000 |        0.0061 |        0.0058 |        0.0067 |               0.0000 |               0.0000 |               0.0000 |
| sent.sentdug.concordance |          0.4352 |          0.4127 |          0.4637 |       16.4780 |       16.1782 |       17.2698 |               0.3764 |               0.0062 |               0.6353 |
| sent.sentdug.merge       |          0.0068 |          0.0067 |          0.0070 |        0.4027 |        0.3947 |        0.4113 |               2.1757 |               2.1102 |               2.2122 |

> `python bin/calc_benchmark_stats.py docs/jem_reports/ultima_benchmarks.tsv`

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
One of the spots which processed the Ultima CRAMs:
- m7i.48xlarge
- 196vcpu
- $2.73/hr
- Ran for ~38min (cost of $1.365)

### Ultima Input crams
I ran DNAscope on the following subset of Ultima crams, I wanted at least one each for every GIAB sample (I'll be re-running the entire set once the initial data review of whats here checks out).

```bash
ls -1 /fsx/data/cram_data/ug_init/*cram
/fsx/data/cram_data/ug_init/408622-NA12878-Z0025-CTCGAGATTGATGAT-1pctDS.cram
/fsx/data/cram_data/ug_init/408622-NA12878-Z0025-CTCGAGATTGATGAT.cram
/fsx/data/cram_data/ug_init/408622-NA24143-Z0022-CGGAGTGCATGCATGAT.cram
/fsx/data/cram_data/ug_init/408622-NA24143-Z0113-CAGTTCATCTGTGAT.cram
/fsx/data/cram_data/ug_init/408622-NA24149-Z0008-CACATCCTGCATGTGAT.cram
/fsx/data/cram_data/ug_init/408622-NA24385-Z0032-CTCTGTATTGCAGAT.cram
/fsx/data/cram_data/ug_init/408622-NA24631-Z0016-CATCCTGTGCGCATGAT.cram
/fsx/data/cram_data/ug_init/408622-NA24631-Z0115-CGGCTAGATGCAGAT.cram
/fsx/data/cram_data/ug_init/408622-NA24694-Z0116-CAGTTATGTGCTGAT.cram
/fsx/data/cram_data/ug_init/408622-NA24695-Z0120-CTGCTGCGGAGCATGAT.cram
```

### CRAM Coverage + Other QC Stats 

> **running**

### DNAscope Commands

```bash
timestamp=$(date +%Y%m%d%H%M%S);
export TMPDIR=/fsx/scratch/sentdug_tmp_$timestamp;
mkdir -p $TMPDIR;

trap "rm -rf "$TMPDIR" || echo '$TMPDIR rm fails' >> results/day/hg38/UG1_ANA1-HG004_DBC0_0/align/sent/snv/sentdug/log/vcfs/UG1_ANA1-HG004_DBC0_0.sent.sentdug.1-24.snv.log 2>&1" EXIT;

ulimit -n 65536;


# Find the jemalloc library in the active conda environment
jemalloc_path=$(find "$CONDA_PREFIX" -name "libjemalloc*" | grep -E '\.so|\.dylib' | head -n 1);

# Check if jemalloc was found and set LD_PRELOAD accordingly
if [[ -n "$jemalloc_path" ]]; then
    LD_PRELOAD="$jemalloc_path";
    echo "LD_PRELOAD set to: $LD_PRELOAD" >> results/day/hg38/UG1_ANA1-HG004_DBC0_0/align/sent/snv/sentdug/log/vcfs/UG1_ANA1-HG004_DBC0_0.sent.sentdug.1-24.snv.log;
else
    echo "libjemalloc not found in the active conda environment $CONDA_PREFIX.";
    exit 3;
fi

LD_PRELOAD=$LD_PRELOAD /fsx/data/cached_envs/sentieon-genomics-202503/bin/sentieon driver \
        -t 182 \            
        -r /fsx/data/genomic_data/organism_references/H_sapiens/hg38/broad_hg38/Homo_sapiens_assembly38.fasta \
        -i results/day/hg38/UG1_ANA1-HG004_DBC0_0/align/sent/UG1_ANA1-HG004_DBC0_0.cram \            
        --interval chr1:,chr2:,chr3:,chr4:,chr5:,chr6:,chr7:,chr8:,chr9:,chr10:,chr11:,chr12:,chr13:,chr14:,chr15:,chr16:,chr17:,chr18:,chr19:,chr20:,chr21:,chr22:,chrX:,chrY: \             
        --read_filter UltimaReadFilter \            
        --algo DNAscope \
        --model "/fsx/data/cached_envs/sentieon-genomics-202503/bundles/SentieonUltima1.0/dnascope.model" \
        --pcr_indel_model none \            
        --emit_mode variant \            
        results/day/hg38/UG1_ANA1-HG004_DBC0_0/align/sent/snv/sentdug/vcfs/1-24/UG1_ANA1-HG004_DBC0_0.sent.sentdug.1-24.snv.gvcf \
            >> results/day/hg38/UG1_ANA1-HG004_DBC0_0/align/sent/snv/sentdug/log/vcfs/UG1_ANA1-HG004_DBC0_0.sent.sentdug.1-24.snv.log 2>&1;

LD_PRELOAD=$LD_PRELOAD /fsx/data/cached_envs/sentieon-genomics-202503/bin/sentieon driver \
        -t 182 \            
        -r /fsx/data/genomic_data/organism_references/H_sapiens/hg38/broad_hg38/Homo_sapiens_assembly38.fasta \
        --algo DNAModelApply \            
        --model /fsx/data/cached_envs/sentieon-genomics-202503/bundles/SentieonUltima1.0/dnascope.model \            
        -v results/day/hg38/UG1_ANA1-HG004_DBC0_0/align/sent/snv/sentdug/vcfs/1-24/UG1_ANA1-HG004_DBC0_0.sent.sentdug.1-24.snv.gvcf \
        results/day/hg38/UG1_ANA1-HG004_DBC0_0/align/sent/snv/sentdug/vcfs/1-24/UG1_ANA1-HG004_DBC0_0.sent.sentdug.1-24.snv.vcf \
            >> results/day/hg38/UG1_ANA1-HG004_DBC0_0/align/sent/snv/sentdug/log/vcfs/UG1_ANA1-HG004_DBC0_0.sent.sentdug.1-24.snv.log 2>&1;

```

