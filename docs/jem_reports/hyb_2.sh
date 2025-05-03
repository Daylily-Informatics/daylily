
export SENTIEON_LICENSE=/fsx/data/cached_envs/Life_Sciences_Manufacturing_Corporation_eval.lic

export PATH=$PATH:/fsx/data/cached_envs/sentieon-genomics-202503/bin/

jemalloc_path=$(find "$CONDA_PREFIX" -name "libjemalloc*" | grep -E '\.so|\.dylib' | head -n 1); 
echo "JEMALLOC=$jemalloc_path"

LD_PRELOAD="$jemalloc_path";

ulimit -n 65536 

echo "starting hybrid test"

timestamp=$(date +%Y%m%d%H%M%S);
export TMPDIR=/fsx/scratch/sentdontr_tmp_$timestamp;
mkdir -p $TMPDIR;


trap "rm -rf \"$TMPDIR\" || echo '$TMPDIR rm fails' >> {log} 2>&1" EXIT;

LD_PRELOAD=$LD_PRELOAD sentieon-cli -v dnascope-hybrid \
	     -t 186 \
	     -r /fsx/data/genomic_data/organism_references/H_sapiens/hg38/broad_hg38/Homo_sapiens_assembly38.fasta \
	     --sr_r1_fastq /fsx/data/genomic_data/organism_reads/H_sapiens/giab/google-brain/novaseq/pcr_free/30x/30x/HG001.novaseq.pcr-free.30x.R1.fastq.gz \
	     --sr_r2_fastq /fsx/data/genomic_data/organism_reads/H_sapiens/giab/google-brain/novaseq/pcr_free/30x/30x/HG001.novaseq.pcr-free.30x.R2.fastq.gz \
	     --sr_readgroups "@RG\tID:hg001-1\tSM:hg001\tLB:hg001-LB-1\tPL:ILLUMINA" \
	     --lr_aln /fsx/data/cram_data/ont_crams/HG001-sup-PAW81754.haplotagged.cram \
             --lr_align_input \
             --lr_input_ref /fsx/data/genomic_data/organism_references/H_sapiens/hg38/broad_hg38/Homo_sapiens_assembly38.fasta \
	     -m /fsx/data/cached_envs/sentieon-genomics-202503/bundles/HybridIlluminaONT1.0.bundle \
	     -b /fsx/data/genomic_data/organism_annotations/H_sapiens/hg38/controls/giab/snv/v4.2.1/HG001/wgsHC/HG001.bed \
	     --longread_tech ONT \
	     HG001_hybrid.vcf.gz

echo "hybrid test complete"
