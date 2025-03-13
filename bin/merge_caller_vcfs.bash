#!/usr/bin/env bash

set -euo pipefail

# Define input files
ALIGNERS=("bwa2a" "strobe" "sent")
SAMPLE="RIH0_ANA0-HG002-19_DBC0_0"
CALLERS=("clair3" "deep" "lfq2" "sentd" "oct")
CALLER_NAMES=("clair3" "deepvariant" "lofreq2" "sentieon" "octopus")

# Annotate each caller's VCF

for ai in "${!ALIGNERS[@]}"; do

  for i in "${!CALLERS[@]}"; do
    CALLER=${CALLERS[$i]}
    ALIGNER=${ALIGNERS[$ai]}
    CALLER_FILE="${SAMPLE}.${ALIGNER}.${CALLERS[$i]}.snv.sort.rh.vcf.gz"
    CALLER_LABEL="$ALIGNER"-"${CALLER_NAMES[$i]}"
    OUTPUT_FILE="${CALLER_LABEL}.annotated.vcf.gz"

    echo "Annotating $CALLER_FILE with caller: $CALLER_LABEL"

    #    -x INFO/AF,FORMAT/AF \

    echo "A" $CALLER .. $CALLER_FILE
    bcftools annotate \
      -h <(echo '##INFO=<ID=aligner_caller,Number=.,Type=String,Description="Source aligner-caller">') \
      --set-id +'%CHROM:%POS:%REF:%ALT' \
      --set "INFO/aligner_caller=$ALIGNER\-$CALLER" \
      -Oz \
      -o ${OUTPUT_FILE} $CALLER_FILE

    tabix $OUTPUT_FILE
    bcftools index $OUTPUT_FILE

  done
done

# Merge annotated VCFs into one
CONCAT_VCF="${SAMPLE}.concat_callers.vcf.gz"



bcftools concat -a --threads 8 -Oz \
  $(ls *anno*.gz) \
  > "$CONCAT_VCF"

# Sort and index merged VCF
tabix $CONCAT_VCF
bcftools index $CONCAT_VCF

SORTED_VCF="${SAMPLE}.merged_callers.sorted.vcf.gz"
bcftools sort -Oz -o "$SORTED_VCF" $CONCAT_VCF
tabix $SORTED_VCF
bcftools index $SORTED_VCF

bcftools norm -Oz -m+any \
  $SORTED_VCF \
  -o all_callers.merged.vcf.gz

tabix all_callers.merged.vcf.gz
bcftools index all_callers.merged.vcf.gz

bcftools sort -Oz -o all_callers.merged.sort.vcf.gz all_callers.merged.vcf.gz
bcftools index all_callers.merged.sort.vcf.gz
tabix all_callers.merged.sort.vcf.gz

# Optional stats for QC
echo "Generating stats for merged VCF"
#bcftools stats "$SORTED_VCF" > "${SORTED_VCF%.vcf.gz}.stats.txt"

echo "Merged VCF created successfully: $SORTED_VCF"

# Reminder for RTG vcfeval usage
echo "Now you can run RTG vcfeval with:"
echo "rtg vcfeval -b truth.vcf.gz -c $SORTED_VCF -t reference.sdf -o eval_output"
