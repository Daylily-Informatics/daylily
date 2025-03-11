#!/usr/bin/env bash

set -euo pipefail

# Define input files
SAMPLE="RIH0_ANA0-HG002-19_DBC0_0.bwa2a"
CALLERS=("clair3" "deep" "lfq2" "oct" "sentd")
CALLER_NAMES=("clair3" "deepvariant" "lofreq2" "octopus" "sentieon")

# Annotate each caller's VCF
for i in "${!CALLERS[@]}"; do
  CALLER_FILE="${SAMPLE}.${CALLERS[$i]}.snv.sort.vcf.gz"
  CALLER_LABEL="${CALLER_NAMES[$i]}"
  OUTPUT_FILE="${CALLER_LABEL}.annotated.vcf.gz"

  echo "Annotating $CALLER_FILE with caller: $CALLER_LABEL"
  bcftools annotate --threads 4 -Oz \
    --set-id '+' \
    --annotations "INFO/caller=${CALLER_LABEL}" \
    "$CALLER_FILE" \
    > "$OUTPUT_FILE"

  bcftools index "$OUTPUT_FILE"
done

# Merge annotated VCFs into one
MERGED_VCF="${SAMPLE}.merged_callers.vcf.gz"

bcftools merge --threads 8 -Oz \
  --info-rules caller:join \
  clair3.annotated.vcf.gz \
  deepvariant.annotated.vcf.gz \
  lofreq2.annotated.vcf.gz \
  octopus.annotated.vcf.gz \
  sentieon.annotated.vcf.gz \
  > "$MERGED_VCF"

# Sort and index merged VCF
SORTED_VCF="${SAMPLE}.merged_callers.sorted.vcf.gz"
bcftools sort -Oz -o "$SORTED_VCF" "$MERGED_VCF"
bcftools index "$SORTED_VCF"

# Optional stats for QC
echo "Generating stats for merged VCF"
bcftools stats "$SORTED_VCF" > "${SORTED_VCF%.vcf.gz}.stats.txt"

echo "Merged VCF created successfully: $SORTED_VCF"

# Reminder for RTG vcfeval usage
echo "Now you can run RTG vcfeval with:"
echo "rtg vcfeval -b truth.vcf.gz -c $SORTED_VCF -t reference.sdf -o eval_output"
