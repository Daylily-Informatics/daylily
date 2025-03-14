
# Define input files                                                                                                                                                                                                                                                                                                    
ALIGNERS=("bwa2a")
SAMPLE="R0_H37013-Tumor_Tumor_0-S1"
CALLERS=("clair3" "deep" "lfq2" "oct")
CALLER_NAMES=("clair3" "deepvariant" "lofreq2" "octopus")

# Annotate each caller's VCF                                                                                                                                                                                                                                                                                            

for ai in "${!ALIGNERS[@]}"; do

  for i in "${!CALLERS[@]}"; do
    CALLER=${CALLERS[$i]}
    ALIGNER=${ALIGNERS[$ai]}
    CALLER_FILE="align/${ALIGNER}/snv/${CALLER}/${SAMPLE}.${ALIGNER}.${CALLERS[$i]}.snv.sort.vcf.gz"
    CALLER_LABEL="$ALIGNER"-"${CALLER_NAMES[$i]}"
    OUTPUT_FILE="${CALLER_LABEL}.annotated.vcf.gz"

    echo "Annotating $CALLER_FILE with caller: $CALLER_LABEL"

    #    -x INFO/AF,FORMAT/AF \                                                                                                                                                                                                                                                                                         

    echo "A" $CALLER .. $CALLER_FILE
    bcftools annotate \
      -h <(echo '##INFO=<ID=aligner_caller,Number=.,Type=String,Description="Source aligner-caller">') \
      --set-id +'%CHROM:%POS:%REF:%ALT' \
      --set "INFO/aligner_caller=$ALIGNER\-$CALLER" \
      -x INFO/AF,FORMAT/AF \
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

#### If combining across different aligners (strobe aligner specifically), the AF fields get wacky and need to be completely stripped                                                                                                                                                                                   
# bcftools annotate -x INFO/AF RIH0_ANA0-HG002-19_DBC0_0.merged_callers.sorted.vcf.gz -Oz -o no_AF.vcf.gz                                                                                                                                                                                                               

bcftools norm -Oz -m+any \
  $SORTED_VCF \
  -o all_callers.merged.vcf.gz

tabix all_callers.merged.vcf.gz
bcftools index all_callers.merged.vcf.gz


FIN=$SAMPLE.merged_allcallers.sort.vcf.gz
bcftools sort -Oz -o $FIN all_callers.merged.vcf.gz
bcftools index $FIN
tabix $FIN

bcftools stats "$FIN" > "${FIN}.stats.txt"
