version 1.0

task annotateKnownVariants {
  input {
    File vcf
    File vcf_tbi
    File? validated_variants
    File? validated_variants_tbi
  }

  Int space_needed_gb = 10 + round(size([vcf, vcf_tbi, validated_variants, validated_variants_tbi], "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/bcftools-cwl:1.12"
    memory: "8GB"
    disks:  "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = "validated_annotated_pipeline_variants.vcf.gz"
  command <<<
    set -eou pipefail

    PIPELINE_VCF=~{vcf}
    VALIDATED_VCF=~{if defined(validated_variants) then validated_variants else ""}

    if [[ -n $VALIDATED_VCF ]]; then
        /opt/bcftools/bin/bcftools view -f PASS -Oz -o pass_filtered_validated_variants.vcf.gz $VALIDATED_VCF
        /opt/bcftools/bin/bcftools index -t pass_filtered_validated_variants.vcf.gz
        /opt/bcftools/bin/bcftools annotate -Oz -o ~{outfile} -a pass_filtered_validated_variants.vcf.gz -m 'VALIDATED' $PIPELINE_VCF
        /opt/bcftools/bin/bcftools index -t ~{outfile}
    else
        cp $PIPELINE_VCF ~{outfile}
        cp $PIPELINE_VCF.tbi ~{outfile}.tbi
    fi
  >>>

  output {
    File validated_annotated_vcf = outfile
    File validated_annotated_vcf_tbi = outfile + ".tbi"
  }
}
