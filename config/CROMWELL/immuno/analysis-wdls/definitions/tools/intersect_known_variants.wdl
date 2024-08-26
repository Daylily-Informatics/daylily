version 1.0

task intersectKnownVariants {
  input {
    File vcf
    File vcf_tbi
    File? validated_variants
    File? validated_variants_tbi
  }

  Int space_needed_gb = 10 + round(2*size([vcf, vcf_tbi, validated_variants, validated_variants_tbi], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "8GB"
    docker: "mgibio/bcftools-cwl:1.12"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = "pass_filtered_validated_variants.vcf.gz"
  command <<<
    set -eou

    if [ ~{defined(validated_variants)} = true ]; then
        #filter the validated vcf to ensure there are only passing variants, then re-index
        /opt/bcftools/bin/bcftools view -f PASS -Oz -o "~{outfile}" "~{validated_variants}"
        /opt/bcftools/bin/bcftools index -t "~{outfile}"
        #intersect the two vcfs; output will contain only passing variants
        #-n specifies that the output should contain only variants found in both files
        #-w results in a single output vcf containing the intersection
        #-p specifies the directory that will contain output files (vcf, index, and summary files)
        #-Oz specifies the output format as compressed
        /opt/bcftools/bin/bcftools isec -f PASS -n=2 -w1 -p validated -Oz "~{vcf}" "~{outfile}"
    else
        mkdir validated
        cp "~{vcf}" validated/0000.vcf.gz
        cp "~{vcf_tbi}" validated/0000.vcf.gz.tbi
    fi
  >>>

  output {
    File validated_and_pipeline_vcf = "validated/0000.vcf.gz"
    File validated_and_pipeline_vcf_tbi = "validated/0000.vcf.gz.tbi"
  }
}

workflow wf {
  input {
    File vcf
    File vcf_tbi
    File? validated_variants
    File? validated_variants_tbi
  }

  call intersectKnownVariants {
    input:
    vcf=vcf,
    vcf_tbi=vcf_tbi,
    validated_variants=validated_variants,
    validated_variants_tbi=validated_variants_tbi
  }
}
