version 1.0

task filterVcfCodingVariant {
  input {
    File vcf
  }

  Int space_needed_gb = 10 + round(2*size(vcf, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "4GB"
    docker: "mgibio/vep_helper-cwl:vep_105.0_v1"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = "annotated.coding_variant_filtered"
  command <<<
    /usr/bin/perl /usr/bin/vcf_check.pl ~{vcf} ~{outfile} \
    /usr/bin/perl /opt/vep/src/ensembl-vep/filter_vep \
    --format vcf \
    -o ~{outfile} \
    --ontology \
    --filter "Consequence is coding_sequence_variant" \
    -i ~{vcf}
  >>>

  output {
    File filtered_vcf = outfile
  }
}
