version 1.0

task filterVcfCustomAlleleFreq {
  input {
    File vcf
    Float maximum_population_allele_frequency
    String field_name
  }

  Int space_needed_gb = 10 + round(size(vcf, "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/vep_helper-cwl:vep_105.0_v1"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = "annotated.af_filtered.vcf"
  command <<<
    /usr/bin/perl /usr/bin/vcf_check.pl ~{vcf} ~{outfile} \
    /usr/bin/perl /opt/vep/src/ensembl-vep/filter_vep --format vcf -o ~{outfile} -i ~{vcf} \
    --filter "~{field_name} < ~{maximum_population_allele_frequency} or not ~{field_name}"
  >>>

  output {
    File filtered_vcf = outfile
  }
}

workflow wf {
  input {
    File vcf
    Float maximum_population_allele_frequency
    String field_name
  }

  call filterVcfCustomAlleleFreq {
    input:
    vcf=vcf,
    maximum_population_allele_frequency=maximum_population_allele_frequency,
    field_name=field_name
  }
}
