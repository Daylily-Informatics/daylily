version 1.0

task removeEndTags {
  input {
    File vcf
    File vcf_tbi
  }

  Int space_needed_gb = 10 + round(size(vcf, "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "4GB"
    docker: "mgibio/bcftools-cwl:1.12"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = "pindel.noend.vcf.gz"
  command <<<
    /opt/bcftools/bin/bcftools annotate -x INFO/END -Oz -o ~{outfile} ~{vcf}
  >>>

  output {
    File processed_vcf = outfile
  }
}

workflow wf {
  input {
    File vcf
    File vcf_tbi
  }

  call removeEndTags {
    input:
    vcf=vcf,
    vcf_tbi=vcf_tbi
  }
}
