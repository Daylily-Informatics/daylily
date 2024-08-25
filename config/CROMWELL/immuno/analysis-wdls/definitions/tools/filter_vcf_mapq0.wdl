version 1.0

task filterVcfMapq0 {
  input {
    File vcf
    File tumor_bam
    File tumor_bam_bai
    String sample_name
    Float threshold
  }

  Float bam_size = size([tumor_bam, tumor_bam_bai], "GB")
  Int space_needed_gb = 10 + round(bam_size + 2*size(vcf, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/mapq0-filter:v0.5.4"
    memory: "8GB"
    bootDiskSizeGb: 10
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = "mapq_filtered.vcf.gz"
  command <<<
    /bin/bash /usr/bin/mapq0_vcf_filter.sh `pwd` ~{vcf} ~{tumor_bam} ~{threshold} ~{sample_name}
  >>>

  output {
    File mapq0_filtered_vcf = outfile
  }
}

workflow wf {
  input {
    File vcf
    File tumor_bam
    File tumor_bam_bai
    String sample_name
    Float threshold
  }
  call filterVcfMapq0 {
    input:
    vcf=vcf,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    sample_name=sample_name,
    threshold=threshold
  }
}
