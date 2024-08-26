version 1.0

task biscuitPileup {
  input {
    File bam
    File reference
    File reference_fai
  }
  Int cores = 4

  Int space_needed_gb = 10 + round(2*size([bam, reference], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "48GB"
    cpu: cores
    docker: "mgibio/biscuit:0.3.8"
    disks: "local-disk ~{space_needed_gb} HDD"
    bootDiskSizeGb: space_needed_gb
  }

  command <<<
    set -eo pipefail
    /usr/bin/biscuit pileup -q ~{cores} -w pileup_stats.txt  "~{reference}" "~{bam}" | /opt/htslib/bin/bgzip > out.vcf
  >>>

  output {
    File vcf = "output.vcf"
  }
}
