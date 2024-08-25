version 1.0

task biscuitMarkdup {
  input {
    File bam
    File bam_bai
  }

  Int cores = 4
  Int space_needed_gb = 10 + round(2*size(bam, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    cpu: cores
    memory: "24GB"
    docker: "mgibio/biscuit:0.3.8"
    disks: "local-disk ~{space_needed_gb} HDD"
    bootDiskSizeGb: space_needed_gb
  }

  command <<<
    set -eou pipefail
    /usr/bin/biscuit markdup "~{bam}" /dev/stdout | /usr/bin/sambamba sort -t ~{cores} -m 15G -o "markdup.bam" /dev/stdin
  >>>

  output {
    File markdup_bam = "markdup.bam"
  }
}
