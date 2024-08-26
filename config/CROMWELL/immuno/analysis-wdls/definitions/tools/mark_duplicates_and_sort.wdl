version 1.0

task markDuplicatesAndSort {
  input {
    File bam
    String output_name = "MarkedSorted.bam"
  }
  String metrics_file_name = sub(output_name, "\.bam$", ".mark_dups_metrics.txt")
  Int space_needed_gb = 20 + round(5*size(bam, "GB"))
  #estimate 15M reads per Gb size of bam
  #markdup is listed as 2Gb per 100M reads
  Int mem_needed_gb = round(((size(bam, "GB")*15)/100)*2)+32
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/sambamba:0.8.2--h98b6b92_2"
    memory: "~{mem_needed_gb}GB"
    cpu: 16
    # add space to shift bam around via stdin/stdout and a bit more
    bootDiskSizeGb: space_needed_gb
    # add space for input bam, output bam, and a bit more
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -o pipefail
    set -o errexit
    sambamba markdup -t 8 ~{bam} /dev/stdout 2>~{metrics_file_name} \
        | sambamba sort -t 8 -m 16G -o ~{output_name} /dev/stdin
  >>>

  output {
    File sorted_bam = output_name
    File sorted_bam_bai = output_name + ".bai"
    File metrics_file = metrics_file_name
  }
}

workflow wf {
  input {
    File bam
    String? output_name
  }
  call markDuplicatesAndSort {
    input:
    bam=bam,
    output_name=output_name
  }
}
