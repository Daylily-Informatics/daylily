version 1.0

task biscuitAlign {
  input {
    File reference_index
    File fastq1
    File fastq2
    String read_group_id
  }

  Int cores = 12
  Int space_needed_gb = 10 + round(2*size([reference_index, fastq1, fastq2], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "32GB"
    cpu: cores
    docker: "mgibio/biscuit:0.3.8"
    disks: "local-disk ~{space_needed_gb} HDD"
    bootDiskSizeGb: space_needed_gb
  }

  command <<<
    set -eou pipefail
    /usr/bin/biscuit align -t ~{cores} -M -R "~{read_group_id}" "~{reference_index}" "~{fastq1}" "~{fastq2}" \
    | /usr/bin/sambamba view -S -f bam -l 0 /dev/stdin \
    | /usr/bin/sambamba sort -t ~{cores} -m 8G -o "aligned.bam" /dev/stdin
  >>>

  output {
    File aligned_bam = "aligned.bam"
  }
}
