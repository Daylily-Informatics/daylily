version 1.0

task addStringAtLineBgzipped {
  input {
    File input_file
    Int line_number
    String some_text
    String output_name = basename(input_file, ".gz") + ".commented.gz"
  }

  Int space_needed_gb = 10 + round(2*size(input_file, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/samtools:1.11-h6270b1f_0"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} HDD"
    bootDiskSizeGb: space_needed_gb
  }

  command <<<
    zcat ~{input_file} | awk -v n=~{line_number} -v s="~{some_text}" 'NR == n {print s} {print}' | /opt/htslib/bin/bgzip
  >>>

  output {
    File output_file = output_name
  }
}
