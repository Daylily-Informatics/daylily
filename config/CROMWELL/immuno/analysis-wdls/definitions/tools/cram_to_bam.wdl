version 1.0

task cramToBam {
  input {
    File cram
    File cram_index
    File reference
    File reference_index
    File reference_dict
  }

  Int space_needed_gb = 10 + round(size([cram, cram_index, reference, reference_index, reference_dict], "GB") * 3)
  runtime {
    memory: "4GB"
    docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = basename(cram, ".cram") + ".bam"

  command <<<
    /usr/local/bin/samtools view -b -o "~{outfile}" -T "~{reference}" "~{cram}"
  >>>

  output {
    File bam = outfile
  }
}
