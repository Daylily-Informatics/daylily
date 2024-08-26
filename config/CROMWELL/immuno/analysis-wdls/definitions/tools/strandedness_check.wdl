version 1.0

task strandednessCheck {
  input {
    File reference_annotation
    File kallisto_index
    File cdna_fasta
    File reads1
    File reads2
  }

  Int space_needed_gb = 10 + round(2*size([reference_annotation, kallisto_index, cdna_fasta, reads1, reads2], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "16GB"
    bootDiskSizeGb: space_needed_gb  # default
    cpu: 1
    docker: "mgibio/checkstrandedness:v1"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = basename(reads1, ".fastq") + "strandness_check.txt"
  command <<<
    check_strandedness --print_commands \
        --gtf ~{reference_annotation} --kallisto_index ~{kallisto_index} --transcripts ~{cdna_fasta} \
        --reads_1 ~{reads1} --reads_2 ~{reads2} -n 100000 > ~{outfile}
  >>>

  output {
    File strandedness_check = outfile
  }
}

workflow wf {
  call strandednessCheck { input: }
}
