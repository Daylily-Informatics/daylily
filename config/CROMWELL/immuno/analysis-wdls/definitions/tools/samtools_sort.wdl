version 1.0

task samtoolsSort {
  input {
    String output_filename = "sorted.bam"
    File input_bam
  }

  Int cores = 1
  Int space_needed_gb = 10 + round(3*size(input_bam, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "4GB"
    cpu: cores
    docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    /usr/local/bin/samtools sort -o ~{output_filename} -@ ~{cores} ~{input_bam}
  >>>

  output {
    File sorted_bam = output_filename
  }
}

workflow wf {
  input {
    String? output_filename
    File input_bam
  }

  call samtoolsSort {
    input:
    input_bam=input_bam,
    output_filename=output_filename
  }
}
