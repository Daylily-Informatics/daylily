version 1.0

task nameSort {
  input {
    File bam
  }
  Int cores = 8
  String outfile = basename(bam, ".bam") + ".NameSorted.bam"

  Int input_size_gb = 10 + 5*round(size(bam, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/sambamba-cwl:0.6.4"
    memory: "26GB"
    cpu: cores
    disks: "local-disk ~{input_size_gb} HDD"
  }

  command <<<
    /usr/bin/sambamba sort "~{bam}" -t ~{cores} -m "22GB" -n -o "~{outfile}"
  >>>

  output {
    File name_sorted_bam = outfile
  }
}

workflow wf {
  input {
    File bam
  }
  call nameSort {
    input: bam=bam
  }
}
