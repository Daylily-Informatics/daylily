version 1.0

task kallisto {
  input {
    File kallisto_index
    Array[Array[File]] fastqs
    String strand = "unstranded"  # [first, second, unstranded]
  }

  Int cores = 8
  Int space_needed_gb = 10 + round(size(flatten(fastqs), "GB") + size(kallisto_index, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "32GB"
    cpu: cores
    docker: "quay.io/biocontainers/kallisto:0.46.1--h4f7b962_0"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  Map[String, String] strandness = {
    "first": "--rf-stranded",
    "second": "--fr-stranded",
    "unstranded": ""
  }
  command <<<
    kallisto quant -t ~{cores} -b 100 --fusion -o kallisto \
    ~{strandness[strand]} -i ~{kallisto_index} ~{sep=" " flatten(fastqs)}

  >>>

  output {
    File expression_transcript_table = "kallisto/abundance.tsv"
    File expression_transcript_h5 = "kallisto/abundance.h5"
    File fusion_evidence = "kallisto/fusion.txt"
  }
}

workflow wf {
  input {
    File kallisto_index
    Array[Array[File]] fastqs
    String strand = "unstranded"  # [first, second, unstranded]
  }

  call kallisto {
    input:
    kallisto_index=kallisto_index,
    fastqs=fastqs,
    strand=strand
  }
}
