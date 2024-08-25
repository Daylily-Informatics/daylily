version 1.0

task stringtie {
  input {
    String strand = "unstranded"  # [first, second, unstranded]
    File reference_annotation
    String sample_name
    File bam
  }

  Int cores = 12
  Int space_needed_gb = 10 + round(size([bam, reference_annotation], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "16GB"
    cpu: cores
    docker: "quay.io/biocontainers/stringtie:2.1.4--h7e0af3c_0"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String transcripts = "stringtie_transcripts.gtf"
  String expression = "stringtie_gene_expression.tsv"
  Map[String, String] strandness = {
    "first": "--rf", "second": "--fr", "unstranded": ""
  }
  command <<<
    /usr/local/bin/stringtie -o ~{transcripts} -A ~{expression} \
    -p ~{cores} -e ~{strandness[strand]} -G ~{reference_annotation} \
    -l ~{sample_name} ~{bam}
  >>>

  output {
    File transcript_gtf = transcripts
    File gene_expression_tsv = expression
  }
}

workflow wf {
  input {
    String? strand
    File reference_annotation
    String sample_name
    File bam
  }

  call stringtie {
    input:
    strand=strand,
    reference_annotation=reference_annotation,
    sample_name=sample_name,
    bam=bam
  }
}
