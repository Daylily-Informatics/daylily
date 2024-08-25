version 1.0

task generateQcMetrics {
  input {
    File bam
    File refFlat
    File? ribosomal_intervals
    String strand = "unstranded"  # [first, second, unstranded]
  }

  Int space_needed_gb = 10 + round(size([bam, refFlat, ribosomal_intervals], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "18GB"
    docker: "mgibio/rnaseq:1.0.0"
    cpu: 1
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  # the mismatch between first and second here is intentional (and the nomenclature clash is stupid)
  # see https://github.com/griffithlab/rnaseq_tutorial/blob/master/manuscript/supplementary_tables/supplementary_table_5.md
  Map[String, String] strandness = {
    "first": "SECOND_READ_TRANSCRIPTION_STRAND",
    "second": "FIRST_READ_TRANSCRIPTION_STRAND",
    "unstranded": "NONE"
  }
  command <<<
    /usr/bin/java -Xmx16g -jar /opt/picard/picard.jar CollectRnaSeqMetrics \
    O=rna_metrics.txt CHART=rna_coverage_by_transcript_position.pdf REF_FLAT=~{refFlat} \
    ~{if (defined(ribosomal_intervals)) then "RIBOSOMAL_INTERVALS=~{select_first([ribosomal_intervals])}" else ""} \
    STRAND=~{strandness[strand]} I=~{bam}
  >>>

  output {
    File metrics= "rna_metrics.txt"
    File? chart = "rna_coverage_by_transcript_position.pdf"
  }
}
