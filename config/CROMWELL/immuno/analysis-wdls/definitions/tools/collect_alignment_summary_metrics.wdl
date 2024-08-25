version 1.0

task collectAlignmentSummaryMetrics {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    String metric_accumulation_level
  }

  Int space_needed_gb = 10 + round(size([bam, bam_bai, reference, reference_fai, reference_dict],"GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "48GB"
    docker: "broadinstitute/picard:2.23.6"
    disks: "local-disk ~{space_needed_gb} HDD"
  }
  String bamroot = basename(bam, ".bam")
  String summary_metrics = "~{bamroot}.AlignmentSummaryMetrics.txt"
  command <<<
    /usr/bin/java -Xmx32g -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics INPUT=~{bam} OUTPUT=~{summary_metrics} REFERENCE_SEQUENCE=~{reference} METRIC_ACCUMULATION_LEVEL=~{metric_accumulation_level}
  >>>

  output {
    File alignment_summary_metrics = summary_metrics
  }
}

workflow wf {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    String metric_accumulation_level
  }

  call collectAlignmentSummaryMetrics {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    metric_accumulation_level=metric_accumulation_level
  }
}
