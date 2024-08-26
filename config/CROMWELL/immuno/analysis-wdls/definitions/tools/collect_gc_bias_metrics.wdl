version 1.0

task collectGcBiasMetrics {
  input {
    String sample_name = "final"
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    String metric_accumulation_level
  }


  Float bam_size_gb = size([bam, bam_bai], "GB")
  Float reference_size_gb = size([reference, reference_fai, reference_dict], "GB")
  Int space_needed_gb = 10 + round(bam_size_gb + reference_size_gb)
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "48GB"
    docker: "broadinstitute/picard:2.23.6"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String bias_metrics = sample_name + ".GcBiasMetrics.txt"
  String bias_metrics_summary = sample_name + ".GcBiasMetricsSummary.txt"
  String bias_metrics_chart = sample_name + ".GcBiasMetricsChart.pdf"
  command <<<
    /usr/bin/java -Xmx32g -jar /usr/picard/picard.jar CollectGcBiasMetrics \
    O=~{bias_metrics} \
    CHART=~{bias_metrics_chart} \
    S=~{bias_metrics_summary} \
    I=~{bam} R=~{reference} \
    METRIC_ACCUMULATION_LEVEL=~{metric_accumulation_level}
  >>>

  output {
    File gc_bias_metrics = bias_metrics
    File gc_bias_metrics_chart = bias_metrics_chart
    File gc_bias_metrics_summary = bias_metrics_summary
  }
}
