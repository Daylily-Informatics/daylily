version 1.0

task collectHsMetrics {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    String metric_accumulation_level

    File bait_intervals
    File target_intervals
    Boolean per_target_coverage = false
    Boolean per_base_coverage = false
    Int? minimum_base_quality
    Int? minimum_mapping_quality

    String output_prefix = "out"
  }

  Int space_needed_gb = 10 + round(size([bam, bam_bai, reference, reference_fai, reference_dict, bait_intervals, target_intervals], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "60GB"
    docker: "broadinstitute/picard:2.23.6"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String bamroot = basename(bam, ".bam")
  String hs_txt = "~{bamroot}.~{output_prefix}-HsMetrics.txt"
  String per_target_txt = "~{bamroot}.~{output_prefix}-PerTargetCoverage.txt"
  String per_base_txt = "~{bamroot}.~{output_prefix}-PerBaseCoverage.txt"
  command <<<
    /usr/bin/java -Xmx48g -jar /usr/picard/picard.jar CollectHsMetrics \
    O=~{hs_txt} \
    I=~{bam} \
    R=~{reference} \
    TARGET_INTERVALS=~{target_intervals} \
    METRIC_ACCUMULATION_LEVEL=~{metric_accumulation_level} \
    BAIT_INTERVALS=~{bait_intervals} \
    ~{if per_target_coverage then "PER_TARGET_COVERAGE=~{per_target_txt}" else ""} \
    ~{if per_base_coverage then "PER_BASE_COVERAGE=~{per_base_txt}" else ""} \
    ~{if defined(minimum_mapping_quality) then "MINIMUM_MAPPING_QUALITY=~{minimum_mapping_quality}" else ""} \
    ~{if defined(minimum_base_quality) then "MINIMUM_BASE_QUALITY=~{minimum_base_quality}" else ""}
  >>>

  output {
    File hs_metrics = hs_txt
    File? per_target_coverage_metrics = per_target_txt
    File? per_base_coverage_metrics = per_base_txt
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

    File bait_intervals
    File target_intervals
    Boolean per_target_coverage = false
    Boolean per_base_coverage = false
    Int? minimum_base_quality
    Int? minimum_mapping_quality

    String output_prefix = "out"
  }

  call collectHsMetrics {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    metric_accumulation_level=metric_accumulation_level,
    bait_intervals=bait_intervals,
    target_intervals=target_intervals,
    per_target_coverage=per_target_coverage,
    per_base_coverage=per_base_coverage,
    minimum_base_quality=minimum_base_quality,
    minimum_mapping_quality=minimum_mapping_quality,
    output_prefix=output_prefix
  }
}
