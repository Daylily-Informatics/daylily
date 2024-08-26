version 1.0

import "../types.wdl"  # !UnusedImport

import "../tools/collect_hs_metrics.wdl" as hs

workflow hsMetrics {
  input {
    File bam
    File bam_bai

    Int? minimum_mapping_quality
    Int? minimum_base_quality

    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals

    File reference
    File reference_fai
    File reference_dict
  }

  scatter(interval in summary_intervals) {
    call hs.collectHsMetrics as collectSummaryHsMetrics{
      input:
      bam=bam,
      bam_bai=bam_bai,
      reference=reference,
      reference_fai=reference_fai,
      reference_dict=reference_dict,
      bait_intervals=interval.file,
      target_intervals=interval.file,
      output_prefix="summary-~{interval.label}",
      minimum_mapping_quality=minimum_mapping_quality,
      minimum_base_quality=minimum_base_quality,
      metric_accumulation_level="ALL_READS",
      per_target_coverage=false,
      per_base_coverage=false
    }
  }

  scatter(interval in per_base_intervals) {
    call hs.collectHsMetrics as collectPerBaseHsMetrics {
      input:
      bam=bam,
      bam_bai=bam_bai,
      reference=reference,
      reference_fai=reference_fai,
      reference_dict=reference_dict,
      bait_intervals=interval.file,
      target_intervals=interval.file,
      output_prefix="base-~{interval.label}",
      minimum_mapping_quality=minimum_mapping_quality,
      minimum_base_quality=minimum_base_quality,
      metric_accumulation_level="ALL_READS",
      per_target_coverage=false,
      per_base_coverage=true
    }
  }

  scatter(interval in per_target_intervals) {
    call hs.collectHsMetrics as collectPerTargetHsMetrics{
      input:
      bam=bam,
      bam_bai=bam_bai,
      reference=reference,
      reference_fai=reference_fai,
      reference_dict=reference_dict,
      bait_intervals=interval.file,
      target_intervals=interval.file,
      output_prefix="target-~{interval.label}",
      minimum_mapping_quality=minimum_mapping_quality,
      minimum_base_quality=minimum_base_quality,
      metric_accumulation_level="ALL_READS",
      per_target_coverage=true,
      per_base_coverage=false
    }
  }

  output {
    Array[File] per_base_coverage_metrics = select_all(collectPerBaseHsMetrics.per_base_coverage_metrics)
    Array[File] per_base_hs_metrics = collectPerBaseHsMetrics.hs_metrics
    Array[File] per_target_coverage_metrics = select_all(collectPerTargetHsMetrics.per_target_coverage_metrics)
    Array[File] per_target_hs_metrics = collectPerTargetHsMetrics.hs_metrics
    Array[File] summary_hs_metrics = collectSummaryHsMetrics.hs_metrics
  }
}
