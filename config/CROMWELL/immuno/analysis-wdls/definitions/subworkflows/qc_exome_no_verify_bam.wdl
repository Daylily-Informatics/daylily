version 1.0

import "../types.wdl"  # !UnusedImport

import "../subworkflows/hs_metrics.wdl" as hm
import "../tools/collect_alignment_summary_metrics.wdl" as casm
import "../tools/collect_hs_metrics.wdl" as chm
import "../tools/collect_insert_size_metrics.wdl" as cism
import "../tools/samtools_flagstat.wdl" as sf


workflow qcExomeNoVerifyBam {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bai
    String picard_metric_accumulation_level
    Int? minimum_mapping_quality
    Int? minimum_base_quality
    File bait_intervals
    File target_intervals
    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals
  }

  call cism.collectInsertSizeMetrics {
    input:
    bam=bam,
    bam_bai=bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    metric_accumulation_level=picard_metric_accumulation_level
  }

  call casm.collectAlignmentSummaryMetrics {
    input:
    bam=bam,
    bam_bai=bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    metric_accumulation_level=picard_metric_accumulation_level
  }

  call chm.collectHsMetrics as collectRoiHsMetrics{
    input:
    bam=bam,
    bam_bai=bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    metric_accumulation_level="ALL_READS",
    bait_intervals=bait_intervals,
    target_intervals=target_intervals,
    per_target_coverage=false,
    per_base_coverage=false,
    output_prefix="roi",
    minimum_mapping_quality=minimum_mapping_quality,
    minimum_base_quality=minimum_base_quality
  }


  call hm.hsMetrics as collectDetailedHsMetrics {
    input:
    bam=bam,
    bam_bai=bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    minimum_mapping_quality=minimum_mapping_quality,
    minimum_base_quality=minimum_base_quality,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals
  }


  call sf.samtoolsFlagstat {
    input:
    bam=bam,
    bam_bai=bai
  }

  output {
    QCMetrics metrics = object {
      insert_size_metrics: collectInsertSizeMetrics.insert_size_metrics,
      insert_size_histogram: collectInsertSizeMetrics.insert_size_histogram,
      alignment_summary_metrics: collectAlignmentSummaryMetrics.alignment_summary_metrics,
      hs_metrics: collectRoiHsMetrics.hs_metrics,
      per_target_coverage_metrics: collectDetailedHsMetrics.per_target_coverage_metrics,
      per_target_hs_metrics: collectDetailedHsMetrics.per_target_hs_metrics,
      per_base_coverage_metrics: collectDetailedHsMetrics.per_base_coverage_metrics,
      per_base_hs_metrics: collectDetailedHsMetrics.per_base_hs_metrics,
      summary_hs_metrics: collectDetailedHsMetrics.summary_hs_metrics,
      flagstats: samtoolsFlagstat.flagstats
    }
  }
}
