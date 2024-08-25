version 1.0

import "../subworkflows/hs_metrics.wdl" as hm
import "../tools/collect_alignment_summary_metrics.wdl" as casm
import "../tools/collect_gc_bias_metrics.wdl" as cgbm
import "../tools/collect_insert_size_metrics.wdl" as cism
import "../tools/collect_wgs_metrics.wdl" as cwm
import "../tools/samtools_flagstat.wdl" as sf
import "../tools/verify_bam_id.wdl" as vbi
import "../types.wdl"  # !UnusedImport

workflow qcWgs {
  input {
    String sample_name = "final"
    File bam
    File bam_bai

    File reference
    File reference_fai
    File reference_dict

    File intervals

    File omni_vcf
    File omni_vcf_tbi

    String picard_metric_accumulation_level = "ALL_READS"
    Int? minimum_mapping_quality
    Int? minimum_base_quality

    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals
  }

  call cism.collectInsertSizeMetrics {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    metric_accumulation_level=picard_metric_accumulation_level
  }

  call casm.collectAlignmentSummaryMetrics {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    metric_accumulation_level=picard_metric_accumulation_level
  }

  call cgbm.collectGcBiasMetrics {
    input:
    sample_name=sample_name,
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    metric_accumulation_level=picard_metric_accumulation_level
  }

  call cwm.collectWgsMetrics {
    input:
    sample_name=sample_name,
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    intervals=intervals
  }

  call sf.samtoolsFlagstat {
    input:
    bam=bam,
    bam_bai=bam_bai
  }

  call vbi.verifyBamId {
    input:
    bam=bam,
    bam_bai=bam_bai,
    vcf=omni_vcf
  }

  call hm.hsMetrics as collectHsMetrics {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals,
    minimum_mapping_quality=minimum_mapping_quality,
    minimum_base_quality=minimum_base_quality
  }

  output {
    QCMetrics metrics = object {
      insert_size_metrics: collectInsertSizeMetrics.insert_size_metrics,
      insert_size_histogram: collectInsertSizeMetrics.insert_size_histogram,
      alignment_summary_metrics: collectAlignmentSummaryMetrics.alignment_summary_metrics,
      gc_bias_metrics: collectGcBiasMetrics.gc_bias_metrics,
      gc_bias_metrics_chart: collectGcBiasMetrics.gc_bias_metrics_chart,
      gc_bias_metrics_summary: collectGcBiasMetrics.gc_bias_metrics_summary,
      wgs_metrics: collectWgsMetrics.wgs_metrics,
      flagstats: samtoolsFlagstat.flagstats,
      verify_bam_id_metrics: verifyBamId.verify_bam_id_metrics,
      verify_bam_id_depth: verifyBamId.verify_bam_id_depth,
      per_base_coverage_metrics: collectHsMetrics.per_base_coverage_metrics,
      per_base_hs_metrics: collectHsMetrics.per_base_hs_metrics,
      per_target_coverage_metrics: collectHsMetrics.per_target_coverage_metrics,
      per_target_hs_metrics: collectHsMetrics.per_target_hs_metrics,
      summary_hs_metrics: collectHsMetrics.summary_hs_metrics
    }
  }
}
