version 1.0

import "../types.wdl"  # !UnusedImport

import "../tools/collect_insert_size_metrics.wdl" as cism
import "../tools/collect_alignment_summary_metrics.wdl" as casm
import "../tools/collect_hs_metrics.wdl" as chm
import "../tools/samtools_flagstat.wdl" as sf
import "../tools/select_variants.wdl" as sv
import "../tools/verify_bam_id.wdl" as vbi
import "../subworkflows/hs_metrics.wdl" as hm

workflow qcExome {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    File bait_intervals
    File target_intervals
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

  call chm.collectHsMetrics as collectRoiHsMetrics {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    metric_accumulation_level="ALL_READS",
    bait_intervals=bait_intervals,
    target_intervals=target_intervals,
    per_target_coverage = false,
    per_base_coverage = false,
    output_prefix = "roi",
    minimum_mapping_quality=minimum_mapping_quality,
    minimum_base_quality=minimum_base_quality
  }

  call hm.hsMetrics as collectDetailedHsMetrics {
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

  call sf.samtoolsFlagstat {
    input:
    bam=bam,
    bam_bai=bam_bai
  }

  call sv.selectVariants {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    vcf=omni_vcf,
    vcf_tbi=omni_vcf_tbi,
    interval_list=target_intervals
  }

  call vbi.verifyBamId {
    input:
    bam=bam,
    bam_bai=bam_bai,
    vcf=selectVariants.filtered_vcf
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
      flagstats: samtoolsFlagstat.flagstats,
      verify_bam_id_metrics: verifyBamId.verify_bam_id_metrics,
      verify_bam_id_depth: verifyBamId.verify_bam_id_depth
    }
  }
}
