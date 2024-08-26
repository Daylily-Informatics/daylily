version 1.0

import "subworkflows/qc_wgs.wdl" as qw
import "subworkflows/sequence_to_bqsr.wdl" as stb
import "types.wdl"  # !UnusedImport

workflow alignmentWgs {
  input {
    File reference
    File reference_fai
    File reference_dict
    File reference_alt
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_0123

    Array[SequenceData] sequence
    TrimmingOptions? trimming

    File omni_vcf
    File omni_vcf_tbi

    File intervals
    String picard_metric_accumulation_level
    Array[File] bqsr_known_sites
    Array[File] bqsr_known_sites_tbi

    Int? minimum_mapping_quality
    Int? minimum_base_quality
    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals
    String? sample_name
  }

  call stb.sequenceToBqsr as alignment {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    reference_alt=reference_alt,
    reference_amb=reference_amb,
    reference_ann=reference_ann,
    reference_bwt=reference_bwt,
    reference_pac=reference_pac,
    reference_0123=reference_0123,
    unaligned=sequence,
    trimming=trimming,
    bqsr_known_sites=bqsr_known_sites,
    bqsr_known_sites_tbi=bqsr_known_sites_tbi,
    final_name=sample_name
  }

  call qw.qcWgs as qc {
    input:
    sample_name=sample_name,
    bam=alignment.final_bam,
    bam_bai=alignment.final_bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    omni_vcf=omni_vcf,
    omni_vcf_tbi=omni_vcf_tbi,
    intervals=intervals,
    picard_metric_accumulation_level=picard_metric_accumulation_level,
    minimum_mapping_quality=minimum_mapping_quality,
    minimum_base_quality=minimum_base_quality,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals
  }

  output {
    File bam = alignment.final_bam
    File bam_bai = alignment.final_bam_bai
    File bai = alignment.final_bai
    File mark_duplicates_metrics = alignment.mark_duplicates_metrics_file
    QCMetrics qc_metrics = qc.metrics
  }
}
