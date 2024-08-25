version 1.0

import "types.wdl"

import "subworkflows/sequence_to_bam_nonhuman.wdl" as stbn
import "subworkflows/qc_wgs_nonhuman.wdl" as qwn

workflow alignmentWgsNonhuman {
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
    String? final_name
    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals

    String picard_metric_accumulation_level
    Int? minimum_mapping_quality
    Int? minimum_base_quality
  }

  call stbn.sequenceToBamNonhuman as alignment {
    input:
    reference=reference,
    reference_alt=reference_alt,
    reference_amb=reference_amb,
    reference_ann=reference_ann,
    reference_bwt=reference_bwt,
    reference_pac=reference_pac,
    reference_0123=reference_0123,
    unaligned=sequence,
    trimming=trimming,
    final_name=final_name
  }

  call qwn.qcWgsNonhuman as qc {
    input:
    bam=alignment.final_bam,
    bam_bai=alignment.final_bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    picard_metric_accumulation_level=picard_metric_accumulation_level,
    minimum_mapping_quality=minimum_mapping_quality,
    minimum_base_quality=minimum_base_quality,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals
  }

  output {
    File bam = alignment.final_bam
    File mark_duplicates_metrics = alignment.mark_duplicates_metrics_file
    QCMetrics qc_metrics = qc.metrics
  }
}
