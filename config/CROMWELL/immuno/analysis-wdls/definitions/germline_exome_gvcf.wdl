version 1.0

import "types.wdl"
import "alignment_exome.wdl" as ae
import "tools/index_cram.wdl" as ic
import "tools/bam_to_cram.wdl" as btc
import "tools/freemix.wdl" as f
import "subworkflows/gatk_haplotypecaller_iterator.wdl" as ghi

workflow germlineExomeGvcf {
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
    Array[File] bqsr_known_sites
    Array[File] bqsr_known_sites_tbi
    File bait_intervals
    File target_intervals
    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals
    File omni_vcf
    File omni_vcf_tbi
    String picard_metric_accumulation_level
    String emit_reference_confidence  # enum ['NONE', 'BP_RESOLUTION', 'GVCF']
    Array[String] gvcf_gq_bands
    Array[Array[String]] intervals
    Int? ploidy
    Int? qc_minimum_mapping_quality
    Int? qc_minimum_base_quality
  }

  call ae.alignmentExome as alignmentAndQc {
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
    sequence=sequence,
    trimming=trimming,
    bqsr_known_sites=bqsr_known_sites,
    bqsr_known_sites_tbi=bqsr_known_sites_tbi,
    bait_intervals=bait_intervals,
    target_intervals=target_intervals,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals,
    omni_vcf=omni_vcf,
    omni_vcf_tbi=omni_vcf_tbi,
    picard_metric_accumulation_level=picard_metric_accumulation_level,
    qc_minimum_mapping_quality=qc_minimum_mapping_quality,
    qc_minimum_base_quality=qc_minimum_base_quality
  }

  call f.freemix {
    input:
    verify_bam_id_metrics=select_first([alignmentAndQc.qc_metrics.verify_bam_id_metrics])
  }

  call ghi.gatkHaplotypecallerIterator as generateGvcfs {
    input:
    bam=alignmentAndQc.bam,
    bai=alignmentAndQc.bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    emit_reference_confidence=emit_reference_confidence,
    gvcf_gq_bands=gvcf_gq_bands,
    intervals=intervals,
    contamination_fraction=freemix.out,
    ploidy=ploidy
  }

  call btc.bamToCram {
    input:
    bam=alignmentAndQc.bam,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict
  }

  call ic.indexCram {
    input: cram=bamToCram.cram
  }

  output {
    File cram = indexCram.indexed_cram
    File cram_crai = indexCram.indexed_cram_crai
    File mark_duplicates_metrics = alignmentAndQc.mark_duplicates_metrics
    QCMetrics qc_metrics = alignmentAndQc.qc_metrics
    Array[File] gvcf = generateGvcfs.gvcf
  }
}
