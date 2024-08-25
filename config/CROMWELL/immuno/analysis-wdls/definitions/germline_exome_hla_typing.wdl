version 1.0

import "germline_exome.wdl" as ge
import "types.wdl"  # !UnusedImport
import "tools/optitype_dna.wdl" as od
import "tools/phlat.wdl" as ph

workflow germlineExomeHlaTyping {
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
    Array[File] bqsr_known_sites
    Array[File] bqsr_known_sites_tbi
    File bait_intervals
    File target_intervals
    Int target_interval_padding = 100
    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals
    File omni_vcf
    File omni_vcf_tbi
    String picard_metric_accumulation_level
    Array[String] gvcf_gq_bands
    Array[Array[String]] intervals
    Int? ploidy
    File vep_cache_dir_zip
    String vep_ensembl_assembly
    String vep_ensembl_version
    String vep_ensembl_species
    File? synonyms_file
    Boolean? annotate_coding_only
    Int? qc_minimum_mapping_quality
    Int? qc_minimum_base_quality
    Array[VepCustomAnnotation] vep_custom_annotations
    String? optitype_name
    Float germline_filter_gnomAD_maximum_population_allele_frequency
  }

  call ge.germlineExome {
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
    bqsr_known_sites=bqsr_known_sites,
    bqsr_known_sites_tbi=bqsr_known_sites_tbi,
    bait_intervals=bait_intervals,
    target_intervals=target_intervals,
    target_interval_padding=target_interval_padding,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals,
    omni_vcf=omni_vcf,
    omni_vcf_tbi=omni_vcf_tbi,
    picard_metric_accumulation_level=picard_metric_accumulation_level,
    gvcf_gq_bands=gvcf_gq_bands,
    intervals=intervals,
    ploidy=ploidy,
    vep_cache_dir_zip=vep_cache_dir_zip,
    vep_ensembl_assembly=vep_ensembl_assembly,
    vep_ensembl_version=vep_ensembl_version,
    vep_ensembl_species=vep_ensembl_species,
    synonyms_file=synonyms_file,
    annotate_coding_only=annotate_coding_only,
    vep_custom_annotations=vep_custom_annotations,
    qc_minimum_mapping_quality=qc_minimum_mapping_quality,
    qc_minimum_base_quality=qc_minimum_base_quality,
    germline_filter_gnomAD_maximum_population_allele_frequency=germline_filter_gnomAD_maximum_population_allele_frequency
  }

  call od.optitypeDna as optitype {
    input:
    reference=reference,
    reference_fai=reference_fai,
    cram=germlineExome.cram,
    cram_crai=germlineExome.cram_crai,
    optitype_name=optitype_name
  }

  call ph.phlat {
    input:
    phlat_name="phlat_normal",
    reference=reference,
    reference_fai=reference_fai,
    cram=germlineExome.cram,
    cram_crai=germlineExome.cram_crai
  }

  output {
    File cram = germlineExome.cram
    File mark_duplicates_metrics = germlineExome.mark_duplicates_metrics
    QCMetrics qc_metrics = germlineExome.qc_metrics
    File raw_vcf = germlineExome.raw_vcf
    File raw_vcf_tbi = germlineExome.raw_vcf_tbi
    File final_vcf = germlineExome.final_vcf
    File final_vcf_tbi = germlineExome.final_vcf_tbi
    File filtered_vcf = germlineExome.filtered_vcf
    File filtered_vcf_tbi = germlineExome.filtered_vcf_tbi
    File vep_summary = germlineExome.vep_summary
    File optitype_tsv = optitype.optitype_tsv
    File optitype_plot = optitype.optitype_plot
    File phlat_summary = phlat.phlat_summary 
  }
}
