version 1.0

import "types.wdl"

import "alignment_wgs.wdl" as aw
import "tumor_only_detect_variants.wdl" as todv

import "tools/bam_to_cram.wdl" as btc
import "tools/index_cram.wdl" as ic

workflow tumorOnlyWgs {
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

    File target_intervals
    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals

    File omni_vcf
    File omni_vcf_tbi

    String picard_metric_accumulation_level

    File roi_intervals

    File vep_cache_dir_zip
    String vep_ensembl_assembly
    String vep_ensembl_version
    String vep_ensembl_species
    String? vep_pick  # enum ["pick", "flag_pick", "pick_allele", "per_gene", "pick_allele_gene", "flag_pick_allele", "flag_pick_allele_gene"]
    Array[VepCustomAnnotation] vep_custom_annotations

    File? synonyms_file
    String sample_name

    File docm_vcf
    File docm_vcf_tbi

    Int? readcount_minimum_mapping_quality
    Int? readcount_minimum_base_quality
    Int qc_minimum_mapping_quality = 0
    Int qc_minimum_base_quality = 0

    Array[String] variants_to_table_fields = ["CHROM", "POS", "ID", "REF", "ALT", "set", "AC", "AF"]
    Array[String] variants_to_table_genotype_fields = ["GT", "AD"]
    Array[String] vep_to_table_fields = ["HGVSc", "HGVSp"]
    Int? varscan_min_coverage
    Float? varscan_min_var_freq
    Int? varscan_min_reads
    Float maximum_population_allele_frequency = 0.001
}

  call aw.alignmentWgs as alignmentAndQc {
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
    omni_vcf=omni_vcf,
    omni_vcf_tbi=omni_vcf_tbi,
    intervals=target_intervals,
    picard_metric_accumulation_level=picard_metric_accumulation_level,
    bqsr_known_sites=bqsr_known_sites,
    bqsr_known_sites_tbi=bqsr_known_sites_tbi,
    minimum_mapping_quality=qc_minimum_mapping_quality,
    minimum_base_quality=qc_minimum_base_quality,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals,
    sample_name=sample_name,
  }

  call todv.tumorOnlyDetectVariants as detectVariants {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=alignmentAndQc.bam,
    bam_bai=alignmentAndQc.bam_bai,
    roi_intervals=roi_intervals,
    varscan_min_coverage=varscan_min_coverage,
    varscan_min_var_freq=varscan_min_var_freq,
    varscan_min_reads=varscan_min_reads,
    maximum_population_allele_frequency=maximum_population_allele_frequency,
    vep_cache_dir_zip=vep_cache_dir_zip,
    synonyms_file=synonyms_file,
    vep_pick=vep_pick,
    vep_ensembl_assembly=vep_ensembl_assembly,
    vep_ensembl_version=vep_ensembl_version,
    vep_ensembl_species=vep_ensembl_species,
    variants_to_table_fields=variants_to_table_fields,
    variants_to_table_genotype_fields=variants_to_table_genotype_fields,
    vep_to_table_fields=vep_to_table_fields,
    sample_name=sample_name,
    docm_vcf=docm_vcf,
    docm_vcf_tbi=docm_vcf_tbi,
    vep_custom_annotations=vep_custom_annotations,
    readcount_minimum_mapping_quality=readcount_minimum_mapping_quality,
    readcount_minimum_base_quality=readcount_minimum_base_quality
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
    File mark_duplicates_metrics = alignmentAndQc.mark_duplicates_metrics
    QCMetrics qc_metrics = alignmentAndQc.qc_metrics
    File varscan_vcf = detectVariants.varscan_vcf
    File varscan_vcf_tbi = detectVariants.varscan_vcf_tbi
    File docm_gatk_vcf = detectVariants.docm_gatk_vcf
    File annotated_vcf = detectVariants.annotated_vcf
    File annotated_vcf_tbi = detectVariants.annotated_vcf_tbi
    File final_vcf = detectVariants.final_vcf
    File final_vcf_tbi = detectVariants.final_vcf_tbi
    File final_tsv = detectVariants.final_tsv
    File vep_summary = detectVariants.vep_summary
    File tumor_snv_bam_readcount_tsv = detectVariants.tumor_snv_bam_readcount_tsv
    File tumor_indel_bam_readcount_tsv = detectVariants.tumor_indel_bam_readcount_tsv
  }
}
