version 1.0

import "types.wdl"

import "alignment_exome_nonhuman.wdl" as aen
import "detect_variants_nonhuman.wdl" as dvn

import "tools/bam_to_cram.wdl" as btc
import "tools/index_cram.wdl" as ic
import "tools/interval_list_expand.wdl" as ile

workflow somaticExomeNonhuman {
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

    String tumor_name = "tumor"
    Array[SequenceData] tumor_sequence

    String normal_name = "normal"
    Array[SequenceData] normal_sequence

    TrimmingOptions? trimming

    File bait_intervals
    File target_intervals
    Int target_interval_padding = 100
    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals

    String picard_metric_accumulation_level
    Int qc_minimum_mapping_quality = 0
    Int qc_minimum_base_quality = 0

    Int strelka_cpu_reserved = 8
    Int scatter_count = 50

    Int? varscan_strand_filter
    Int? varscan_min_coverage
    Float? varscan_min_var_freq
    Float? varscan_p_value
    Float? varscan_max_normal_freq

    Float? fp_min_var_freq

    File vep_cache_dir_zip
    String vep_ensembl_assembly
    String vep_ensembl_version
    String vep_ensembl_species
    File? synonyms_file
    Boolean annotate_coding_only = false
    # one of [pick, flag_pick, pick-allele, per_gene, pick_allele_gene, flag_pick_allele, flag_pick_allele_gene]
    String? vep_pick
    Boolean cle_vcf_filter = false

    Float? filter_somatic_llr_threshold
    Float? filter_somatic_llr_tumor_purity
    Float? filter_somatic_llr_normal_contamination_rate

    Array[String] vep_to_table_fields = ["Consequence", "SYMBOL", "Feature"]
    Array[String] variants_to_table_genotype_fields = ["GT", "AD"]
    Array[String] variants_to_table_fields = ["CHROM", "POS", "ID", "REF", "ALT", "set", "AC", "AF"]

    String tumor_sample_name
    String normal_sample_name
  }

  call aen.alignmentExomeNonhuman as tumorAlignmentAndQc {
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
    sequence=tumor_sequence,
    trimming=trimming,
    bait_intervals=bait_intervals,
    target_intervals=target_intervals,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals,
    picard_metric_accumulation_level=picard_metric_accumulation_level,
    qc_minimum_mapping_quality=qc_minimum_mapping_quality,
    qc_minimum_base_quality=qc_minimum_base_quality,
    final_name="~{tumor_name}.bam"
  }

  call aen.alignmentExomeNonhuman as normalAlignmentAndQc {
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
    sequence=normal_sequence,
    trimming=trimming,
    bait_intervals=bait_intervals,
    target_intervals=target_intervals,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals,
    picard_metric_accumulation_level=picard_metric_accumulation_level,
    qc_minimum_mapping_quality=qc_minimum_mapping_quality,
    qc_minimum_base_quality=qc_minimum_base_quality,
    final_name="~{normal_name}.bam"
  }

  call ile.intervalListExpand as padTargetIntervals {
    input:
    interval_list=target_intervals,
    roi_padding=target_interval_padding
  }

  call dvn.detectVariantsNonhuman as detectVariants {
    input:
    tumor_bam=tumorAlignmentAndQc.bam,
    tumor_bam_bai=tumorAlignmentAndQc.bam_bai,
    normal_bam=normalAlignmentAndQc.bam,
    normal_bam_bai=normalAlignmentAndQc.bam_bai,
    roi_intervals=padTargetIntervals.expanded_interval_list,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    strelka_cpu_reserved=strelka_cpu_reserved,
    scatter_count=scatter_count,
    varscan_strand_filter=varscan_strand_filter,
    varscan_min_coverage=varscan_min_coverage,
    varscan_min_var_freq=varscan_min_var_freq,
    varscan_p_value=varscan_p_value,
    varscan_max_normal_freq=varscan_max_normal_freq,
    fp_min_var_freq=fp_min_var_freq,
    filter_somatic_llr_threshold=filter_somatic_llr_threshold,
    filter_somatic_llr_tumor_purity=filter_somatic_llr_tumor_purity,
    filter_somatic_llr_normal_contamination_rate=filter_somatic_llr_normal_contamination_rate,
    vep_cache_dir_zip=vep_cache_dir_zip,
    vep_ensembl_assembly=vep_ensembl_assembly,
    vep_ensembl_version=vep_ensembl_version,
    vep_ensembl_species=vep_ensembl_species,
    synonyms_file=synonyms_file,
    annotate_coding_only=annotate_coding_only,
    vep_pick=vep_pick,
    cle_vcf_filter=cle_vcf_filter,
    variants_to_table_fields=variants_to_table_fields,
    variants_to_table_genotype_fields=variants_to_table_genotype_fields,
    vep_to_table_fields=vep_to_table_fields,
    tumor_sample_name=tumor_sample_name,
    normal_sample_name=normal_sample_name,
    strelka_exome_mode=true
  }

  call btc.bamToCram as tumorBamToCram {
    input:
    bam=tumorAlignmentAndQc.bam,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict
  }

  call ic.indexCram as tumorIndexCram {
    input: cram=tumorBamToCram.cram
  }

  call btc.bamToCram as normalBamToCram {
    input:
    bam=normalAlignmentAndQc.bam,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict
  }

  call ic.indexCram as normalIndexCram {
    input: cram=normalBamToCram.cram
  }

  output {
    File tumor_cram = tumorIndexCram.indexed_cram
    File tumor_cram_crai = tumorIndexCram.indexed_cram_crai
    File tumor_mark_duplicates_metrics = tumorAlignmentAndQc.mark_duplicates_metrics
    QCMetrics tumor_qc_metrics = tumorAlignmentAndQc.qc_metrics

    File normal_cram = normalIndexCram.indexed_cram
    File normal_cram_crai = normalIndexCram.indexed_cram_crai
    File normal_mark_duplicates_metrics = normalAlignmentAndQc.mark_duplicates_metrics
    QCMetrics normal_qc_metrics = normalAlignmentAndQc.qc_metrics

    File mutect_unfiltered_vcf = detectVariants.mutect_unfiltered_vcf
    File mutect_unfiltered_vcf_tbi = detectVariants.mutect_unfiltered_vcf_tbi
    File mutect_filtered_vcf = detectVariants.mutect_filtered_vcf
    File mutect_filtered_vcf_tbi = detectVariants.mutect_filtered_vcf_tbi

    File strelka_unfiltered_vcf = detectVariants.strelka_unfiltered_vcf
    File strelka_unfiltered_vcf_tbi = detectVariants.strelka_unfiltered_vcf_tbi
    File strelka_filtered_vcf = detectVariants.strelka_filtered_vcf
    File strelka_filtered_vcf_tbi = detectVariants.strelka_filtered_vcf_tbi

    File varscan_unfiltered_vcf = detectVariants.varscan_unfiltered_vcf
    File varscan_unfiltered_vcf_tbi = detectVariants.varscan_unfiltered_vcf_tbi
    File varscan_filtered_vcf = detectVariants.varscan_filtered_vcf
    File varscan_filtered_vcf_tbi = detectVariants.varscan_filtered_vcf_tbi

    File final_vcf = detectVariants.final_vcf
    File final_vcf_tbi = detectVariants.final_vcf_tbi
    File final_filtered_vcf = detectVariants.final_filtered_vcf
    File final_filtered_vcf_tbi = detectVariants.final_filtered_vcf_tbi

    File final_tsv = detectVariants.final_tsv
    File vep_summary = detectVariants.vep_summary
    File tumor_snv_bam_readcount_tsv = detectVariants.tumor_snv_bam_readcount_tsv
    File tumor_indel_bam_readcount_tsv = detectVariants.tumor_indel_bam_readcount_tsv
    File normal_snv_bam_readcount_tsv = detectVariants.normal_snv_bam_readcount_tsv
    File normal_indel_bam_readcount_tsv = detectVariants.normal_indel_bam_readcount_tsv
  }
}
