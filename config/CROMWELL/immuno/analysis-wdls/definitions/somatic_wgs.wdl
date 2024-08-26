version 1.0

import "subworkflows/sequence_to_bqsr.wdl" as s2b
import "subworkflows/qc_wgs.wdl" as qe

import "detect_variants_wgs.wdl" as dvw

import "tools/bam_to_cram.wdl" as btc
import "tools/concordance.wdl" as c
import "tools/index_cram.wdl" as ic
import "tools/manta_somatic.wdl" as ms
import "tools/cnvkit_batch.wdl" as cb

import "types.wdl"

workflow somaticWgs {
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

    Array[File] bqsr_known_sites
    Array[File] bqsr_known_sites_tbi

    File target_intervals
    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals

    File omni_vcf
    File omni_vcf_tbi

    String picard_metric_accumulation_level

    Int? qc_minimum_mapping_quality
    Int? qc_minimum_base_quality

    Int strelka_cpu_reserved = 8
    File? strelka_call_regions
    File? strelka_call_regions_tbi
    Int scatter_count
    Int? varscan_strand_filter
    Int? varscan_min_coverage
    Float? varscan_min_var_freq
    Float? varscan_p_value
    Float? varscan_max_normal_freq
    Float? fp_min_var_freq
    File docm_vcf
    File docm_vcf_tbi
    Boolean filter_docm_variants = true

    Float? filter_somatic_llr_threshold
    Float? filter_somatic_llr_tumor_purity
    Float? filter_somatic_llr_normal_contamination_rate
    File vep_cache_dir_zip
    String vep_ensembl_assembly
    String vep_ensembl_version
    String vep_ensembl_species
    File? synonyms_file
    Boolean? annotate_coding_only
    String? vep_pick  # enum ["pick", "flag_pick", "pick_allele", "per_gene", "pick_allele_gene", "flag_pick_allele", "flag_pick_allele_gene"]
    Boolean cle_vcf_filter = false
    Array[String] variants_to_table_fields = ["CHROM", "POS", "ID", "REF", "ALT", "set", "AC", "AF"]
    Array[String] variants_to_table_genotype_fields = ["GT", "AD"]
    Array[String] vep_to_table_fields = ["HGVSc", "HGVSp"]
    Array[VepCustomAnnotation] vep_custom_annotations
    Boolean manta_non_wgs = false
    Boolean? manta_output_contigs
    File somalier_vcf
    String tumor_sample_name
    String normal_sample_name
    File? validated_variants
    File? validated_variants_tbi

    #approximate size of split target bins for CNVkit; if not set a suitable window size will be set by CNVkit automatically
    Int? cnvkit_target_average_size
  }


  call s2b.sequenceToBqsr as tumorAlignment {
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
    unaligned=tumor_sequence,
    trimming=trimming,
    bqsr_known_sites=bqsr_known_sites,
    bqsr_known_sites_tbi=bqsr_known_sites_tbi,
    final_name=tumor_name
  }
  call qe.qcWgs as tumorQc {
    input:
    bam=tumorAlignment.final_bam,
    bam_bai=tumorAlignment.final_bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    intervals=target_intervals,
    summary_intervals=summary_intervals,
    omni_vcf=omni_vcf,
    omni_vcf_tbi=omni_vcf_tbi,
    picard_metric_accumulation_level=picard_metric_accumulation_level,
    minimum_mapping_quality=qc_minimum_mapping_quality,
    minimum_base_quality=qc_minimum_base_quality
  }

  call s2b.sequenceToBqsr as normalAlignment {
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
    unaligned=normal_sequence,
    trimming=trimming,
    bqsr_known_sites=bqsr_known_sites,
    bqsr_known_sites_tbi=bqsr_known_sites_tbi,
    final_name=normal_name
  }

  call qe.qcWgs as normalQc {
    input:
    bam=normalAlignment.final_bam,
    bam_bai=normalAlignment.final_bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    intervals=target_intervals,
    summary_intervals=summary_intervals,
    omni_vcf=omni_vcf,
    omni_vcf_tbi=omni_vcf_tbi,
    picard_metric_accumulation_level=picard_metric_accumulation_level,
    minimum_mapping_quality=qc_minimum_mapping_quality,
    minimum_base_quality=qc_minimum_base_quality
  }

  call c.concordance {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam_1=tumorAlignment.final_bam,
    bam_1_bai=tumorAlignment.final_bam_bai,
    bam_2=normalAlignment.final_bam,
    bam_2_bai=normalAlignment.final_bam_bai,
    vcf=somalier_vcf
  }

  call dvw.detectVariantsWgs as detectVariants {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_bam=tumorAlignment.final_bam,
    tumor_bam_bai=tumorAlignment.final_bam_bai,
    normal_bam=normalAlignment.final_bam,
    normal_bam_bai=normalAlignment.final_bam_bai,
    roi_intervals=target_intervals,
    strelka_exome_mode=false,
    strelka_cpu_reserved=strelka_cpu_reserved,
    strelka_call_regions=strelka_call_regions,
    strelka_call_regions_tbi=strelka_call_regions_tbi,
    scatter_count=scatter_count,
    varscan_strand_filter=varscan_strand_filter,
    varscan_min_coverage=varscan_min_coverage,
    varscan_min_var_freq=varscan_min_var_freq,
    varscan_p_value=varscan_p_value,
    varscan_max_normal_freq=varscan_max_normal_freq,
    fp_min_var_freq=fp_min_var_freq,
    docm_vcf=docm_vcf,
    docm_vcf_tbi=docm_vcf_tbi,
    filter_docm_variants=filter_docm_variants,
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
    vep_custom_annotations=vep_custom_annotations,
    validated_variants=validated_variants,
    validated_variants_tbi=validated_variants_tbi
  }

  call ms.mantaSomatic as manta {
    input:
    normal_bam=normalAlignment.final_bam,
    normal_bam_bai=normalAlignment.final_bam_bai,
    tumor_bam=tumorAlignment.final_bam,
    tumor_bam_bai=tumorAlignment.final_bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    non_wgs=manta_non_wgs,
    output_contigs=manta_output_contigs
  }

  call cb.cnvkitBatch as cnvkit {
    input:
    tumor_bam=tumorAlignment.final_bam,
    tumor_bam_bai=tumorAlignment.final_bam_bai,
    normal_bam=normalAlignment.final_bam,
    normal_bam_bai=normalAlignment.final_bam_bai,
    reference_fasta=reference,
    target_average_size=cnvkit_target_average_size,
    method="wgs"
  }

  call btc.bamToCram as tumorBamToCram {
    input:
    bam=tumorAlignment.final_bam,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict
  }

  call ic.indexCram as tumorIndexCram {
    input: cram=tumorBamToCram.cram
  }

  call btc.bamToCram as normalBamToCram {
    input:
    bam=normalAlignment.final_bam,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict
  }

  call ic.indexCram as normalIndexCram {
    input: cram=normalBamToCram.cram
  }


  output {
    File tumor_cram = tumorIndexCram.indexed_cram
    File tumor_mark_duplicates_metrics = tumorAlignment.mark_duplicates_metrics_file
    QCMetrics tumor_qc_metrics = tumorQc.metrics
##normal alignment and qc
    File normal_cram = normalIndexCram.indexed_cram
    File normal_mark_duplicates_metrics = normalAlignment.mark_duplicates_metrics_file
    QCMetrics normal_qc_metrics = normalQc.metrics
##variant calling
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
    File docm_filtered_vcf = detectVariants.docm_filtered_vcf
    File docm_filtered_vcf_tbi = detectVariants.docm_filtered_vcf_tbi
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
    #manta
    File? diploid_variants = manta.diploid_variants
    File? diploid_variants_tbi = manta.diploid_variants_tbi
    File? somatic_variants = manta.somatic_variants
    File? somatic_variants_tbi = manta.somatic_variants_tbi
    File all_candidates = manta.all_candidates
    File all_candidates_tbi = manta.all_candidates_tbi
    File small_candidates = manta.small_candidates
    File small_candidates_tbi = manta.small_candidates_tbi
    File? tumor_only_variants = manta.tumor_only_variants
    File? tumor_only_variants_tbi = manta.tumor_only_variants_tbi
    #cnvkit
    File? intervals_antitarget = cnvkit.intervals_antitarget
    File? intervals_target = cnvkit.intervals_target
    File? normal_antitarget_coverage = cnvkit.normal_antitarget_coverage
    File? normal_target_coverage = cnvkit.normal_target_coverage
    File? reference_coverage = cnvkit.reference_coverage
    File tumor_antitarget_coverage = cnvkit.tumor_antitarget_coverage
    File tumor_target_coverage = cnvkit.tumor_target_coverage
    File tumor_bin_level_ratios = cnvkit.tumor_bin_level_ratios
    File tumor_segmented_ratios = cnvkit.tumor_segmented_ratios
    File? cn_diagram = cnvkit.cn_diagram
    File? cn_scatter_plot = cnvkit.cn_scatter_plot

##sample concordance check
    File somalier_concordance_metrics = concordance.somalier_pairs
    File somalier_concordance_statistics = concordance.somalier_samples
  }
}
