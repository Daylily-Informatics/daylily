version 1.0

import "alignment_exome.wdl" as ae
import "detect_variants.wdl" as dv

import "tools/add_string_at_line.wdl" as asal
import "tools/add_string_at_line_bgzipped.wdl" as asalb
import "tools/concordance.wdl" as c
import "tools/interval_list_expand.wdl" as ile
import "tools/index_vcf.wdl" as iv
import "tools/index_cram.wdl" as ic
import "tools/bam_to_cram.wdl" as btc

import "types.wdl"  # !UnusedImport

workflow somaticExomeCle {
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

    Array[SequenceData] tumor_sequence
    String tumor_name = "tumor"
    Array[SequenceData] normal_sequence
    String normal_name = "normal"

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
    Int qc_minimum_mapping_quality = 0
    Int qc_minimum_base_quality = 0
    Int strelka_cpu_reserved = 8
    Int scatter_count
    Int? varscan_strand_filter
    Int? varscan_min_coverage
    Float? varscan_min_var_freq
    Float? varscan_p_value
    Float? varscan_max_normal_freq
    Float? fp_min_var_freq
    File docm_vcf
    File docm_vcf_tbi
    String? gnomad_field_name
    Float? filter_gnomADe_maximum_population_allele_frequency
    Boolean filter_docm_variants = true
    Int filter_minimum_depth = 20
    Float? filter_somatic_llr_threshold
    Float? filter_somatic_llr_tumor_purity
    Float? filter_somatic_llr_normal_contamination_rate
    File vep_cache_dir_zip
    String vep_ensembl_assembly
    String vep_ensembl_version
    String vep_ensembl_species
    File? synonyms_file
    Boolean? annotate_coding_only
    String? vep_pick
    Boolean cle_vcf_filter = false
    Array[String] variants_to_table_fields = ["CHROM", "POS", "ID", "REF", "ALT", "set", "AC", "AF"]
    Array[String] variants_to_table_genotype_fields = ["GT", "AD"]
    Array[String] vep_to_table_fields = ["HGVSc", "HGVSp"]
    Array[VepCustomAnnotation] vep_custom_annotations
    File somalier_vcf
    String disclaimer_text = "This laboratory developed test (LDT) was developed and its performance characteristics determined by the CLIA Licensed Environment laboratory at the McDonnell Genome Institute at Washington University (MGI-CLE, CLIA #26D2092546, CAP #9047655), Dr. David H. Spencer MD, PhD, FCAP, Medical Director. 4444 Forest Park Avenue, Rm 4127 St. Louis, Missouri 63108 (314) 286-1460 Fax: (314) 286-1810. The MGI-CLE laboratory is regulated under CLIA as certified to perform high-complexity testing. This test has not been cleared or approved by the FDA."
    String disclaimer_version
    String tumor_sample_name
    String normal_sample_name
  }

  call ae.alignmentExome as tumorAlignmentAndQc {
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
    qc_minimum_base_quality=qc_minimum_base_quality,
    final_name="~{tumor_name}.bam"
  }

  call ae.alignmentExome as normalAlignmentAndQc {
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
    qc_minimum_base_quality=qc_minimum_base_quality,
    final_name="~{normal_name}.bam"
  }

  call c.concordance {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam_1=tumorAlignmentAndQc.bam,
    bam_1_bai=tumorAlignmentAndQc.bam_bai,
    bam_2=normalAlignmentAndQc.bam,
    bam_2_bai=normalAlignmentAndQc.bam_bai,
    vcf=somalier_vcf
  }

  call ile.intervalListExpand as padTargetIntervals {
    input:
    interval_list=target_intervals,
    roi_padding=target_interval_padding
  }

  call dv.detectVariants {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_bam=tumorAlignmentAndQc.bam,
    tumor_bam_bai=tumorAlignmentAndQc.bam_bai,
    normal_bam=normalAlignmentAndQc.bam,
    normal_bam_bai=normalAlignmentAndQc.bam_bai,
    roi_intervals=padTargetIntervals.expanded_interval_list,
    strelka_exome_mode=true,
    strelka_cpu_reserved=strelka_cpu_reserved,
    scatter_count=scatter_count,
    varscan_strand_filter=varscan_strand_filter,
    varscan_min_coverage=varscan_min_coverage,
    varscan_min_var_freq=varscan_min_var_freq,
    varscan_p_value=varscan_p_value,
    varscan_max_normal_freq=varscan_max_normal_freq,
    fp_min_var_freq=fp_min_var_freq,
    docm_vcf=docm_vcf,
    docm_vcf_tbi=docm_vcf_tbi,
    gnomad_field_name=gnomad_field_name,
    filter_gnomADe_maximum_population_allele_frequency=filter_gnomADe_maximum_population_allele_frequency,
    filter_docm_variants=filter_docm_variants,
    filter_minimum_depth=filter_minimum_depth,
    vep_cache_dir_zip=vep_cache_dir_zip,
    vep_ensembl_assembly=vep_ensembl_assembly,
    vep_ensembl_version=vep_ensembl_version,
    vep_ensembl_species=vep_ensembl_species,
    synonyms_file=synonyms_file,
    annotate_coding_only=annotate_coding_only,
    vep_pick=vep_pick,
    cle_vcf_filter=cle_vcf_filter,
    filter_somatic_llr_threshold=filter_somatic_llr_threshold,
    filter_somatic_llr_tumor_purity=filter_somatic_llr_tumor_purity,
    filter_somatic_llr_normal_contamination_rate=filter_somatic_llr_normal_contamination_rate,
    variants_to_table_fields=variants_to_table_fields,
    variants_to_table_genotype_fields=variants_to_table_genotype_fields,
    vep_to_table_fields=vep_to_table_fields,
    tumor_sample_name=tumor_sample_name,
    normal_sample_name=normal_sample_name,
    vep_custom_annotations=vep_custom_annotations,
  }

  call asal.addStringAtLine as addDisclaimerToFinalTsv {
    input:
    input_file=detectVariants.final_tsv,
    line_number=1,
    some_text="#~{disclaimer_text}",
    output_name=basename(detectVariants.final_tsv)
  }

  call asal.addStringAtLine as addDisclaimerVersionToFinalTsv {
    input:
    input_file=addDisclaimerToFinalTsv.output_file,
    line_number=2,
    some_text="#The software version is ~{disclaimer_version}",
    output_name=basename(addDisclaimerToFinalTsv.output_file)
  }

  call asalb.addStringAtLineBgzipped as addDisclaimerToFinalFilteredVcf {
    input:
    input_file=detectVariants.final_filtered_vcf,
    line_number=2,
    some_text="##DisclaimerText=~{disclaimer_text}",
    output_name=basename(detectVariants.final_filtered_vcf)
  }

  call asalb.addStringAtLineBgzipped as addDisclaimerVersionToFinalFilteredVcf {
    input:
    input_file=addDisclaimerToFinalFilteredVcf.output_file,
    line_number=3,
    some_text="##CLESoftwareVersion=~{disclaimer_version}",
    output_name=basename(addDisclaimerToFinalFilteredVcf.output_file)
  }

  call iv.indexVcf as annotatedFilterVcfIndex {
    input: vcf=addDisclaimerVersionToFinalFilteredVcf.output_file
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
    QCMetrics normal_qc_metrics = tumorAlignmentAndQc.qc_metrics
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
    File final_filtered_vcf = annotatedFilterVcfIndex.indexed_vcf
    File final_filtered_vcf_tbi = annotatedFilterVcfIndex.indexed_vcf_tbi
    File final_tsv = addDisclaimerVersionToFinalTsv.output_file
    File vep_summary = detectVariants.vep_summary
    File tumor_snv_bam_readcount_tsv = detectVariants.tumor_snv_bam_readcount_tsv
    File tumor_indel_bam_readcount_tsv = detectVariants.tumor_indel_bam_readcount_tsv
    File normal_snv_bam_readcount_tsv = detectVariants.normal_snv_bam_readcount_tsv
    File normal_indel_bam_readcount_tsv = detectVariants.normal_indel_bam_readcount_tsv
    File somalier_concordance_metrics = concordance.somalier_pairs
    File somalier_concordance_statistics = concordance.somalier_samples
  }
}
