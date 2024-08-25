version 1.0

import "types.wdl"

import "alignment_wgs.wdl" as aw
import "subworkflows/germline_detect_variants.wdl" as gdv
import "subworkflows/single_sample_sv_callers.wdl" as sssc
import "tools/add_string_at_line.wdl" as asal
import "tools/add_string_at_line_bgzipped.wdl" as asalb
import "tools/bam_to_cram.wdl" as btc
import "tools/index_cram.wdl" as ic
import "tools/index_vcf.wdl" as iv

workflow germlineWgs {
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
    String picard_metric_accumulation_level
    Array[String] gvcf_gq_bands
    Array[Array[String]] intervals
    Int? ploidy
    File qc_intervals
    File variant_reporting_intervals
    File vep_cache_dir_zip
    String vep_ensembl_assembly
    String vep_ensembl_version
    String vep_ensembl_species
    Array[String]? vep_plugins
    File? synonyms_file
    Boolean? annotate_coding_only
    Array[File] bqsr_known_sites
    Array[File] bqsr_known_sites_tbi
    Int? minimum_mapping_quality
    Int? minimum_base_quality
    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals
    Array[VepCustomAnnotation] vep_custom_annotations
    Boolean? cnvkit_diagram
    Boolean? cnvkit_drop_low_coverage
    String cnvkit_method = "wgs"  # enum ["hybrid", "amplicon", "wgs"]
    File? cnvkit_reference_cnn
    Boolean? cnvkit_scatter_plot
    Boolean? cnvkit_male_reference
    String? cnvkit_vcf_name
    File? manta_call_regions
    File? manta_call_regions_tbi
    Boolean? manta_non_wgs
    Boolean? manta_output_contigs
    File? smoove_exclude_regions
    Int merge_max_distance
    Int merge_min_svs
    Boolean merge_same_type
    Boolean merge_same_strand
    Boolean merge_estimate_sv_distance
    Int merge_min_sv_size
    Float? sv_filter_alt_abundance_percentage
    Int? sv_filter_paired_count
    Int? sv_filter_split_count
    Float? cnv_filter_deletion_depth
    Float? cnv_filter_duplication_depth
    Array[String]? variants_to_table_fields
    Array[String]? variants_to_table_genotype_fields
    Array[String]? vep_to_table_fields
    Int? cnv_filter_min_size
    File? blocklist_bedpe
    String disclaimer_text = "Workflow source can be found at https://github.com/genome/analysis-workflows"
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
    intervals=qc_intervals,
    picard_metric_accumulation_level=picard_metric_accumulation_level,
    bqsr_known_sites=bqsr_known_sites,
    bqsr_known_sites_tbi=bqsr_known_sites_tbi,
    minimum_mapping_quality=minimum_mapping_quality,
    minimum_base_quality=minimum_base_quality,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals
  }

  call gdv.germlineDetectVariants as detectVariants {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=alignmentAndQc.bam,
    bai=alignmentAndQc.bai,
    verify_bam_id_metrics=select_first([alignmentAndQc.qc_metrics.verify_bam_id_metrics]),
    gvcf_gq_bands=gvcf_gq_bands,
    intervals=intervals,
    ploidy=ploidy,
    vep_cache_dir_zip=vep_cache_dir_zip,
    synonyms_file=synonyms_file,
    annotate_coding_only=annotate_coding_only,
    limit_variant_intervals=variant_reporting_intervals,
    vep_custom_annotations=vep_custom_annotations,
    vep_ensembl_assembly=vep_ensembl_assembly,
    vep_ensembl_version=vep_ensembl_version,
    vep_ensembl_species=vep_ensembl_species,
    vep_plugins=vep_plugins,
    vep_to_table_fields=vep_to_table_fields,
    variants_to_table_fields=variants_to_table_fields,
    variants_to_table_genotype_fields=variants_to_table_genotype_fields
  }

  call asalb.addStringAtLineBgzipped as addDisclaimerFilteredVcf {
    input:
    input_file=detectVariants.filtered_vcf,
    line_number=2,
    some_text="##disclaimer=~{disclaimer_text}",
    output_name=basename(detectVariants.filtered_vcf)
  }

  call iv.indexVcf as indexDisclaimerFilteredVcf {
    input: vcf=addDisclaimerFilteredVcf.output_file
  }

  call asalb.addStringAtLineBgzipped as addDisclaimerFinalVcf {
    input:
    input_file=detectVariants.final_vcf,
    line_number=2,
    some_text="##disclaimer=~{disclaimer_text}",
    output_name=basename(detectVariants.final_vcf)
  }

  call iv.indexVcf as indexDisclaimerFinalVcf {
    input: vcf=addDisclaimerFinalVcf.output_file
  }

  call asal.addStringAtLine as addDisclaimerFilteredTsv {
    input:
    input_file=detectVariants.filtered_tsv,
    line_number=1,
    some_text="#~{disclaimer_text}",
    output_name=basename(detectVariants.filtered_tsv)
  }

  call asal.addStringAtLine as addDisclaimerFinalTsv {
    input:
    input_file=detectVariants.final_tsv,
    line_number=1,
    some_text="#~{disclaimer_text}",
    output_name=basename(detectVariants.final_tsv)
  }

  call sssc.singleSampleSvCallers as svDetectVariants {
    input:
    bam=alignmentAndQc.bam,
    bam_bai=alignmentAndQc.bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    cnvkit_diagram=cnvkit_diagram,
    cnvkit_drop_low_coverage=cnvkit_drop_low_coverage,
    cnvkit_method=cnvkit_method,
    cnvkit_reference_cnn=cnvkit_reference_cnn,
    cnvkit_scatter_plot=cnvkit_scatter_plot,
    cnvkit_male_reference=cnvkit_male_reference,
    cnvkit_vcf_name=cnvkit_vcf_name,
    cnv_deletion_depth=cnv_filter_deletion_depth,
    cnv_duplication_depth=cnv_filter_duplication_depth,
    cnv_filter_min_size=cnv_filter_min_size,
    manta_call_regions=manta_call_regions,
    manta_non_wgs=manta_non_wgs,
    manta_output_contigs=manta_output_contigs,
    smoove_exclude_regions=smoove_exclude_regions,
    merge_max_distance=merge_max_distance,
    merge_min_svs=merge_min_svs,
    merge_same_type=merge_same_type,
    merge_same_strand=merge_same_strand,
    merge_estimate_sv_distance=merge_estimate_sv_distance,
    merge_min_sv_size=merge_min_sv_size,
    snps_vcf=detectVariants.final_vcf,
    sv_alt_abundance_percentage=sv_filter_alt_abundance_percentage,
    sv_paired_count=sv_filter_paired_count,
    sv_split_count=sv_filter_split_count,
    genome_build=vep_ensembl_assembly,
    blocklist_bedpe=blocklist_bedpe
  }

  call asalb.addStringAtLineBgzipped as addDisclaimerSurvivorSvVcf {
    input:
    input_file=svDetectVariants.survivor_merged_vcf,
    line_number=2,
    some_text="##disclaimer=~{disclaimer_text}",
    output_name=basename(svDetectVariants.survivor_merged_vcf)
  }

  call asalb.addStringAtLineBgzipped as addDisclaimerBcftoolsSvVcf {
    input:
    input_file=svDetectVariants.bcftools_merged_vcf,
    line_number=2,
    some_text="##disclaimer=~{disclaimer_text}",
    output_name=basename(svDetectVariants.bcftools_merged_vcf)
  }

  call asal.addStringAtLine as addDisclaimerSurvivorSvTsv {
    input:
    input_file=svDetectVariants.survivor_merged_annotated_tsv,
    line_number=1,
    some_text="#~{disclaimer_text}",
    output_name=basename(svDetectVariants.survivor_merged_annotated_tsv)
  }

  call asal.addStringAtLine as addDisclaimerBcftoolsSvTsv {
    input:
    input_file=svDetectVariants.bcftools_merged_annotated_tsv,
    line_number=1,
    some_text="#~{disclaimer_text}",
    output_name=basename(svDetectVariants.bcftools_merged_annotated_tsv)
  }

  call asal.addStringAtLine as addDisclaimerBcftoolsFilteredSvTsv {
    input:
    input_file=svDetectVariants.bcftools_merged_filtered_annotated_tsv,
    line_number=1,
    some_text="#~{disclaimer_text}",
    output_name=basename(svDetectVariants.bcftools_merged_filtered_annotated_tsv)
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
    File raw_vcf = detectVariants.raw_vcf
    File raw_vcf_tbi = detectVariants.raw_vcf_tbi
    File final_vcf = indexDisclaimerFinalVcf.indexed_vcf
    File final_vcf_tbi = indexDisclaimerFinalVcf.indexed_vcf_tbi
    File filtered_vcf = indexDisclaimerFilteredVcf.indexed_vcf
    File filtered_vcf_tbi = indexDisclaimerFilteredVcf.indexed_vcf_tbi
    File vep_summary = detectVariants.vep_summary
    File? cn_diagram = svDetectVariants.cn_diagram
    File? cn_scatter_plot = svDetectVariants.cn_scatter_plot
    File tumor_antitarget_coverage = svDetectVariants.tumor_antitarget_coverage
    File tumor_target_coverage = svDetectVariants.tumor_target_coverage
    File tumor_bin_level_ratios = svDetectVariants.tumor_bin_level_ratios
    File tumor_segmented_ratios = svDetectVariants.tumor_segmented_ratios
    File cnvkit_vcf = svDetectVariants.cnvkit_vcf
    File cnvnator_cn_file = svDetectVariants.cnvnator_cn_file
    File cnvnator_root = svDetectVariants.cnvnator_root
    File cnvnator_vcf = svDetectVariants.cnvnator_vcf
    File? manta_diploid_variants = svDetectVariants.manta_diploid_variants
    File? manta_somatic_variants = svDetectVariants.manta_somatic_variants
    File manta_all_candidates = svDetectVariants.manta_all_candidates
    File manta_small_candidates = svDetectVariants.manta_small_candidates
    File? manta_tumor_only_variants = svDetectVariants.manta_tumor_only_variants
    File smoove_output_variants = svDetectVariants.smoove_output_variants
    File final_tsv = addDisclaimerFinalTsv.output_file
    File filtered_tsv = addDisclaimerFilteredTsv.output_file
    File cnvkit_filtered_vcf = svDetectVariants.cnvkit_filtered_vcf
    File cnvnator_filtered_vcf = svDetectVariants.cnvnator_filtered_vcf
    File manta_filtered_vcf = svDetectVariants.manta_filtered_vcf
    File smoove_filtered_vcf = svDetectVariants.smoove_filtered_vcf
    File survivor_merged_vcf = addDisclaimerSurvivorSvVcf.output_file
    File survivor_merged_annotated_tsv = addDisclaimerSurvivorSvTsv.output_file
    File bcftools_merged_vcf = addDisclaimerBcftoolsSvVcf.output_file
    File bcftools_merged_annotated_tsv = addDisclaimerBcftoolsSvTsv.output_file
    File bcftools_merged_filtered_annotated_tsv = addDisclaimerBcftoolsFilteredSvTsv.output_file
  }
}
