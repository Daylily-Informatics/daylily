version 1.0

import "../subworkflows/cnvkit_single_sample.wdl" as css
import "../subworkflows/sv_depth_caller_filter.wdl" as sdcf
import "../subworkflows/sv_paired_read_caller_filter.wdl" as sprcf
import "../subworkflows/merge_svs.wdl" as mss
import "../tools/bgzip.wdl" as b
import "../tools/cnvnator.wdl" as c
import "../tools/index_vcf.wdl" as iv
import "../tools/manta_somatic.wdl" as ms
import "../tools/smoove.wdl" as s


workflow singleSampleSvCallers {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    Boolean? cnvkit_diagram
    Boolean? cnvkit_drop_low_coverage
    String cnvkit_method  # enum ["hybrid", "amplicon", "wgs"]
    File? cnvkit_reference_cnn
    Boolean? cnvkit_scatter_plot
    Boolean? cnvkit_male_reference
    String? cnvkit_vcf_name
    File? manta_call_regions
    File? manta_call_regions_tbi
    Boolean? manta_non_wgs
    Boolean? manta_output_contigs
    Int merge_max_distance
    Int merge_min_svs
    Boolean merge_same_type
    Boolean merge_same_strand
    Boolean merge_estimate_sv_distance
    Int merge_min_sv_size
    File? smoove_exclude_regions
    File? snps_vcf
    String genome_build
    Float? sv_alt_abundance_percentage
    Int? sv_paired_count
    Int? sv_split_count
    Float? cnv_deletion_depth
    Float? cnv_duplication_depth
    Int? cnv_filter_min_size
    File? blocklist_bedpe
  }

  call css.cnvkitSingleSample as runCnvkit {
    input:
    diagram=cnvkit_diagram,
    drop_low_coverage=cnvkit_drop_low_coverage,
    method=cnvkit_method,
    reference_cnn=cnvkit_reference_cnn,
    tumor_bam=bam,
    tumor_bam_bai=bam_bai,
    scatter_plot=cnvkit_scatter_plot,
    male_reference=cnvkit_male_reference,
    cnvkit_vcf_name=cnvkit_vcf_name,
    segment_filter="cn",
    reference=reference,
    reference_fai=reference_fai
  }

  call b.bgzip as runCnvkitRawBgzip {
    input: file=runCnvkit.cnvkit_vcf
  }

  call iv.indexVcf as runCnvkitRawIndex {
    input: vcf=runCnvkitRawBgzip.bgzipped_file
  }

  call sdcf.svDepthCallerFilter as runCnvkitFilter {
    input:
    deletion_depth=cnv_deletion_depth,
    duplication_depth=cnv_duplication_depth,
    min_sv_size=cnv_filter_min_size,
    output_vcf_name="filtered_cnvkit.vcf",
    sv_vcf=runCnvkit.cnvkit_vcf,
    vcf_source="cnvkit"
  }

  call c.cnvnator as runCnvnator {
    input:
    bam=bam,
    reference=reference,
    reference_fai=reference_fai,
    sample_name="cnvnator"
  }

  call b.bgzip as runCnvnatorRawBgzip {
    input: file=runCnvnator.vcf
  }

  call iv.indexVcf as runCnvnatorRawIndex {
    input: vcf=runCnvnatorRawBgzip.bgzipped_file
  }

  call sdcf.svDepthCallerFilter as runCnvnatorFilter {
    input:
    deletion_depth=cnv_deletion_depth,
    duplication_depth=cnv_duplication_depth,
    min_sv_size=cnv_filter_min_size,
    output_vcf_name="filtered_cnvnator.vcf",
    sv_vcf=runCnvnator.vcf,
    vcf_source="cnvnator"
  }

  call ms.mantaSomatic as runManta {
    input:
    call_regions=manta_call_regions,
    non_wgs=manta_non_wgs,
    output_contigs=manta_output_contigs,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_bam=bam,
    tumor_bam_bai=bam_bai
  }

  call sprcf.svPairedReadCallerFilter as runMantaFilter {
    input:
    abundance_percentage=sv_alt_abundance_percentage,
    bam=bam,
    deletion_depth=cnv_deletion_depth,
    duplication_depth=cnv_duplication_depth,
    output_vcf_name="filtered_manta.vcf",
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    snps_vcf=snps_vcf,
    sv_paired_count=sv_paired_count,
    sv_split_count=sv_split_count,
    sv_vcf=select_first([runManta.tumor_only_variants]),
    vcf_source="manta"
  }

  call s.smoove as runSmoove {
    input:
    bams=[bam],
    exclude_regions=smoove_exclude_regions,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict
  }

  call sprcf.svPairedReadCallerFilter as runSmooveFilter {
    input:
    abundance_percentage=sv_alt_abundance_percentage,
    bam=bam,
    deletion_depth=cnv_deletion_depth,
    duplication_depth=cnv_duplication_depth,
    output_vcf_name="filtered_smoove.vcf",
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    snps_vcf=snps_vcf,
    sv_paired_count=sv_paired_count,
    sv_split_count=sv_split_count,
    sv_vcf=runSmoove.output_vcf,
    vcf_source="smoove"
  }

  call mss.mergeSvs as runMerge {
    input:
    estimate_sv_distance=merge_estimate_sv_distance,
    genome_build=genome_build,
    max_distance_to_merge=merge_max_distance,
    minimum_sv_calls=merge_min_svs,
    minimum_sv_size=merge_min_sv_size,
    same_strand=merge_same_strand,
    same_type=merge_same_type,
    snps_vcf=snps_vcf,
    blocklist_bedpe=blocklist_bedpe,
    sv_vcfs=[runCnvkitFilter.filtered_vcf, runCnvnatorFilter.filtered_vcf, runMantaFilter.filtered_vcf, runSmooveFilter.filtered_vcf]
  }

  output {
    File? cn_diagram = runCnvkit.cn_diagram
    File? cn_scatter_plot = runCnvkit.cn_scatter_plot
    File tumor_antitarget_coverage = runCnvkit.tumor_antitarget_coverage
    File tumor_target_coverage = runCnvkit.tumor_target_coverage
    File tumor_bin_level_ratios = runCnvkit.tumor_bin_level_ratios
    File tumor_segmented_ratios = runCnvkit.tumor_segmented_ratios
    File cnvkit_filtered_vcf = runCnvkitFilter.filtered_vcf
    File cnvkit_vcf = runCnvkitRawIndex.indexed_vcf
    File cnvnator_cn_file = runCnvnator.cn_file
    File cnvnator_filtered_vcf = runCnvnatorFilter.filtered_vcf
    File cnvnator_root = runCnvnator.root_file
    File cnvnator_vcf = runCnvnatorRawIndex.indexed_vcf
    File? manta_diploid_variants = runManta.diploid_variants
    File manta_filtered_vcf = runMantaFilter.filtered_vcf
    File? manta_somatic_variants = runManta.somatic_variants
    File manta_all_candidates = runManta.all_candidates
    File manta_small_candidates = runManta.small_candidates
    File? manta_tumor_only_variants = runManta.tumor_only_variants
    File smoove_output_variants = runSmoove.output_vcf
    File smoove_filtered_vcf = runSmooveFilter.filtered_vcf
    File survivor_merged_vcf = runMerge.survivor_merged_sv_vcf
    File survivor_merged_annotated_tsv = runMerge.survivor_merged_annotated_tsv
    File bcftools_merged_vcf = runMerge.bcftools_merged_sv_vcf
    File bcftools_merged_annotated_tsv = runMerge.bcftools_merged_annotated_tsv
    File bcftools_merged_filtered_annotated_tsv = runMerge.bcftools_merged_filtered_annotated_tsv
  }
}
