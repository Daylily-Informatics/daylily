version 1.0

import "subworkflows/vcf_readcount.wdl" as vr
import "subworkflows/filter_vcf_nonhuman.wdl" as fvn
import "subworkflows/mutect.wdl" as m
import "subworkflows/strelka_and_post_processing.wdl" as sapp
import "subworkflows/varscan_pre_and_post_processing.wdl" as vpapp
import "tools/add_vep_fields_to_table.wdl" as avftt
import "tools/bgzip.wdl" as b
import "tools/combine_variants.wdl" as cv
import "tools/index_vcf.wdl" as iv
import "tools/variants_to_table.wdl" as vtt
import "tools/vep.wdl" as v
import "tools/vt_decompose.wdl" as vd

workflow detectVariantsNonhuman {
  input {
    File reference
    File reference_fai
    File reference_dict

    String tumor_sample_name
    File tumor_bam
    File tumor_bam_bai

    String normal_sample_name
    File normal_bam
    File normal_bam_bai

    File roi_intervals
    Int scatter_count

    Boolean strelka_exome_mode
    Int strelka_cpu_reserved = 8

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
    Boolean? annotate_coding_only
    File? synonyms_file
    # one of [pick, flag_pick, pick-allele, per_gene, pick_allele_gene, flag_pick_allele, flag_pick_allele_gene]
    String? vep_pick
    Array[String] vep_plugins = ["Frameshift", "Wildtype"]

    Int? readcount_minimum_base_quality
    Int? readcount_minimum_mapping_quality

    Float filter_mapq0_threshold = 0.15
    Float? filter_somatic_llr_threshold
    Float? filter_somatic_llr_tumor_purity
    Float? filter_somatic_llr_normal_contamination_rate
    Int filter_minimum_depth = 1
    Boolean cle_vcf_filter = false
    Array[String] variants_to_table_fields = ["CHROM", "POS", "ID", "REF", "ALT", "set", "AC", "AF"]
    Array[String] variants_to_table_genotype_fields = ["GT", "AD"]
    Array[String] vep_to_table_fields = []
  }

  call m.mutect {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_sample_name=tumor_sample_name,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    interval_list=roi_intervals,
    scatter_count=scatter_count,
    fp_min_var_freq=fp_min_var_freq
  }

  call sapp.strelkaAndPostProcessing as strelka {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    interval_list=roi_intervals,
    exome_mode=strelka_exome_mode,
    cpu_reserved=strelka_cpu_reserved,
    normal_sample_name=normal_sample_name,
    tumor_sample_name=normal_sample_name,
    fp_min_var_freq=fp_min_var_freq
  }

  call vpapp.varscanPreAndPostProcessing as varscan {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    interval_list=roi_intervals,
    scatter_count=scatter_count,
    strand_filter=varscan_strand_filter,
    min_coverage=varscan_min_coverage,
    varscan_min_var_freq=varscan_min_var_freq,
    p_value=varscan_p_value,
    max_normal_freq=varscan_max_normal_freq,
    normal_sample_name=normal_sample_name,
    tumor_sample_name=tumor_sample_name
  }

  call cv.combineVariants as combine {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,

    mutect_vcf=mutect.filtered_vcf,
    mutect_vcf_tbi=mutect.filtered_vcf_tbi,

    strelka_vcf=strelka.filtered_vcf,
    strelka_vcf_tbi=strelka.filtered_vcf_tbi,

    varscan_vcf=varscan.filtered_vcf,
    varscan_vcf_tbi=varscan.filtered_vcf_tbi,
  }

  call vd.vtDecompose as decompose {
    input:
    vcf=combine.combined_vcf,
    vcf_tbi=combine.combined_vcf_tbi
  }

  call iv.indexVcf as decomposeIndex {
    input: vcf=decompose.decomposed_vcf
  }

  call v.vep as annotateVariants {
    input:
    vcf=decomposeIndex.indexed_vcf,
    reference=reference,
    reference_dict=reference_dict,
    reference_fai=reference_fai,
    cache_dir_zip=vep_cache_dir_zip,
    coding_only=annotate_coding_only,
    ensembl_assembly=vep_ensembl_assembly,
    ensembl_species=vep_ensembl_species,
    ensembl_version=vep_ensembl_version,
    pick=vep_pick,
    plugins=vep_plugins,
    synonyms_file=synonyms_file
  }

  call vr.vcfReadcount {
    input:
    vcf=annotateVariants.annotated_vcf,
    tumor_sample_name=tumor_sample_name,
    normal_sample_name=normal_sample_name,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    minimum_base_quality=readcount_minimum_base_quality,
    minimum_mapping_quality=readcount_minimum_mapping_quality
  }

  call fvn.filterVcfNonhuman as filterVcf {
    input:
    vcf=vcfReadcount.readcounted_vcf,
    filter_mapq0_threshold=filter_mapq0_threshold,
    filter_somatic_llr_threshold=filter_somatic_llr_threshold,
    filter_somatic_llr_tumor_purity=filter_somatic_llr_tumor_purity,
    filter_somatic_llr_normal_contamination_rate=filter_somatic_llr_normal_contamination_rate,
    filter_minimum_depth=filter_minimum_depth,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    do_cle_vcf_filter=cle_vcf_filter,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    normal_sample_name=normal_sample_name,
    tumor_sample_name=tumor_sample_name
  }

  call b.bgzip as annotatedFilterBgzip {
    input: file=filterVcf.filtered_vcf
  }

  call iv.indexVcf as annotatedFilterIndex {
    input: vcf=annotatedFilterBgzip.bgzipped_file
  }

  call vtt.variantsToTable {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    vcf=annotatedFilterIndex.indexed_vcf,
    vcf_tbi=annotatedFilterIndex.indexed_vcf_tbi,
    fields=variants_to_table_fields,
    genotype_fields=variants_to_table_genotype_fields
  }

  call avftt.addVepFieldsToTable {
    input:
    vcf=annotatedFilterIndex.indexed_vcf,
    vep_fields=vep_to_table_fields,
    tsv=variantsToTable.variants_tsv
  }

  output {
    File mutect_unfiltered_vcf = mutect.unfiltered_vcf
    File mutect_unfiltered_vcf_tbi = mutect.unfiltered_vcf_tbi
    File mutect_filtered_vcf = mutect.filtered_vcf
    File mutect_filtered_vcf_tbi = mutect.filtered_vcf_tbi

    File strelka_unfiltered_vcf = strelka.unfiltered_vcf
    File strelka_unfiltered_vcf_tbi = strelka.unfiltered_vcf_tbi
    File strelka_filtered_vcf = strelka.filtered_vcf
    File strelka_filtered_vcf_tbi = strelka.filtered_vcf_tbi

    File varscan_unfiltered_vcf = varscan.unfiltered_vcf
    File varscan_unfiltered_vcf_tbi = varscan.unfiltered_vcf_tbi
    File varscan_filtered_vcf = varscan.filtered_vcf
    File varscan_filtered_vcf_tbi = varscan.filtered_vcf_tbi

    File vep_summary = annotateVariants.vep_summary

    File final_vcf = vcfReadcount.readcounted_vcf
    File final_vcf_tbi = vcfReadcount.readcounted_vcf_tbi

    File final_filtered_vcf = annotatedFilterIndex.indexed_vcf
    File final_filtered_vcf_tbi = annotatedFilterIndex.indexed_vcf_tbi

    File final_tsv = addVepFieldsToTable.annotated_variants_tsv

    File tumor_snv_bam_readcount_tsv = vcfReadcount.tumor_snv_bam_readcount_tsv
    File tumor_indel_bam_readcount_tsv = vcfReadcount.tumor_indel_bam_readcount_tsv

    File normal_snv_bam_readcount_tsv = vcfReadcount.normal_snv_bam_readcount_tsv
    File normal_indel_bam_readcount_tsv = vcfReadcount.normal_indel_bam_readcount_tsv
  }
}
