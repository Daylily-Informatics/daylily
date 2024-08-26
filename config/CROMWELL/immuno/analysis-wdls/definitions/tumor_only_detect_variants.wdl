version 1.0

import "types.wdl"

import "subworkflows/docm_germline.wdl" as dg
import "subworkflows/varscan_germline.wdl" as vg
import "subworkflows/vcf_readcount_annotator.wdl" as vra

import "tools/add_vep_fields_to_table.wdl" as avftt
import "tools/bam_readcount.wdl" as br
import "tools/bgzip.wdl" as b
import "tools/filter_vcf_coding_variant.wdl" as fvcv
import "tools/filter_vcf_custom_allele_freq.wdl" as fvcaf
import "tools/germline_combine_variants.wdl" as gcv
import "tools/index_vcf.wdl" as iv
import "tools/select_variants.wdl" as sv
import "tools/variants_to_table.wdl" as vtt
import "tools/vep.wdl" as v
import "tools/vt_decompose.wdl" as vd

workflow tumorOnlyDetectVariants {
  input {
    File reference
    File reference_fai
    File reference_dict

    String sample_name
    File bam
    File bam_bai

    String gnomad_field_name = "gnomADe_AF"  # only change with gnomad_filter_annotation
    File roi_intervals

    Int? varscan_strand_filter
    Int? varscan_min_coverage
    Int? varscan_min_reads
    Float? varscan_min_var_freq
    Float? varscan_p_value
    Float maximum_population_allele_frequency = 0.001

    File vep_cache_dir_zip
    String vep_ensembl_assembly
    String vep_ensembl_version
    String vep_ensembl_species
    # one of [pick, flag_pick, pick-allele, per_gene, pick_allele_gene, flag_pick_allele, flag_pick_allele_gene]
    String? vep_pick
    Array[String] vep_plugins = ["Frameshift", "Wildtype"]
    Array[VepCustomAnnotation] vep_custom_annotations

    Boolean annotate_coding_only = false
    File? synonyms_file

    Array[String] variants_to_table_fields = ["CHROM", "POS", "ID", "REF", "ALT", "set", "AC", "AF"]
    Array[String] variants_to_table_genotype_fields = ["GT", "AD"]
    Array[String] vep_to_table_fields = ["HGVSc", "HGVSp"]

    File docm_vcf
    File docm_vcf_tbi

    Int? readcount_minimum_base_quality
    Int? readcount_minimum_mapping_quality
  }

  call vg.varscanGermline as varscan {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    interval_list=roi_intervals,
    strand_filter=varscan_strand_filter,
    min_coverage=varscan_min_coverage,
    varscan_min_var_freq=varscan_min_var_freq,
    min_reads=varscan_min_reads,
    p_value=varscan_p_value,
    sample_name=sample_name
  }

  call dg.docmGermline as docm {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=bam,
    bam_bai=bam_bai,
    interval_list=roi_intervals,
    docm_vcf=docm_vcf,
    docm_vcf_tbi=docm_vcf_tbi
  }

  call gcv.germlineCombineVariants as combineVariants {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    varscan_vcf=varscan.filtered_vcf,
    varscan_vcf_tbi=varscan.filtered_vcf_tbi,
    docm_vcf=docm.filtered_vcf,
    docm_vcf_tbi=docm.filtered_vcf_tbi
  }

  call vd.vtDecompose as decompose {
    input:
    vcf=combineVariants.combined_vcf,
    vcf_tbi=combineVariants.combined_vcf_tbi
  }

  call v.vep as annotateVariants {
    input:
    vcf=decompose.decomposed_vcf,
    cache_dir_zip=vep_cache_dir_zip,
    synonyms_file=synonyms_file,
    coding_only=annotate_coding_only,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    custom_annotations=vep_custom_annotations,
    pick=vep_pick,
    ensembl_assembly=vep_ensembl_assembly,
    ensembl_version=vep_ensembl_version,
    ensembl_species=vep_ensembl_species,
    plugins=vep_plugins
  }

  call br.bamReadcount {
    input:
    vcf=annotateVariants.annotated_vcf,
    sample=sample_name,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=bam,
    bam_bai=bam_bai,
    min_mapping_quality=readcount_minimum_mapping_quality,
    min_base_quality=readcount_minimum_base_quality
  }

  call vra.vcfReadcountAnnotator as addBamReadcountToVcf {
    input:
    vcf=annotateVariants.annotated_vcf,
    snv_bam_readcount_tsv=bamReadcount.snv_bam_readcount_tsv,
    indel_bam_readcount_tsv=bamReadcount.indel_bam_readcount_tsv,
    data_type="DNA",
    sample_name=sample_name
  }

  call iv.indexVcf as index {
    input:
    vcf=addBamReadcountToVcf.annotated_bam_readcount_vcf
  }

  call sv.selectVariants as hardFilter {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    vcf=index.indexed_vcf,
    vcf_tbi=index.indexed_vcf_tbi,
    interval_list=roi_intervals,
    exclude_filtered=true
  }

  call fvcaf.filterVcfCustomAlleleFreq as afFilter {
    input:
    field_name=gnomad_field_name,
    vcf=hardFilter.filtered_vcf,
    maximum_population_allele_frequency=maximum_population_allele_frequency
  }

  call fvcv.filterVcfCodingVariant as codingVariantFilter {
    input:
    vcf=afFilter.filtered_vcf
  }

  call b.bgzip as bgzipFiltered {
    input:
    file=codingVariantFilter.filtered_vcf
  }

  call iv.indexVcf as indexFiltered {
    input:
    vcf=bgzipFiltered.bgzipped_file
  }

  call vtt.variantsToTable {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    vcf=indexFiltered.indexed_vcf,
    vcf_tbi=indexFiltered.indexed_vcf_tbi,
    fields=variants_to_table_fields,
    genotype_fields=variants_to_table_genotype_fields
  }

  call avftt.addVepFieldsToTable {
    input:
    vcf=indexFiltered.indexed_vcf,
    vep_fields=vep_to_table_fields,
    tsv=variantsToTable.variants_tsv
  }

  output {
    File varscan_vcf = varscan.unfiltered_vcf
    File varscan_vcf_tbi = varscan.unfiltered_vcf_tbi

    File docm_gatk_vcf = docm.unfiltered_vcf

    File annotated_vcf = index.indexed_vcf
    File annotated_vcf_tbi = index.indexed_vcf_tbi

    File final_vcf = indexFiltered.indexed_vcf
    File final_vcf_tbi = indexFiltered.indexed_vcf_tbi

    File final_tsv = addVepFieldsToTable.annotated_variants_tsv
    File vep_summary = annotateVariants.vep_summary

    File tumor_snv_bam_readcount_tsv = bamReadcount.snv_bam_readcount_tsv
    File tumor_indel_bam_readcount_tsv = bamReadcount.indel_bam_readcount_tsv
  }
}
