version 1.0

import "../types.wdl"  # !UnusedImport

import "../subworkflows/gatk_haplotypecaller_iterator.wdl" as ghi
import "../subworkflows/germline_filter_vcf.wdl" as gfv
import "../tools/add_vep_fields_to_table.wdl" as avftt
import "../tools/freemix.wdl" as f
import "../tools/picard_merge_vcfs.wdl" as pmv
import "../tools/staged_rename.wdl" as sr
import "../tools/variants_to_table.wdl" as vtt
import "../tools/vep.wdl" as v


workflow germlineDetectVariants {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bai
    Array[String] gvcf_gq_bands
    Array[Array[String]] intervals
    File verify_bam_id_metrics
    Int? ploidy
    File vep_cache_dir_zip
    String vep_ensembl_assembly
    String vep_ensembl_version
    String vep_ensembl_species
    Array[String] vep_plugins = ["Frameshift", "Wildtype"]
    File? synonyms_file
    Boolean? annotate_coding_only
    Array[VepCustomAnnotation] vep_custom_annotations
    File limit_variant_intervals
    Array[String] variants_to_table_fields = ["CHROM","POS","ID","REF","ALT"]
    Array[String]? variants_to_table_genotype_fields
    Array[String]? vep_to_table_fields
    String final_tsv_prefix = "variants"
    String gnomad_field_name = "gnomADe_AF"  # only change with gnomad_filter annotation
    Float germline_filter_gnomAD_maximum_population_allele_frequency
  }

  call f.freemix {
    input:
    verify_bam_id_metrics=verify_bam_id_metrics
  }

  call ghi.gatkHaplotypecallerIterator as haplotypeCaller {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=bam,
    bai=bai,
    emit_reference_confidence="NONE",
    gvcf_gq_bands=gvcf_gq_bands,
    intervals=intervals,
    contamination_fraction=freemix.out,
    ploidy=ploidy
  }

  call pmv.picardMergeVcfs as mergeVcfs {
    input: vcfs=haplotypeCaller.gvcf
  }

  call v.vep as annotateVariants {
    input:
    vcf=mergeVcfs.merged_vcf,
    cache_dir_zip=vep_cache_dir_zip,
    ensembl_assembly=vep_ensembl_assembly,
    ensembl_version=vep_ensembl_version,
    ensembl_species=vep_ensembl_species,
    synonyms_file=synonyms_file,
    coding_only=annotate_coding_only,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    custom_annotations=vep_custom_annotations,
    plugins=vep_plugins
  }


  call gfv.germlineFilterVcf as filterVcf {
    input:
    annotated_vcf=annotateVariants.annotated_vcf,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    gnomad_field_name=gnomad_field_name,
    germline_filter_gnomAD_maximum_population_allele_frequency=germline_filter_gnomAD_maximum_population_allele_frequency,
    limit_variant_intervals=limit_variant_intervals
  }

  call vtt.variantsToTable as  filteredVariantsToTable {
    input:
    vcf=filterVcf.filtered_vcf,
    vcf_tbi=filterVcf.filtered_vcf_tbi,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    fields=variants_to_table_fields,
    genotype_fields=variants_to_table_genotype_fields
  }

  call avftt.addVepFieldsToTable as filteredAddVepFieldsToTable {
    input:
    vcf=filterVcf.filtered_vcf,
    tsv=filteredVariantsToTable.variants_tsv,
    vep_fields=vep_to_table_fields,
    prefix=final_tsv_prefix
  }

  call sr.stagedRename as setFilteredTsvName {
    input:
    original=filteredAddVepFieldsToTable.annotated_variants_tsv,
    name="annotated.filtered.tsv"
  }

  call vtt.variantsToTable as finalVariantsToTable {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    vcf=filterVcf.final_vcf,
    vcf_tbi=filterVcf.final_vcf_tbi,
    fields=variants_to_table_fields,
    genotype_fields=variants_to_table_genotype_fields
}

  call avftt.addVepFieldsToTable as finalAddVepFieldsToTable {
    input:
    vcf=filterVcf.final_vcf,
    vep_fields=vep_to_table_fields,
    tsv=finalVariantsToTable.variants_tsv,
    prefix=final_tsv_prefix
  }

  call sr.stagedRename as setFinalTsvName {
    input:
    original=finalAddVepFieldsToTable.annotated_variants_tsv,
    name="annotated.filtered.final.tsv"
  }

  output {
    File raw_vcf = mergeVcfs.merged_vcf
    File raw_vcf_tbi = mergeVcfs.merged_vcf_tbi
    File final_vcf = filterVcf.final_vcf
    File final_vcf_tbi = filterVcf.final_vcf_tbi
    File filtered_vcf = filterVcf.filtered_vcf
    File filtered_vcf_tbi = filterVcf.filtered_vcf_tbi
    File vep_summary = annotateVariants.vep_summary
    File final_tsv = setFinalTsvName.replacement
    File filtered_tsv = setFilteredTsvName.replacement
  }
}
