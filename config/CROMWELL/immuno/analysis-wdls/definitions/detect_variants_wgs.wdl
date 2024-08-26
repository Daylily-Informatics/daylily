version 1.0

import "subworkflows/docm_cle.wdl" as dc
import "subworkflows/filter_vcf.wdl" as fv
import "subworkflows/mutect.wdl" as m
import "subworkflows/strelka_and_post_processing.wdl" as sapp
import "subworkflows/varscan_pre_and_post_processing.wdl" as vpapp
import "subworkflows/vcf_readcount_annotator.wdl" as vra
import "tools/add_vep_fields_to_table.wdl" as avftt
import "tools/bam_readcount.wdl" as br
import "tools/bgzip.wdl" as b
import "tools/combine_variants_wgs.wdl" as cvw
import "tools/docm_add_variants.wdl" as dav
import "tools/index_vcf.wdl" as iv
import "tools/variants_to_table.wdl" as vtt
import "tools/vep.wdl" as v
import "tools/vt_decompose.wdl" as vd
import "types.wdl"  # !UnusedImport

workflow detectVariantsWgs {
  input {
    File reference
    File reference_fai
    File reference_dict
    File tumor_bam
    File tumor_bam_bai
    File normal_bam
    File normal_bam_bai
    File roi_intervals
    Boolean strelka_exome_mode
    Int strelka_cpu_reserved = 8
    File? strelka_call_regions
    File? strelka_call_regions_tbi
    Int? readcount_minimum_base_quality
    Int? readcount_minimum_mapping_quality
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
    File vep_cache_dir_zip
    String vep_ensembl_assembly
    String vep_ensembl_version
    String vep_ensembl_species
    File? synonyms_file
    Boolean? annotate_coding_only
    String? vep_pick  # enum ["pick", "flag_pick", "pick_allele", "per_gene", "pick_allele_gene", "flag_pick_allele", "flag_pick_allele_gene"]
    Array[String] vep_plugins = ["Frameshift", "Wildtype"]

    String gnomad_field_name = "gnomADe_AF"  # only change with gnomad_filter_annotation
    Float filter_gnomADe_maximum_population_allele_frequency = 0.001
    Float filter_mapq0_threshold = 0.15
    Int filter_minimum_depth = 1
    Float? filter_somatic_llr_threshold
    Float? filter_somatic_llr_tumor_purity
    Float? filter_somatic_llr_normal_contamination_rate
    Boolean cle_vcf_filter = false
    Array[String] variants_to_table_fields = ["CHROM", "POS", "ID", "REF", "ALT", "set", "AC", "AF"]
    Array[String] variants_to_table_genotype_fields = ["GT", "AD"]
    Array[String] vep_to_table_fields = ["HGVSc", "HGVSp"]
    String tumor_sample_name
    String normal_sample_name
    Array[VepCustomAnnotation] vep_custom_annotations
    File? validated_variants
    File? validated_variants_tbi
  }

  call m.mutect {
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
    tumor_sample_name=tumor_sample_name,
    fp_min_var_freq=fp_min_var_freq
  }

  call sapp.strelkaAndPostProcessing as strelka {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,

    tumor_sample_name=tumor_sample_name,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,

    normal_sample_name=normal_sample_name,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,

    interval_list=roi_intervals,
    exome_mode=strelka_exome_mode,
    cpu_reserved=strelka_cpu_reserved,
    call_regions=strelka_call_regions,
    call_regions_tbi=strelka_call_regions_tbi,
    fp_min_var_freq=fp_min_var_freq
  }

  call vpapp.varscanPreAndPostProcessing as varscan {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,

    tumor_sample_name=tumor_sample_name,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,

    normal_sample_name=normal_sample_name,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,

    interval_list=roi_intervals,
    scatter_count=scatter_count,
    strand_filter=varscan_strand_filter,
    min_coverage=varscan_min_coverage,
    varscan_min_var_freq=varscan_min_var_freq,
    p_value=varscan_p_value,
    max_normal_freq=varscan_max_normal_freq
  }

  call dc.docmCle as docm {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,

    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,

    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,

    interval_list=roi_intervals,
    docm_vcf=docm_vcf,
    docm_vcf_tbi=docm_vcf_tbi,
    filter_docm_variants=filter_docm_variants
  }

  call cvw.combineVariantsWgs as combine {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    mutect_vcf=mutect.filtered_vcf,
    mutect_vcf_tbi=mutect.filtered_vcf_tbi,
    strelka_vcf=strelka.filtered_vcf,
    strelka_vcf_tbi=strelka.filtered_vcf_tbi,
    varscan_vcf=varscan.filtered_vcf,
    varscan_vcf_tbi=varscan.filtered_vcf_tbi
  }

  call dav.docmAddVariants {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    callers_vcf=combine.combined_vcf,
    callers_vcf_tbi=combine.combined_vcf_tbi,
    docm_vcf=docm.docm_variants_vcf,
    docm_vcf_tbi=docm.docm_variants_vcf_tbi
  }

  call vd.vtDecompose as decompose {
    input:
    vcf=docmAddVariants.merged_vcf,
    vcf_tbi=docmAddVariants.merged_vcf_tbi
  }

  call iv.indexVcf as decomposeIndex {
    input:
    vcf=decompose.decomposed_vcf
  }

  call v.vep as annotateVariants {
    input:
    vcf=decomposeIndex.indexed_vcf,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    cache_dir_zip=vep_cache_dir_zip,
    ensembl_assembly=vep_ensembl_assembly,
    ensembl_version=vep_ensembl_version,
    ensembl_species=vep_ensembl_species,
    synonyms_file=synonyms_file,
    coding_only=annotate_coding_only,
    custom_annotations=vep_custom_annotations,
    pick=vep_pick,
    plugins=vep_plugins
  }

  call br.bamReadcount as tumorBamReadcount {
    input:
    vcf=annotateVariants.annotated_vcf,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    sample=tumor_sample_name,
    bam=tumor_bam,
    bam_bai=tumor_bam_bai,
    min_base_quality=readcount_minimum_base_quality,
    min_mapping_quality=readcount_minimum_mapping_quality
  }

  call br.bamReadcount as normalBamReadcount {
    input:
    vcf=annotateVariants.annotated_vcf,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    sample=normal_sample_name,
    bam=normal_bam,
    bam_bai=normal_bam_bai,
    min_base_quality=readcount_minimum_base_quality,
    min_mapping_quality=readcount_minimum_mapping_quality
  }

  call vra.vcfReadcountAnnotator as addTumorBamReadcountToVcf {
    input:
    vcf=annotateVariants.annotated_vcf,
    snv_bam_readcount_tsv=tumorBamReadcount.snv_bam_readcount_tsv,
    indel_bam_readcount_tsv=tumorBamReadcount.indel_bam_readcount_tsv,
    data_type="DNA",
    sample_name=tumor_sample_name
  }

  call vra.vcfReadcountAnnotator as addNormalBamReadcountToVcf {
    input:
    vcf=addTumorBamReadcountToVcf.annotated_bam_readcount_vcf,
    snv_bam_readcount_tsv=normalBamReadcount.snv_bam_readcount_tsv,
    indel_bam_readcount_tsv=normalBamReadcount.indel_bam_readcount_tsv,
    data_type="DNA",
    sample_name=normal_sample_name
  }

  call iv.indexVcf as index {
    input: vcf=addNormalBamReadcountToVcf.annotated_bam_readcount_vcf
  }

  call fv.filterVcf {
    input:
    vcf=index.indexed_vcf,
    vcf_tbi=index.indexed_vcf_tbi,
    filter_gnomADe_maximum_population_allele_frequency=filter_gnomADe_maximum_population_allele_frequency,
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
    tumor_sample_name=tumor_sample_name,
    gnomad_field_name=gnomad_field_name,
    validated_variants=validated_variants,
    validated_variants_tbi=validated_variants_tbi
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
    tsv=variantsToTable.variants_tsv,
    vep_fields=vep_to_table_fields
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
    File docm_filtered_vcf = docm.docm_variants_vcf
    File docm_filtered_vcf_tbi = docm.docm_variants_vcf_tbi
    File final_vcf = index.indexed_vcf
    File final_vcf_tbi = index.indexed_vcf_tbi
    File final_filtered_vcf = annotatedFilterIndex.indexed_vcf
    File final_filtered_vcf_tbi = annotatedFilterIndex.indexed_vcf_tbi
    File final_tsv = addVepFieldsToTable.annotated_variants_tsv
    File vep_summary = annotateVariants.vep_summary
    File tumor_snv_bam_readcount_tsv = tumorBamReadcount.snv_bam_readcount_tsv
    File tumor_indel_bam_readcount_tsv = tumorBamReadcount.indel_bam_readcount_tsv
    File normal_snv_bam_readcount_tsv = normalBamReadcount.snv_bam_readcount_tsv
    File normal_indel_bam_readcount_tsv = normalBamReadcount.indel_bam_readcount_tsv
  }
}
