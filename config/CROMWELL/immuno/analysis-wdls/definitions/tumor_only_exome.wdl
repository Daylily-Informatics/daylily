version 1.0

import "types.wdl"

import "alignment_exome.wdl" as ae
import "tumor_only_detect_variants.wdl" as todv

import "tools/bam_to_cram.wdl" as btc
import "tools/index_cram.wdl" as ic
import "tools/interval_list_expand.wdl" as ile


workflow tumorOnlyExome {
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
    File bait_intervals
    File target_intervals
    Array[LabelledFile] per_base_intervals
    Array[LabelledFile] per_target_intervals
    Array[LabelledFile] summary_intervals
    File omni_vcf
    File omni_vcf_tbi
    String picard_metric_accumulation_level
    Int target_interval_padding = 100
    Int? varscan_strand_filter
    Int? varscan_min_coverage
    Float? varscan_min_var_freq
    Float? varscan_p_value
    Int? varscan_min_reads
    Float? fp_min_var_freq
    Float maximum_population_allele_frequency = 0.001

    File vep_cache_dir_zip
    String vep_ensembl_assembly
    String vep_ensembl_version
    String vep_ensembl_species
    File? synonyms_file
    String vep_pick  # enum ["pick", "flag_pick", "pick_allele", "per_gene", "pick_allele_gene", "flag_pick_allele", "flag_pick_allele_gene"]
    Array[String] variants_to_table_fields = ["CHROM", "POS", "REF", "ALT", "set"]
    Array[String] variants_to_table_genotype_fields = ["GT", "AD", "AF", "DP"]
    Array[String] vep_to_table_fields = ["Consequence", "SYMBOL", "Feature_type", "Feature", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "HGNC_ID", "Existing_variation", "gnomADe_AF", "CLIN_SIG", "SOMATIC", "PHENO"]

    String sample_name
    File docm_vcf
    File docm_vcf_tbi
    Array[VepCustomAnnotation] vep_custom_annotations
    Int? qc_minimum_mapping_quality
    Int? qc_minimum_base_quality
    Int? readcount_minimum_mapping_quality
    Int? readcount_minimum_base_quality
  }

  call ae.alignmentExome as alignmentAndQc {
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
    qc_minimum_base_quality=qc_minimum_base_quality
  }

  call ile.intervalListExpand as padTargetIntervals {
    input:
    interval_list=target_intervals,
    roi_padding=target_interval_padding
  }

  call todv.tumorOnlyDetectVariants as detectVariants {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=alignmentAndQc.bam,
    bam_bai=alignmentAndQc.bam_bai,
    roi_intervals=padTargetIntervals.expanded_interval_list,
    varscan_strand_filter=varscan_strand_filter,
    varscan_min_coverage=varscan_min_coverage,
    varscan_min_var_freq=varscan_min_var_freq,
    varscan_p_value=varscan_p_value,
    varscan_min_reads=varscan_min_reads,
    maximum_population_allele_frequency=maximum_population_allele_frequency,
    vep_cache_dir_zip=vep_cache_dir_zip,
    vep_ensembl_assembly=vep_ensembl_assembly,
    vep_ensembl_version=vep_ensembl_version,
    vep_ensembl_species=vep_ensembl_species,
    synonyms_file=synonyms_file,
    vep_pick=vep_pick,
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
    input:
    cram=bamToCram.cram
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
