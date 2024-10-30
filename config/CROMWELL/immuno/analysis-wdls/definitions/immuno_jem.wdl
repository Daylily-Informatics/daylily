version 1.0


# pipelines
import "somatic_exome.wdl" as se
# others
import "types.wdl"  # !UnusedImport


#
# These structs are needed only because MiniWDL, used by some of our
# scripts, has issues with Object types. We avoid this problem by
# enforcing struct boundaries to help the parser understand what we're
# doing.
#
# The reason to use these instead of primitives in outputs is
# to encode intended output directory for each file.
#




struct Variants {
  Array[File?] mutect
  Array[File?] strelka
  Array[File?] varscan
  Array[File?] docm
}

struct Cnv {
  Array[File?] cnvkit
}

struct Sv {
  Array[File?] manta
}

struct Somatic {
  Variants variants
  Array[File] final
  Cnv cnv
  Sv sv
}

struct Germline {
  Array[File?] variants
}

workflow immuno {
  input {

    # --------- RNAseq Inputs ------------------------------------------

    # File reference_annotation
    # Array[SequenceData] rna_sequence
    # String sample_name

    # File trimming_adapters
    # String trimming_adapter_trim_end
    # Int trimming_adapter_min_overlap
    # Int trimming_max_uncalled
    # Int trimming_min_readlength

    # File kallisto_index
    # File gene_transcript_lookup_table
    # String? strand  # [first, second, unstranded]
    # File refFlat
    # File ribosomal_intervals
    # File star_fusion_genome_dir_zip
    # Boolean examine_coding_effect = true
    # String? fusioninspector_mode  # enum [inspect validate]
    # File cdna_fasta
    # File agfusion_database
    # Boolean agfusion_annotate_noncanonical = true
    # Float? min_ffpm_level = 0.05

    # --------- Somatic Exome Inputs -----------------------------------

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
    String tumor_sample_name
    Array[SequenceData] tumor_sequence

    String normal_name = "normal"
    String normal_sample_name
    Array[SequenceData] normal_sequence

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
    Int scatter_count = 50

    Int? varscan_strand_filter
    Int? varscan_min_coverage
    Float? varscan_min_var_freq
    Float? varscan_p_value
    Float? varscan_max_normal_freq

    Float? fp_min_var_freq

    File docm_vcf
    File docm_vcf_tbi

    Boolean filter_docm_variants = true

    String? gnomad_field_name
    Float? filter_gnomADe_maximum_population_allele_frequency

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

    Array[String] vep_to_table_fields = ["HGVSc", "HGVSp"]
    Array[String] variants_to_table_genotype_fields = ["GT", "AD"]
    Array[String] variants_to_table_fields = ["CHROM", "POS", "ID", "REF", "ALT", "set", "AC", "AF"]
    Array[VepCustomAnnotation] vep_custom_annotations

    File? manta_call_regions
    File? manta_call_regions_tbi
    Boolean manta_non_wgs = true
    Boolean? manta_output_contigs

    File somalier_vcf
    File? validated_variants
    File? validated_variants_tbi

    # --------- Germline Inputs ----------------------------------------

    Array[String] gvcf_gq_bands
    Array[Array[String]] gatk_haplotypecaller_intervals
    Int? ploidy
    String? optitype_name
    Float germline_filter_gnomAD_maximum_population_allele_frequency = 1.1
  }

  call geht.germlineExomeHlaTyping as germlineExome {
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
    target_interval_padding=target_interval_padding,
    per_base_intervals=per_base_intervals,
    per_target_intervals=per_target_intervals,
    summary_intervals=summary_intervals,
    omni_vcf=omni_vcf,
    omni_vcf_tbi=omni_vcf_tbi,
    picard_metric_accumulation_level=picard_metric_accumulation_level,
    gvcf_gq_bands=gvcf_gq_bands,
    intervals=gatk_haplotypecaller_intervals,
    ploidy=ploidy,
    vep_cache_dir_zip=vep_cache_dir_zip,
    vep_ensembl_assembly=vep_ensembl_assembly,
    vep_ensembl_version=vep_ensembl_version,
    vep_ensembl_species=vep_ensembl_species,
    vep_custom_annotations=vep_custom_annotations,
    synonyms_file=synonyms_file,
    annotate_coding_only=annotate_coding_only,
    qc_minimum_mapping_quality=qc_minimum_mapping_quality,
    qc_minimum_base_quality=qc_minimum_base_quality,
    optitype_name="optitype_normal",
    germline_filter_gnomAD_maximum_population_allele_frequency=germline_filter_gnomAD_maximum_population_allele_frequency
  }


  output {

    # ---------- Germline Outputs --------------------------------------

    Germline germline = object {
      variants: [
        germlineExome.final_vcf,
        germlineExome.final_vcf_tbi,
        germlineExome.filtered_vcf,
        germlineExome.filtered_vcf_tbi,
        germlineExome.vep_summary
      ]
    }

    Array[File] hla_typing = flatten([
      [germlineExome.optitype_tsv,
       germlineExome.optitype_plot,
       optitype.optitype_tsv,
       optitype.optitype_plot,
       germlineExome.phlat_summary,
       phlat.phlat_summary,
       extractAlleles.allele_file,
       hlaConsensus.consensus_alleles_file],
      hlaConsensus.hla_call_files
    ])


    # --------- Other Outputs ------------------------------------------


  }
}
