version 1.0

import "../tools/annotate_known_variants.wdl" as akv
import "../tools/filter_vcf_custom_allele_freq.wdl" as fvcaf
import "../tools/filter_vcf_mapq0.wdl" as fvm
import "../tools/filter_vcf_cle.wdl" as fvc
import "../tools/filter_vcf_depth.wdl" as fvd
import "../tools/filter_vcf_somatic_llr.wdl" as fvsl
import "../tools/staged_rename.wdl" as sr

workflow filterVcf {
  input {
    File vcf
    File vcf_tbi
    Float filter_mapq0_threshold
    Float filter_gnomADe_maximum_population_allele_frequency
    String gnomad_field_name
    File tumor_bam
    File tumor_bam_bai
    Boolean do_cle_vcf_filter
    Float? filter_somatic_llr_threshold
    Float? filter_somatic_llr_tumor_purity
    Float? filter_somatic_llr_normal_contamination_rate
    File reference
    File reference_fai
    File reference_dict
    Int filter_minimum_depth
    String tumor_sample_name
    String normal_sample_name
    # both or neither
    File? validated_variants
    File? validated_variants_tbi
  }

  if (defined(validated_variants)) {
    call akv.annotateKnownVariants {
      input:
      vcf=vcf,
      vcf_tbi=vcf_tbi,
      validated_variants=validated_variants,
      validated_variants_tbi=validated_variants_tbi
    }
  }

  call fvcaf.filterVcfCustomAlleleFreq as filterVcfGnomadeAlleleFreq {
    input:
    vcf=select_first([annotateKnownVariants.validated_annotated_vcf, vcf]),
    maximum_population_allele_frequency=filter_gnomADe_maximum_population_allele_frequency,
    field_name=gnomad_field_name
  }

  call fvm.filterVcfMapq0 {
    input:
    vcf=filterVcfGnomadeAlleleFreq.filtered_vcf,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    sample_name=tumor_sample_name,
    threshold=filter_mapq0_threshold
  }

  if (do_cle_vcf_filter) {
    call fvc.filterVcfCle {
      input:
      vcf=filterVcfMapq0.mapq0_filtered_vcf,
      filter=do_cle_vcf_filter
    }
  }

  call fvd.filterVcfDepth {
    input:
    vcf=select_first([filterVcfCle.cle_filtered_vcf, filterVcfMapq0.mapq0_filtered_vcf]),
    minimum_depth=filter_minimum_depth,
    sample_names=[normal_sample_name, tumor_sample_name]
  }

  call fvsl.filterVcfSomaticLlr {
    input:
    vcf=filterVcfDepth.depth_filtered_vcf,
    threshold=filter_somatic_llr_threshold,
    tumor_purity=filter_somatic_llr_tumor_purity,
    normal_contamination_rate=filter_somatic_llr_normal_contamination_rate,
    tumor_sample_name=tumor_sample_name,
    normal_sample_name=normal_sample_name
  }

  call sr.stagedRename as setFinalVcfName {
    input:
    original=filterVcfSomaticLlr.somatic_llr_filtered_vcf,
    name="annotated_filtered.vcf"
  }

  output {
    File filtered_vcf = setFinalVcfName.replacement
  }
}
