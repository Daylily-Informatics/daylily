version 1.0

import "../tools/filter_vcf_mapq0.wdl" as fvm
import "../tools/filter_vcf_cle.wdl" as fvc
import "../tools/filter_vcf_depth.wdl" as fvd
import "../tools/filter_vcf_somatic_llr.wdl" as fvsl
import "../tools/staged_rename.wdl" as sr

workflow filterVcfNonhuman {
  input {
    File vcf
    Float filter_mapq0_threshold
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
  }

  call fvm.filterVcfMapq0 {
    input:
    vcf=vcf,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    sample_name=tumor_sample_name,
    threshold=filter_mapq0_threshold
  }

  call fvc.filterVcfCle {
    input:
    vcf=filterVcfMapq0.mapq0_filtered_vcf,
    filter=do_cle_vcf_filter
  }

  call fvd.filterVcfDepth {
    input:
    vcf=filterVcfCle.cle_filtered_vcf,
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
