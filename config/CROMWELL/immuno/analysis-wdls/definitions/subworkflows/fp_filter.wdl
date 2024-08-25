version 1.0

import "../tools/bgzip.wdl" as b
import "../tools/fp_filter.wdl" as ff
import "../tools/index_vcf.wdl" as iv
import "../tools/normalize_variants.wdl" as nv
import "../tools/vcf_sanitize.wdl" as vs
import "../tools/vt_decompose.wdl" as vd
import "../tools/select_variants.wdl" as sv

workflow fpFilter {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    File vcf
    File vcf_tbi
    String variant_caller
    String? sample_name
    Float? fp_min_var_freq
  }

  call vs.vcfSanitize as sanitizeVcf {
    input: vcf=vcf
  }

  call nv.normalizeVariants {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    vcf=sanitizeVcf.sanitized_vcf,
    vcf_tbi=sanitizeVcf.sanitized_vcf_tbi
  }

  call vd.vtDecompose as decomposeVariants {
    input:
    vcf=normalizeVariants.normalized_vcf,
    vcf_tbi=normalizeVariants.normalized_vcf_tbi
  }

  call iv.indexVcf as index {
    input: vcf=decomposeVariants.decomposed_vcf
  }

  call ff.fpFilter as firstFilter {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=bam,
    vcf=index.indexed_vcf,
    sample_name=sample_name,
    fp_min_var_freq=fp_min_var_freq,
    output_vcf_basename = variant_caller + "_full"
  }

  call b.bgzip as fpBgzip {
    input: file=firstFilter.filtered_vcf
  }

  call iv.indexVcf as fpIndex {
    input: vcf=fpBgzip.bgzipped_file
  }

  call sv.selectVariants as hardFilter {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    vcf=fpIndex.indexed_vcf,
    vcf_tbi=fpIndex.indexed_vcf_tbi,
    exclude_filtered=true,
    output_vcf_basename = variant_caller + "_filtered"
  }

  output {
    File unfiltered_vcf = fpIndex.indexed_vcf
    File unfiltered_vcf_tbi = fpIndex.indexed_vcf_tbi
    File filtered_vcf = hardFilter.filtered_vcf
    File filtered_vcf_tbi = hardFilter.filtered_vcf_tbi
  }
}
