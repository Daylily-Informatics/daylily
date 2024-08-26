version 1.0

import "../tools/filter_vcf_coding_variant.wdl" as fvcv
import "../tools/filter_vcf_custom_allele_freq.wdl" as fvcaf
import "../tools/staged_rename.wdl" as sr
import "../tools/bgzip.wdl" as b
import "../tools/index_vcf.wdl" as iv
import "../tools/select_variants.wdl" as sv

workflow germlineFilterVcf {
  input {
    File annotated_vcf
    Float germline_filter_gnomAD_maximum_population_allele_frequency
    String gnomad_field_name
    File limit_variant_intervals
    File reference
    File reference_fai
    File reference_dict
  }

  call fvcv.filterVcfCodingVariant as codingVariantFilter {
    input: vcf=annotated_vcf
  }

  call fvcaf.filterVcfCustomAlleleFreq as gnomadFrequencyFilter {
    input:
    vcf=codingVariantFilter.filtered_vcf,
    maximum_population_allele_frequency=germline_filter_gnomAD_maximum_population_allele_frequency,
    field_name=gnomad_field_name
  }

  call sr.stagedRename as setFilteredVcfName {
    input:
    original=gnomadFrequencyFilter.filtered_vcf,
    name="annotated.filtered.vcf"
  }

  call b.bgzip as bgzipFilteredVcf {
    input: file=setFilteredVcfName.replacement
  }

  call iv.indexVcf as indexFilteredVcf {
    input: vcf=bgzipFilteredVcf.bgzipped_file
  }

  call sv.selectVariants as limitVariants {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    vcf=indexFilteredVcf.indexed_vcf,
    vcf_tbi=indexFilteredVcf.indexed_vcf_tbi,
    interval_list=limit_variant_intervals,
    exclude_filtered=true,
    output_vcf_basename="annotated.filtered.final"
  }

  output {
    File filtered_vcf = indexFilteredVcf.indexed_vcf
    File filtered_vcf_tbi = indexFilteredVcf.indexed_vcf_tbi
    File final_vcf = limitVariants.filtered_vcf
    File final_vcf_tbi = limitVariants.filtered_vcf_tbi
  }
}
