version 1.0

import "../subworkflows/bgzip_and_index.wdl" as bi
import "../subworkflows/fp_filter.wdl" as ff
import "../tools/intervals_to_bed.wdl" as itb
import "../tools/varscan_germline.wdl" as vg

workflow varscanGermline {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    File interval_list
    Int? strand_filter
    Int? min_coverage
    Float? varscan_min_var_freq
    Int? min_reads
    Float? p_value
    String sample_name
  }

  call itb.intervalsToBed {
    input:
    interval_list=interval_list
  }

  call vg.varscanGermline as varscan {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=bam,
    bam_bai=bam_bai,
    roi_bed=intervalsToBed.interval_bed,
    strand_filter=strand_filter,
    min_coverage=min_coverage,
    varscan_min_var_freq=varscan_min_var_freq,
    min_reads=min_reads,
    p_value=p_value,
    sample_name=sample_name
  }

  call bi.bgzipAndIndex {
    input:
    vcf=varscan.variants
  }

  call ff.fpFilter as filter {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=bam,
    bam_bai=bam_bai,
    vcf=bgzipAndIndex.indexed_vcf,
    vcf_tbi=bgzipAndIndex.indexed_vcf_tbi,
    variant_caller="varscan",
    sample_name=sample_name,
    fp_min_var_freq=varscan_min_var_freq
  }

  output {
    File unfiltered_vcf = filter.unfiltered_vcf
    File unfiltered_vcf_tbi = filter.unfiltered_vcf_tbi
    File filtered_vcf = filter.filtered_vcf
    File filtered_vcf_tbi = filter.filtered_vcf_tbi
  }
}
