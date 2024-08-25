version 1.0

import "../tools/index_vcf.wdl" as iv
import "../tools/bgzip.wdl" as b
import "../tools/single_sample_docm_filter.wdl" as ssdf
import "../tools/docm_gatk_haplotype_caller.wdl" as dghc

workflow docmGermline {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    File docm_vcf
    File docm_vcf_tbi
    File interval_list
  }

  call dghc.docmGatkHaplotypeCaller as gatkHaplotypecaller {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=bam,
    bam_bai=bam_bai,
    docm_vcf=docm_vcf,
    docm_vcf_tbi=docm_vcf_tbi,
    interval_list=interval_list
  }

  call ssdf.singleSampleDocmFilter as docmFilter {
    input:
    docm_out=gatkHaplotypecaller.docm_raw_variants
  }

  call b.bgzip {
    input:
    file=docmFilter.docm_filter_out
  }

  call iv.indexVcf as index {
    input:
    vcf=bgzip.bgzipped_file
  }

  output {
    File unfiltered_vcf = gatkHaplotypecaller.docm_raw_variants
    File filtered_vcf = index.indexed_vcf
    File filtered_vcf_tbi = index.indexed_vcf_tbi
  }
}
