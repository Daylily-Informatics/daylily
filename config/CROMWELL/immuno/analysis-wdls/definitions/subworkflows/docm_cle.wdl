version 1.0

import "../tools/bgzip.wdl" as b
import "../tools/docm_gatk_haplotype_caller.wdl" as dghc
import "../tools/filter_vcf_docm.wdl" as fvd
import "../tools/index_vcf.wdl" as iv
import "../tools/vt_decompose.wdl" as vd

workflow docmCle {
  input {
    File reference
    File reference_fai
    File reference_dict

    File tumor_bam
    File tumor_bam_bai

    File normal_bam
    File normal_bam_bai

    File docm_vcf
    File docm_vcf_tbi

    File interval_list

    Boolean filter_docm_variants
  }

  call dghc.docmGatkHaplotypeCaller as gatkHaplotypeCaller {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=tumor_bam,
    bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    interval_list=interval_list,
    docm_vcf=docm_vcf,
    docm_vcf_tbi=docm_vcf_tbi
  }

  call b.bgzip {
    input: file=gatkHaplotypeCaller.docm_raw_variants
  }

  call iv.indexVcf as index {
    input: vcf=bgzip.bgzipped_file
  }

  call vd.vtDecompose as decompose {
    input:
    vcf=index.indexed_vcf,
    vcf_tbi=index.indexed_vcf_tbi
  }

  call fvd.filterVcfDocm as docmFilter {
    input:
    docm_raw_variants=decompose.decomposed_vcf,
    normal_bam=normal_bam,
    tumor_bam=tumor_bam,
    filter_docm_variants=filter_docm_variants
  }

  call b.bgzip as bgzip2 {
    input: file=docmFilter.docm_filtered_variants
  }

  call iv.indexVcf as index2 {
    input: vcf=bgzip2.bgzipped_file
  }

  output {
    File docm_variants_vcf = index2.indexed_vcf
    File docm_variants_vcf_tbi = index2.indexed_vcf_tbi
  }
}
