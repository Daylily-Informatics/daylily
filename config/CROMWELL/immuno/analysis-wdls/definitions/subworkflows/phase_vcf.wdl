version 1.0

import "../subworkflows/bgzip_and_index.wdl" as bi
import "../tools/index_vcf.wdl" as iv
import "../tools/pvacseq_combine_variants.wdl" as pcv
import "../tools/read_backed_phasing.wdl" as rbp
import "../tools/replace_vcf_sample_name.wdl" as rvsn
import "../tools/select_variants.wdl" as sev
import "../tools/sort_vcf.wdl" as sov

workflow phaseVcf {
  input {
    File somatic_vcf
    File somatic_vcf_tbi
    File germline_vcf
    File reference
    File reference_fai
    File reference_dict
    File bam      # could be .cram
    File bam_bai  # could be .crai
    String normal_sample_name
    String tumor_sample_name
  }

  call rvsn.replaceVcfSampleName as renameGermlineVcf {
    input:
    input_vcf=germline_vcf,
    sample_to_replace=normal_sample_name,
    new_sample_name=tumor_sample_name
  }

  call iv.indexVcf as indexRenamedGermline {
    input: vcf=renameGermlineVcf.renamed_vcf
  }

  call sev.selectVariants as selectSomaticTumorSample {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    vcf=somatic_vcf,
    vcf_tbi=somatic_vcf_tbi,
    output_vcf_basename="somatic_tumor_only",
    samples_to_include=[tumor_sample_name]
  }

  call iv.indexVcf as indexFilteredSomatic {
    input: vcf=selectSomaticTumorSample.filtered_vcf
  }

  call pcv.pvacseqCombineVariants as combineVariants {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    germline_vcf=indexRenamedGermline.indexed_vcf,
    germline_vcf_tbi=indexRenamedGermline.indexed_vcf_tbi,
    somatic_vcf=indexFilteredSomatic.indexed_vcf,
    somatic_vcf_tbi=indexFilteredSomatic.indexed_vcf_tbi
  }

  call sov.sortVcf as sort {
    input:
    vcf=combineVariants.combined_vcf,
    reference_dict=reference_dict
  }

  call bi.bgzipAndIndex {
    input: vcf=sort.sorted_vcf
  }

  call rbp.readBackedPhasing  {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=bam,
    bam_index=bam_bai,
    vcf=bgzipAndIndex.indexed_vcf,
    vcf_tbi=bgzipAndIndex.indexed_vcf_tbi
  }

  call bi.bgzipAndIndex as bgzipAndIndexPhasedVcf {
    input: vcf=readBackedPhasing.phased_vcf
  }

  output {
    File phased_vcf = bgzipAndIndexPhasedVcf.indexed_vcf
    File phased_vcf_tbi = bgzipAndIndexPhasedVcf.indexed_vcf_tbi
  }
}
