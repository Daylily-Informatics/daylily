version 1.0

import "../subworkflows/fp_filter.wdl" as ff
import "../subworkflows/strelka_process_vcf.wdl" as spv
import "../tools/index_vcf.wdl" as iv
import "../tools/merge_vcf.wdl" as mv
import "../tools/replace_vcf_sample_name.wdl" as rvsn
import "../tools/select_variants.wdl" as sv
import "../tools/strelka.wdl" as s

workflow strelkaAndPostProcessing {
  input {
    File reference
    File reference_fai
    File reference_dict

    String tumor_sample_name
    File tumor_bam
    File tumor_bam_bai

    String normal_sample_name
    File normal_bam
    File normal_bam_bai

    File interval_list

    Boolean exome_mode
    Int cpu_reserved = 8

    File? call_regions
    File? call_regions_tbi

    Float? fp_min_var_freq
  }

  call s.strelka {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    exome_mode=exome_mode,
    cpu_reserved=cpu_reserved,
    call_regions=call_regions,
    call_regions_tbi=call_regions_tbi
  }

  scatter(vcf in [strelka.snvs, strelka.indels]) {
    call spv.strelkaProcessVcf as process {
      input: vcf=vcf
    }
  }

  call mv.mergeVcf as merge {
    input:
    vcfs=process.processed_vcf,
    vcf_tbis=process.processed_vcf_tbi
  }

  call rvsn.replaceVcfSampleName as renameTumorSample {
    input:
    input_vcf=merge.merged_vcf,
    sample_to_replace="TUMOR",
    new_sample_name=tumor_sample_name
  }

  call rvsn.replaceVcfSampleName as renameNormalSample {
    input:
    input_vcf=renameTumorSample.renamed_vcf,
    sample_to_replace="NORMAL",
    new_sample_name=normal_sample_name
  }

  call iv.indexVcf as indexFull {
    input: vcf=renameNormalSample.renamed_vcf
  }

  call sv.selectVariants as regionFilter {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    vcf=indexFull.indexed_vcf,
    vcf_tbi=indexFull.indexed_vcf_tbi,
    interval_list=interval_list
  }

  call ff.fpFilter as filter {
    input:
    reference=reference,
    reference_dict=reference_dict,
    reference_fai=reference_fai,
    bam=tumor_bam,
    bam_bai=tumor_bam_bai,
    vcf=regionFilter.filtered_vcf,
    vcf_tbi=regionFilter.filtered_vcf_tbi,
    sample_name=tumor_sample_name,
    variant_caller="strelka",
    fp_min_var_freq=fp_min_var_freq
  }

  output {
    File unfiltered_vcf = filter.unfiltered_vcf
    File unfiltered_vcf_tbi = filter.unfiltered_vcf_tbi

    File filtered_vcf = filter.filtered_vcf
    File filtered_vcf_tbi = filter.filtered_vcf_tbi
  }
}
