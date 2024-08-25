version 1.0

import "../tools/index_vcf.wdl" as iv
import "../tools/bam_readcount.wdl" as br
import "../subworkflows/vcf_readcount_annotator.wdl" as vra

workflow vcfReadcount {
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

    Int? minimum_base_quality
    Int? minimum_mapping_quality

    File vcf
  }

  call br.bamReadcount as tumorBamReadcount {
    input:
    vcf=vcf,
    sample=tumor_sample_name,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=tumor_bam,
    bam_bai=tumor_bam_bai,
    min_base_quality=minimum_base_quality,
    min_mapping_quality=minimum_mapping_quality
  }

  call br.bamReadcount as normalBamReadcount {
    input:
    vcf=vcf,
    sample=normal_sample_name,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=normal_bam,
    bam_bai=normal_bam_bai,
    min_base_quality=minimum_base_quality,
    min_mapping_quality=minimum_mapping_quality
  }

  call vra.vcfReadcountAnnotator as addTumorBamReadcountToVcf {
    input:
    vcf=vcf,
    snv_bam_readcount_tsv=tumorBamReadcount.snv_bam_readcount_tsv,
    indel_bam_readcount_tsv=tumorBamReadcount.indel_bam_readcount_tsv,
    data_type="DNA",
    sample_name=tumor_sample_name
  }

  call vra.vcfReadcountAnnotator as addNormalBamReadcountToVcf {
    input:
    vcf=addTumorBamReadcountToVcf.annotated_bam_readcount_vcf,
    snv_bam_readcount_tsv=normalBamReadcount.snv_bam_readcount_tsv,
    indel_bam_readcount_tsv=normalBamReadcount.indel_bam_readcount_tsv,
    data_type="DNA",
    sample_name=normal_sample_name
  }

  call iv.indexVcf as index {
    input: vcf=addNormalBamReadcountToVcf.annotated_bam_readcount_vcf
  }

  output {
    File readcounted_vcf = index.indexed_vcf
    File readcounted_vcf_tbi = index.indexed_vcf_tbi

    File tumor_snv_bam_readcount_tsv = tumorBamReadcount.snv_bam_readcount_tsv
    File tumor_indel_bam_readcount_tsv = tumorBamReadcount.indel_bam_readcount_tsv

    File normal_snv_bam_readcount_tsv = normalBamReadcount.snv_bam_readcount_tsv
    File normal_indel_bam_readcount_tsv = normalBamReadcount.indel_bam_readcount_tsv
  }
}
