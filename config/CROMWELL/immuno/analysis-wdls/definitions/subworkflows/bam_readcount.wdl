version 1.0

import "../tools/normalize_variants.wdl" as nv
import "../tools/vt_decompose.wdl" as vd
import "../tools/bam_readcount.wdl" as br

workflow bamReadcount {
  input {
    File vcf
    File vcf_tbi
    String sample
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bam_bai
    Int? min_base_quality
    Int? min_mapping_quality
  }

  call nv.normalizeVariants {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    vcf=vcf,
    vcf_tbi=vcf_tbi
  }

  call vd.vtDecompose as decomposeVariants {
    input:
    vcf=normalizeVariants.normalized_vcf,
    vcf_tbi=normalizeVariants.normalized_vcf_tbi
  }

  call br.bamReadcount as readcount {
    input:
    vcf=decomposeVariants.decomposed_vcf,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=bam,
    bam_bai=bam_bai,
    sample=sample,
    min_mapping_quality=min_mapping_quality,
    min_base_quality=min_base_quality
  }

  output {
    File normalized_vcf = decomposeVariants.decomposed_vcf
    File snv_bam_readcount_tsv = readcount.snv_bam_readcount_tsv
    File indel_bam_readcount_tsv = readcount.indel_bam_readcount_tsv
  }
}
