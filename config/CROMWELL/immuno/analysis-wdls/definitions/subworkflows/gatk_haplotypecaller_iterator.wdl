version 1.0

import "../tools/gatk_haplotype_caller.wdl" as ghc

workflow gatkHaplotypecallerIterator {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bai
    String emit_reference_confidence  # enum [NONE, BP_RESOLUTION, GVCF]
    Array[String] gvcf_gq_bands
    Array[Array[String]] intervals
    File? dbsnp_vcf
    File? dbsnp_vcf_tbi
    File? contamination_fraction
    Int? max_alternate_alleles
    Int? ploidy
    String? read_filter
    String output_prefix = ""
  }


  scatter(interval_sublist in intervals) {
    # TODO(john): if only alphanumeric
    String base = if length(interval_sublist) == 1 then interval_sublist[0] else "output"
    call ghc.gatkHaplotypeCaller as haplotypeCaller {
      input:
      reference=reference,
      reference_fai=reference_fai,
      reference_dict=reference_dict,
      bam=bam,
      bai=bai,
      emit_reference_confidence=emit_reference_confidence,
      gvcf_gq_bands=gvcf_gq_bands,
      intervals=interval_sublist,
      dbsnp_vcf=dbsnp_vcf,
      dbsnp_vcf_tbi=dbsnp_vcf_tbi,
      contamination_fraction=contamination_fraction,
      max_alternate_alleles=max_alternate_alleles,
      ploidy=ploidy,
      read_filter=read_filter,
      output_file_name = output_prefix + base + ".g.vcf.gz"
    }
  }

  output {
    Array[File] gvcf = haplotypeCaller.gvcf
    Array[File] gvcf_tbi = haplotypeCaller.gvcf_tbi
  }
}
