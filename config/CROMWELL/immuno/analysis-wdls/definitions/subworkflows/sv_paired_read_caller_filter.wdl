version 1.0

import "../tools/filter_sv_vcf_read_support.wdl" as fsvrs
import "../tools/duphold.wdl" as d
import "../tools/filter_sv_vcf_depth.wdl" as fsvd
import "../tools/bgzip.wdl" as b
import "../tools/index_vcf.wdl" as iv

workflow svPairedReadCallerFilter {
  input {
    Float? abundance_percentage
    File bam
    Float? deletion_depth
    Float? duplication_depth
    String? output_vcf_name
    File reference
    File reference_fai
    File reference_dict
    File? snps_vcf
    Int? sv_paired_count
    Int? sv_split_count
    File sv_vcf
    String vcf_source  # enum ["manta", "smoove"]
  }

  call fsvrs.filterSvVcfReadSupport as readSupportFilter {
    input:
    abundance_percentage=abundance_percentage,
    input_vcf=sv_vcf,
    paired_count=sv_paired_count,
    split_count=sv_split_count,
    vcf_source=vcf_source
  }

  call d.duphold as dupholdAnnotate {
    input:
    bam=bam,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    snps_vcf=snps_vcf,
    sv_vcf=readSupportFilter.filtered_sv_vcf
  }

  call fsvd.filterSvVcfDepth as depthFilter {
    input:
    input_vcf=dupholdAnnotate.annotated_sv_vcf,
    deletion_depth=deletion_depth,
    duplication_depth=duplication_depth,
    output_vcf_name=output_vcf_name,
    vcf_source="duphold"
  }

  call b.bgzip as filteredVcfBgzip {
    input: file=depthFilter.filtered_sv_vcf
  }

  call iv.indexVcf as filteredVcfIndex {
    input: vcf=filteredVcfBgzip.bgzipped_file
  }

  output {
    File filtered_vcf = filteredVcfIndex.indexed_vcf
    File filtered_vcf_tbi = filteredVcfIndex.indexed_vcf_tbi
  }
}
