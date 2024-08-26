version 1.0

import "../tools/filter_sv_vcf_size.wdl" as fsvs
import "../tools/filter_sv_vcf_depth.wdl" as fsvd
import "../tools/bgzip.wdl" as b
import "../tools/index_vcf.wdl" as iv

workflow svDepthCallerFilter {
  input {
    Float? deletion_depth
    Float? duplication_depth
    Int? min_sv_size
    String? output_vcf_name
    File sv_vcf
    String vcf_source  # enum ["cnvkit" "cnvnator"]
  }

  call fsvs.filterSvVcfSize as sizeFilter {
    input:
    input_vcf=sv_vcf,
    size_method="min_len",
    sv_size=min_sv_size
  }

  call fsvd.filterSvVcfDepth as depthFilter {
    input:
    input_vcf=sizeFilter.filtered_sv_vcf,
    deletion_depth=deletion_depth,
    duplication_depth=duplication_depth,
    output_vcf_name=output_vcf_name,
    vcf_source=vcf_source
  }

  call b.bgzip as filteredVcfBgzip {
    input: file=depthFilter.filtered_sv_vcf
  }

  call iv.indexVcf as filteredVcfIndex {
    input: vcf=filteredVcfBgzip.bgzipped_file
  }

  output {
    File filtered_vcf=filteredVcfIndex.indexed_vcf
  }
}
