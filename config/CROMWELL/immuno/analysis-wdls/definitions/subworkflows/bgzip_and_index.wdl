version 1.0

import "../tools/bgzip.wdl" as b
import "../tools/index_vcf.wdl" as i

workflow bgzipAndIndex {
  input {
    File vcf
  }

  call b.bgzip {
    input: file=vcf
  }

  call i.indexVcf as index {
    input: vcf=bgzip.bgzipped_file
  }

  output {
    File indexed_vcf = index.indexed_vcf
    File indexed_vcf_tbi = index.indexed_vcf_tbi
  }
}
