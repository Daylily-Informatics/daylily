version 1.0

import "../tools/add_strelka_gt.wdl" as asg
import "../tools/bgzip.wdl" as b
import "../tools/index_vcf.wdl" as iv

workflow strelkaProcessVcf {
  input {
    File vcf
  }

  call asg.addStrelkaGt as addGt {
    input: vcf=vcf
  }

  call b.bgzip {
    input: file=addGt.processed_vcf
  }

  call iv.indexVcf as index {
    input: vcf=bgzip.bgzipped_file
  }

  output {
    File processed_vcf = index.indexed_vcf
    File processed_vcf_tbi = index.indexed_vcf_tbi
  }
}
