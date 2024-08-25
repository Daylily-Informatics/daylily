version 1.0

import "../tools/pindel.wdl" as p
import "../tools/cat_out.wdl" as co

workflow pindelCat {
  input {
    File reference
    File reference_fai
    File reference_dict
    File tumor_bam
    File tumor_bam_bai
    File normal_bam
    File normal_bam_bai
    File region_file
    Int insert_size = 400
    String tumor_sample_name
    String normal_sample_name
  }

  call p.pindel {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    region_file=region_file,
    insert_size=insert_size,
    tumor_sample_name=tumor_sample_name,
    normal_sample_name=normal_sample_name
  }

  call co.catOut as cat {
    input: pindel_outs=[pindel.deletions, pindel.insertions, pindel.tandems, pindel.long_insertions, pindel.inversions]
  }

  output {
    File per_region_pindel_out = cat.pindel_out
  }
}
