version 1.0

import "../tools/bam_to_fastq.wdl" as btf
import "../tools/trim_fastq.wdl" as tf
import "../tools/biscuit_align.wdl" as ba
import "../tools/index_bam.wdl" as ib
import "../tools/biscuit_markdup.wdl" as bm

workflow bamToTrimmedFastqAndBiscuitAlignments {
  input {
    File bam
    File adapters
    String adapter_trim_end
    Int adapter_min_overlap
    Int max_uncalled
    Int min_readlength
    String read_group_id
    File reference_index
  }

  call btf.bamToFastq {
    input: bam=bam
  }

  call tf.trimFastq {
    input:
    reads1=bamToFastq.fastq1,
    reads2=bamToFastq.fastq2,
    adapters=adapters,
    adapter_trim_end=adapter_trim_end,
    adapter_min_overlap=adapter_min_overlap,
    max_uncalled=max_uncalled,
    min_readlength=min_readlength
  }

  call ba.biscuitAlign {
    input:
    reference_index=reference_index,
    fastq1=trimFastq.fastq1,
    fastq2=trimFastq.fastq2,
    read_group_id=read_group_id
  }

  call ib.indexBam {
    input: bam=biscuitAlign.aligned_bam
  }

  call bm.biscuitMarkdup {
    input:
    bam=indexBam.indexed_bam,
    bam_bai=indexBam.indexed_bam_bai
  }

  output {
    File aligned_bam = biscuitMarkdup.markdup_bam
  }
}
