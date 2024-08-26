version 1.0

import "../tools/trim_fastq.wdl" as tf
import "../tools/sequence_to_fastq.wdl" as stf
import "../tools/hisat2_align.wdl" as ha
import "../types.wdl"  # !UnusedImport

workflow sequenceToTrimmedFastqAndHisatAlignments {
  input {
    SequenceData unaligned
    File adapters
    String adapter_trim_end
    Int adapter_min_overlap
    Int max_uncalled
    Int min_readlength
    String read_group_id
    Array[String] read_group_fields
    File reference_index
    File reference_index_1ht2
    File reference_index_2ht2
    File reference_index_3ht2
    File reference_index_4ht2
    File reference_index_5ht2
    File reference_index_6ht2
    File reference_index_7ht2
    File reference_index_8ht2
    String? strand  # [first, second, unstranded]
    Boolean? unzip_fastqs
  }

  call stf.sequenceToFastq {
    input:
    bam=unaligned.sequence.bam,
    fastq1=unaligned.sequence.fastq1,
    fastq2=unaligned.sequence.fastq2,
    unzip_fastqs=unzip_fastqs
  }

  call tf.trimFastq {
    input:
    reads1=sequenceToFastq.read1_fastq,
    reads2=sequenceToFastq.read2_fastq,
    adapters=adapters,
    adapter_trim_end=adapter_trim_end,
    adapter_min_overlap=adapter_min_overlap,
    max_uncalled=max_uncalled,
    min_readlength=min_readlength
  }

  call ha.hisat2Align {
    input:
    reference_index=reference_index,
    reference_index_1ht2=reference_index_1ht2,
    reference_index_2ht2=reference_index_2ht2,
    reference_index_3ht2=reference_index_3ht2,
    reference_index_4ht2=reference_index_4ht2,
    reference_index_5ht2=reference_index_5ht2,
    reference_index_6ht2=reference_index_6ht2,
    reference_index_7ht2=reference_index_7ht2,
    reference_index_8ht2=reference_index_8ht2,
    fastq1=trimFastq.fastqs[0],
    fastq2=trimFastq.fastqs[1],
    read_group_id=read_group_id,
    read_group_fields=read_group_fields,
    strand=strand
  }

  output {
    Array[File] fastqs = trimFastq.fastqs
    File aligned_bam = hisat2Align.aligned_bam
  }
}
