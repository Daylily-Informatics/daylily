version 1.0

import "../types.wdl"  # !UnusedImport

import "../tools/sequence_align_and_tag.wdl" as saat
import "../tools/merge_bams.wdl" as mb
import "../tools/name_sort.wdl" as ns
import "../tools/mark_duplicates_and_sort.wdl" as mdas
import "../tools/index_bam.wdl" as ib

workflow sequenceToBamNonhuman {
  input {
    Array[SequenceData] unaligned
    File reference
    File reference_alt
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_0123
    TrimmingOptions? trimming
    String final_name = "final.bam"
  }

  scatter(sequence in unaligned) {
    call saat.sequenceAlignAndTag as align {
      input:
      unaligned=sequence,
      reference=reference,
      reference_alt=reference_alt,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_0123=reference_0123,
      trimming=trimming
    }
  }

  call mb.mergeBams as merge {
    input:
    bams=align.aligned_bam,
    name=final_name
  }

  call ns.nameSort {
    input:
    bam=merge.merged_bam
  }

  call mdas.markDuplicatesAndSort {
    input:
    bam=nameSort.name_sorted_bam,
    output_name=final_name
  }

  call ib.indexBam {
    input:
    bam=markDuplicatesAndSort.sorted_bam
  }

  output {
    File final_bam = indexBam.indexed_bam
    File final_bam_bai = indexBam.indexed_bam_bai
    File mark_duplicates_metrics_file = markDuplicatesAndSort.metrics_file
  }
}
