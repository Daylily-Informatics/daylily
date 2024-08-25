version 1.0

import "../types.wdl"

import "../tools/aligned_seq_fda_stats.wdl" as aligned_seq_fda_stats
import "../tools/fastqc.wdl" as fastqc_tool
import "../tools/generate_fda_tables.wdl" as generate_fda_tables
import "../tools/md5sum.wdl" as md5sum_tool
import "../tools/unaligned_seq_fda_stats.wdl" as unaligned_seq_fda_stats

struct FdaMetrics {
  Array[File] fastqc_data
  Array[File] table_metrics
  File md5sums
  File table
}

workflow generateFdaMetricsForBamOrFastqs {
  input {
    Array[SequenceData] unaligned_sequence
    Array[File] aligned_data

    File reference
    File reference_index

    File? alignment_summary_metrics
    File? duplication_metrics
    File? insert_size_metrics
    File? hs_metrics
    File? rna_metrics
    File? flagstat
    File? unaligned_rna_table

    String? sequencing_platform
    String? sequencing_instrument
    String? sequencing_kit
    String? sequencing_type
    String? single_or_paired_end
    String? spike_in_error_rate
    String? total_dna
    String? reference_genome
    String? total_rna
    String? rin_score
    String? freq_normalization_method
    String? annotation_file
    String? sample_name

    String table_file_name
    String table_num
    String? table_source
    String stats_output_name
    String md5sum_output_name

    Boolean is_aligned
  }

  scatter(seq in unaligned_sequence) {
    #yields an Array[Array[File]] outside of the scatter
    Array[File] unaligned_data = select_all([seq.sequence.bam, seq.sequence.fastq1, seq.sequence.fastq2])
  }

  Array[File] bam_or_fastqs = flatten(select_all([flatten(unaligned_data), aligned_data]))
  call fastqc_tool.fastqc {
    input:
      files = bam_or_fastqs
  }
 
  if(is_aligned) {
    call aligned_seq_fda_stats.alignedSeqFdaStats {
      input:
        reference = reference,
        reference_index = reference_index,
        files = aligned_data,
        output_name = stats_output_name
    }
  }
  if(!is_aligned) { #else, but no else in WDL 1.0
    scatter(idx in range(length(unaligned_data))) {
      call unaligned_seq_fda_stats.unalignedSeqFdaStats {
        input:
          files = unaligned_data[idx],
          output_name = stats_output_name,
          suffix = idx + 1
      }
    }
  }

  call md5sum_tool.md5sum {
    input:
      files = bam_or_fastqs,
      output_name = md5sum_output_name
  }
 
  call generate_fda_tables.generateFdaTables {
    input:
      table_file_name = table_file_name,
      table_num = table_num,
      source = table_source,
      md5sum_file = md5sum.md5_file,
      fastqc_zips = fastqc.fastqc_data,
      unaligned_metrics = unalignedSeqFdaStats.unaligned_stats,
      aligned_metrics = alignedSeqFdaStats.aligned_stats,
      sample_name = sample_name,
      sequencing_platform = sequencing_platform,
      sequencing_instrument = sequencing_instrument,
      sequencing_kit = sequencing_kit,
      single_or_paired_end = single_or_paired_end,
      sequencing_type = sequencing_type,
      spike_in_error_rate = spike_in_error_rate,
      alignment_summary_metrics = alignment_summary_metrics,
      duplication_metrics = duplication_metrics,
      insert_size_metrics = insert_size_metrics,
      hs_metrics = hs_metrics,
      flagstat = flagstat,
      total_dna = total_dna,
      reference_genome = reference_genome,
      total_rna = total_rna,
      rin_score = rin_score,
      freq_normalization_method = freq_normalization_method,
      annotation_file = annotation_file,
      rna_metrics = rna_metrics,
      unaligned_rna_table = unaligned_rna_table
  }

  output {
    FdaMetrics metrics = object {
      fastqc_data: fastqc.fastqc_data,
      table_metrics: select_all(select_first([unalignedSeqFdaStats.unaligned_stats, [alignedSeqFdaStats.aligned_stats]])),
      md5sums: md5sum.md5_file,
      table: generateFdaTables.table
    }
  }
}
